#include "thimble_lattice.hpp"
using std::abs;
using std::exp;
using std::log;


//simple dirac delta function
template <class T>
double dd(T i, T j)
{
  double Delta(0);
  if (i == j)
  {
    Delta = 1;
  }
  return Delta;
}

//*************************************thimble_system***************************

thimble_system::thimble_system(int x_dim, int t_dim, double flow_time, long unsigned int seed) : 
Nx(x_dim), 
Nt(t_dim), 
tau(flow_time), 
Npath(2*Nt), 
Nrpath(Npath - 4), 
Ntot(Nrpath*Nx),
rng_seed(seed),
J(0,0),
J_conj(0,0),
rel_path(""),
file_name("0"),
j(0,1),
dx(1.),
dt(0.75),
generator(rng_seed),
uniform_double(0, 1),
abcd(0, 1),
acceptance_rate(0.)
{ 
  //determining the number of timesteps for the ODE solvers from the flow time  
  h = 0.1; //sets the base size 
  number_of_timesteps = int(ceil(tau/h)); //calculates how many steps this corresponds to (overshooting in the case of it not being exact)
  h = tau/number_of_timesteps; //readusting the size of h to prevent overshooting
}

thimble_system::~thimble_system()
{
  //empty destructor
}
void thimble_system::add_scalar_field()
{
  //adds a scalar field to the vector of scalar fields
  scalars.emplace_back(scalar_field(Nx, Nt, dt, dx, this));
}

void thimble_system::add_scalar_field(double mass)
{
  //adds a scalar field and a mass at the same time, otherwise same as above
  scalars.emplace_back(scalar_field(Nx, Nt, dt, dx, this));
  scalars[scalars.size() - 1].set_mass(mass);
}

void thimble_system::add_interaction(double coupling, std::vector<int> powers)
{
  //adds an interaction to the list of all possible interactions
  interactions.emplace_back(interaction(coupling, powers));
}

void thimble_system::add_interaction(double coupling, int powers)
{
  std::vector<int> powers_vector(scalars.size(), powers);
  add_interaction(coupling, powers_vector);
}

void thimble_system::set_path(std::string new_path)
{
  //this lets you set the path for saving the data to a different directory.
  std::string boost_path = new_path;
  boost_path = boost_path.erase(boost_path.size() - 1); //removes the trailing / from the path before passing to the boost library
  boost::filesystem::path p(boost_path);
  if (!exists(p))
  {
    create_directory(p); //creates the directory if it doesn't yet exist
  }
  rel_path = new_path;
}

dcomp thimble_system::calc_dS(uint site, uint field, uint field_type)
{
  dcomp dS = 0; //derivative of the action
  dcomp interaction = 0;
  dS = scalars[field].free_action_derivative(site, field_type); //free field kinetic contribution
  for (int i = 0; i < interactions.size(); ++i)
  {
    //looping through the first derivatives of all the interactions (derivatives with respect to this field)
    interaction += interactions[i].first_derivative(site, field, *this, field_type); 
  }
  interaction *= dx*j*(scalars[field].path[scalars[field].calc_n(site)] + scalars[field].path_offset[scalars[field].calc_n(site)])/2.; //(delta_n + delta_n-1) factor
  dS += interaction;
  return dS;
}

dcomp thimble_system::calc_dS(uint site, uint field_type)
{
  int field = 0;
  int internal_site = site; //this essentially takes the busy work out of calculating which field the Jacobian is dealing with
  if (site > Ntot - 1)
  {
    while (internal_site > Ntot - 1)
    {
      internal_site -= Ntot;
      ++field;
    }
  }
  return calc_dS(internal_site, field, field_type);
}

dcomp thimble_system::calc_ddS(uint site_1, uint site_2, uint field_1, uint field_2, uint field_type)
{
  dcomp ddS = 0;
  dcomp interaction = 0;
  if (field_1 == field_2)
  {
    ddS = scalars[field_1].free_action_second_derivative(site_1, site_2); //free component
  }
  
  if (site_1 == site_2) //only add the interactions which are on equal sites
  {
    for (int i = 0; i < interactions.size(); ++i)
    {
      interaction += interactions[i].second_derivative(site_1, field_1, field_2, *this, field_type);
    }
    interaction *= -1.*j*dx*(scalars[field_1].path[scalars[field_1].calc_n(site_1)] + scalars[field_1].path_offset[scalars[field_1].calc_n(site_1)])/2.;
  }
  ddS += interaction;
  return ddS;
}

dcomp thimble_system::calc_ddS(uint site_1, uint site_2, uint field_type)
{
  uint field_1 = 0;
  uint field_2 = 0;
  dcomp ddS;

  while(site_1 > Ntot - 1)
  {
    site_1 -= Ntot;
    ++field_1;
  }
  
  while(site_2 > Ntot - 1)
  {
    site_2 -= Ntot;
    ++field_2;
  }

  ddS = calc_ddS(site_1, site_2, field_1, field_2, field_type);
  return ddS;
}

field_id_return thimble_system::calc_field(int master_site)
{ //returns a struct with the field ID and site as members
  field_id_return field_setup;
  field_setup.field_number = 0;
  while (master_site > Ntot)
  {
    master_site -= Ntot;
    ++field_setup.field_number;
  }
  field_setup.site_number = master_site;
  return field_setup;
}

void thimble_system::sync_ajustment(dcomp ajustment[])
{ //forces the intermediate fields in the ode solver back into the relevant fields
  for(int i = 0; i < scalars.size(); ++i)
  {
    for(int r = 0; r < Ntot; ++r)
    {
      scalars[i].fields[4][r] = ajustment[i*Ntot + r];
    }
  }
}

matrix<dcomp> thimble_system::calc_jacobian(bool proposal)
{
  //function calculates the Jacobian and it's determinant from either the proposed or orignal fields
  dcomp* working_scalar = new dcomp[Njac];
  dcomp* ajustment_scalar = new dcomp[Njac];
  dcomp* k1_scalar = new dcomp[Njac]; //RK45 variables for the scalar fields
  dcomp* k2_scalar = new dcomp[Njac];
  dcomp* k3_scalar = new dcomp[Njac];
  dcomp* k4_scalar = new dcomp[Njac];

  dcomp* ajustment_jac = new dcomp[NjacSquared];
  dcomp* k1_jac = new dcomp[NjacSquared];
  dcomp* k2_jac = new dcomp[NjacSquared];
  dcomp* k3_jac = new dcomp[NjacSquared];
  dcomp* k4_jac = new dcomp[NjacSquared];
  int ajustment = 4;

  matrix<dcomp> Jac(Njac, Njac);
  dcomp jac_element;

  //identifying if it's the proposal or exising fields we wish to flow
  int proposal_or = 2;
  if (proposal)
  {
    proposal_or = 3;
  }
  for (uint i = 0; i < scalars.size(); ++i)
  {
    for (uint k = 0; k < Ntot; ++k)
    {
      //resetting the flowed field to the base field. At this point we're at tau = 0 so the flowed field and base field should be identical
      scalars[i].fields[proposal_or - 2][k] = scalars[i].fields[proposal_or][k];
      working_scalar[i*Ntot + k] = scalars[i].fields[proposal_or][k];
    }
  }

  //setting up an identity matrix
  for(uint r = 0; r < Njac; ++r)
  {
    for(uint c = 0; c < Njac; ++c)
    {
      if (r == c)
      {
        Jac.set_element(r, c, 1);
      }
      else
      {
        Jac.set_element(r, c, 0);
      }
    }
  }
  //standard implementation of RK45 for an autonomous system
  for (uint i = 0; i < number_of_timesteps; ++i)
  {
    for (uint r = 0; r < Njac; ++r)
    {
      k1_scalar[r] = h*conj(calc_dS(r, proposal_or - 2));
      ajustment_scalar[r] = working_scalar[r] + k1_scalar[r]/2.;
      for (int c = 0; c < Njac; ++c)
      {
        k1_jac[r + Njac*c] = 0.;
        for (int s = 0; s < Njac; ++s)
        {
          k1_jac[r + Njac*c] += h*conj(calc_ddS(r, s, proposal_or - 2)*Jac.get_element(s, c));
        }
        ajustment_jac[r + Njac*c] = Jac.get_element(r, c) + k1_jac[r + Njac*c]/2.;
      }
    }
    sync_ajustment(ajustment_scalar);
    //ajustment scalar does the job of holding field values that are calculated intermittenlty in the RK45 method
    //sync_ajustment pushes the values stored in the array in this function back out to the scalarfields, which then use the value to calculate the ds and dds functions

    for (uint r = 0; r < Njac; ++r)
    {
      k2_scalar[r] = h*conj(calc_dS(r, ajustment));
      for (uint c = 0; c < Njac; ++c)
      {
        k2_jac[r + Njac*c] = 0;
        for(uint s = 0; s < Njac; ++s)
        {
          k2_jac[r + Njac*c] += h*conj(calc_ddS(r, s, ajustment)*ajustment_jac[s + Njac*c]);
        }
      }
    }

    for (uint r = 0; r < Njac; ++r)
    {
      ajustment_scalar[r] = working_scalar[r] + k2_scalar[r]/2.;
      for (uint c = 0; c < Njac; ++c)
      {
        ajustment_jac[r + Njac*c] = Jac.get_element(r,c) + k2_jac[r + Njac*c]/2.;
      }
    }
    sync_ajustment(ajustment_scalar);

    for (uint r = 0; r < Njac; ++r)
    {
      k3_scalar[r] = h*conj(calc_dS(r, ajustment));
      for (uint c = 0; c < Njac; ++c)
      {
        k3_jac[r + Njac*c] = 0;
        for (uint s = 0; s < Njac; ++s)
        {
          k3_jac[r + Njac*c] += h*conj(calc_ddS(r, s, ajustment)*ajustment_jac[s + Njac*c]);
        }
      }
    }

    for (uint r = 0; r < Njac; ++r)
    {
      ajustment_scalar[r] = working_scalar[r] + k3_scalar[r];
      for (uint c = 0; c < Njac; ++c)
      {
        ajustment_jac[r + Njac*c] = Jac.get_element(r, c) + k3_jac[r + Njac*c];
      }
    }
    sync_ajustment(ajustment_scalar);

    for (uint r = 0; r < Njac; ++r)
    {
      k4_scalar[r] = h*conj(calc_dS(r, ajustment));
      for (uint c = 0; c < Njac; ++c)
      {
        k4_jac[r + Njac*c] = 0;
        for(uint s = 0; s < Njac; ++s)
        {
          k4_jac[r + Njac*c] += h*conj(calc_ddS(r, s, ajustment)*ajustment_jac[s + Njac*c]);
        }
      }
    }
    for (uint r = 0; r < Njac; ++r)
    {
      working_scalar[r] += (k1_scalar[r] + 2.*k2_scalar[r] + 2.*k3_scalar[r] + k4_scalar[r])/6.;
      for (uint c = 0; c < Njac; ++c)
      {
        jac_element = Jac.get_element(r, c);
        jac_element += (k1_jac[r + Njac*c] + 2.*k2_jac[r + Njac*c] + 2.*k3_jac[r + Njac*c] + k4_jac[r + Njac*c])/6.;
        Jac.set_element(r, c, jac_element);
      }
    }
    //returning the flowed fields to the scalar fields object
    for (uint i = 0; i < scalars.size(); ++i)
    {
      for (uint k = 0; k < Ntot; ++k)
      {
        scalars[i].fields[proposal_or - 2][k] = working_scalar[i*Ntot + k];
      }
    }
  }
  
  delete[] working_scalar;
  delete[] ajustment_scalar;
  delete[] k1_scalar;
  delete[] k2_scalar;
  delete[] k3_scalar;
  delete[] k4_scalar;
  delete[] ajustment_jac;
  delete[] k1_jac;
  delete[] k2_jac;
  delete[] k3_jac;
  delete[] k4_jac;

  return Jac;
}

dcomp thimble_system::calc_S(uint field_type)
{
  dcomp S = 0;
  int n;
  //looping through all the fields and sites to add their free free field contributions
  for (uint i = 0; i < scalars.size(); ++i)
  {
    for (uint k = 0; k < Ntot; ++k)
    {
      S += scalars[i].free_action(k, field_type);
    }
  }
  //looping through all the interactions to add them to the action too
  for (uint i = 0; i < interactions.size(); ++i)
  {
    for (uint k = 0; k < Ntot; ++k)
    {
      n = scalars[0].calc_n(k);
      S += -1.*j*dx*(scalars[0].path[n] + scalars[0].path_offset[n])*interactions[i].base(k, *this, field_type)/2.;
    }
  }
  return S;
}

int thimble_system::update()
{
  bool proposal = true;
  dcomp proposed_action, proposed_detJ;
  mydouble log_proposal;
  mydouble matrix_exponenet, proposed_matrix_exponenet, exponenet, check;
  int output = 0; //this is the return value
  
  int field_update_id = field_choice(generator);;
  matrix<dcomp> Delta = sweep_field_proposal(field_update_id);

  //creating new basefield condtions
  for(uint i = 0; i < scalars.size(); ++i)
  {
    for(uint k = 0; k < Ntot; ++k)
    {
      scalars[i].fields[3][k] = scalars[i].fields[2][k] + Delta.get_element(k + i*Ntot, 0);
    }
  }
  
  //calculating the Jacobian, it's determinant, conjugate, and the action of the proposed field state
  matrix<dcomp> proposed_J = calc_jacobian(proposal);
  matrix<dcomp> proposed_J_conj = proposed_J.conjugate();
  proposed_action = calc_S(1);
  log_proposal = (mydouble) log(abs(proposed_J.get_det()));

  //matrix multiplication required to calculate the accpetance exponenet
  matrix<dcomp> Delta_transpose = Delta.transpose();

  matrix_exponenet = (mydouble) real((Delta_transpose*J*J_conj*Delta).get_element(0, 0));
  proposed_matrix_exponenet = (mydouble) real((Delta_transpose*proposed_J*proposed_J_conj*Delta).get_element(0,0)); //hacky solution
  //exponenet for the MC test
  exponenet = ((mydouble) real(S - proposed_action)) + 2.*log_proposal - 2.*((mydouble) log(real(J.get_det()))) + matrix_exponenet/pow(scalars[field_update_id].delta, 2) - proposed_matrix_exponenet/pow(scalars[field_update_id].delta, 2);
  check = uniform_double(generator);
  if (exp(exponenet) > check)
  {
    //proposal was accepted, transfering all the proposed parameters to the storage
    S = proposed_action;
    J = proposed_J;
    J_conj = proposed_J_conj;
    //transfering over all the scalar fields
    for (uint i = 0; i < scalars.size(); ++i)
    {
      for (uint k = 0; k < Ntot; ++k)
      {
        scalars[i].fields[0][k] = scalars[i].fields[1][k];
        scalars[i].fields[2][k] = scalars[i].fields[3][k];
      }
    }
    output = 1;
  }
  return output;
}

matrix<dcomp> thimble_system::sweep_proposal()
{
  //this function generates a proposal vector Delta, based on sweep updating, such that the entire lattice is updated by roughly the same amount.
  dcomp* eta = new dcomp[Njac];
  uint field_id;
  for(uint i = 0; i < Njac; ++i)
  {
    //setting up the proposal on the complex manifold
    field_id = i % Ntot;
    eta[i] = scalars[field_id].proposal(generator) + j*scalars[field_id].proposal(generator);
  }

  //returning the proposal to the real manifold
  matrix<dcomp> Delta = J.solve(eta);
  for(uint i = 0; i < Njac; ++i)
  {
    //taking only the elements that fit in the reduced space
    Delta.set_element(i, 0, real(Delta.get_element(i, 0)));
  }
  delete[] eta;

  return Delta;
}

matrix<dcomp> thimble_system::sweep_field_proposal(int field_choice)
{
  //this function generates a proposal vector Delta, based on sweep updating, such that the entire lattice is updated by roughly the same amount.
  dcomp* eta = new dcomp[Njac];
  for(uint i = 0; i < Njac; ++i)
  {
    //setting up the proposal on the complex manifold
    eta[i] = 0;
  }

  for(uint i = (Ntot*field_choice); i < (Ntot*(field_choice + 1)); ++i)
  {
    eta[i] = scalars[field_choice].proposal(generator) + j*scalars[field_choice].proposal(generator);
  }

  //returning the proposal to the real manifold
  matrix<dcomp> Delta = J.solve(eta);
  for(uint i = 0; i < Njac; ++i)
  {
    //taking only the elements that fit in the reduced space
    Delta.set_element(i, 0, real(Delta.get_element(i, 0)));
  }
  delete[] eta;
  return Delta;
}

matrix<dcomp> thimble_system::site_proposal()
{
  int target_site = uniform_int(generator);
  int target_field = target_site % Ntot;
  dcomp proposal = scalars[target_field].proposal(generator) + j*scalars[target_field].proposal(generator);
  dcomp* eta = new dcomp[Njac];
  for (uint i = 0; i < Njac; ++i)
  {
    eta[i] = 0;
  }
  eta[target_site] = proposal;
  matrix<dcomp> Delta = J.solve(eta);
  for(uint i = 0; i < Njac; ++i)
  {
    Delta.set_element(i, 0, real(Delta.get_element(i, 0)));
  }

  delete[] eta;
  return Delta;
}

void thimble_system::simulate(uint n_burn_in, uint n_simulation, uint n_existing)
{
  pre_simulation_check();
  Njac = scalars.size()*Ntot; //setting the size of the Jacobian "matricies"
  NjacSquared = pow(Njac, 2);
  Nsys = Ntot*scalars.size(); //total number of sites in the system

  for(uint i = 0; i < scalars.size(); ++i)
  {
    //this gives each field it's own gaussian to draw values from, with the appropreate sigma.
    std::normal_distribution<double> redo(0, scalars[i].sigma);
    scalars[i].proposal = redo;
  }
  
  std::uniform_int_distribution<int> redo_2(0, Nsys - 1);
  //this resets the site selection system, informing it of the existence of all the fields
  uniform_int = redo_2;

  std::uniform_int_distribution<int> redo_3(0, scalars.size() - 1);
  field_choice = redo_3;

  dcomp* state_storage = new dcomp[(Njac + 2)*n_simulation]; //This stores the data between updates and will be saved to a file
  std::ofstream data_storage;
  std::ofstream aux_storage;

  double *a = new double[Nx];
  double *b = new double[Nx];
  double *c = new double[Nx];
  double *d = new double[Nx];
  if(n_existing == 0)
  {
    //initialising the fields
    for (uint i = 0; i < scalars.size(); ++i)
    {
      for (uint k = 0; k < Nx; ++k)
      {
        a[k] = abcd(generator);
        b[k] = abcd(generator);
        c[k] = abcd(generator);
        d[k] = abcd(generator);
      }
      scalars[i].initialise(a, b, c, d);
    }
  }
  J.resize(Njac, Njac);
  J_conj.resize(Njac, Njac);
  J = calc_jacobian();
  J_conj = J.conjugate();
  S = calc_S(0);
  //setup is now complete, the Jacobian, it's conjugate, and it's determinant have been calculated, and the scalars are primed.
  for (uint i = 0; i < n_burn_in; ++i)
  {
    update();
  }

  for(uint i = 0; i < n_simulation; ++i)
  {
    acceptance_rate += update(); //calculating the acceptance rate from the return value of the update
    //storing the values of the fields, the action, and the Jacobian in a unified vector to be written to file later
    for (int k = 0; k < scalars.size(); ++k)
    {
      for (int r = 0; r < Ntot; ++r)
      {
        //writing each scalar into the array
        state_storage[i*(Njac + 2) + k*Ntot + r] = scalars[k].fields[0][r];
      }
    }
    state_storage[i*(Njac + 2) + scalars.size()*Ntot] = S;
    state_storage[i*(Njac + 2) + scalars.size()*Ntot + 1] = log(real(J.get_det())) +j*arg(J.get_det());
    //adding auxiliary parameters, action and the log determinant
  }
  acceptance_rate /= (n_simulation + n_existing);
  if (n_existing == 0)
  {
    //opens the file for header writing if there aren't any existing records
    data_storage.open(rel_path + file_name);
    //saving initial conditions, random seed, and simulation parameters
    data_storage << rng_seed << ",";
    for(int i = 0; i < scalars.size(); ++i)
    {
      for(int k = 0; k < Nx; ++k)
      {
        data_storage << real(scalars[i].field_0[k]) << "," << imag(scalars[i].field_0[k]) << ",";
      }
    }
    //saving the data previously stored
    for (uint i = 0; i < scalars.size(); ++i)
    {
      for (uint k = 0; k < Nx; ++k)
      {
        data_storage << real(scalars[i].field_1[k]) << "," << imag(scalars[i].field_1[k]) << ",";
      }
    }
    for (int i = 0; i < scalars.size(); ++i)
    {
      data_storage << scalars[i].delta << ",";
    }
    data_storage << tau << "," << acceptance_rate << std::endl;
  }
  else
  {
    //opens the file for appending if there are existing records 
    data_storage.open(rel_path + file_name, std::ios_base::app);
  }
  for(uint i = 0; i < n_simulation; ++i)
  {
    for (uint k = 0; k < Njac + 1; ++k)
    {
      //split the complex data into real and imaginary components
      data_storage << real(state_storage[i*(Njac + 2) + k]) << "," << imag(state_storage[i*(Njac + 2) + k]) << ",";
    }
    data_storage << real(state_storage[i*(Njac + 2) + Njac + 1]) << "," << imag(state_storage[i*(Njac + 2) + Njac + 1]) << std::endl;
  }
  data_storage.close();

  //this bit stores the conditions on the real manifold for the last iteration before the simulation closes. This can then be used to recover the state 
  aux_storage.open(rel_path + file_name + "_aux");
  for (uint i = 0; i < scalars.size(); ++i)
  {
    for (uint k = 0; k < Ntot; ++k)
    {
      aux_storage << real(scalars[i].fields[2][k]) << ",";
    }
  }
  aux_storage << n_simulation + n_existing;
  aux_storage.close();
  delete[] state_storage;
  delete[] a;
  delete[] b;
  delete[] c;
  delete[] d;
}

void thimble_system::set_field_mass(int field_number, double new_mass)
{
  //external interface for the user to set a new mass for the field
  scalars[field_number].set_mass(new_mass);
}

void thimble_system::set_dt(double new_dt)
{
  //external interface to let the user resize the lattice
  for(uint i = 0; i < scalars.size(); ++i)
  {
    scalars[i].set_dt(new_dt);
  }
  dt = new_dt;
}

void thimble_system::set_dx(double new_dx)
{
  //external interface to let the user resize the lattice (space direction this time)
  for(uint i = 0; i < scalars.size(); ++i)
  {
    scalars[i].set_dx(new_dx);
  }
  dx = new_dx;
}

void thimble_system::set_occupation_number(int field_number, int new_occupation_number)
{
  //Allows the user to set a uniform occupation number for a field
  scalars[field_number].set_occupation_number(new_occupation_number);
}

void thimble_system::set_occupation_number(int field_number, int new_occupation_number[])
{
  //Allows the user to set a non-uniform occupation number for a given field
  scalars[field_number].set_occupation_number(new_occupation_number);
}

void thimble_system::set_occupation_number(int field_number, std::vector<int> new_occupation_number)
{
  //thin wrapper around the above function, allowing it to take either a native C style array or a more C++ style vector.
  int* internal = new int[Nx];
  for(int i = 0; i < Nx; ++i)
  {
    internal[i] = new_occupation_number[i];
  }
  scalars[field_number].set_occupation_number(internal);
  delete[] internal;
}

void thimble_system::set_proposal_size(int field_number, double new_delta)
{
  //Alows the user to specify the proposal step size for the simulation.
  scalars[field_number].delta = new_delta;
  scalars[field_number].sigma = new_delta/pow(2, 0.5);
}

void thimble_system::pre_simulation_check()
{ //this checks all the masses are in the correct range for the simulation

  std::vector<uint> test_failed; //vector of all fields that failed the check
  double test_condition, new_dt, possible_new_dt;
  new_dt = pow(10, 100); //stupidly large number to start with
  for(uint i = 0; i < scalars.size(); ++i)
  {
    if(Nx == 1)
    {
      test_condition = scalars[i].m;
      //this takes into account that if it's a 1D simulation, the momentum plays no role
    }
    else
    {
      test_condition = sqrt(4/pow(dx,2) + scalars[i].squareMass);
      //this is the condition where the momentum does have a role, the first term is the maximum lattice momentum
    }
    if (test_condition > (2/dt))
    {
      test_failed.push_back(i);
      possible_new_dt = 2./test_condition;
      if(new_dt > possible_new_dt)
      {
        new_dt = possible_new_dt; 
      }
    }
  }

  if(test_failed.size() != 0)
  {
    int number_of_new_sites = int(ceil((dt/new_dt)*Nt));
    printf("Resolution insufficient for the mass scales involved. The following fields have too high mass: \n");
    for(uint i = 0; i < test_failed.size(); ++i)
    {
      printf("Scalar field %i \n", test_failed[i]);
    }
    printf("This can be solved by reducing the temporal link size to at least %f. Doing this will require at least %i time sites to maintain the existing range. \n", new_dt, number_of_new_sites);
    exit(11);
  }
}

void thimble_system::restart(std::string data_path, std::string aux_path, int n_new_simulation)
{
  std::ifstream file(data_path);
  std::vector<std::string> header;
  std::string line;
  getline(file, line);
  boost::algorithm::split(header, line, boost::is_any_of(","));
  for (uint i = 0; i < scalars.size(); ++i)
  {
    for (uint k = 0; k < Nx; ++k)
    {
      //loading in the initial time data from the header
      scalars[i].field_0[k] = std::stod(header[2*k + i*Nx + 1]) + j*std::stod(header[2*k + i*Nx + 2]);
      scalars[i].field_1[k] = std::stod(header[2*k + i*Nx + 1 + scalars.size()*Nx]) + j*std::stod(header[2*k + i*Nx + 2 + scalars.size()*Nx]);
      scalars[i].field_2[k] = -1.*pow(dt, 2) + 2.*scalars[i].field_1[k] - scalars[i].field_0[k];
    }
    //loading in the flow constants
    scalars[i].calculate_C();
  }
  
  std::ifstream aux_file(aux_path);
  std::vector<std::string> field_state;
  getline(aux_file, line);
  boost::algorithm::split(field_state, line, boost::is_any_of(","));
  acceptance_rate = std::stod(header[0])*std::stod(field_state[scalars.size()*Ntot]);
  for(uint i = 0; i < scalars.size(); ++i)
  {
    for(uint k = 0; k < Ntot; ++k)
    {
      scalars[i].fields[2][k] = std::stod(field_state[i*Ntot + k]);
    }
  }
  simulate(0, n_new_simulation, int(std::stod(field_state[scalars.size()*Ntot])));
}

void thimble_system::test()
{
  
  for (int i = 0; i < Ntot; ++i)
  {
    printf("C[4][%i] = %f%+fi \n", i, std::real(scalars[0].C[4][i]), std::imag(scalars[0].C[4][i]));
  }
  
  
  for (int i = 0; i < Ntot; ++i)
  {
    std::cout << i << "\t" << scalars[0].path[i] << std::endl;
  }

  for (int i = 0; i < Nx; ++i)
  {
    printf("field_0[%i] = %f%+fi \n", i, std::real(scalars[0].field_0[i]), std::imag(scalars[0].field_0[i]));
  }

  for (int i = 0; i < Nx; ++i)
  {
    printf("field_1[%i] = %f%+fi \n", i, std::real(scalars[0].field_1[i]), std::imag(scalars[0].field_1[i]));
  }
  for (int i = 0; i < Nx; ++i)
  {
    printf("field_2[%i] = %f%+fi \n", i, std::real(scalars[0].field_2[i]), std::imag(scalars[0].field_2[i]));
  }
  
}