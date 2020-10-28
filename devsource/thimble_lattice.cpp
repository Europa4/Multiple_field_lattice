#include "thimble_lattice.hpp"
using std::abs;
using std::exp;
using std::log;
using namespace boost::numeric::odeint;


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

void ode_handler::operator() (const std::vector<dcomp> &x, std::vector<dcomp> &dx, const double t)
{
  for(uint i = 0; i < sys.Njac + sys.NjacSquared; ++i)
  {
    dx[i] = 0;
  }

  for (uint i = 0; i < sys.Njac; ++i)
  {
    dx[i] = conj(sys.calc_dS(i, x));
  }
  for (uint r = 0; r < sys.Njac; ++r)
  {
    for (uint c = 0; c < sys.Njac; ++c)
    {
      for (uint s = 0; s < sys.Njac; ++s)
      {
        dx[sys.Njac + r + c*sys.Njac] += conj(sys.calc_ddS(r, s, x)*x[s + c*sys.Njac + sys.Njac]);
      }
    }
  }
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
acceptance_rate(0.),
rand_n(0)
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
  interaction *= -1.*dx*j*(scalars[field].path[scalars[field].calc_n(site)] + scalars[field].path_offset[scalars[field].calc_n(site)])/2.; //(delta_n + delta_n-1) factor
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

dcomp thimble_system::calc_dS(uint site, const std::vector<dcomp>& field)
{
  //used only by the ODE solver
  for (int i = 0; i < scalars.size(); ++i)
  {
    for (int k = 0; k < Ntot; ++k)
    {
      scalars[i].fields[4][k] = field[i*Ntot + k];
    }
  }

  return calc_dS(site, 4);
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

dcomp thimble_system::calc_ddS(uint site_1, uint site_2, const std::vector<dcomp>& field)
{
  //used only by the ODE solver
  for (int i = 0; i < scalars.size(); ++i)
  {
    for (int k = 0; k < Ntot; ++k)
    {
      scalars[i].fields[4][k] = field[i*Ntot + k];
    }
  }
  return calc_ddS(site_1, site_2, 4);
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

matrix<dcomp> thimble_system::calc_jacobian(bool proposal, int n_iteration)
{
  //function calculates the Jacobian and it's determinant from either the proposed or orignal fields

  int ajustment = 4;

  matrix<dcomp> Jac(Njac, Njac);
  dcomp jac_element;

  std::vector<dcomp> vec(Njac + NjacSquared);
  dcomp test;
  //identifying if it's the proposal or exising fields we wish to flow
  proposal_or = 2;
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
      vec[i*Ntot + k] = scalars[i].fields[proposal_or][k];
    }
  }

  //setting up an identity matrix
  for(uint r = 0; r < Njac; ++r)
  {
    for(uint c = 0; c < Njac; ++c)
    {
      if (r == c)
      {
        vec[Njac + r + Njac*c] = 1.;
      }
      else
      {
        vec[Njac + r + Njac*c] = 0.;
      }
    }
  }
  //standard implementation of RK45 for an autonomous system
  ode_handler ode_sys(*this);
  //integrate_adaptive(make_controlled<runge_kutta_cash_karp54<std::vector<dcomp>>>(1.e-4, 1.e-4), ode_sys, vec, 0.0, tau, h);
  //integrate_const(stepper, ode_sys, vec, 0., tau, h);
  flow(vec);
  for (int r = 0; r < Njac; ++r)
  {
    for (int c = 0; c < Njac; ++c)
    {
      Jac.set_element(r, c, vec[r + c*Njac + Njac]);
    }
  }

  for (int i = 0; i < scalars.size(); ++i)
  {
    for (int k = 0; k < Ntot; ++k)
    {
      scalars[i].fields[proposal_or - 2][k] = vec[k + i*Ntot];
    }
  }

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

  
  int field_update_id = field_choice(generator);
  matrix<dcomp> Delta = sweep_field_proposal(field_update_id);

  //creating new basefield condtions
  for(uint k = 0; k < Ntot; ++k)
  {
    scalars[field_update_id].fields[3][k] = scalars[field_update_id].fields[2][k] + Delta.get_element(k + field_update_id*Ntot, 0);
  }

  
  //calculating the Jacobian, it's determinant, conjugate, and the action of the proposed field state
  flow_check = true;
  matrix<dcomp> proposed_J = calc_jacobian(proposal);
  if(flow_check)
  {
    matrix<dcomp> proposed_J_conj = proposed_J.conjugate();
    proposed_action = calc_S(1);
    log_proposal = (mydouble) log(abs(proposed_J.get_det()));


    //matrix multiplication required to calculate the accpetance exponenet
    matrix<dcomp> Delta_transpose = Delta.transpose();

    matrix_exponenet = (mydouble) real((Delta_transpose*J*J_conj*Delta).get_element(0, 0));
    proposed_matrix_exponenet = (mydouble) real((Delta_transpose*proposed_J*proposed_J_conj*Delta).get_element(0,0)); //hacky solution
    //exponenet for the MC test
    exponenet = ((mydouble) real(S - proposed_action)) + 2.*log_proposal - 2.*((mydouble) log(abs(J.get_det()))) + matrix_exponenet/pow(scalars[field_update_id].delta, 2) - proposed_matrix_exponenet/pow(scalars[field_update_id].delta, 2);
    
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
  for (int i = 0; i < NjacSquared; ++i)
  {
    //printf("J[%i] = %f%+fi\n", i, std::real(J.storage[i]), std::imag(J.storage[i]));
  }
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
  propogate();
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
    ;
    //acceptance_rate += update(); //calculating the acceptance rate from the return value of the update
    acceptance_rate += update();
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
    state_storage[i*(Njac + 2) + scalars.size()*Ntot + 1] = log(abs(J.get_det())) +j*arg(J.get_det());
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

double thimble_system::p_rand()
{
  double r_2 = 1.414213562;
  ++rand_n;
  return std::cos(r_2*rand_n);
}

void thimble_system::print_field(int field_id, int field_type)
{
  for(int i = 0; i < Ntot; ++i)
  {
    printf("Field[%i] = %f%+fi \n", i, std::real(scalars[field_id].fields[field_type][i]), std::imag(scalars[field_id].fields[field_type][i]));
  }
}

void thimble_system::propogate()
{
  //Alright this code is grotty as hell, but since it's called at the start of the simulation and not during, it doesn't have to be as well optimised as the stuff in the simulation routines.
  //It propogates the field's initial values as per the "classical" solution, taking into account interactions
  for (int k = 0; k < Nx; ++k)
  {
    for (int f = 0; f < scalars.size(); ++f)
    {
      scalars[f].fields[2][0 + Nrpath*k] = scalars[f].field_1[k];
      scalars[f].fields[2][1 + Nrpath*k] = -1.0*dt*dt*(scalars[f].squareMass*scalars[f].fields[2][0 + Nrpath*k]) + 2.0*scalars[f].fields[2][0 + Nrpath*k] - scalars[f].field_0[k];
    }
    for (int i = 1; i < (int) (Nrpath/2); ++i)
    {
      for (int f = 0; f < scalars.size(); ++f)
      {
        scalars[f].fields[2][i + Nrpath*k + 1] = -dt*dt*scalars[f].squareMass*scalars[f].fields[2][i + Nrpath*k] + 2.0*scalars[f].fields[2][i + Nrpath*k] - scalars[f].fields[2][i + Nrpath*k - 1];
        for (int I = 0; I < interactions.size(); ++I)
        {
          scalars[f].fields[2][i + Nrpath*k + 1] -= dt*dt*interactions[I].first_derivative(i + Nrpath*k, f, *this, 2); 
        }
      }
    }
  }
/*
//sets up the return leg of the contour
  for(uint f = 0; f < scalars.size(); ++f)
  {
    for (uint k = 0; k < Nx; ++k)
    {
      for (uint n = 0; n < Nt - 2; ++n)
      {
        scalars[f].fields[2][Nrpath*(k + 1) - n - 1] = scalars[f].fields[2][Nrpath*k + n + 1];
      }
    }
  }
  */

  //now that the propogation is complete, we need to clear the classical data from the first site, and calculate C for each field.
  for (int f = 0; f < scalars.size(); ++f)
  {
    for (int k = 0; k < Nx; ++k)
    {
      scalars[f].fields[2][k*Nrpath] = 0;
    }

    for (int i = 0; i < Ntot; ++i)
    {
      //copying the field from the base into the full
      scalars[f].fields[0][i] = scalars[f].fields[2][i];
      scalars[f].fields[1][i] = scalars[f].fields[2][i];
      scalars[f].fields[3][i] = scalars[f].fields[2][i];
    }

    for (int i = 0; i < Nx; ++i)
    {
      scalars[f].field_2[i] = scalars[f].fields[2][i*Nrpath + 1];
    }

    scalars[f].calculate_C();
  }
}

void thimble_system::flow_rhs(const std::vector<dcomp> &x, std::vector<dcomp> &dx, const double t)
{
  //DEPRECIATED
  for (int i = 0; i < Njac; ++i)
  {
    dx[i] = conj(calc_dS(i, x));
  }
  for (int r = 0; r < Njac; ++r)
  {
    for (int c = 0; c < Njac; ++c)
    {
      for (int s = 0; s < Njac; ++s)
      {
        dx[Njac + r + c*Njac] = conj(calc_ddS(r, s, x)*x[r + c*Njac]);
      }
    }
  }
}

void thimble_system::action_output()
{
  dcomp action = calc_S(0);
  printf("action = %f%+fi\n", std::real(action), std::imag(action));
}

bool thimble_system::flow(std::vector<dcomp>& vec)
{
  double t(0), timestep(h);
  ode_handler ode_sys(*this);

  //setting up the boost ODE stepper
  runge_kutta_cash_karp54<std::vector<dcomp>> stepper;
  auto adaptive_stepper = make_controlled(1.e-4, 1.e-4, stepper);

  while (t < tau)
  {
    //if the simulation slows to a crawl, it's probably tending to infinity because the non-linearlity has screwed things up.
    //Given those steps should be rejected anyway, it's faster to scrub them early. This also prevents bugs due to the values overflowing.
    if(timestep > h/1000.)
    {
      //ensuring there's no overstepping
      if (t + timestep > tau)
      {
        timestep = tau - t;
      }
      adaptive_stepper.try_step(ode_sys, vec, t, timestep);
    }
    else
    {
      flow_check = false;
      return flow_check;
    }
  }
  return flow_check;
}

void thimble_system::test()
{
  for (int i = 0; i < Ntot; ++i)
  {
    printf("site %i, \t positive space site = %i \n", i, scalars[0].positive_space_site[i]);
  }
}

