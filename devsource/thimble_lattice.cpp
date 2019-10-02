#include "thimble_lattice.h"

//simple dirac delta function
template <class T>
T dd(T i, T j)
{
  T Delta(0);
  if (i == j)
  {
    Delta = 1;
  }
  return Delta;
}


//interaction class for mediating between fields
interaction::interaction(double Coupling, std::vector<int> Powers) : coupling(Coupling), powers(Powers) 
{

}

dcomp interaction::base(int site, thimble_system* current_system, bool ajustment)
{
  dcomp interaction_contribution = coupling;
  if (ajustment) //checks if this is to be done on the final field, or on the fields used as part of the ODE solver
  {
    for(int i = 0; i < powers.size(); ++i)
    {
      //applies all the fields raised to the relevant power
      interaction_contribution *= pow(current_system->scalars[i].ajustment_field[site], powers[i]); 
    }
  }
  else
  {
    for(int i = 0; i < powers.size(); ++i)
    {
      interaction_contribution *= pow(current_system->scalars[i].flowed_field[site], powers[i]);
    }
  }
  
  return interaction_contribution;
}

dcomp interaction::first_derivative(int site, int field, thimble_system* current_system, bool ajustment)
{
  dcomp interaction_contribution = coupling;
  if (ajustment)
  {
    //all the non-derivative fields up to the derivative
    for (int i = 0; i < field; ++i)
    {
      interaction_contribution *= pow(current_system->scalars[i].ajustment_field[site], powers[i]);
    }
    //contribution of the derivative field
    interaction_contribution *= double(powers[field])*pow(current_system->scalars[field].ajustment_field[site], powers[field] - 1);

    //contribution of all non-derivative fields from beyond the derivative field value
    for (int i = field + 1; i < powers.size(); ++i)
    {
      interaction_contribution *= pow(current_system->scalars[i].ajustment_field[site], powers[i]);
    }
  }
  else{
    //all the non-derivative fields up to the derivative
    for (int i = 0; i < field; ++i)
    {
      interaction_contribution *= pow(current_system->scalars[i].flowed_field[site], powers[i]);
    }
    //contribution of the derivative field
    interaction_contribution *= double(powers[field])*pow(current_system->scalars[field].flowed_field[site], powers[field] - 1);

    //contribution of all non-derivative fields from beyond the derivative field value
    for (int i = field + 1; i < powers.size(); ++i)
    {
      interaction_contribution *= pow(current_system->scalars[i].flowed_field[site], powers[i]);
    }
  }
  return interaction_contribution;
}

dcomp interaction::second_derivative(int site, int field_1, int field_2, thimble_system* current_system, bool ajustment)
{
  dcomp interaction_contribution = coupling;
  if (field_1 == field_2)
  {
    if (ajustment)
    {
      for (int i = 0; i < field_1; ++i)
      {
        interaction_contribution *= pow(current_system->scalars[i].ajustment_field[site], powers[i]);
      }
      //contribution of the derivative field
      interaction_contribution *= double(powers[field_1])*double(powers[field_1] - 1)*pow(current_system->scalars[field_1].ajustment_field[site], powers[field_1] - 2);

      //contribution of all non-derivative fields from beyond the derivative field value
      for (int i = field_1 + 1; i < powers.size(); ++i)
      {
        interaction_contribution *= pow(current_system->scalars[i].ajustment_field[site], powers[i]);
      }
    }
    else
    {
      for (int i = 0; i < field_1; ++i)
      {
        interaction_contribution *= pow(current_system->scalars[i].flowed_field[site], powers[i]);
      }
      //contribution of the derivative field
      interaction_contribution *= double(powers[field_1])*double(powers[field_1] - 1)*pow(current_system->scalars[field_1].flowed_field[site], powers[field_1] - 2);

      //contribution of all non-derivative fields from beyond the derivative field value
      for (int i = field_1 + 1; i < powers.size(); ++i)
      {
        interaction_contribution *= pow(current_system->scalars[i].flowed_field[site], powers[i]);
      }
    }
  }
  else
  {
    interaction_contribution = 0;
  }
  return interaction_contribution;
}

//scalar field class to be incorporated into a system class (class composition)*******************************************************
scalar_field::scalar_field(int x_dim, int t_dim, gsl_rng * rngPointer) : occupation_number(new int[x_dim]),
  Nx(x_dim), //member initialiser list 
  Nt(t_dim), 
  Npath(2*Nt), 
  Nrpath(Npath - 4), //parameters of the lattice
  Ntot(Nrpath*Nx), 
  m(0),
  squareMass(0), 
  dt(0.75),
  dx(1),
  path(new double[Nrpath]),
  path_offset(new double[Nrpath]),
  is_flowed(false),
  field_0(new dcomp[Nx]),
  field_1(new dcomp[Nx]),
  field_2(new dcomp[Nx]),
  base_field(new dcomp[Ntot]), //configuring the lattice arrays
  flowed_field(new dcomp[Ntot]),
  proposed_base_field(new dcomp[Ntot]),
  proposed_flowed_field(new dcomp[Ntot]),
  ajustment_field(new dcomp[Ntot]),
  positive_time_site(new int[Ntot]), //offset arrays to speed up computation
  positive_space_site(new int[Ntot]),
  negative_time_site(new int[Ntot]),
  negative_space_site(new int[Ntot]),
  my_rngPointer(rngPointer),
  j(0,1)
{
  //class constructor
  
  for (int i = 0; i < Ntot; ++i)
  {
    //zero initialising the fields
    base_field[i] = 0;
    flowed_field[i] = 0;
  }
  
  for (int i = 0; i < Nx; ++i)
  {
    //zero initialising the occupation number
    occupation_number[i] = 0;
  }

  for (int i = 0; i < (int) (Nrpath/2); ++i)
  {
    path[i] = dt;
    path[i + (int)(Nrpath/2)] = -1.*dt;
  }

  path_offset[0] = -1.*dt;
  for (int i = 0; i < Nrpath; ++i)
  {
    path_offset[i + 1] = path[i];
  }
}

scalar_field::~scalar_field()
{
    //class destructor
    delete[] base_field;
    delete[] flowed_field;
    delete[] proposed_base_field;
    delete[] proposed_flowed_field;
    delete[] ajustment_field;
    delete[] occupation_number;
    delete[] field_0;
    delete[] field_1;
    delete[] field_2;
    delete[] positive_time_site;
    delete[] positive_space_site;
    delete[] negative_time_site;
    delete[] negative_space_site;
    delete[] path;
    delete[] path_offset;
}

scalar_field::scalar_field(const scalar_field &obj) : occupation_number(new int[obj.Nx]), //object copy constructor
Nx(obj.Nx), //*large* member initalisation list 
Nt(obj.Nt),
Npath(obj.Npath),
Nrpath(obj.Nrpath),
Ntot(obj.Ntot),
m(obj.m),
squareMass(obj.squareMass),
dt(obj.dt),
dx(obj.dx),
path(new double[obj.Nrpath]),
path_offset(new double[obj.Nrpath]),
is_flowed(obj.is_flowed),
field_0(new dcomp[obj.Nx]),
field_1(new dcomp[obj.Nx]),
field_2(new dcomp[obj.Nx]),
positive_time_site(new int[obj.Ntot]),
positive_space_site(new int[obj.Ntot]),
negative_time_site(new int[obj.Ntot]),
negative_space_site(new int[obj.Ntot]),
my_rngPointer(obj.my_rngPointer),
j(obj.j),
base_field(new dcomp[obj.Ntot]),
flowed_field(new dcomp[obj.Ntot]),
proposed_base_field(new dcomp[obj.Ntot]),
proposed_flowed_field(new dcomp[obj.Ntot]),
ajustment_field(new dcomp[obj.Ntot])
{
  for (int i = 0; i < Nx; ++i)
  {
    occupation_number[i] = obj.occupation_number[i]; //setting values for the arrays that are copied over from the original object
    field_0[i] = obj.field_0[i];
    field_1[i] = obj.field_1[i];
  }
  
  for(int i = 0; i < Ntot; ++i)
  {
    positive_time_site[i] = obj.positive_time_site[i];
    positive_space_site[i] = obj.positive_space_site[i];
    negative_time_site[i] = obj.negative_time_site[i];
    negative_space_site[i] = obj.negative_space_site[i];
    base_field[i] = obj.base_field[i];
    flowed_field[i] = obj.flowed_field[i];
    proposed_base_field[i] = obj.proposed_base_field[i];
    proposed_flowed_field[i] = obj.proposed_flowed_field[i];
    ajustment_field[i] = obj.ajustment_field[i];
  }

  for (int i = 0; i < Nrpath; ++i)
  {
    path[i] = obj.path[i];
    path_offset[i] = obj.path_offset[i];
  }
}

void scalar_field::set_occupation_number(int new_occupation_number[])
{
  for (int i = 0; i < Nx; ++i)
  {
    occupation_number[i] = new_occupation_number[i];
  }
}

void scalar_field::set_occupation_number(int new_occupation_number)
{
  for(int i = 0; i < Nx; ++i)
  {
    occupation_number[i] = new_occupation_number;
  }
}

void scalar_field::set_mass(double new_mass)
{
  //allows the user to set a new mass for the scalar field
  if(new_mass < 0)
  {
    printf("Error: Mass must be greather than or equal to 0 \n");
    exit(5);
  }
  m = new_mass;
  squareMass = pow(new_mass,2);
}

void scalar_field::initialise()
{
    //this is called *after* the constructor, because we need to specify the mass
  double a[Nx], b[Nx], c[Nx], d[Nx];
  double p, omega_p, omega_tilde, Omega_p, V;
  
  for (int i = 0; i < Nx; ++i)
  {
      a[i] = gsl_ran_gaussian(my_rngPointer, 1); //Yeah naming a variable a, b, c, d isn't helpful, but it's what we call them in the maths. Besides, they're local to this function.
      b[i] = gsl_ran_gaussian(my_rngPointer, 1);
      c[i] = gsl_ran_gaussian(my_rngPointer, 1);
      d[i] = gsl_ran_gaussian(my_rngPointer, 1);
      field_0[i] = 0;
      field_1[i] = 0;
  } //random number arrays for Mou's initial conditions, and the initial states of the phi field
  
  V = Nx*dx;

  for (int i = 0; i < Nx; ++i)
	{
    p = i*2.*pi/(Nx*dx);
    omega_p = pow(pow(p,2) + pow(m,2),0.5);
    omega_tilde = acos(1 - pow(omega_p*dt,2)/2)/dt;
    Omega_p = sin(omega_p*dt)/dt; //configuring variables for this momentum
    if ((Nx - i)%Nx == 0)
    {
        //corner mode case
        field_0[i] += a[i]*pow(e,j*p*(i*dx))*pow(occupation_number[i] + 0.5,0.5)/pow(Omega_p,0.5);
        field_1[i] += pow(e,j*p*(i*dx))*(a[i]*cos(omega_tilde*dt) + c[i]*Omega_p*dt)*pow(occupation_number[i] + 0.5,0.5)/pow(Omega_p,0.5);
    }
    else
    {
        //bulk mode case
        field_0[i] += ((a[i] + j*b[i])*pow(e,j*p*(i*dx))/pow(2*Omega_p,0.5) + (c[i] - j*b[i])*pow(e,-1.*j*p*(i*dx))/pow(2*Omega_p,0.5))*pow(occupation_number[i] + 0.5,0.5);
        field_1[i] += pow(e,j*p*(i*dx))*(a[i]*cos(omega_tilde*dt) + c[i]*Omega_p*dt)*pow(occupation_number[i] + 0.5,0.5)/pow(Omega_p,0.5);
    }

    field_0[i] = field_0[i]/V; //rescaling for the volume. hbar is taken to be one.
    field_1[i] = field_1[i]/V;

    //manual force to check values
    //field_0[i] = 0.8;
    //field_1[i] = 1.0;
  }
  
  for(int k = 0; k < Nx; ++k)
  {
      base_field[0 + Nrpath*k] = field_1[k];
      base_field[1 + Nrpath*k] = -1.0*dt*dt*(squareMass*base_field[0 + Nrpath*k]) + 2.0*base_field[0 + Nrpath*k] - field_0[k];
      for (int i = 1; i < (int) (Nrpath/2); ++i)
      {
        base_field[i + Nrpath*k + 1] = -dt*dt*squareMass*base_field[i + Nrpath*k] + 2.0*base_field[i + Nrpath*k] - base_field[i + Nrpath*k - 1];
      }
      
      for(int i = 0; i < (int) (Nrpath/2); ++i)
      {
        base_field[Nrpath*(k + 1) - i - 1] = base_field[Nrpath*k + i + 1]; //sets up the return leg of the contour
      }

  }

  //clearing the classical data from the first site, note the initial condition data is saved in field_0 and field_1
  for(int i = 0; i < Nx; ++i)
  {
    base_field[i] = 0;
  }

  //setting up the co-ordinate shifted arrays. Could combine them, but this won't be called much, so the legibility is prioritised over speed
  //I imagine the O3 flag for g++ automatically tidies them up anyway
  for (int k = 0; k < Nx; ++k)
  {
    for (int i = 0; i < Nrpath - 1; ++i)
    {
      positive_time_site[k*Nrpath + i] = k*Nrpath + i + 1;
    }
    positive_time_site[(k + 1)*Nrpath - 1] = k*Nrpath;
    
    for (int i = 1; i < Nrpath; ++i)
    {
      negative_time_site[k*Nrpath + i] = k*Nrpath + i - 1;
    }
    negative_time_site[k*Nrpath] = (k + 1)*Nrpath - 1;
  }
  
  for (int i = 0; i < Nrpath; ++i)
  {
    for (int k = 0; k < Nx - 1; ++k)
    {
      positive_space_site[k*Nrpath + i] = (k + 1)*Nrpath + i;
    }
    positive_space_site[(Nx - 1)*Nrpath + i] = i;
    
    for (int k = 1; k < Nx; ++k)
    {
      negative_space_site[k*Nrpath + i] = (k - 1)*Nrpath + i;
    }
    negative_space_site[i] = (Nx - 1)*Nrpath + i;
  }

  for (int i = 0; i < Ntot; ++i)
  {
    flowed_field[i] = base_field[i];
  }
}

dcomp scalar_field::free_action(int site, bool ajustment)
{
  //Standard P^2 - m^2 action
  int n = calc_n(site);
  dcomp* work_field;
  if (ajustment)
  {
    work_field = ajustment_field;
  }
  else
  {
    work_field = flowed_field;
  }
  
  dcomp S = pow(work_field[positive_time_site[site]] - work_field[site], 2)/(2.*path[n]) 
  - (((path[n] + path_offset[n])/2)*(pow(work_field[positive_space_site[site]] - work_field[site], 2))/(2*pow(dx,2)) + squareMass*pow(work_field[site], 2));
  return S;
}

dcomp scalar_field::free_action_derivative(int site, bool ajustment)
{
  //derivative of the above action
  int n = calc_n(site);
  dcomp* work_field;
  if(ajustment_field)
  {
    work_field = ajustment_field;
  }
  else
  {
    work_field = flowed_field;
  }
  
  dcomp dS = (work_field[site] - work_field[positive_time_site[site]])/path[n] + (work_field[site] - work_field[negative_time_site[site]])/path_offset[n]
  - ((path[n] + path_offset[n])/2)*((2.*work_field[site] - work_field[positive_space_site[site]] - work_field[negative_space_site[site]])/pow(dx,2) + squareMass*work_field[site]);
  return dS;
}

dcomp scalar_field::free_action_second_derivative(int site_1, int site_2)
{
  //second derivative calculation
  int n = calc_n(site_1);
  dcomp ddS = (dd(site_1, site_2) + dd(positive_time_site[site_1], site_2))/path[n] + (dd(site_1, site_2) - dd(negative_time_site[site_1], site_2))/path_offset[n]
    - ((path[n] + path_offset[n])/2.)*((2*dd(site_1, site_2) - dd(positive_space_site[site_1], site_2) - dd(negative_space_site[site_1], site_2))/pow(dx, 2) - squareMass*dd(site_1, site_2));
  return ddS;
}

void scalar_field::set_dt(double new_dt)
{
  for (int i = 0; i < Nrpath; ++i)
  {
    path[i] *= new_dt/dt;
  }
  dt = new_dt;
}

dcomp scalar_field::edge_effects(int site)
{
  dcomp effect = 0;
  int n = calc_n(site);
  int x = calc_x(site);
  if (n == 0)
  {
    effect = 2.*field_2[x]/dt;
  }
  else if (n == 1)
  {
    effect = -1.*field_1[x]/dt;
  }
  else if (n == Nrpath - 1)
  {
    effect = 1.*field_1[x]/dt;
  }
  effect *= flowed_field[site];
  return effect;
}

dcomp scalar_field::edge_effects_derivative(int site)
{
  dcomp effect = 0;
  int n = calc_n(site);
  int x = calc_x(site);
  if (n == 0)
  {
    effect = 2.*field_2[x]/dt;
  }
  else if (n == 1)
  {
    effect = -1.*field_1[x]/dt;
  }
  else if (n == Nrpath - 1)
  {
    effect = field_1[x]/dt;
  }
  return effect;
}

//decomposing the total lattice position into the single timeslice position (for the dt array)
int scalar_field::calc_n(int site)
{
  int n = site%Nrpath;
  return n;
}

int scalar_field::calc_x(int site)
{
  int x = int((site - calc_n(site))/Nrpath);
  return x;
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
jac_defined(false),
rel_path(""),
j(0,1)
{
    //thimble system constructor
    const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    my_rngPointer = gsl_rng_alloc (T);
    gsl_rng_set(my_rngPointer, seed);
    
    //determining the number of timesteps for the ODE solvers from the flow time  
    h = 0.02; //sets the base size 
    number_of_timesteps = int(ceil(tau/h)); //calculates how many steps this corresponds to (overshooting in the case of it not being exact)
    h = tau/number_of_timesteps; //readusting the size of h to prevent overshooting
}

thimble_system::~thimble_system()
{
  gsl_rng_free(my_rngPointer);
  if (jac_defined)
  {
    delete[] J;
    delete[] proposed_J;
    delete[] invJ;
    delete[] proposed_invJ;
  }
}

void thimble_system::add_scalar_field()
{
  scalars.emplace_back(scalar_field(Nx, Nt, my_rngPointer));
}

void thimble_system::add_scalar_field(double mass)
{
  scalars.emplace_back(scalar_field(Nx, Nt, my_rngPointer));
  scalars[scalars.size() - 1].set_mass(mass);
}

void thimble_system::add_interaction(double coupling, std::vector<int> powers)
{
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

dcomp thimble_system::calc_dS(int site, int field, bool ajustment)
{
  dcomp dS; //derivative of the action
  dcomp interaction = 0;
  dS = scalars[field].free_action_derivative(site, ajustment); //free field kinetic contribution
  for (int i = 0; i < interactions.size(); ++i)
  {
    //looping through the first derivatives of all the interactions (derivatives with respect to this field)
    interaction += interactions[i].first_derivative(site, field, this); 
  }
  interaction *= (scalars[field].path[scalars[field].calc_n(site)] + scalars[field].path_offset[scalars[field].calc_n(site)])/2.; //(delta_n + delta_n-1) factor
  dS += interaction;
  return dS;
}

dcomp thimble_system::calc_dS(int site, bool ajustment)
{
  int field = 0;
  int internal_site = site; //this essentially takes the busy work out of calculating which field the Jacobian is dealing with
  if (site > Ntot)
  {
    while (internal_site > Ntot)
    {
      internal_site -= Ntot;
      ++field;
    }
  }
  return calc_dS(internal_site, field, ajustment);
}

dcomp thimble_system::calc_ddS(int site_1, int site_2, int field_1, int field_2, bool ajustment)
{
  dcomp ddS;
  dcomp interaction = 0;
  if (field_1 == field_2)
  {
    ddS = scalars[field_1].free_action_second_derivative(site_1, site_2); //free component
  }
  
  if (site_1 == site_2) //only add the interactions which are on equal sites
  {
    for (int i = 0; i < interactions.size(); ++i)
    {
      interaction *= interactions[i].second_derivative(site_1, field_1, field_2, this, ajustment);
    }
    interaction *= (scalars[field_1].path[scalars[field_1].calc_n(site_1)] + scalars[field_1].path_offset[scalars[field_1].calc_n(site_1)])/2.;
  }
  ddS += interaction;
  return ddS;
}

dcomp thimble_system::calc_ddS(int site_1, int site_2, bool ajustment)
{
  int field_1 = 0;
  int field_2 = 0;
  dcomp ddS;

  while(site_1 > Ntot)
  {
    site_1 -= Ntot;
    ++field_1;
  }
  
  while(site_2 > Ntot)
  {
    site_2 -= Ntot;
    ++field_2;
  }

  ddS = calc_ddS(site_1, site_2, field_1, field_2, ajustment);
  return ddS;
}

field_id_return thimble_system::calc_field(int master_site)
{
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
{
  for(int i = 0; i < scalars.size(); ++i)
  {
    for(int r = 0; r < Ntot; ++r)
    {
      scalars[i].ajustment_field[r] = ajustment[i*Ntot + r];
    }
  }
}

dcomp thimble_system::calc_jacobian(dcomp Jac[], bool proposal)
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
  bool ajustment = true;

  int s;
  gsl_permutation* p = gsl_permutation_alloc(Njac);
  gsl_matrix_complex* mJ = gsl_matrix_complex_alloc(Njac, Njac);
  gsl_complex det_gsl;

  //identifying if it's the proposal or exising fields we wish to flow
  if (proposal)
  {
    for (int i = 0; i < scalars.size(); ++i)
    {
      for (int k = 0; k < Ntot; ++i)
      {
        working_scalar[i*Ntot + k] = scalars[i].proposed_base_field[k];
      }
    }
  }
  else
  {
    for(int i = 0; i < scalars.size(); ++i)
    {
      for(int k = 0; k < Ntot; ++k)
      {
        working_scalar[i*Ntot + k] = scalars[i].base_field[k];
      }
    }
  }
  //setting up an identity matrix
  for(int r = 0; r < Njac; ++r)
  {
    for(int c = 0; c < Njac; ++c)
    {
      if (r == c)
      {
        Jac[r + Njac*c] = 1.;
      }
      else
      {
        Jac[r + Njac*c] = 0.;
      }
    }
  }
  //standard implementation of RK45
  for (int i = 0; i < number_of_timesteps; ++i)
  {
    for (int r = 0; r < Njac; ++r)
    {
      k1_scalar[r] = h*std::conj(-1.*j*calc_dS(r));
      ajustment_scalar[r] = working_scalar[r] + k1_scalar[r]/2.;
      for (int c = 0; c < Njac; ++c)
      {
        k1_jac[r + Ntot*c] = h*std::conj(-1.*j*calc_ddS(r, c)*Jac[r + Ntot*c]);
        ajustment_jac[r + Ntot*c] = Jac[r + Ntot*c] + k1_jac[r + Ntot*c]/2.;
      }
    }
    sync_ajustment(ajustment_scalar);
    //ajustment scalar does the job of holding field values that are calculated intermittenlty in the RK45 method
    //sync_ajustment pushes the values stored in the array in this function back out to the scalarfields, which then use the value to calculate the ds and dds functions

    for (int r = 0; r < Njac; ++r)
    {
      k2_scalar[r] = h*std::conj(-1.*j*calc_dS(r, ajustment));
      for (int c = 0; c < Njac; ++c)
      {
        k2_jac[r + Ntot*c] = h*std::conj(calc_ddS(r, c, ajustment));
      }
    }

    for (int r = 0; r < Njac; ++r)
    {
      ajustment_scalar[r] = working_scalar[r] + k2_scalar[r]/2.;
      for (int c = 0; c < Njac; ++c)
      {
        ajustment_jac[r + Ntot*c] = Jac[r + Ntot*c] + k2_jac[r + Ntot*c]/2.;
      }
    }
    sync_ajustment(ajustment_scalar);

    for (int r = 0; r < Njac; ++r)
    {
      k3_scalar[r] = h*std::conj(-1.*j*calc_dS(r, ajustment));
      for (int c = 0; c < Njac; ++c)
      {
        k3_jac[r + Ntot*c] = h*std::conj(-1.*j*calc_ddS(r, c, ajustment));
      }
    }

    for (int r = 0; r < Njac; ++r)
    {
      ajustment_scalar[r] = working_scalar[r] + k3_scalar[r];
      for (int c = 0; c < Njac; ++c)
      {
        ajustment_jac[r + Ntot*c] = Jac[r + Ntot*c] + k3_jac[r + Ntot*c];
      }
    }
    sync_ajustment(ajustment_scalar);

    for (int r = 0; r < Njac; ++r)
    {
      k4_scalar[r] = h*std::conj(-1.*j*calc_dS(r, ajustment));
      for (int c = 0; c < Njac; ++c)
      {
        k4_jac[r + Ntot*c] = h*std::conj(-1.*j*calc_ddS(r, c, ajustment));
      }
    }

    for (int r = 0; r < Njac; ++r)
    {
      working_scalar[r] = (k1_scalar[r] + 2.*k2_scalar[r] + 2.*k3_scalar[r] + k4_scalar[r])/6.;
      for (int c = 0; c < Njac; ++c)
      {
        Jac[r + Ntot*c] = (k1_jac[r] + 2.*k2_jac[r] + 2.*k3_jac[r] + k4_jac[r])/6.;
      }
    }
  }
  //returning the flowed fields to the scalar fields class
  if(proposal)
  {
    for (int i = 0; i < scalars.size(); ++i)
    {
      for (int k = 0; k < Ntot; ++i)
      {
        scalars[i].proposed_flowed_field[k] = working_scalar[i*Ntot + k];
      }
    }
  }
  else
  {
    for (int i = 0; i < scalars.size(); ++i)
    {
      for (int k = 0; k < Ntot; ++k)
      {
        scalars[i].flowed_field[k] = working_scalar[i*Ntot +k];
      }
    }
  }

  //casting from our complex array to a GSL matrix
  for(int r = 0; r < Njac; ++r)
  {
    for(int c = 0; c < Njac; ++c)
    {
      gsl_matrix_complex_set(mJ, r, c, gsl_complex_rect(std::real(Jac[r + c*Ntot]), std::imag(Jac [r + c*Ntot])));
    }
  }
  gsl_linalg_complex_LU_decomp(mJ, p, &s);
  det_gsl = gsl_linalg_complex_LU_det(mJ, s); //this actually calculates the determinant
  
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
  gsl_permutation_free(p);
  gsl_matrix_complex_free(mJ);

  return GSL_REAL(det_gsl) + j*GSL_IMAG(det_gsl);
}

void thimble_system::simulate(int n_burn_in, int n_simulation)
{
  Njac = scalars.size()*Ntot; //setting the size of the Jacobian "matricies"
  NjacSquared = pow(Njac, 2);

  J = new dcomp[NjacSquared]; //setting up the arrays that hhold the Jacobian matricies
  proposed_J = new dcomp[NjacSquared];
  invJ = new dcomp[NjacSquared];
  proposed_invJ = new dcomp[NjacSquared];
  jac_defined = true; //this ensures the correct memory management happens
}