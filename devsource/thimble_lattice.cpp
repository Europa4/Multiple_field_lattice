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

dcomp interaction::base(int site, thimble_system* current_system)
{
  dcomp interaction_contribution = coupling;
  for(int i = 0; i < powers.size(); ++i)
  {
    interaction_contribution *= pow(current_system->scalars[i].flowed_field[site], powers[i]);
  }
  return interaction_contribution;
}

dcomp interaction::first_derivative(int site, int field, thimble_system* current_system)
{
  dcomp interaction_contribution = coupling;
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
  return interaction_contribution;
}

dcomp interaction::second_derivative(int site, int field_1, int field_2, thimble_system* current_system)
{
  dcomp interaction_contribution = coupling;
  if (field_1 == field_2)
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
  base_field(new dcomp[Ntot]), //configuring the lattice arrays
  flowed_field(new dcomp[Ntot]),
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
    delete[] occupation_number;
    delete[] field_0;
    delete[] field_1;
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
positive_time_site(new int[obj.Ntot]),
positive_space_site(new int[obj.Ntot]),
negative_time_site(new int[obj.Ntot]),
negative_space_site(new int[obj.Ntot]),
my_rngPointer(obj.my_rngPointer),
j(obj.j),
base_field(new dcomp[obj.Ntot]),
flowed_field(new dcomp[obj.Ntot])
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

dcomp scalar_field::free_action(int site)
{
  //Standard P^2 - m^2 action
  int n = calc_n(site);
  dcomp S = pow(flowed_field[positive_time_site[site]] - flowed_field[site], 2)/(2.*path[n]) 
  - (((path[n] + path_offset[n])/2)*(pow(flowed_field[positive_space_site[site]] - flowed_field[site], 2))/(2*pow(dx,2)) + squareMass*pow(flowed_field[site], 2));
  return S;
}

dcomp scalar_field::free_action_derivative(int site)
{
  //derivative of the above action
  int n = calc_n(site);
  dcomp dS = (flowed_field[site] - flowed_field[positive_time_site[site]])/path[n] + (flowed_field[site] - flowed_field[negative_time_site[site]])/path_offset[n]
  - ((path[n] + path_offset[n])/2)*((2.*flowed_field[site] - flowed_field[positive_space_site[site]] - flowed_field[negative_space_site[site]])/pow(dx,2) + squareMass*flowed_field[site]);
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

//decomposing the total lattice position into the single timeslice position (for the dt array)
int scalar_field::calc_n(int site)
{
  int n = site%Nrpath;
  return n;
}
//*************************************thimble_system***************************

thimble_system::thimble_system(int x_dim, int t_dim, double flow_time, long unsigned int seed) : Nx(x_dim), Nt(t_dim), tau(flow_time), Npath(2*Nt), Nrpath(Npath - 4), Ntot(Nrpath*Nx), rng_seed(seed)
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