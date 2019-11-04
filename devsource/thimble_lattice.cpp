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


//interaction class for mediating between fields
interaction::interaction(double Coupling, std::vector<int> Powers) : coupling(Coupling), powers(Powers) 
{

}

dcomp interaction::base(int site, thimble_system* current_system, int field_type)
{
  dcomp interaction_contribution = coupling;
  for(int i = 0; i < powers.size(); ++i)
  {
    //applies all the fields raised to the relevant power
    interaction_contribution *= pow(current_system->scalars[i].fields[field_type][site], powers[i]); 
  }

  return interaction_contribution;
}

dcomp interaction::first_derivative(int site, int field, thimble_system* current_system, int field_type)
{
  dcomp interaction_contribution = coupling;
  //all the non-derivative fields up to the derivative
  for (int i = 0; i < field; ++i)
  {
    interaction_contribution *= pow(current_system->scalars[i].fields[field_type][site], powers[i]);
  }
  //contribution of the derivative field
  interaction_contribution *= double(powers[field])*pow(current_system->scalars[field].fields[field_type][site], powers[field] - 1);

  //contribution of all non-derivative fields from beyond the derivative field value
  for (int i = field + 1; i < powers.size(); ++i)
  {
    interaction_contribution *= pow(current_system->scalars[i].fields[field_type][site], powers[i]);
  }
  return interaction_contribution;
}

dcomp interaction::second_derivative(int site, int field_1, int field_2, thimble_system* current_system, int field_type)
{
  dcomp interaction_contribution = coupling;
  if (field_1 == field_2)
  {
    for (int i = 0; i < field_1; ++i)
    {
      interaction_contribution *= pow(current_system->scalars[i].fields[field_type][site], powers[i]);
    }
    //contribution of the derivative field
    interaction_contribution *= double(powers[field_1])*double(powers[field_1] - 1)*pow(current_system->scalars[field_1].fields[field_type][site], powers[field_1] - 2);

    //contribution of all non-derivative fields from beyond the derivative field value
    for (int i = field_1 + 1; i < powers.size(); ++i)
    {
      interaction_contribution *= pow(current_system->scalars[i].fields[field_type][site], powers[i]);
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
  positive_time_site(new int[Ntot]), //offset arrays to speed up computation
  positive_space_site(new int[Ntot]),
  negative_time_site(new int[Ntot]),
  negative_space_site(new int[Ntot]),
  my_rngPointer(rngPointer),
  j(0,1)
{
  //class constructor
  fields[0] = new dcomp[Ntot]; //configuring the lattice arrays
  fields[1] = new dcomp[Ntot];
  fields[2] = new dcomp[Ntot];
  fields[3] = new dcomp[Ntot];
  fields[4] = new dcomp[Ntot];

  C[0] = new dcomp[Ntot]; //configuring Mou's constant arrays
  C[1] = new dcomp[Ntot];
  C[2] = new dcomp[Ntot];
  C[3] = new dcomp[Ntot];
  C[4] = new dcomp[Ntot];

  int n = 0;
  
  for (int i = 0; i < Ntot; ++i)
  {
    //zero initialising the fields
    for(int k = 0; k < 5; ++k)
    {
      fields[k][i] = 0;
    }
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
  for (int i = 0; i < Nrpath - 1; ++i)
  {
    path_offset[i + 1] = path[i];
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

}

scalar_field::~scalar_field()
{
    //class destructor
    for (int i = 0; i < 5; ++i)
    {
      delete[] fields[i];
      delete[] C[i];
    }
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
j(obj.j)
{
  for (int i = 0; i < 5; ++i)
  {
    fields[i] = new dcomp[obj.Ntot];
    C[i] = new dcomp[obj.Ntot];
  }

  for (int i = 0; i < Nx; ++i)
  {
    occupation_number[i] = obj.occupation_number[i]; //setting values for the arrays that are copied over from the original object
    field_0[i] = obj.field_0[i];
    field_1[i] = obj.field_1[i];
    field_2[i] = obj.field_2[i];
  }
  
  for(int i = 0; i < Ntot; ++i)
  {
    positive_time_site[i] = obj.positive_time_site[i];
    positive_space_site[i] = obj.positive_space_site[i];
    negative_time_site[i] = obj.negative_time_site[i];
    negative_space_site[i] = obj.negative_space_site[i];
    for (int k = 0; k < 5; ++k)
    {
      fields[k][i] = obj.fields[k][i];
      C[k][i] = obj.C[k][i];
    }
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
      fields[2][0 + Nrpath*k] = field_1[k];
      fields[2][1 + Nrpath*k] = -1.0*dt*dt*(squareMass*fields[2][0 + Nrpath*k]) + 2.0*fields[2][0 + Nrpath*k] - field_0[k];
      for (int i = 1; i < (int) (Nrpath/2); ++i)
      {
        fields[2][i + Nrpath*k + 1] = -dt*dt*squareMass*fields[2][i + Nrpath*k] + 2.0*fields[2][i + Nrpath*k] - fields[2][i + Nrpath*k - 1];
      }
      
      for(int i = 0; i < (int) (Nrpath/2); ++i)
      {
        fields[2][Nrpath*(k + 1) - i - 1] = fields[2][Nrpath*k + i + 1]; //sets up the return leg of the contour
      }

  }

  //clearing the classical data from the first site, note the initial condition data is saved in field_0 and field_1
  for(int i = 0; i < Nx; ++i)
  {
    fields[2][i] = 0;
  }


  for (int i = 0; i < Ntot; ++i)
  {
    fields[0][i] = fields[2][i];
  }

  //setting up Mou's constant arrays
  for (int i = 0; i < Ntot; ++i)
  {
    int n = calc_n(i);
    C[0][i] = -1.*dx*j*(1/path[n] + 1/path_offset[n] + (path[n] + path_offset[n]/2.)*(-2./pow(dx, 2) - squareMass));
    C[1][i] = j*dx/path[n];
    C[2][i] = j*dx/path_offset[n];
    C[3][i] = -1.*dx*j*(path[n] + path_offset[n])/(2.*pow(dx, 2));
    C[4][i] = 0.;
  }
  //edge terms for the edge effects
  for (int i = 0; i < Nx; ++i)
  {
    C[4][i*Nrpath] = -2.*j*field_2[i]/dt;
    C[4][i*Nrpath + 1] = j*field_1[i]/dt;
    C[4][(i + 1)*Nrpath - 1] = -1.*j*field_1[i]/dt;
  }

  //setting up Mou's constant arrays
  for (int i = 0; i < Ntot; ++i)
  {
    int n = calc_n(i);
    C[0][i] = -1.*dx*j*(1/path[n] + 1/path_offset[n] + ((path[n] + path_offset[n])/2.)*(-2./pow(dx, 2) - squareMass));
    C[1][i] = j*dx/path[n];
    C[2][i] = j*dx/path_offset[n];
    C[3][i] = -1.*dx*j*(path[n] + path_offset[n])/(2.*pow(dx, 2));
    C[4][i] = 0.;
  }
  //edge terms for the edge effects
  for (int i = 0; i < Nx; ++i)
  {
    C[4][i*Nrpath] = -2.*j*field_2[i]/dt;
    C[4][i*Nrpath + 1] = j*field_1[i]/dt;
    C[4][(i + 1)*Nrpath - 1] = -1.*j*field_1[i]/dt;
  }

  //applying the anti-periodic bounday terms
  for(int i = 0; i < Nx; ++i)
  {
    C[1][(i + 1)*Nrpath - 1] *= -1.;
    C[2][i*Nrpath] *= -1.;
  }

  for(int i = 0; i < Nx; ++i)
  {
    field_2[i] = fields[2][i*Nrpath + 1];
  }
}

dcomp scalar_field::free_action(int site, int field_type)
{
  //Standard P^2 - m^2 action
  dcomp S = 0;
  int n = calc_n(site);
  if(n != Nrpath - 1)
  {
    S = -1.*(C[1][site]/2.)*pow(fields[field_type][positive_time_site[site]] - fields[field_type][site], 2)
      - pow(dx, 2)*C[3][site]*(pow(fields[field_type][positive_space_site[site]] - fields[field_type][site], 2)/(2.*pow(dx, 2)) + squareMass*pow(fields[field_type][site], 2)/2.)
      + C[4][site]*fields[field_type][site];
  }
  else
  {
    //This is the anti periodic boundary term
    S = (C[1][site]/2.)*pow(fields[field_type][positive_time_site[site]] + fields[field_type][site], 2)
      - pow(dx, 2)*C[3][site]*(pow(fields[field_type][positive_space_site[site]] - fields[field_type][site], 2)/(2.*pow(dx, 2)) + squareMass*pow(fields[field_type][site], 2)/2.)
      + C[4][site]*fields[field_type][site];
  }
  return S;
}

dcomp scalar_field::free_action_derivative(int site, int field_type)
{
  //derivative of the above action
  dcomp dS = C[0][site]*fields[field_type][site] + C[1][site]*fields[field_type][positive_time_site[site]] + C[2][site]*fields[field_type][negative_time_site[site]]
    + C[3][site]*(fields[field_type][positive_space_site[site]] + fields[field_type][negative_space_site[site]]) + C[4][site];
  return dS;
}

dcomp scalar_field::free_action_second_derivative(int site_1, int site_2)
{
  //second derivative calculation
  dcomp ddS = C[0][site_1]*dd(site_1, site_2) + C[1][site_1]*dd(positive_time_site[site_1], site_2) + C[2][site_1]*dd(negative_time_site[site_1], site_2)
    + C[3][site_1]*(dd(positive_space_site[site_1], site_2) + dd(negative_space_site[site_1], site_2));
  return ddS;
}

void scalar_field::set_dt(double new_dt)
{
  for (int i = 0; i < Nrpath; ++i)
  {
    path[i] *= new_dt/dt; //rescales the path to account for a new lattice spacing
  }
  dt = new_dt;
}

dcomp scalar_field::edge_effects(int site, int field_type)
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
  effect *= fields[field_type][site];
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
J(0,0),
J_conj(0,0),
rel_path(""),
file_name("0"),
j(0,1),
dx(1.),
dt(0.75),
sigma(0.07071067812),
delta(0.1)
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

dcomp thimble_system::calc_dS(int site, int field, int field_type)
{
  dcomp dS = 0; //derivative of the action
  dcomp interaction = 0;
  dS = scalars[field].free_action_derivative(site, field_type); //free field kinetic contribution
  for (int i = 0; i < interactions.size(); ++i)
  {
    //looping through the first derivatives of all the interactions (derivatives with respect to this field)
    interaction += interactions[i].first_derivative(site, field, this, field_type); 
  }
  interaction *= -1.*dx*j*(scalars[field].path[scalars[field].calc_n(site)] + scalars[field].path_offset[scalars[field].calc_n(site)])/2.; //(delta_n + delta_n-1) factor
  dS += interaction;
  return dS;
}

dcomp thimble_system::calc_dS(int site, int field_type)
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
  return calc_dS(internal_site, field, field_type);
}

dcomp thimble_system::calc_ddS(int site_1, int site_2, int field_1, int field_2, int field_type)
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
      interaction *= interactions[i].second_derivative(site_1, field_1, field_2, this, field_type);
    }
    interaction *= -1.*j*dx*(scalars[field_1].path[scalars[field_1].calc_n(site_1)] + scalars[field_1].path_offset[scalars[field_1].calc_n(site_1)])/2.;
  }
  //printf("site = %i %i, free action = %f%+fi, interaction = %f%+fi \n", site_1, site_2, real(ddS), imag(ddS), real(interaction), imag(interaction));
  ddS += interaction;
  return ddS;
}

dcomp thimble_system::calc_ddS(int site_1, int site_2, int field_type)
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
  for (int i = 0; i < scalars.size(); ++i)
  {
    for (int k = 0; k < Ntot; ++k)
    {
      //resetting the flowed field to the base field. At this point we're at tau = 0 so the flowed field and base field should be identical
      scalars[i].fields[proposal_or - 2][k] = scalars[i].fields[proposal_or][k];
      working_scalar[i*Ntot + k] = scalars[i].fields[proposal_or][k];
    }
  }

  //setting up an identity matrix
  for(int r = 0; r < Njac; ++r)
  {
    for(int c = 0; c < Njac; ++c)
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
  for (int i = 0; i < number_of_timesteps; ++i)
  {
    for (int r = 0; r < Njac; ++r)
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

    for (int r = 0; r < Njac; ++r)
    {
      k2_scalar[r] = h*conj(calc_dS(r, ajustment));
      for (int c = 0; c < Njac; ++c)
      {
        k2_jac[r + Njac*c] = 0;
        for(int s = 0; s < Njac; ++s)
        {
          k2_jac[r + Njac*c] += h*conj(calc_ddS(r, s, ajustment)*ajustment_jac[s + Njac*c]);
        }
      }
    }

    for (int r = 0; r < Njac; ++r)
    {
      ajustment_scalar[r] = working_scalar[r] + k2_scalar[r]/2.;
      for (int c = 0; c < Njac; ++c)
      {
        ajustment_jac[r + Njac*c] = Jac.get_element(r,c) + k2_jac[r + Njac*c]/2.;
      }
    }
    sync_ajustment(ajustment_scalar);

    for (int r = 0; r < Njac; ++r)
    {
      k3_scalar[r] = h*conj(calc_dS(r, ajustment));
      for (int c = 0; c < Njac; ++c)
      {
        k3_jac[r + Njac*c] = 0;
        for (int s = 0; s < Njac; ++s)
        {
          k3_jac[r + Njac*c] += h*conj(calc_ddS(r, s, ajustment)*ajustment_jac[s + Njac*c]);
        }
      }
    }

    for (int r = 0; r < Njac; ++r)
    {
      ajustment_scalar[r] = working_scalar[r] + k3_scalar[r];
      for (int c = 0; c < Njac; ++c)
      {
        ajustment_jac[r + Njac*c] = Jac.get_element(r, c) + k3_jac[r + Njac*c];
      }
    }
    sync_ajustment(ajustment_scalar);

    for (int r = 0; r < Njac; ++r)
    {
      k4_scalar[r] = h*conj(calc_dS(r, ajustment));
      for (int c = 0; c < Njac; ++c)
      {
        k4_jac[r + Njac*c] = 0;
        for(int s = 0; s < Njac; ++s)
        {
          k4_jac[r + Njac*c] += h*conj(calc_ddS(r, s, ajustment)*ajustment_jac[s + Njac*c]);
        }
      }
    }
    for (int r = 0; r < Njac; ++r)
    {
      working_scalar[r] += (k1_scalar[r] + 2.*k2_scalar[r] + 2.*k3_scalar[r] + k4_scalar[r])/6.;
      for (int c = 0; c < Njac; ++c)
      {
        jac_element = Jac.get_element(r, c);
        jac_element += (k1_jac[r + Njac*c] + 2.*k2_jac[r + Njac*c] + 2.*k3_jac[r + Njac*c] + k4_jac[r + Njac*c])/6.;
        Jac.set_element(r, c, jac_element);
      }
    }
    //returning the flowed fields to the scalar fields object
    for (int i = 0; i < scalars.size(); ++i)
    {
      for (int k = 0; k < Ntot; ++k)
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

dcomp thimble_system::calc_S(int field_type)
{
  dcomp S = 0;
  int n;
  //looping through all the fields and sites to add their free free field contributions
  for (int i = 0; i < scalars.size(); ++i)
  {
    for (int k = 0; k < Ntot; ++k)
    {
      S += scalars[i].free_action(k, field_type);
    }
  }
  //looping through all the interactions to add them to the action too
  for (int i = 0; i < interactions.size(); ++i)
  {
    for (int k = 0; k < Ntot; ++k)
    {
      n = scalars[0].calc_n(k);
      S += -1.*j*dx*(scalars[0].path[n] + scalars[0].path_offset[n])*interactions[i].base(k, this, field_type)/2.;
    }
  }
  return S;
}

int thimble_system::update()
{
  bool proposal = true;
  dcomp proposed_action, proposed_detJ;
  mydouble log_proposal;
  dcomp* eta = new dcomp[Njac];
  mydouble matrix_exponenet, proposed_matrix_exponenet, exponenet, check;
  int output = 0; //this is the return value

  for(int i = 0; i < Njac; ++i)
  {
    //setting up the proposal on the complex manifold
    eta[i] = gsl_ran_gaussian(my_rngPointer, sigma) + j*gsl_ran_gaussian(my_rngPointer, sigma);
  }

  //returning the proposal to the real manifold
  matrix<dcomp> Delta = J.solve(eta);
  for(int i = 0; i < Njac; ++i)
  {
    //taking only the elements that fit in the reduced space
    Delta.set_element(i, 0, real(Delta.get_element(i, 0)));
  }

  //creating new basefield condtions
  for(int i = 0; i < scalars.size(); ++i)
  {
    for(int k = 0; k < Ntot; ++k)
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
  exponenet = ((mydouble) real(S - proposed_action)) + 2.*log_proposal - 2.*((mydouble) log(real(J.get_det()))) + matrix_exponenet/pow(delta, 2) - proposed_matrix_exponenet/pow(delta, 2);
  check = gsl_rng_uniform(my_rngPointer);
  if (exp(exponenet) > check)
  {
    //proposal was accepted, transfering all the proposed parameters to the storage
    S = proposed_action;
    J = proposed_J;
    J_conj = proposed_J_conj;
    //transfering over all the scalar fields
    for (int i = 0; i < scalars.size(); ++i)
    {
      for (int k = 0; k < Ntot; ++k)
      {
        scalars[i].fields[0][k] = scalars[i].fields[1][k];
        scalars[i].fields[2][k] = scalars[i].fields[3][k];
      }
    }
    output = 1;
  }
  delete[] eta;
  return output;
}

void thimble_system::simulate(int n_burn_in, int n_simulation)
{
  Njac = scalars.size()*Ntot; //setting the size of the Jacobian "matricies"
  NjacSquared = pow(Njac, 2);

  dcomp* state_storage = new dcomp[(Njac + 2)*n_simulation]; //This stores the data between updates and will be saved to a file
  std::ofstream data_storage;
  

  //initialising the fields
  for (int i = 0; i < scalars.size(); ++i)
  {
    scalars[i].initialise();
  }
  J.resize(Njac, Njac);
  J_conj.resize(Njac, Njac);
  J = calc_jacobian();
  J_conj = J.conjugate();
  S = calc_S(0);
  //setup is now complete, the Jacobian, it's conjugate, and it's determinant have been calculated, and the scalars are primed.
  for (int i = 0; i < n_burn_in; ++i)
  {
    update();
  }
  acceptance_rate = 0.;
  for(int i = 0; i < n_simulation; ++i)
  {
    acceptance_rate += update(); //calculating the acceptance rate from the return value of the update
    //storing the values of the fields, the action, and the Jacobian in a unified vector to be written to file later
    for (int k = 0; k < scalars.size(); ++k)
    {
      for (int r = 0; r < Ntot; ++r)
      {
        state_storage[i*(Njac + 2) + k*Ntot + r] = scalars[k].fields[0][r];
      }
    }
    state_storage[i*(Njac + 2) + scalars.size()*Ntot] = S;
    state_storage[i*(Njac + 2) + scalars.size()*Ntot + 1] = log(real(J.get_det())) +j*arg(J.get_det());
  }
  acceptance_rate /= n_simulation;

  data_storage.open(rel_path + file_name);
  data_storage << rng_seed << ",";
  for(int i = 0; i < scalars.size(); ++i)
  {
    for(int k = 0; k < Nx; ++k)
    {
      data_storage << real(scalars[i].field_0[k]) << "," << imag(scalars[i].field_0[k]) << ",";
    }
  }

  for (int i = 0; i < scalars.size(); ++i)
  {
    for (int k = 0; k < Nx; ++k)
    {
      data_storage << real(scalars[i].field_1[k]) << "," << imag(scalars[i].field_1[k]) << ",";
    }
  }
  data_storage << delta << "," << tau << "," << acceptance_rate << std::endl;
  for(int i = 0; i < n_simulation; ++i)
  {
    for (int k = 0; k < Njac + 1; ++k)
    {
      data_storage << real(state_storage[i*(Njac + 2) + k]) << "," << imag(state_storage[i*(Njac + 2) + k]) << ",";
    }
    data_storage << real(state_storage[i*(Njac + 2) + Njac + 1]) << "," << imag(state_storage[i*(Njac + 2) + Njac + 1]) << std::endl;
  }
  data_storage.close();
  delete[] state_storage;
}

void thimble_system::test()
{
  add_scalar_field(1.);
  scalars[0].initialise();
  matrix<dcomp> Jac = calc_jacobian();
}