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
  //empty constructor
}

//basic interaction, assuming no derivatives
dcomp interaction::base(int site, thimble_system &current_system, int field_type)
{
  dcomp interaction_contribution = coupling;

  for(int i = 0; i < powers.size(); ++i)
  {
    //applies all the fields raised to the relevant power
    interaction_contribution *= pow(current_system.scalars[i].fields[field_type][site], powers[i]); 
  }
  return interaction_contribution;
}

//differentiated with respect to the field (field) at site (site)
dcomp interaction::first_derivative(int site, int field, thimble_system &current_system, int field_type)
{
  dcomp interaction_contribution = coupling;
  //all the non-derivative fields up to the derivative
  for (int i = 0; i < field; ++i)
  {
    interaction_contribution *= pow(current_system.scalars[i].fields[field_type][site], powers[i]);
  }

  //contribution of the derivative field
  interaction_contribution *= double(powers[field])*pow(current_system.scalars[field].fields[field_type][site], powers[field] - 1);

  //contribution of all non-derivative fields from beyond the derivative field value
  for (int i = field + 1; i < powers.size(); ++i)
  {
    interaction_contribution *= pow(current_system.scalars[i].fields[field_type][site], powers[i]);
  }
  return interaction_contribution;
}

//second derivative set up similarly to the one above. Takes into account that unless the derivates and fields are the same then it should return zero.
dcomp interaction::second_derivative(int site, int field_1, int field_2, thimble_system &current_system, int field_type)
{
  dcomp interaction_contribution = coupling;
  if (field_1 == field_2)
  {
    for (int i = 0; i < field_1; ++i)
    {
      interaction_contribution *= pow(current_system.scalars[i].fields[field_type][site], powers[i]);
    }
    //contribution of the derivative field
    interaction_contribution *= double(powers[field_1])*double(powers[field_1] - 1)*pow(current_system.scalars[field_1].fields[field_type][site], powers[field_1] - 2);

    //contribution of all non-derivative fields from beyond the derivative field value
    for (int i = field_1 + 1; i < powers.size(); ++i)
    {
      interaction_contribution *= pow(current_system.scalars[i].fields[field_type][site], powers[i]);
    }
  }
  else
  {
    interaction_contribution = 0;
  }
  return interaction_contribution;
}

//scalar field class to be incorporated into a system class (class composition)*******************************************************
scalar_field::scalar_field(int x_dim, int t_dim, double system_dt, double system_dx, thimble_system* current_host) : occupation_number(new int[x_dim]),
  //member initialiser list 
  m(0),
  squareMass(0), 
  dt(system_dt),
  dx(system_dx),
  path(new double[current_host->Nrpath]),
  path_offset(new double[current_host->Nrpath]),
  is_flowed(false),
  field_0(new dcomp[current_host->Nx]),
  field_1(new dcomp[current_host->Nx]),
  field_2(new dcomp[current_host->Nx]),
  positive_time_site(new uint[current_host->Ntot]), //offset arrays to speed up computation
  positive_space_site(new uint[current_host->Ntot]),
  negative_time_site(new uint[current_host->Ntot]),
  negative_space_site(new uint[current_host->Ntot]),
  j(0,1),
  sigma(0.07071067812),
  delta(0.1),
  host(current_host)
{
  //class constructor
  fields[0] = new dcomp[current_host->Ntot]; //configuring the lattice arrays
  fields[1] = new dcomp[current_host->Ntot];
  fields[2] = new dcomp[current_host->Ntot];
  fields[3] = new dcomp[current_host->Ntot];
  fields[4] = new dcomp[current_host->Ntot];

  C[0] = new dcomp[current_host->Ntot]; //configuring Mou's constant arrays
  C[1] = new dcomp[current_host->Ntot];
  C[2] = new dcomp[current_host->Ntot];
  C[3] = new dcomp[current_host->Ntot];
  C[4] = new dcomp[current_host->Ntot];

  int n = 0;
  
  for (int i = 0; i < current_host->Ntot; ++i)
  {
    //zero initialising the fields
    for(int k = 0; k < 5; ++k)
    {
      fields[k][i] = 0;
    }
  }
  
  for (int i = 0; i < current_host->Nx; ++i)
  {
    //zero initialising the occupation number
    occupation_number[i] = 0;
  }

  for (int i = 0; i < (int) (current_host->Nrpath/2); ++i)
  {
    path[i] = dt;
    path[i + (int)(current_host->Nrpath/2)] = -1.*dt;
  }

  path_offset[0] = -1.*dt;
  for (int i = 0; i < current_host->Nrpath - 1; ++i)
  {
    path_offset[i + 1] = path[i];
  }

  //setting up the co-ordinate shifted arrays. Could combine them, but this won't be called much, so the legibility is prioritised over speed
  //I imagine the O3 flag for g++ automatically tidies them up anyway
  for (int k = 0; k < current_host->Nx; ++k)
  {
    for (int i = 0; i < current_host->Nrpath - 1; ++i)
    {
      positive_time_site[k*current_host->Nrpath + i] = k*current_host->Nrpath + i + 1;
    }
    positive_time_site[(k + 1)*current_host->Nrpath - 1] = k*current_host->Nrpath;
    
    for (int i = 1; i < current_host->Nrpath; ++i)
    {
      negative_time_site[k*current_host->Nrpath + i] = k*current_host->Nrpath + i - 1;
    }
    negative_time_site[k*current_host->Nrpath] = (k + 1)*current_host->Nrpath - 1;
  }
  
  for (int i = 0; i < current_host->Nrpath; ++i)
  {
    for (int k = 0; k < current_host->Nx - 1; ++k)
    {
      positive_space_site[k*current_host->Nrpath + i] = (k + 1)*current_host->Nrpath + i;
    }
    positive_space_site[(current_host->Nx - 1)*current_host->Nrpath + i] = i;
    
    for (int k = 1; k < current_host->Nx; ++k)
    {
      negative_space_site[k*current_host->Nrpath + i] = (k - 1)*current_host->Nrpath + i;
    }
    negative_space_site[i] = (current_host->Nx - 1)*current_host->Nrpath + i;
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

scalar_field::scalar_field(const scalar_field &obj) : occupation_number(new int[obj.host->Nx]), //object copy constructor
//*large* member initalisation list 
m(obj.m),
squareMass(obj.squareMass),
dt(obj.dt),
dx(obj.dx),
path(new double[obj.host->Nrpath]),
path_offset(new double[obj.host->Nrpath]),
is_flowed(obj.is_flowed),
field_0(new dcomp[obj.host->Nx]),
field_1(new dcomp[obj.host->Nx]),
field_2(new dcomp[obj.host->Nx]),
positive_time_site(new uint[obj.host->Ntot]),
positive_space_site(new uint[obj.host->Ntot]),
negative_time_site(new uint[obj.host->Ntot]),
negative_space_site(new uint[obj.host->Ntot]),
j(obj.j),
sigma(obj.sigma),
delta(obj.delta),
host(obj.host)
{
  for (int i = 0; i < 5; ++i)
  {
    fields[i] = new dcomp[obj.host->Ntot];
    C[i] = new dcomp[obj.host->Ntot];
  }

  for (int i = 0; i < host->Nx; ++i)
  {
    occupation_number[i] = obj.occupation_number[i]; //setting values for the arrays that are copied over from the original object
    field_0[i] = obj.field_0[i];
    field_1[i] = obj.field_1[i];
    field_2[i] = obj.field_2[i];
  }
  
  for(int i = 0; i < host->Ntot; ++i)
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

  for (int i = 0; i < host->Nrpath; ++i)
  {
    path[i] = obj.path[i];
    path_offset[i] = obj.path_offset[i];
  }
}

//sets the occupation number of the momentum modes of a distribution
void scalar_field::set_occupation_number(int new_occupation_number[])
{
  for (int i = 0; i < host->Nx; ++i)
  {
    occupation_number[i] = new_occupation_number[i];
  }
}

//same as above, but for a constant value
void scalar_field::set_occupation_number(int new_occupation_number)
{
  for(int i = 0; i < host->Nx; ++i)
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

void scalar_field::initialise(double a[], double b[], double c[], double d[])
{
  //this is called *after* the constructor, because we need to specify the mass
  double p, omega_p, omega_tilde, Omega_p, V;

  
  for (int i = 0; i < host->Nx; ++i)
  {
    field_0[i] = 0;
    field_1[i] = 0;
  } 

  V = host->Nx*dx;

  for (int i = 0; i < host->Nx; ++i)
	{
    //p = i*2.*pi/(Nx*dx);
    for(uint q = 0; q < host->Nx; ++q)
    {
      p = pow(2.*(1. - cos(q*2.*pi/(dx*host->Nx)))/pow(dx, 2), 0.5);
      omega_p = pow(pow(p,2) + pow(m,2),0.5);
      omega_tilde = acos(1 - pow(omega_p*dt,2)/2)/dt;
      Omega_p = sin(omega_p*dt)/dt; //configuring variables for this momentum
      if ((host->Nx - q)%host->Nx == 0)
      {
        //corner mode case
        field_0[i] += a[q]*pow(e,j*p*(i*dx))*pow(occupation_number[q] + 0.5,0.5)/pow(Omega_p,0.5);
        field_1[i] += pow(e,j*p*(i*dx))*(a[q]*cos(omega_tilde*dt) + c[q]*Omega_p*dt)*pow(occupation_number[q] + 0.5,0.5)/pow(Omega_p,0.5);
      }
      else
      {
        //bulk mode case
        field_0[i] += ((a[q] + j*b[q])*pow(e,j*p*(i*dx))/pow(2*Omega_p,0.5) + (a[q] - j*b[q])*pow(e,-1.*j*p*(i*dx))/pow(2*Omega_p,0.5))*pow(occupation_number[q] + 0.5,0.5);
        field_1[i] += pow(e, j*p*(i*dx))*((a[q] + j*b[q])*cos(omega_tilde*dt) + (c[q] + j*d[q])*Omega_p*dt)*pow(occupation_number[q] + 0.5, 0.5)/pow(2*Omega_p, 2) + conj(pow(e, j*p*(i*dx))*((a[q] + j*b[q])*cos(omega_tilde*dt) + (c[q] + j*d[q])*Omega_p*dt)*pow(occupation_number[q] + 0.5, 0.5)/pow(2*Omega_p, 2));
      }
    }
    field_0[i] = field_0[i]/V; //rescaling for the volume. hbar is taken to be one.
    field_1[i] = field_1[i]/V;

    //manual force to check values
    //field_0[i] = 0.8;
    //field_1[i] = 1.0;
  }

  for(int k = 0; k < host->Nx; ++k)
  {
      fields[2][0 + host->Nrpath*k] = field_1[k];
      fields[2][1 + host->Nrpath*k] = -1.0*dt*dt*(squareMass*fields[2][0 + host->Nrpath*k]) + 2.0*fields[2][0 + host->Nrpath*k] - field_0[k];
      for (int i = 1; i < (int) (host->Nrpath/2); ++i)
      {
        fields[2][i + host->Nrpath*k + 1] = -dt*dt*squareMass*fields[2][i + host->Nrpath*k] + 2.0*fields[2][i + host->Nrpath*k] - fields[2][i + host->Nrpath*k - 1];
      }
      
      for(int i = 0; i < (int) (host->Nrpath/2); ++i)
      {
        fields[2][host->Nrpath*(k + 1) - i - 1] = fields[2][host->Nrpath*k + i + 1]; //sets up the return leg of the contour
      }
  }

  //clearing the classical data from the first site, note the initial condition data is saved in field_0 and field_1
  for(int i = 0; i < host->Nx; ++i)
  {
    fields[2][i] = 0;
  }

  for (int i = 0; i < host->Ntot; ++i)
  {
    fields[0][i] = fields[2][i];
  }
  
  for(int i = 0; i < host->Nx; ++i)
  {
    field_2[i] = fields[2][i*host->Nrpath + 1];
  }
  calculate_C();
  /*
  for (int i = 0; i < host -> Nrpath; ++i)
  {
    printf("%f%+fi \n", real(fields[2][i]), imag(fields[2][i]));
  }
  */
}

void scalar_field::calculate_C()
{
  //setting up Mou's constant arrays
  for (int i = 0; i < host->Ntot; ++i)
  {
    int n = calc_n(i);
    C[0][i] = -1.*dx*j*(1/path[n] + 1/path_offset[n] + ((path[n] + path_offset[n])/2.)*(-2./pow(dx, 2) - squareMass));
    C[1][i] = j*dx/path[n];
    C[2][i] = j*dx/path_offset[n];
    C[3][i] = -1.*dx*j*(path[n] + path_offset[n])/(2.*pow(dx, 2));
    C[4][i] = 0.;
  }
  //edge terms for the edge effects
  for (int i = 0; i < host->Nx; ++i)
  {
    C[4][i*host->Nrpath] = -2.*j*dx*field_2[i]/dt;
    C[4][i*host->Nrpath + 1] = j*dx*field_1[i]/dt;
    C[4][(i + 1)*host->Nrpath - 1] = -1.*j*dx*field_1[i]/dt;
  }

  //applying the anti-periodic bounday terms
  for(int i = 0; i < host->Nx; ++i)
  {
    C[1][(i + 1)*host->Nrpath - 1] *= -1.;
    C[2][i*host->Nrpath] *= -1.;
  }
}

dcomp scalar_field::free_action(uint site, uint field_type)
{
  //Standard P^2 - m^2 action
  dcomp S = 0;
  int n = calc_n(site);
  if(n != host->Nrpath - 1)
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

dcomp scalar_field::free_action_derivative(uint site, uint field_type)
{
  //derivative of the above action
  dcomp dS = C[0][site]*fields[field_type][site] + C[1][site]*fields[field_type][positive_time_site[site]] + C[2][site]*fields[field_type][negative_time_site[site]]
    + C[3][site]*(fields[field_type][positive_space_site[site]] + fields[field_type][negative_space_site[site]]) + C[4][site];
  return dS;
}

dcomp scalar_field::free_action_second_derivative(uint site_1, uint site_2)
{
  //second derivative calculation
  dcomp ddS = C[0][site_1]*dd(site_1, site_2) + C[1][site_1]*dd(positive_time_site[site_1], site_2) + C[2][site_1]*dd(negative_time_site[site_1], site_2)
    + C[3][site_1]*(dd(positive_space_site[site_1], site_2) + dd(negative_space_site[site_1], site_2));
  return ddS;
}

void scalar_field::set_dt(double new_dt)
{
  for (int i = 0; i < host->Nrpath; ++i)
  {
    path[i] *= new_dt/dt; //rescales the path to account for a new lattice spacing
    path_offset[i] *= new_dt/dt;
  }
  dt = new_dt;
}

//decomposing the total lattice position into the single timeslice position (for the dt array)
uint scalar_field::calc_n(uint site)
{
  uint n = site%host->Nrpath;
  return n;
}

//decomposing the composite lattice position to calculate the spacesclice position
int scalar_field::calc_x(int site)
{
  int x = int((site - calc_n(site))/host->Nrpath);
  return x;
}