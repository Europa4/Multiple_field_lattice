#include "thimble_lattice.h"
//scalar field class to be incorporated into a system class (class composition)

scalar_field::scalar_field(int x_dim, int t_dim, gsl_rng * rngPointer) : occupation_number(new int[x_dim]),
  Nx(x_dim), //member initialiser list 
  Nt(t_dim), 
  Npath(2*Nt), 
  Nrpath(Npath - 4), //parameters of the lattice
  Ntot(Nrpath*Nx), 
  m(0),
  squareMass(0), 
  dt(0),
  dx(0),
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
}

void scalar_field::set_occupation_number(int new_occupation_number[])
{
  for (int i = 0; i < Nx; ++i)
  {
    occupation_number[i] = new_occupation_number[i];
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
        field_0[i] += a[i]*pow(e,j*p*(i*dx))*pow(np[i] + 0.5,0.5)/pow(Omega_p,0.5);
        field_1[i] += pow(e,j*p*(i*dx))*(a[i]*cos(omega_tilde*dt) + c[i]*Omega_p*dt)*pow(np[i] + 0.5,0.5)/pow(Omega_p,0.5);
    }
    else
    {
        //bulk mode case
        field_0[i] += ((a[i] + j*b[i])*pow(e,j*p*(i*dx))/pow(2*Omega_p,0.5) + (c[i] - j*b[i])*pow(e,-1.*j*p*(i*dx))/pow(2*Omega_p,0.5))*pow(np[i] + 0.5,0.5);
        field_1[i] += pow(e,j*p*(i*dx))*(a[i]*cos(omega_tilde*dt) + c[i]*Omega_p*dt)*pow(np[i] + 0.5,0.5)/pow(Omega_p,0.5);
    }

    field_0[i] = field_0[i]/V; //rescaling for the volume. hbar is taken to be one.
    field_1[i] = field_1[i]/V;

    //manual force to check values
    field_0[i] = 0.8;
    field_1[i] = 1.0;
  }
  
    for(int k = 0; k < Nx; ++k)
    {
        base_field[0 + Nrpath*k] = field_1[k];
        base_field[1 + Nrpath*k] = -1.0*dt*dt*(squareMass*base_field[0 + Nrpath*k] + 2.0*base_field[0 + Nrpath*k] - field_0[k]);
        for (int i = 1; i < (int) (Nrpath/2); ++i)
        {
          base_field[i + Nrpath*k + 1] = -dt*dt*(squareMass*base_field[i + Nrpath*k] + 2.0*base_field[i + Nrpath*k] - base_field[i + Nrpath*k - 1]);
        }
        
        for(int i = 0; i < (int) (Nrpath/2); ++i)
	{
	  base_field[Nrpath*(k + 1) - i - 1] = base_field[Nrpath*k + i + 1]; //sets up the return leg of the contour
	}
  

  //clearing the classical data from the first site
  for(int i = 0; i < Nx; ++i)
  {
    base_field[i] = 0;
  }

  //setting up the co-ordinate shifted arrays. Could combine them, but this won't be called much, so the legibility is prioritised over speed
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

//*************************************thimble_system***************************

thimble_system::thimble_system(int x_dim, int t_dim, double flow_time, long unsigned int seed) : Nx(x_dim), Nt(t_dim), tau(flow_time), Npath(2*Nt), Nrpath(Npath - 4), Ntot(Nrpath*Nx), phi(Nx, Nt), rng_seed(seed)
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
