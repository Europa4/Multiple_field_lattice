#include "Prot.h"
#include <unistd.h>
#include <random>
//c_phi constructor
c_phi::c_phi(double field_mass, double tau, double Lambda, double Delta, unsigned long int seed)
{
  call_constructor(field_mass, tau, Lambda, Delta, seed);
}

void c_phi::call_constructor(double field_mass, double tau, double Lambda, double Delta, unsigned long int seed)
{
  //function called once from the c_phi constructor
    mass = field_mass; //initialising members of the class
    squareMass = pow(mass, 2);
    flowTime = tau;
    lambda = Lambda;
    isFlowed = false;
    delta = Delta;
    sigma_squared = 0.5*pow(delta, 2);
    sigma = pow(0.5, 0.5)*delta;
    delta_squared = pow(delta, 2);
    rng_seed = seed;
    rel_path = "";
    
    
    const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    my_rngPointer = gsl_rng_alloc (T);
    gsl_rng_set(my_rngPointer, seed);
    
    double a[Nx],b[Nx],c[Nx],d[Nx],np[Nx];
    double p, omega_p, omega_tilde, Omega_p, V;
    
    
    //determining the number of timesteps for the ODE solvers from the flow time  
    h = 0.02; //sets the base size 
    number_of_timesteps = int(ceil(flowTime/h)); //calculates how many steps this corresponds to (overshooting in the case of it not being exact)
    h = flowTime/number_of_timesteps; //readusting the size of h to prevent overshooting
    
    for (int i = 0; i < Nx; ++i)
    {
        a[i] = gsl_ran_gaussian(my_rngPointer, 1); //Yeah naming a variable a, b, c, d isn't helpful, but it's what we call them in the maths. Besides, they're local to this function.
        b[i] = gsl_ran_gaussian(my_rngPointer, 1);
        c[i] = gsl_ran_gaussian(my_rngPointer, 1);
        d[i] = gsl_ran_gaussian(my_rngPointer, 1);
        np[i] = 0;
        phi0[i] = 0;
        phi1[i] = 0;
    } //random number arrays for Mou's initial conditions, the population number, and the initial states of the phi field


    V = Nx*dx; //lattice volume


    for (int i = 0; i < Nx; ++i)
    {
        p = i*2.*pi/(Nx*dx);
        omega_p = pow(pow(p,2) + pow(mass,2),0.5);
        omega_tilde = acos(1 - pow(omega_p*deltaT,2)/2)/deltaT;
        Omega_p = sin(omega_p*deltaT)/deltaT; //configuring variables for this momentum
        if ((Nx - i)%Nx == 0)
        {
            //corner mode case
            phi0[i] += a[i]*pow(e,j*p*(i*dx))*pow(np[i] + 0.5,0.5)/pow(Omega_p,0.5);
            phi1[i] += pow(e,j*p*(i*dx))*(a[i]*cos(omega_tilde*deltaT) + c[i]*Omega_p*deltaT)*pow(np[i] + 0.5,0.5)/pow(Omega_p,0.5);
        }
        else
        {
            //bulk mode case
            phi0[i] += ((a[i] + j*b[i])*pow(e,j*p*(i*dx))/pow(2*Omega_p,0.5) + (c[i] - j*b[i])*pow(e,-1.*j*p*(i*dx))/pow(2*Omega_p,0.5))*pow(np[i] + 0.5,0.5);
            phi1[i] += pow(e,j*p*(i*dx))*(a[i]*cos(omega_tilde*deltaT) + c[i]*Omega_p*deltaT)*pow(np[i] + 0.5,0.5)/pow(Omega_p,0.5);
        }

        phi0[i] = phi0[i]/V; //rescaling for the volume. hbar is taken to be one.
        phi1[i] = phi1[i]/V;

        //manual force to check values
        //phi0[i] = 0.8;
        //phi1[i] = 1.0;
    }
    
    
    //Now that we have our initial conditions we can propogate the classical solution through the lattice
    for(int k = 0; k < Nx; ++k)
    {
        baseField[0 + Nrpath*k] = phi1[k];
        baseField[1 + Nrpath*k] = -1.0*deltaT*deltaT*(squareMass*baseField[0 + Nrpath*k] + lambda/6*pow(baseField[0 + Nrpath*k],3)) + 2.0*baseField[0 + Nrpath*k] - phi0[k];
        for (int i = 1; i < (int) (Nrpath/2); ++i)
        {
            baseField[i + Nrpath*k + 1] = -deltaT*deltaT*(squareMass*baseField[i + Nrpath*k] + lambda/6*pow(baseField[i + Nrpath*k],3)) + 2.0*baseField[i + Nrpath*k] - baseField[i + Nrpath*k - 1];
        }
        
        for(int i = 0; i < (int) (Nrpath/2); ++i)
	{
	  baseField[Nrpath*(k + 1) - i - 1] = baseField[Nrpath*k + i + 1]; //sets up the return leg of the contour
	}
    }
    //Calculating the C array
    for (int k = 0; k < Nx; ++k)
    {
      for (int i = 0; i < Nrpath; ++i)
      {
	C[i + k*Nrpath][0] = -1.*j*dx*(1./dt[i] + 1./dt_minus_one[i] + ((dt[i] + dt_minus_one[i])/2.)*(-2./pow(dx,2) - squareMass));
	C[i + k*Nrpath][1] = j*dx/dt[i];
	C[i + k*Nrpath][2] = j*dx/dt_minus_one[i];
	C[i + k*Nrpath][3] = -1.*j*dx*(dt[i] + dt_minus_one[i])/(2.*pow(dx,2));
	C[i + k*Nrpath][4] = j*dx*(dt[i] + dt_minus_one[i])*lambda/12.;
	C[i + k*Nrpath][5] = 0;
	C[i + k*Nrpath][6] = 0;
      }
      //Special case C terms
      C[k*Nrpath][5] = -2.0*j*dx*baseField[1 + k*Nrpath]/deltaT;
      C[k*Nrpath][6] = j*dx*lambda*deltaT*baseField[k*Nrpath];
      
      C[1 + k*Nrpath][5] = j*dx*baseField[k*Nrpath]/deltaT;
      
      C[(k + 1)*Nrpath - 1][5] = -1.0*j*dx*baseField[k*Nrpath]/deltaT;
    }

    //Taking into account the antiperiodic boundary conditions
    for (int i = 0; i < Nx; ++i)
    {
        C[(i+1)*Nrpath - 1][1] = -1.*C[(i+1)*Nrpath - 1][1];
        C[i*Nrpath][2] = -1.*C[i*Nrpath][2];
    }
    
    //clearing the classical data from the first phi site to be used in the simulation
    for (int i = 0; i < Nx; ++i)
    {
      baseField[i*Nrpath] = 0.0;
    }
    
    //populating the shifted co-ordinate arrays
    //It's possible, easy even, to optimise these a little better by combining the loops. However, this should be called once for each field which shouldn't be very many, and this keeps it more legible.
    //with the -O3 flag I imagine g++ auto optimises that out anyway
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

void c_phi::calc_S(bool proposal)
{
  int site; //short hand variable for current site. Equivilent to r in the notes.
  dcomp working_field[Ntot];
  dcomp internal_S = 0.;
  if (proposal == false)
  {
    for (int i = 0; i < Ntot; ++i)
    {
      working_field[i] = flowedField[i];
    }
  }
  else
  {
    for (int i = 0; i < Ntot; ++i)
    {
      working_field[i] = proposed_flowedField[i];
    }
  }
  for (int k = 0; k < Nx; ++k)
  {
    for (int i = 0; i < Nrpath - 1; ++i)
    {
      //applying the sum from eqn 18 in Mou's notes
      site = k*Nrpath + i;
      internal_S += -1.*C[site][1]*pow(working_field[positive_time_site[site]] - working_field[site], 2)/2.
      + pow(dx,2)*C[site][3]*(-1.*pow(working_field[positive_space_site[site]] - working_field[site], 2)/(2.*pow(dx,2)) -1.*squareMass*pow(working_field[site],2)/2. - lambda*pow(working_field[site],4)/24.) 
      + C[site][5]*working_field[site] 
      + C[site][6]*pow(working_field[site],3)/3.;
    }
    //taking into account the anti periodic boundary conditions on the last time site of each loop
    site = (k + 1)*Nrpath - 1;
    internal_S += C[site][1]*pow(working_field[positive_time_site[site]] + working_field[site], 2)/2. + pow(dx,2)*C[site][3]*(-1.*pow(working_field[positive_space_site[site]] - working_field[site], 2)/(2.*pow(dx,2))
	 -1.*squareMass*pow(working_field[site],2)/2. - lambda*pow(working_field[site],4)/24.) + C[site][5]*working_field[site] + C[site][6]*pow(working_field[site],3)/3.;
  }
  
  if (proposal == false)
  {
    S = internal_S;
  }
  else
  {
    proposed_S = internal_S;
  }
}

dcomp c_phi::get_Action()
{
  //returns the action, prints an error if the action hasn't been calculated yet.
  if (isFlowed == true)
  {
    return S;
  }
  else
  {
    printf("Error: Attempting to call action of an unflowed field. Please call scalar_flow() or calc_jacobian() before calling the action. \n" );
    exit(0);
  }
}

dcomp c_phi::calc_dS(int site, dcomp integral_field[Ntot])
{
    //calculates the derivative of the action with respect to the field at "site".
  dcomp dS = 0;
  dS = C[site][0]*integral_field[site] 
  + C[site][1]*integral_field[positive_time_site[site]] 
  + C[site][2]*integral_field[negative_time_site[site]] 
  + C[site][3]*(integral_field[negative_space_site[site]] + integral_field[positive_space_site[site]]) 
  + C[site][4]*pow(integral_field[site],3) 
  + C[site][5] 
  + C[site][6]*pow(integral_field[site],2)
  ;
  return dS;
}

void c_phi::scalar_flow(bool proposal)
{
    dcomp working_field[Ntot]; //this is the field to be used for the actual flowing process
    dcomp ajustment_field[Ntot];
    dcomp k1[Ntot], k2[Ntot], k3[Ntot], k4[Ntot]; //Runge-Kutta constants
    
    if(proposal == false)
    {
      for (int i = 0; i < Ntot; ++i) //initialising the field to be flowed
      {
	working_field[i] = baseField[i];
      }
    }
    else
    {
      for (int i = 0; i < Ntot; ++i)
      {
	working_field[i] = proposed_baseField[i];
      }
    }
    
    for (int i = 0; i < number_of_timesteps; ++i) //evolving the ODEs forward using the RK4 method. 
    {
        for (int k = 0; k < Ntot; ++k)
        {
            k1[k] = h*conj(calc_dS(k,working_field)); //calculating each constant for the RK4 method
	    ajustment_field[k] = working_field[k] + k1[k]/2.; //ajustment_field is used as the updated argument for the RHS function (see the RK method wiki page for details)
	    //printf("k1[%i] = %f%+f \n",k, real(k1[k]), imag(k1[k]));
        }
        
        for (int k = 0; k < Ntot; ++k)
	{
	  k2[k] = h*conj(calc_dS(k,ajustment_field));
	}
	for (int k = 0; k < Ntot; ++k)
	{
	  ajustment_field[k] = working_field[k] + k2[k]/2.;
	}
	
	for (int k = 0; k < Ntot; ++k)
	{
	  k3[k] = h*conj(calc_dS(k,ajustment_field));
	}
	
	for (int k = 0; k < Ntot; ++k)
	{
	  ajustment_field[k] = working_field[k] + k3[k];
	}
	
	for (int k = 0; k < Ntot; ++k)
	{
	  k4[k] = h*conj(calc_dS(k,ajustment_field));
	}
	
	for (int k = 0; k < Ntot; ++k)
	{
	  working_field[k] += (k1[k] + 2.*k2[k] + 2.*k3[k] + k4[k])/6.;
	}
    }
    
    //now that the field has been flowed, it needs to be transfered to the flowedField member
    if (proposal == false)
    {
      for (int i = 0; i < Ntot; ++i)
      {
	flowedField[i] = working_field[i];
      }
    }
    else
    {
      for (int i = 0; i < Ntot; ++i)
      {
	proposed_flowedField[i] = working_field[i];
      }
    }
    isFlowed = true;
}

dcomp c_phi::calc_ddS(int r, int c, dcomp integral_field[Ntot], dcomp integral_J[Ntot][Ntot])
{
  dcomp ddS;
  ddS = (C[r][0] 
    + 3.*C[r][4]*pow(integral_field[r],2) 
    + 2.*C[r][6]*integral_field[r]
    )*integral_J[r][c] 
  + C[r][1]*integral_J[positive_time_site[r]][c] 
  + C[r][2]*integral_J[negative_time_site[r]][c]
  + C[r][3]*(integral_J[positive_space_site[r]][c] + integral_J[negative_space_site[r]][c])
  ;
   return ddS;
}

void c_phi::calc_jacobian(bool proposal)
{
  dcomp working_field[Ntot]; //field used for the actual field flow process
  dcomp ajustment_field[Ntot]; //update term for the RK4 method
  dcomp k1_field[Ntot], k2_field[Ntot], k3_field[Ntot], k4_field[Ntot]; //field RK constants
  
  dcomp working_jacobian[Ntot][Ntot]; //working jacobian term for the flow process
  dcomp ajustment_jacobian[Ntot][Ntot]; 
  dcomp k1_jacobian[Ntot][Ntot], k2_jacobian[Ntot][Ntot], k3_jacobian[Ntot][Ntot], k4_jacobian[Ntot][Ntot];
  
  if (proposal == false)
  {
    for (int i = 0; i < Ntot; ++i)
    {
      working_field[i] = baseField[i];
    }
  }
  else
  {
    for (int i = 0; i < Ntot; ++i)
    {
      working_field[i] = proposed_baseField[i]; //setting up initial conditions based on whether or not this is a proposal or not
    }
  }
  
  
  for (int i = 0; i < Ntot; ++i)
  {
    for (int r = 0; r < Ntot; ++r)
    {
      if (r == i)
      {
	working_jacobian[i][r] = 1.0;
      }
      else
      {
	working_jacobian[i][r] = 0.0;
      }
    }
  }
  
  for (int i = 0; i < number_of_timesteps; ++i)
  {
    for(int r = 0; r < Ntot; ++r)
    {
        k1_field[r] = h*conj(calc_dS(r,working_field));
        ajustment_field[r] = working_field[r] + k1_field[r]/2.;
        
        for (int c = 0; c < Ntot; ++c)
        {
            k1_jacobian[r][c] = h*conj(calc_ddS(r, c, working_field, working_jacobian));
            ajustment_jacobian[r][c] = working_jacobian[r][c] + k1_jacobian[r][c]/2.;
        }
      
    }
    
    for (int r = 0; r < Ntot; ++r)
    {   
        k2_field[r] = h*conj(calc_dS(r,ajustment_field));

        for (int c = 0; c < Ntot; ++c)
        {
            k2_jacobian[r][c] = h*conj(calc_ddS(r, c, ajustment_field, ajustment_jacobian));
        }
    }
    
    for (int r = 0; r < Ntot; ++r)
    {
        ajustment_field[r] = working_field[r] + k2_field[r]/2.;
        for (int c = 0; c < Ntot; ++c)
        {
            ajustment_jacobian[r][c] = working_jacobian[r][c] + k2_jacobian[r][c]/2.;
        }
    }
    
    for (int r = 0; r < Ntot; ++r)
    {
        k3_field[r] = h*conj(calc_dS(r,ajustment_field));
        for (int c = 0; c < Ntot; ++c)
        {
            k3_jacobian[r][c] = h*conj(calc_ddS(r, c, ajustment_field, ajustment_jacobian));
        }
    }
    
    for (int r = 0; r < Ntot; ++r)
    {
        ajustment_field[r] = working_field[r] + k3_field[r];
        for (int c = 0; c < Ntot; ++c)
        {
            ajustment_jacobian[r][c] = working_jacobian[r][c] + k3_jacobian[r][c];
        }
    }
    
    for (int r = 0; r < Ntot; ++r)
    {
        k4_field[r] = h*conj(calc_dS(r, ajustment_field));
        for (int c = 0; c < Ntot; ++c)
        {
            k4_jacobian[r][c] = h*conj(calc_ddS(r, c, ajustment_field, ajustment_jacobian));
        }
    }

    for (int r = 0; r < Ntot; ++r)
    {
        working_field[r] += (k1_field[r] + 2.*k2_field[r] + 2.*k3_field[r] + k4_field[r])/6.;
        for (int c = 0; c < Ntot; ++ c)
        {
            working_jacobian[r][c] += (k1_jacobian[r][c] + 2.*k2_jacobian[r][c] + 2.*k3_jacobian[r][c] + k4_jacobian[r][c])/6.;
        }
    }
  }
  if(proposal == false)
  {
    for (int r = 0; r < Ntot; ++r)
    {
	flowedField[r] = working_field[r];
	for (int c = 0; c < Ntot; ++c)
	{
	    J[r][c] = working_jacobian[r][c];
	}
    }
  }
  else
  {
    for (int r = 0; r < Ntot; ++r)
    {
	proposed_flowedField[r] = working_field[r];
	for (int c = 0; c < Ntot; ++c)
	{
	    proposed_J[r][c] = working_jacobian[r][c];
	}
    }
  }
  isFlowed = true;
}

void c_phi::invert_jacobian(bool proposal)
{
  // inverts the jacobian
  int s;
  gsl_permutation *p = gsl_permutation_alloc(Ntot);
  gsl_matrix_complex * mJ = gsl_matrix_complex_alloc(Ntot, Ntot); //setting up GSL matricies
  gsl_matrix_complex * invmJ = gsl_matrix_complex_alloc(Ntot, Ntot);
  if (proposal == false)
  {
    for (int r = 0; r < Ntot; ++r)
    {
      for (int c = 0; c < Ntot; ++c)
      {
	gsl_matrix_complex_set(mJ, r, c, gsl_complex_rect(real(J[r][c]), imag(J[r][c]))); //populating the GSL matrix with J
      }
    }
    
    gsl_linalg_complex_LU_decomp(mJ, p, &s); //splitting into triangle matrices
    gsl_linalg_complex_LU_invert(mJ, p, invmJ); //inverting the new triangle matricies
    
    for (int r = 0; r < Ntot; ++r)
    {
      for (int c = 0; c < Ntot; ++c)
      {
	invJ[r][c] = GSL_REAL(gsl_matrix_complex_get(invmJ, r, c)) + j*GSL_IMAG(gsl_matrix_complex_get(invmJ, r, c));
      }
    }
  }
  else
  {
    for (int r = 0; r < Ntot; ++r)
    {
      for (int c = 0; c < Ntot; ++c)
      {
	gsl_matrix_complex_set(mJ, r, c, gsl_complex_rect(real(proposed_J[r][c]), imag(proposed_J[r][c]))); //populating the GSL matrix with J
      }
    }
    
    gsl_linalg_complex_LU_decomp(mJ, p, &s); //splitting into triangle matrices
    gsl_linalg_complex_LU_invert(mJ, p, invmJ); //inverting the new triangle matricies
    
    for (int r = 0; r < Ntot; ++r)
    {
      for (int c = 0; c < Ntot; ++c)
      {
	proposed_invJ[r][c] = GSL_REAL(gsl_matrix_complex_get(invmJ, r, c)) + j*GSL_IMAG(gsl_matrix_complex_get(invmJ, r, c));
      }
    }
  }
  gsl_permutation_free(p);
  gsl_matrix_complex_free(mJ);
  gsl_matrix_complex_free(invmJ); //freeing the GSL pointers
}

void c_phi::calc_detJ(bool proposal)
{ //calculates the determinant of J, propsal bool controls whether its the original or proposed scalar field to be calculated from
  int s;
  gsl_permutation * p = gsl_permutation_alloc(Ntot);
  gsl_matrix_complex * mJ = gsl_matrix_complex_alloc(Ntot,Ntot);
  gsl_complex det_gsl;
  if(proposal == false)
  {
    for(int r = 0; r < Ntot; ++r)
    {
      for(int c = 0; c < Ntot; ++c)
      {
	gsl_matrix_complex_set(mJ, r, c, gsl_complex_rect(real(J[r][c]), imag(J[r][c]))); //populating the GSL matrix
      }
    }
  }
  else
  {
    for(int r = 0; r < Ntot; ++r)
    {
      for(int c = 0; c < Ntot; ++c)
      {
	gsl_matrix_complex_set(mJ, r, c, gsl_complex_rect(real(proposed_J[r][c]), imag(proposed_J[r][c])));
      }
    }
  }
  
  gsl_linalg_complex_LU_decomp(mJ, p, &s);
  det_gsl = gsl_linalg_complex_LU_det(mJ, s); //this is the actual inversion
  
  if(proposal == false)
  {
    detJ = GSL_REAL(det_gsl) + j*GSL_IMAG(det_gsl);
  }
  else
  {
    proposed_detJ = GSL_REAL(det_gsl) + j*GSL_IMAG(det_gsl); //returning the values to an array rather than a GSL matrix
  }
  gsl_permutation_free(p);
  gsl_matrix_complex_free(mJ); //freeing the pointers
}

void c_phi::simulate(int n_burn_in, int n_chain, int decorrelation_length, int file_number)
{
  //this controls the monte carlo simulations
  ofstream dataStore;
  if (n_chain % decorrelation_length == 0)
  {
    int n_iterate = round(n_chain/decorrelation_length);
    int counter = 0;
    int Nstore = Ntot + 2;
    dcomp* state_storage = new dcomp[n_iterate*Nstore]; //the first Ntot values are phi, then S, then J
    
    calc_jacobian(); //calcuating initial parameters
    invert_jacobian();
    calc_detJ();
    log_det_J = log(abs(detJ));
    angle_det_J = arg(detJ);
    calc_S();
    
    for (int r = 0; r < Ntot; ++r)
    {
      for (int c = 0; c < Ntot; ++c)
      {
	conj_J[r][c] = conj(J[c][r]);
      }
    }
    
    for (int i = 0; i < n_burn_in; ++i)
    {
      update(); //performing a burn in, no data logging
    }
    
    for (int i = 0; i < n_chain; ++i)
    {
      update(); //simulating. Log data every decorrelation_length steps
      if (i % decorrelation_length == 0)
      {
	for (int r = 0; r < Ntot; ++r)
	{
	  state_storage[counter*Nstore + r] = flowedField[r];
	}
	state_storage[counter*Nstore + Ntot] = S;
	state_storage[counter*Nstore + Ntot + 1] = log_det_J + j*angle_det_J;
	++counter;
      }
    }
    dataStore.open(rel_path + "phi_" + to_string(file_number));
    dataStore << rng_seed << ",";
    for (int i = 0; i < Nx; ++i)
    {
      dataStore << real(phi0[i]) << "," << imag(phi0[i]) << ",";
    }
    for (int i = 0; i < Nx; ++i)
    {
      dataStore << real(phi1[i]) << "," << imag(phi1[i]) << ",";
    }
    dataStore << delta << "," << flowTime << endl;
    for (int i = 0; i < n_iterate; ++i)
    {
      for (int k = 0; k < Nstore - 1; ++k)
      {
	dataStore << real(state_storage[i*Nstore + k]) << "," << imag(state_storage[i*Nstore + k]) << ",";
      }
      dataStore << real(state_storage[i*Nstore + Nstore - 1]) << "," << imag(state_storage[i*Nstore + Nstore - 1]) << endl;
    }
    dataStore.close();
    delete [] state_storage;
  }
  else
  {
    printf("Error: decorrelation length not a factor of chain length");
    exit(0);
  }
}

void c_phi::update()
{
  //this is the actual working subroutine
  dcomp eta[Ntot], raw_delta[Ntot], DeltaJ[Ntot], JDelta[Ntot], proposal_DeltaJ[Ntot], proposal_JDelta[Ntot];
  dcomp Delta[Ntot];
  double log_proposal, deltaS, DeltaJJDelta, proposal_DeltaJJDelta, exponent, check;
  bool proposal = true;

  for (int i = 0; i < Ntot; ++i)
  {
    eta[i] = gsl_ran_gaussian(my_rngPointer, sigma) + j*gsl_ran_gaussian(my_rngPointer, sigma);
  }
  matrix_multiplication(raw_delta, invJ, eta);//using the inverse of J to transport the elements back to the real manifold
  
  for (int i = 0; i < Ntot; ++i)
  {
    Delta[i] = real(raw_delta[i]); //taking only the real parts
    proposed_baseField[i] = baseField[i] + Delta[i]; //creating the proposed base configuration
  }
  calc_jacobian(proposal); //calculating all the new parameters for the proposal
  calc_detJ(proposal); 
  calc_S(proposal);
  log_proposal = log(abs(proposed_detJ));
  
  
  for (int r = 0; r < Ntot; ++r)
  {
    for (int c = 0; c < Ntot; ++c)
    {
      proposed_conj_J[r][c] = conj(proposed_J[c][r]); //calculating the conjugate transpose of the jacobian
    }
  }

  deltaS = real(proposed_S - S);
  matrix_multiplication(DeltaJ, Delta, J); //calculating acceptance parameters, firstly for the existing state
  matrix_multiplication(JDelta, conj_J, Delta);
  DeltaJJDelta = real(dot_product(DeltaJ, JDelta));

  matrix_multiplication(proposal_DeltaJ, Delta, proposed_J); //and now for the new configuration
  matrix_multiplication(proposal_JDelta, proposed_conj_J, Delta); 
  proposal_DeltaJJDelta = real(dot_product(proposal_DeltaJ, proposal_JDelta));

  exponent = -1.*deltaS + 2.*log_proposal - 2.*log_det_J + DeltaJJDelta/delta_squared - proposal_DeltaJJDelta/delta_squared;
  check = gsl_rng_uniform(my_rngPointer);
  
  if (exp(exponent) > check) //checking against the acceptance parameters
  {
    //proposal accepted
    log_det_J = log_proposal;
    angle_det_J = arg(proposed_detJ);
    S = proposed_S;
    detJ = proposed_detJ;
    for (int r = 0; r < Ntot; ++r)
    {
      baseField[r] = proposed_baseField[r];
      flowedField[r] = proposed_flowedField[r];
      for (int c = 0; c < Ntot; ++c)
      {
	J[r][c] = proposed_J[r][c];
	conj_J[r][c] = proposed_conj_J[r][c];
      }
    }
    invert_jacobian();
  }
}


void c_phi::set_path(string new_path)
{
  //this lets you set the path for saving the data to a different directory.
  string boost_path = new_path;
  boost_path = boost_path.erase(boost_path.size() - 1); //removes the trailing / from the path before passing to the boost library
  boost::filesystem::path p(boost_path);
  if (!exists(p))
  {
    create_directory(p); //creates the directory if it doesn't yet exist
  }
  rel_path = new_path;
}

c_phi::~c_phi()
{
  //clears the RNG setup
  gsl_rng_free(my_rngPointer);
}
