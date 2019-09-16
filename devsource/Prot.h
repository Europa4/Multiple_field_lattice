#ifndef PROT_H_INCLUDED
#define PROT_H_INCLUDED

#include <iostream>
#include <fstream>
//#include <ofstream>
#include <complex>
#include <math.h>

#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <boost/filesystem.hpp>
#include <stdio.h>
#include <string>
#include <time.h>
//#include <random>

#include "thimble_lattice.h"

using namespace std;
typedef complex<double> dcomp;

const int Nx = 1;
const int Nt = 10;
const int Npath = 2*Nt;
const int Nrpath = Npath - 4;
const int Ntot = Nrpath*Nx;
const dcomp j(0,1);
//const double pi = 3.14159265359;
//const double e = 2.71828182846;

const double dx = 1.0;
const double deltaT = 0.75;
extern double dt[Nrpath];
extern double dt_minus_one[Nrpath];


//************************Classes********************

class c_phi
{
    //phi class, c_phi constructor is in c_phi.cpp, along with any methods not given here
  private:
    bool isFlowed;
    double mass;
    double squareMass;
    double flowTime;
    double lambda; //interaction parameter
    double h; // timestep for the ode solver
    double delta;//the actual delta used in the maths
    double delta_squared;
    double sigma; //standard deviation for the distributions
    double sigma_squared;//0.5*delta**2. Commonly used constant
    double log_det_J;
    double angle_det_J;
    dcomp S;
    dcomp proposed_S;
    int number_of_timesteps;
    unsigned long int rng_seed;
    gsl_rng * my_rngPointer;
    string rel_path;
    
    
    int positive_time_site[Ntot]; // co-ordinate positions for neighbour sites in a 2-D representation.
    int negative_time_site[Ntot]; 
    int positive_space_site[Ntot]; 
    int negative_space_site[Ntot];
    
    
    dcomp calc_dS(int site, dcomp integral_field[Ntot]);
    dcomp calc_ddS(int r, int c, dcomp integral_field[Ntot], dcomp integral_J[Ntot][Ntot]);
    dcomp conj_J[Ntot][Ntot], proposed_conj_J[Ntot][Ntot];
    
    
    void invert_jacobian(bool proposal = false);
    void calc_detJ(bool proposal = false);
    
  protected:
    void call_constructor(double field_mass, double tau, double Lambda, double delta, unsigned long int seed);
    

  public:
    dcomp C[Ntot][7];//Mou's constant array
    dcomp baseField[Ntot], flowedField[Ntot], proposed_baseField[Ntot], proposed_flowedField[Ntot];
    dcomp phi0[Nx], phi1[Nx], J[Ntot][Ntot], invJ[Ntot][Ntot], proposed_J[Ntot][Ntot], proposed_invJ[Ntot][Ntot], proposed_detJ, detJ;
    
    void calc_S(bool proposal = false);
    bool get_Flowed() {return isFlowed;}
    dcomp get_Action();
    void scalar_flow(bool proposal = false);
    void calc_jacobian(bool proposal = false);
    void simulate(int n_burn_in, int n_simulate, int decorrelation_length = 1, int file_number = 0);
    void update();
    void simulate_lambda_0(int n_burn_in, int n_simulate, int decorrelation_length = 1, int file_number = 0);
    void update_lambda_0();
    void set_path(string new_path);
    
    //constructor and destructor defentions
    c_phi(double field_mass, double tau, double Lambda, double Delta, unsigned long int seed = 0);
    ~c_phi();
};

//************************Predecs********************

void initalise_dt();
void matrix_multiplication(dcomp delta[Ntot], dcomp Mat[Ntot][Ntot], dcomp eta[Ntot]);
void matrix_multiplication(dcomp delta[Ntot], dcomp eta[Ntot], dcomp Mat[Ntot][Ntot]);
dcomp dot_product(dcomp v1[Ntot], dcomp v2[Ntot]);
#endif // PROT_H_INCLUDED
