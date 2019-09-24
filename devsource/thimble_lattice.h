#ifndef THIMBLE_LATTICE_H_INCLUDED
#define THIMBLE_LATTICE_H_INCLUDED

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

typedef std::complex<double> dcomp;

const double pi = 3.14159265359;
const double e = 2.71828182846;

//forward defenitions******************************************************************
class thimble_system;
//*************************************************************************************

class interaction
{
    private:
    double coupling;
    std::vector<int> powers;

    protected:

    public:
    dcomp base(int site, thimble_system* current_system);
    dcomp first_derivative(int site, int field, thimble_system* current_system);
    dcomp second_derivative(int site, int field_1, int field_2, thimble_system* current_system);

    //constructor
    interaction(double Coupling, std::vector<int> Powers);
};

class scalar_field
{
    private:
    int* occupation_number; //initial occupation number during field setup
    int Nx, Nt, Npath, Nrpath, Ntot;
    double m, squareMass; //field mass
    double dt, dx; //lattice spacings
    double* path; //sign of the path around the contour
    double* path_offset;
    bool is_flowed;
    dcomp* field_0;
    dcomp* field_1; //sites 0 and 1, which are integrated out of the full simulation
    int* positive_time_site;
    int* positive_space_site;
    int* negative_time_site;
    int* negative_space_site; //co-ordinate shifted arrays to minimise computation time for neighbour sites
    gsl_rng* my_rngPointer;
    dcomp j;

    //private assigment constructor
    scalar_field& operator = (const scalar_field &t)
    {
        return *this;
    }

    protected:
    dcomp free_action(int site);
    dcomp free_action_derivative(int site);
    dcomp free_action_second_derivative(int site_1, int site_2);
    int calc_n(int site);

    public:
    dcomp* base_field;
    dcomp* flowed_field;

    //bulk of the code is given in scalar_field.cpp
    void initialise();
    void set_occupation_number(int new_occupation_number[]);
    void set_occupation_number(int new_occupation_number);
    void set_mass(double new_mass);
    void set_dx(double new_dx) {dx = new_dx;};
    void set_dt(double new_dt);

    //interfaces
    double get_mass() {return m;};
    double get_square_mass() {return squareMass;};

    //constructor, destructor and copy constructor
    scalar_field(int x_dim, int t_dim, gsl_rng * rngPointer);
    ~scalar_field();
    scalar_field(const scalar_field &obj);
};

class thimble_system
{
    protected:
    int Nx, Nt, Npath, Nrpath, Ntot; //lattice setup parameters
    int number_of_timesteps; //number of iterations for the ode solver
    double tau; //flowtime
    double h; //ode step size
    unsigned long int rng_seed;
    gsl_rng * my_rngPointer; //rng pointer for the system/simulation
    

    public:
    std::vector<scalar_field> scalars; //the scalar fields of the simulation SHOULDN'T BE PUBLIC, ONLY IS FOR TESTING
    std::vector<interaction> interactions; //vector of all field interactions, similarly only public for testing
    
    void add_scalar_field();
    void add_scalar_field(double mass);
    void add_interaction(double coupling, std::vector<int> powers);
    void add_interaction(double coupling, int powers);
    
    //constructor and destructor
    thimble_system(int x_dim, int t_dim, double flow_time, long unsigned int seed);
    ~thimble_system();
    
    //friendship declaration
};

#endif //THIMBLE_LATTICE_H_INCLUDED
