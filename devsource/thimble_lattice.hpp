#ifndef THIMBLE_LATTICE_H_INCLUDED
#define THIMBLE_LATTICE_H_INCLUDED

#include <complex>
#include <math.h>
#include <quadmath.h>

#include <boost/filesystem.hpp>
#include <boost/multiprecision/float128.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#include <stdio.h>
#include <string>
#include <time.h>
#include <fstream>

#include "matrix.hpp"

typedef std::complex<double> dcomp;
//typedef boost::multiprecision::cpp_complex_quad dcomp;
//typedef boost::multiprecision::float128 mydouble;
typedef double mydouble;

const double pi = 3.14159265359;
const double e = 2.71828182846;

//forward defenitions******************************************************************
class thimble_system;
//*************************************************************************************

struct field_id_return
{
    int field_number;
    int site_number;
};

class interaction
{
    private:
    double coupling;
    std::vector<int> powers;

    protected:

    public:
    dcomp base(int site, thimble_system &current_system, int field_type = 0);
    dcomp first_derivative(int site, int field, thimble_system &current_system, int field_type = 0);
    dcomp second_derivative(int site, int field_1, int field_2, thimble_system &current_system, int field_type = 0);

    //constructor
    interaction(double Coupling, std::vector<int> Powers);

    friend class thimble_system;
};

class scalar_field
{
    private:
    int* occupation_number; //initial occupation number during field setup
    int Nx, Nt, Npath, Nrpath, Ntot;
    double m, squareMass; //field mass
    double dt, dx; //lattice spacings

    bool is_flowed;
    dcomp* field_0;
    dcomp* field_1; //sites 0 and 1, which are integrated out of the full simulation
    dcomp* field_2; //used for the edge case calculation, this is the initial state of field site 2 given the initial conditions
    int* positive_time_site;
    int* positive_space_site;
    int* negative_time_site;
    int* negative_space_site; //co-ordinate shifted arrays to minimise computation time for neighbour sites
    dcomp j;

    //private assigment constructor
    scalar_field& operator = (const scalar_field &t)
    {
        return *this;
    }

    protected:

    public:
    dcomp* fields[5]; //contains all the flowed and unflowed fields
    dcomp* C[5];
    double* path; //sign of the path around the contour
    double* path_offset;

    //bulk of the code is given in scalar_field.cpp
    void initialise(double a[], double b[], double c[], double d[]);
    void set_occupation_number(int new_occupation_number[]);
    void set_occupation_number(int new_occupation_number);
    void set_mass(double new_mass);
    void set_dx(double new_dx) {dx = new_dx;};
    void set_dt(double new_dt);    
    dcomp free_action(int site, int field_type = 0);
    dcomp free_action_derivative(int site, int field_type = 0);
    dcomp free_action_second_derivative(int site_1, int site_2);
    dcomp edge_effects(int site, int field_type = 0);
    dcomp edge_effects_derivative(int site);
    int calc_n(int site);
    int calc_x(int site);


    //interfaces
    double get_mass() {return m;};
    double get_square_mass() {return squareMass;};

    //constructor, destructor and copy constructor
    scalar_field(int x_dim, int t_dim, double system_dt, double system_dx);
    ~scalar_field();
    scalar_field(const scalar_field &obj);

    friend class thimble_system;
};

class thimble_system
{
    protected:
    int Nx, Nt, Npath, Nrpath, Ntot, Nsys; //lattice setup parameters
    int number_of_timesteps; //number of iterations for the ode solver
    int Njac;
    int NjacSquared;
    double tau; //flowtime
    double h; //ode step size
    double dx;
    double dt;
    double delta;
    double sigma; //proposal step size when put into the gaussian distribution
    double acceptance_rate;
    unsigned long int rng_seed;
    matrix <dcomp> J;
    matrix <dcomp> J_conj;
    dcomp S;
    std::string rel_path; //path to where the data should be stored
    std::string file_name; //name of the data file generated by this simulation
    dcomp j;
    std::mt19937_64 generator;
    std::uniform_real_distribution<double> uniform_double;
    std::uniform_int_distribution<int> uniform_int;
    std::normal_distribution<double> gaussian;
    std::normal_distribution<double> abcd;

    matrix<dcomp> calc_jacobian(bool proposal = false);
    dcomp calc_dS(int site, int field, int field_type);
    dcomp calc_dS(int site, int field_type);
    dcomp calc_ddS(int site_1, int site_2, int field_1, int field_2, int field_type = 0);
    dcomp calc_ddS(int site_1, int site_2, int field_type = 0);
    field_id_return calc_field(int master_site);
    void sync_ajustment(dcomp ajustment[]);
    int update();
    matrix<dcomp> sweep_proposal();
    matrix<dcomp> site_proposal();
    void pre_simulation_check();
    

    public:
    std::vector<scalar_field> scalars; //the scalar fields of the simulation SHOULDN'T BE PUBLIC, ONLY IS FOR TESTING
    std::vector<interaction> interactions; //vector of all field interactions, similarly only public for testing
    
    void add_scalar_field();
    void add_scalar_field(double mass);
    void add_interaction(double coupling, std::vector<int> powers);
    void add_interaction(double coupling, int powers);
    void set_path(std::string new_path);
    void set_name(std::string new_name) {file_name = new_name;};
    void simulate(int n_burn_in, int n_simulation);
    dcomp calc_S(int field_type = 0);
    void set_field_mass(int field_number, double new_mass);
    void set_dt(double new_dt);
    void set_dx(double new_dx);
    void set_occupation_number(int field_number, int new_occupation_number);
    void set_occupation_number(int field_number, int new_occupation_number[]);
    void set_occupation_number(int field_number, std::vector<int> new_occupation_number);
    void set_proposal_size(double new_delta);
    void test();
    
    //constructor and destructor
    thimble_system(int x_dim, int t_dim, double flow_time, long unsigned int seed);
    ~thimble_system();
    
    //friendship declaration
    friend class scalar_field;
    friend class interaction;
};

#endif //THIMBLE_LATTICE_H_INCLUDED
