#ifndef LATTICE_SCALAR_H_INCLUDED
#define LATTICE_SCALAR_H_H_INCLUDED
class scalar_field
{
    private:
    int* occupation_number; //initial occupation number during field setup
    int Nx, Nt, Npath, Nrpath, Ntot;
    double m, squareMass; //field mass
    double dt, dx; //lattice spacings
    bool is_flowed;
    double* field_0;
    double* field_1; //sites 0 and 1, which are integrated out of the full simulation
    int* positive_time_site;
    int* positive_space_site;
    int* negative_time_site;
    int* negative_space_site; //co-ordinate shifted arrays to minimise computation time for neighbour sites

    protected:

    public:
    
    double* base_field;
    dcomp* flowed_field;

    //bulk of the code is given in scalar_field.cpp
    void initialise();
    void set_occupation_number(int new_occupation_number);
    void set_mass(double new_mass);

    //interfaces
    double get_mass() {return m;};

    //constructor and destructor 
    scalar_field(int x_dim, int t_dim);
    ~scalar_field();
};

#endif //THIMBLE_LATTICE_H_INCLUDED