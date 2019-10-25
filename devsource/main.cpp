#include <iostream>
#include <random>
#include <stdio.h>
#include <unistd.h>
#include "mpi.h"

#include "Prot.h"
using namespace std;
#include <limits>
/*
class tester_class
{
  public:
    tester_class(int size);
    ~tester_class();
    
    int array_size;
    int* test_array;
};

tester_class::tester_class(int size): array_size(size), test_array(new int[array_size])
{
  
}

tester_class::~tester_class()
{
  delete[] test_array;
}

*/
int main(int argc, char **argv)
{
  clock_t t1, t2;
  MPI_Init(NULL, NULL);
  
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size); //setting up the MPI enviroment
  
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);
  
  random_device rd;
  unsigned long int seed;
  uniform_int_distribution<int> dist(0,pow(2,16));
  initalise_dt();
  
  /*
  seed = 5;
  c_phi phi(1.0, 0.5, 0.0, 0.1, seed);
  phi.simulate(pow(10,3), pow(10,5), 1, 50);
  */

 /*
  for(int i = 0; i < 50; ++i)
  {
    seed = dist(rd);
    thimble_system sys(1, 10, 1.5, seed);
    sys.add_scalar_field(1.0);
    sys.set_path("Data/");
    sys.set_name("phi_" + std::to_string(i));
    sys.simulate(pow(10, 3), pow(10, 5));
  }
  */

  matrix<double> A(2,2), B(2,2); //C(2,2);
  
  A.set_element(0,0, 0.);
  A.set_element(0,1, 1.);
  A.set_element(1,0, 2.);
  A.set_element(1,1, 3.);
  
  B.set_element(0,0, 4.);
  B.set_element(0,1, 5.);
  B.set_element(1,0, 6.);
  B.set_element(1,1, 7.);
  
  double b[2], x[2];
  b[0] = 2.;
  b[1] = 3.;

  A.solve(x, b);
  printf("x[0] = %f \n", x[0]);
  printf("x[1] = %f \n", x[1]);
  
  MPI_Finalize(); //closing the MPI enviroment
  return 0;
}
