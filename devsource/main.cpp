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

 
  for(int i = 4; i < 50; ++i)
  {
    seed = dist(rd);
    thimble_system sys(1, 10, 1.5, seed);
    sys.add_scalar_field(1.0);
    sys.set_path("Data/");
    sys.set_name("phi_" + std::to_string(i*world_size + world_rank));
    sys.simulate(pow(10, 3), pow(10, 5));
    printf("simulation %i completed \n", i*world_size + world_rank);
  }
  
 
  
  MPI_Finalize(); //closing the MPI enviroment
  return 0;
}
