
#include <iostream>
#include <random>
#include <stdio.h>
#include <unistd.h>
#include "mpi.h"

#include "thimble_lattice.hpp"
using namespace std;
#include <limits>

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


  std::vector<int> powers = {4};
  /*
  for(int i = 0; i < 50; ++i)
  {
    seed = dist(rd);
    thimble_system sys(1, 10, 1.0, seed);
    sys.add_scalar_field(1.0);
    sys.add_interaction(1., powers);
    sys.set_path("Data_matrix/");
    sys.set_name("phi_" + std::to_string(i*world_size + world_rank));
    sys.simulate(pow(10, 3), pow(10, 5));
    printf("simulation %i completed \n", i*world_size + world_rank);
  }
  */

 while(true)
 {
    seed = dist(rd);
    printf("seed = %i \n", seed);
    thimble_system sys(1, 10, 1.0, seed);
    sys.add_scalar_field(1.0);
    //sys.add_interaction(1., powers);
    sys.set_path("Data_matrix/");
    sys.set_name("phi_" + std::to_string(seed));
    sys.simulate(pow(10, 3), pow(10, 5));
 }
 

  MPI_Finalize(); //closing the MPI enviroment
  return 0;
}



