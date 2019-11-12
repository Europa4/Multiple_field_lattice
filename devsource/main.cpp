
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
  uniform_int_distribution<int> dist(0, pow(2,16));

  std::vector<int> powers = {4};

  for(int i = 0; i < 200; ++i)
  {
    seed = dist(rd);
    printf("simulation %i initiated with seed %i \n", i*world_size + world_rank, seed + 7);
    thimble_system sys(1, 10, .4, seed);
    sys.add_scalar_field(1.0);
    sys.add_scalar_field(3.0);
    //sys.add_interaction(1./24, powers);
    sys.set_path("Data_matrix/");
    sys.set_name("phi_" + std::to_string(i*world_size + world_rank));
    sys.set_dt(0.5);
    sys.simulate(pow(10, 3), pow(10, 5));
    printf("simulation %i completed \n", i*world_size + world_rank + 7);
  }
  
  MPI_Finalize(); //closing the MPI enviroment
  return 0;
}



