
#include <iostream>
#include <random>
#include <stdio.h>
#include <unistd.h>
#include "mpi.h"

#include "thimble_lattice.hpp"
using namespace std;
#include <limits>
#include <chrono>
#include <ctime>

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
  double start_record = 0; //counterintuitively these should be doubles not ints. Integer division causes issues.
  double end_record = 1;
  int cycles;

  cycles = int(ceil((end_record - start_record)/world_size));

  /*
  auto start = std::chrono::system_clock::now();
  std::time_t start_t = std::chrono::system_clock::to_time_t(start);
  std::cout<<"started at " << std::ctime(&start_t) << endl;
  */
  for(int i = 0; i < cycles; ++i)
  {
    seed = dist(rd);
    printf("simulation %i initiated with seed %i \n", i*world_size + world_rank + start_record, seed);
    thimble_system sys(8, 10, 1.5, seed);
    sys.add_scalar_field(1.0);
    //sys.add_scalar_field(0.95);
    //sys.add_interaction(1./24, powers);
    sys.set_path("../Data_1_1_timing/");
    sys.set_name("phi_" + std::to_string(i*world_size + world_rank + start_record));
    sys.set_dt(0.5);
    sys.set_dx(0.75);
    sys.set_occupation_number(0, 1);
    sys.set_proposal_size(2);
    sys.simulate(0, 2);
    printf("simulation %i completed \n", i*world_size + world_rank + start_record);
  }

  /*
  auto end = std::chrono::system_clock::now();
  std::time_t end_t = std::chrono::system_clock::to_time_t(end);
  std::cout<<"ended at " << std::ctime(&end_t) << endl;
  */
  MPI_Finalize(); //closing the MPI enviroment
  return 0;
}



