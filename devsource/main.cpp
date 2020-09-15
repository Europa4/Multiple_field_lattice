
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

void harmonic_oscillator(const std::vector<double> &x, std::vector<double> &dx, const double t)
{
  dx[0] = x[1];
  dx[1] = -1.*x[0];
}


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




  // ode test
  boost::numeric::odeint::runge_kutta4<std::vector<double>> stepper;
  std::vector<double> x = {0., 1.};
  //boost::numeric::odeint::integrate_const(stepper, harmonic_oscillator, x, 0.0, 10.0, 0.01);
  boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::make_controlled<boost::numeric::odeint::runge_kutta_cash_karp54<std::vector<double>>>(1.e-6, 1.e-6), harmonic_oscillator, x, 0.0, 10.0, 0.01);
  //std::cout << x[0] << std::endl;
  //end of test

  std::vector<int> powers = {4};
  //std::vector<int> occ_number = {0, 0, 2, 2, 2, 2, 0, 0};
  std::vector<int> occ_number = {1};
  double start_record = 0; //counterintuitively these should be doubles not ints. Integer division causes issues.
  double end_record = 1;
  int cycles;

  cycles = int(ceil((end_record - start_record)/world_size));
  int start_record_int = int(start_record);
  int end_record_int = int(end_record);
  
  auto start = std::chrono::system_clock::now();
  std::time_t start_t = std::chrono::system_clock::to_time_t(start);
  std::cout<< "started at " << std::ctime(&start_t) << endl;
  
  /*
  seed = 5;
  thimble_system sys(1, 10, 1., seed);
  sys.add_scalar_field(1.0);
  //sys.add_scalar_field(1.0);
  sys.set_path("/run/media/ppxsw1/78fe3857-1897-4617-a65e-83c9aa61be27/self_interaction_test/");
  sys.set_name("phi_" + std::to_string(14099));
  sys.set_dt(0.75);
  sys.set_dx(0.75);
  sys.add_interaction(1./24., {4});
  sys.set_proposal_size(0, 0.25);
  sys.simulate(2 * pow(10, 3), pow(10, 5));
  */
  for(int i = 0; i < cycles; ++i)
  {
    seed = 5;//dist(rd);
    printf("simulation %i initiated with seed %i \n", i*int(world_size) + world_rank + int(start_record), seed);
    thimble_system sys(1, 10, 1.8, seed);
    sys.add_scalar_field(1.0);
    //sys.add_scalar_field(1.0);
    sys.set_path("/run/media/ppxsw1/78fe3857-1897-4617-a65e-83c9aa61be27/boost_free_18/");
    //sys.set_path("test/");
    sys.set_name("phi_" + std::to_string(i*world_size + world_rank + start_record_int));
    sys.set_dt(0.75);
    sys.set_dx(0.75);
    //sys.add_interaction(0.5/24., {4});
    //sys.add_interaction(1./24., {0, 4});
    //sys.add_interaction(1./24., {4, 0});
    //sys.set_occupation_number(0, occ_number);
    sys.set_proposal_size(0, 50.);
    //sys.set_proposal_size(1, 0.25);
    sys.simulate(2 * pow(10, 3), pow(10, 1));
    //sys.simulate(0, 50);
    //sys.simulate(2, 2);
    printf("simulation %i completed \n", i*world_size + world_rank + start_record_int);
  }

  /*
  {
    seed = 100;
    thimble_system sys(1, 4, 1.8, seed);
    sys.add_scalar_field(1.0);
    //sys.add_scalar_field(0.5);
    sys.set_path("/run/media/ppxsw1/78fe3857-1897-4617-a65e-83c9aa61be27/random_test/");
    sys.set_name("1_field");
    sys.set_dt(0.75);
    sys.set_dx(0.75);
    //sys.set_occupation_number(0, occ_number);
    sys.set_proposal_size(0, 0.25);
    //sys.set_proposal_size(1, 0.25);
    sys.simulate(0, 1);
    //sys.simulate(0, 100);
    //sys.simulate(2, 2);
    printf("simulation completed \n");
    printf("\n \n \n");
  }
  
  {
    seed = 100;
    thimble_system sys(1, 4, 1.8, seed);
    sys.add_scalar_field(1.0);
    sys.add_scalar_field(1.0);
    sys.set_path("/run/media/ppxsw1/78fe3857-1897-4617-a65e-83c9aa61be27/random_test/");
    sys.set_name("2_field");
    sys.set_dt(0.75);
    sys.set_dx(0.75);
    //sys.set_occupation_number(0, occ_number);
    sys.set_proposal_size(0, 0.25);
    sys.set_proposal_size(1, 0.25);
    sys.simulate(0, 1);
    //sys.simulate(0, 100);
    //sys.simulate(2, 2);
    printf("simulation completed \n");
  }
  */



  auto end = std::chrono::system_clock::now();
  std::time_t end_t = std::chrono::system_clock::to_time_t(end);
  std::cout << "ended at " << std::ctime(&end_t) << endl;
  
  MPI_Finalize(); //closing the MPI enviroment
  return 0;
}



