//source_function_init.cpp -- code to initialize the source function grid

#include <iostream>
#include <fstream>
#include "nr3.h"
#include "multi_corona_simulator.h" //simulation routines for synthetic atmospheres

//compile: 
//  \/
//mpicxx simulate_coronal_scan.cpp `pkg-config --cflags --libs gsl` -Irequired_procedures/ -Wl,-rpath,/usr/lib64/openmpi/lib/ -o simulate_coronal_scan.x
// add -pg for profiling.
// add -g for line numbers

//call
//./source function init.x 1 3

int main(int argc, char* argv[]) {

  bool forcesim=FALSE;
  bool simulate_IPH=FALSE;
  bool silent=FALSE;
  corona_simulator sim(forcesim,simulate_IPH,silent);
  sim.source_function_init(atof(argv[1]),atof(argv[2]));

  return 0; 
}

