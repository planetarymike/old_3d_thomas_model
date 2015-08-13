//simulate_coronal_scan.cpp --- program to simulate an observation
//provided, given a temperature and a density.

#include <iostream>
#include <fstream>
#include "nr3.h"
#include "corona_simulator.h" //simulation routines for synthetic atmospheres


//compile: 
//  \/
//mpicxx simulate_coronal_scan.cpp `pkg-config --cflags --libs gsl` -Irequired_procedures/ -Wl,-rpath,/usr/lib64/openmpi/lib/ -o simulate_coronal_scan.x
// add -pg for profiling.
// add -g for line numbers


int main(int argc, char* argv[]) {
  // a function to do L-M minimization and return H coronal parameters from spacecraft data 

  //get the observation to perform analysis on from the function call:
  string obsname;
  obsname=argv[1];
  std::cout << "Observation under analysis is: " << obsname << std::endl;
  //get the requested density and temperature from the command line.
  double nexo, Texo;
  nexo=atof(argv[2]);
  Texo=atof(argv[3]);
  double IPHb;
  if (argc > 4) {
    IPHb=atof(argv[4]);
  } else {
    IPHb=0.0;
  }

  corona_simulator sim;
  sim.obs_import(obsname);//load the observation
  sim.get_S(nexo,Texo);//load or simulate the source function
  sim.calc_I(IPHb);//calculate with background
  
  //print out the results:
  std::cout << "\nHere's the result of the calculation:\n"
	    << "\nI_calc = [ ";
  for (int i = 0; i < sim.nobs; i++)
    std::cout << sim.I_calc[i] << ", ";
  std::cout << "\b\b ]\n";
  std::cout << "\n\nYou have reached the end of the program!\nGoodbye.\n\n";

  return 0; 
}

