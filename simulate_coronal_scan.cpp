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
  //a function to simulate the intensities observed for an input
  //geometry and specific coronal parameters.

  //set up the input variables: name of file to fit against, and optional fitting parameters
  string obsname;
  bool userIPH=FALSE;
  double IPHb=0.0;
  bool usercal=FALSE;
  double cal=1.0;
  bool usertemp=FALSE;
  double Texo=1.0;
  bool userdens=FALSE;
  double nexo=1.0;
  bool simulate_IPH=FALSE;//IPH model not yet integrated with this code!
  bool forcesim=FALSE;
  
  for (int i = 1; i < argc; i++)  /* Skip argv[0] (program name). */
    {
      if (strcmp(argv[i], "-b") == 0)  /* Process optional arguments. */
        {
	  userIPH=TRUE;
	  /*
	   * Increment 'i' again (twice total) so that you don't
	   * check these arguments the next time through the loop.
	   */
	  i++;
	  IPHb = atof(argv[i]);  /* Convert string to int. */
	  std::cout << "Fixing background IPH at " << IPHb << " kR.\n";
	  std::cout << std::endl;
	}
      else if (strcmp(argv[i], "-c") == 0)
        {
	  usercal=TRUE;
	  i++;
	  cal = atof(argv[i]);
	  std::cout << "Fixing calibration factor at " << cal << " .\n";
	  std::cout << std::endl;
        }
      else if (strcmp(argv[i], "-T") == 0)
        {
	  usertemp=TRUE;
	  i++;
	  Texo = atof(argv[i]);
	  std::cout << "Fixing temperature at " << Texo << " K.\n";
	  std::cout << std::endl;
        }
      else if (strcmp(argv[i], "-n") == 0)
        {
	  userdens=TRUE;
	  i++;
	  nexo = atof(argv[i]);
	  std::cout << "Fixing density at " << nexo << " / cm3.\n";
	  std::cout << std::endl;
        }
      else if (strcmp(argv[i], "-simIPH") == 0)
        {
	  std::cout << "Warning! Quemerais IPH background model is not yet" 
		    << " integrated with this coronal simulator!";
	  simulate_IPH=TRUE;
        }
      else if (strcmp(argv[i], "-forcesim") == 0)
        {
	  forcesim=TRUE;
        }
      else
	{
	  //get the observation to perform analysis on from the function call:
	  obsname=argv[i];
	  std::cout << "Observation under analysis is: " << obsname << std::endl;	  
        }
    }

  if (!userIPH) {
    IPHb=0.0;
  }
  if (!usercal) {
    cal=1.0;
  }

  corona_simulator sim;
  sim.obs_import(obsname);//load the observation
  sim.get_S(nexo,Texo,forcesim);//load or simulate the source function
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

