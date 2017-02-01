//simulate_coronal_scan.cpp --- program to simulate an observation
//provided, given a temperature and a density.

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
//./simulate_coronal_scan.x -b 0 -T 300 -n 100000 ./mvn_iuv_l1b_outbound-orbit00422-fuv_20141217T160952_v01_r01_obs.dat


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
  bool multidens=FALSE;
  double nexo=1.0;
  double nexo_i=1.0;
  double nexo_f=1.0;
  double nexo_d=1.0;
  bool simulate_IPH=FALSE;
  bool forcesim=FALSE;
  bool nointerp=FALSE;
  bool silent=FALSE;

  //check to see if silent output is requested
  for (int i = 1; i < argc; i++) 
    if (strcmp(argv[i], "-silent") == 0)
      silent=TRUE;
  
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
	}
      else if (strcmp(argv[i], "-c") == 0)
        {
	  usercal=TRUE;
	  i++;
	  cal = atof(argv[i]);
	  if (!silent) {
	    std::cout << "Fixing calibration factor at " << cal << " .\n";
	    std::cout << std::endl;
	  }
        }
      else if (strcmp(argv[i], "-T") == 0)
        {
	  usertemp=TRUE;
	  i++;
	  Texo = atof(argv[i]);
	  if (!silent) {
	    std::cout << "Fixing temperature at " << Texo << " K.\n";
	    std::cout << std::endl;
	  }
        }
      else if (strcmp(argv[i], "-n") == 0)
        {
	  userdens=TRUE;
	  i++;
	  if (strcmp(argv[i],"[") == 0) {
	    //the number density parameter is an array argument
	    multidens=TRUE;
	    i++;
	    nexo_i=atof(argv[i]);
	    i++;
	    nexo_f=atof(argv[i]);
	    i++;
	    nexo_d=atof(argv[i]);
	    i++;//clears the closing brace
	    if (!silent) {
	      std::cout << "Simulating multiple densities.\n"
			<< "From " << nexo_i << "/cc to " << nexo_f << " /cc \n"
			<< "in steps of " << nexo_d << "/cc .\n";
	      std::cout << std::endl;
	    }
	  } else {
	    nexo = atof(argv[i]);
	    nexo_i=nexo;
	    nexo_f=nexo;
	    nexo_d=nexo;
	    if (!silent) {
	      std::cout << "Fixing density at " << nexo << " / cm3.\n";
	      std::cout << std::endl;
	    }
	  }
        }
      else if (strcmp(argv[i], "-simIPH") == 0)
        {
	  simulate_IPH=TRUE;
        }
      else if (strcmp(argv[i], "-forcesim") == 0)
        {
	  forcesim=TRUE;
        }
      else if (strcmp(argv[i], "-nointerp") == 0)
        {
	  nointerp=TRUE;
        }
      else if (strcmp(argv[i], "-silent") == 0)
	silent=TRUE;
      else
	{
	  //get the observation to perform analysis on from the function call:
	  obsname=argv[i];
	  if (!silent) 
	    std::cout << "Observation under analysis is: " << obsname << std::endl;	  
	}
    }

  if (userIPH) {
    if (simulate_IPH) {
      if (!silent) {
	std::cout << "Fixing background IPH multiplier at " << IPHb << " .\n";
	std::cout << std::endl;
      }
    } else {
      if (!silent) {
	std::cout << "Fixing background IPH at " << IPHb << " kR.\n";
	std::cout << std::endl;
      }
    }
  } else {
    if (simulate_IPH) {
      IPHb=1.0;
    } else {
      IPHb=0.0;
    }
  }
  if (!usercal) {
    cal=1.0;
  }

  corona_simulator sim(forcesim,simulate_IPH,silent);
  if (nointerp) {
    //direct simulation for these conditions
    if (!silent)
      std::cout << "Direct simulation proceeding.\n";
    sim.load_obs(obsname);//load the observation
    for (nexo=nexo_i;nexo<=nexo_f;nexo+=nexo_d) {
      sim.get_S(nexo,Texo);//load or simulate the source function
      if (!silent) {
	std::cout << "current_Sinit = " << sim.current_Sinit << std::endl;
      }
      sim.calc_I(IPHb);//calculate with background
      
      //print out the results:
      if (!silent) 
	std::cout << "\nHere's the result of the calculation:\n\n";
      if (multidens)
	std::cout << "nexo = " << nexo << "; ";
      std::cout	<< "I_calc = [ ";
      int i;
      for (i = 0; i < sim.allobsdata[0].nobs-1; i++)
	std::cout << cal*sim.allobsdata[0].current_I_calc[i] << ", ";
      std::cout << cal*sim.allobsdata[0].current_I_calc[i] << " ]\n";
    }
  } else {
    //interpolated simulation
    if (!silent) 
      std::cout << "Proceeding with interpolated computation.\n";
    sim.load_obs(obsname);//load the observation
    VecDoub I_calc;
    for (nexo=nexo_i;nexo<=nexo_f;nexo+=nexo_d) {
      sim.interp_I(0, nexo, Texo, IPHb);//calculate with background
      
      //print out the results:
      if (!silent) 
	std::cout << "\nHere's the result of the calculation:\n\n";
      if (multidens)
	std::cout << "nexo = " << nexo << "; ";
      std::cout	<< "I_calc = [ ";
      int i;
      for (i = 0; i < sim.allobsdata[0].nobs-1; i++)
	std::cout << cal*sim.allobsdata[0].current_I_calc[i] << ", ";
      std::cout << cal*sim.allobsdata[0].current_I_calc[i] << " ]\n";
    }
  }

  if (!silent)
    std::cout << "\n\nYou have reached the end of the program!\nGoodbye.\n\n";

  return 0; 
}

