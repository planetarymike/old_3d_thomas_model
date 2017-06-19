//simulate_coronal_scan.cpp --- program to simulate an observation
//provided, given a temperature and a density.

#include <iostream>
#include <fstream>
#include "nr3.h"
#include "multi_corona_simulator.h" //simulation routines for synthetic atmospheres
#include "eff_interp.h" //table to exchange between temperature and effusion velocity

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
  bool multitemp=FALSE;
  double T_i=1.0;
  double T_f=1.0;
  double T_d=1.0;
  double Texolist[1024];
  int nT;
  int iT;

  bool useresc=FALSE;
  double esc=1.0;
  bool multiesc=FALSE;
  double esc_i=1.0;
  double esc_f=1.0;
  double esc_d=1.0;
  double esclist[1024];

  bool userdens=FALSE;
  bool multidens=FALSE;
  double nexo=1.0;
  double nexo_i=1.0;
  double nexo_f=1.0;
  double nexo_d=1.0;
  double nexolist[1024];
  int nnexo;
  int inexo;

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
	  if (strcmp(argv[i],"[") == 0) {
	    //the temperature parameter is an array argument
	    multitemp=TRUE;
	    i++;
	    T_i=atof(argv[i]);
	    i++;
	    T_f=atof(argv[i]);
	    i++;
	    T_d=atof(argv[i]);
	    i++;//clears the closing brace
	    iT=0;
	    nT=0;
	    for (Texo=T_i;Texo<=T_f;Texo+=T_d) {
	      Texolist[iT]=Texo;
	      nT++;iT++;
	    }
	    std::cout << "nT = " << nT << std::endl;
	    if (!silent) {
	      std::cout << "Simulating multiple temperatures.\n"
			<< "From " << T_i << "K to " << T_f << "K \n"
			<< "in steps of " << T_d << "K .\n";
	      std::cout << std::endl;
	    }
	  } else if (strcmp(argv[i],"{") == 0) {
	    //the temperature parameter is an array argument
	    multitemp=TRUE;
	    i++;
	    nT=0;
	    while (strcmp(argv[i],"}")) {
	      Texolist[nT]=atof(argv[i]);
	      i++;nT++;
	    }
	    //	    i++;//clears the closing brace
	    if (!silent) {
	      std::cout << "Simulating multiple temperatures.\n";
	      std::cout << "{ ";
	      for (iT=0;iT<nT-1;iT++) {
		std::cout << Texolist[iT] << ", ";
	      }
	      std::cout << Texolist[iT] << " }.";
	      std::cout << std::endl;
	    }
	  } else {
	    Texo = atof(argv[i]);
	    nT=1;
	    Texolist[0]=Texo;
	    if (!silent) {
	      std::cout << "Fixing temperature at " << Texo << "K .\n";
	      std::cout << std::endl;
	    }
	  }
        }
      else if (strcmp(argv[i], "-e") == 0)
        {
	  useresc=TRUE;
	  i++;
	  if (strcmp(argv[i],"[") == 0) {
	    //the escape flux parameter is an array argument,
	    //producing multiple temperatures
	    multitemp=TRUE;
	    i++;
	    esc_i=atof(argv[i]);
	    i++;
	    esc_f=atof(argv[i]);
	    i++;
	    esc_d=atof(argv[i]);
	    i++;//clears the closing brace
	    iT=0;
	    nT=0;
	    for (esc=esc_i;esc<=esc_f;esc+=esc_d) {
	      esclist[iT]=esc;
	      nT++;iT++;
	    }
	    std::cout << "nesc = " << nT << std::endl;
	    if (!silent) {
	      std::cout << "Simulating multiple escape fluxes (via temperature).\n"
			<< "From " << esc_i << "/cm2/s to " << esc_f << "/cm2/s \n"
			<< "in steps of " << esc_d << "/cm2/s .\n";
	      std::cout << std::endl;
	    }
	  } else if (strcmp(argv[i],"{") == 0) {
	    //the escape flux parameter is an array argument,
	    //producing multiple temperatures
	    multitemp=TRUE;
	    i++;
	    nT=0;
	    while (strcmp(argv[i],"}")) {
	      esclist[nT]=atof(argv[i]);
	      i++;nT++;
	    }
	    //	    i++;//clears the closing brace
	    if (!silent) {
	      std::cout << "Simulating multiple escape fluxes (via temperature).\n";
	      std::cout << "{ ";
	      for (iT=0;iT<nT-1;iT++) {
		std::cout << esclist[iT] << ", ";
	      }
	      std::cout << esclist[iT] << " }.";
	      std::cout << std::endl;
	    }
	  } else {
	    esc = atof(argv[i]);
	    nT=1;
	    esclist[0]=esc;
	    if (!silent) {
	      std::cout << "Fixing escape flux at " << esc << "/cm2/s .\n";
	      std::cout << std::endl;
	    }
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
	    inexo=0;
	    nnexo=0;
	    for (nexo=nexo_i;nexo<=nexo_f;nexo+=nexo_d) {
	      nexolist[inexo]=nexo;
	      nnexo++;inexo++;
	    }
	    std::cout << "nnexo = " << nnexo << std::endl;
	    if (!silent) {
	      std::cout << "Simulating multiple densities.\n"
			<< "From " << nexo_i << "/cc to " << nexo_f << " /cc \n"
			<< "in steps of " << nexo_d << "/cc .\n";
	      std::cout << std::endl;
	    }
	  } else if (strcmp(argv[i],"{") == 0) {
	    //the number density parameter is an array argument
	    multidens=TRUE;
	    i++;
	    nnexo=0;
	    while (strcmp(argv[i],"}")) {
	      nexolist[nnexo]=atof(argv[i]);
	      i++;nnexo++;
	    }
	    //	    i++;//clears the closing brace
	    if (!silent) {
	      std::cout << "Simulating multiple densities.\n";
	      std::cout << "{ ";
	      for (inexo=0;inexo<nnexo-1;inexo++) {
		std::cout << nexolist[inexo] << ", ";
	      }
	      std::cout << nexolist[inexo] << " }.";
	      std::cout << std::endl;
	    }
	  } else {
	    nexo = atof(argv[i]);
	    nnexo=1;
	    nexolist[0]=nexo;
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

  //make sure the user has not specified both escape flux and temperature
  if (usertemp&&useresc)
    throw("Please specify either temperature or escape flux, but not both.")
  
  corona_simulator sim(forcesim,simulate_IPH,silent);
  eff_interp effinterp(eff_filename,silent);
  if (nointerp) {
    //direct simulation for these conditions
    if (!silent)
      std::cout << "Direct simulation proceeding.\n";
    sim.load_obs(obsname);//load the observation
    for (iT=0;iT<nT;iT++) {
      for (inexo=0;inexo<nnexo;inexo++) {
	nexo=nexolist[inexo];
	if (!useresc)
	  Texo=Texolist[iT];
	if (useresc)
	  Texo=effinterp.interp(esclist[iT]/nexo);
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
	if (multitemp) {
	  if (!useresc)
	    std::cout << "Texo = " << Texo << "; ";
	  if (useresc)
	    std::cout << "esc = " << esclist[iT] << "; Texo = " << Texo << "; ";	      
	}
	std::cout	<< "I_calc = [ ";
	int i;
	for (i = 0; i < sim.allobsdata[0].nobs-1; i++)
	  std::cout << cal*sim.allobsdata[0].current_I_calc[i] << ", ";
	std::cout << cal*sim.allobsdata[0].current_I_calc[i] << " ]\n";
      }
    }
  } else {
    //interpolated simulation
    if (!silent) 
      std::cout << "Proceeding with interpolated computation.\n";
    sim.load_obs(obsname);//load the observation
    VecDoub I_calc;
    for (iT=0;iT<nT;iT++) {
      for (inexo=0;inexo<nnexo;inexo++) {
	nexo=nexolist[inexo];
	if (!useresc)
	  Texo=Texolist[iT];
	if (useresc)
	  Texo=effinterp.interp(esclist[iT]/nexo);
	// std::cout << "inexo = "<<inexo<<std::endl;
	// std::cout << "nexolist[inexo] = "<<nexolist[inexo]<<std::endl;
	sim.interp_I(0, nexo, Texo, IPHb);//calculate with background
	
	//print out the results:
	if (!silent) 
	  std::cout << "\nHere's the result of the calculation:\n\n";
	if (multidens)
	  std::cout << "nexo = " << nexo << "; ";
	if (multitemp) {
	  if (!useresc)
	    std::cout << "Texo = " << Texo << "; ";
	  if (useresc)
	    std::cout << "esc = " << esclist[iT] << "; Texo = " << Texo << "; ";	      
	}
	std::cout	<< "I_calc = [ ";
	int i;
	for (i = 0; i < sim.allobsdata[0].nobs-1; i++)
	  std::cout << cal*sim.allobsdata[0].current_I_calc[i] << ", ";
	std::cout << cal*sim.allobsdata[0].current_I_calc[i] << " ]\n";
      }
    }
  }
  if (!silent)
    std::cout << "\n\nYou have reached the end of the program!\nGoodbye.\n\n";

  return 0; 
}

