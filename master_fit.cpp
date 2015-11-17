//master_fit.cpp -- fitting routines

#include <iostream>
#include <fstream>
#include "nr3.h"
#include "fitmrq_interp.h" //Levenberg-Marquardt chi2 minimization
#include "corona_simulator.h" //simulation routines for synthetic atmospheres


//compile: \/
// SUPERSEDED BY MAKEFILE
// mpicxx master_fit.cpp `pkg-config --cflags --libs gsl` -Irequired_procedures/ 'SRCFNSLOC="./Model_Source_Functions/"' -o master_fit.x
// add -pg for profiling.
// add -g for line numbers


int main(int argc, char* argv[]) {
  // a function to do L-M minimization and return H coronal parameters
  // from spacecraft data

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
  bool printsim=FALSE;
  bool simulate_IPH=FALSE;
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
	  simulate_IPH=TRUE;
        }
      else if (strcmp(argv[i], "-forcesim") == 0)
        {
	  forcesim=TRUE;
        }
      else if (strcmp(argv[i], "-p") == 0)
        {
	  printsim=TRUE;
	}
      else
	{
	  //get the observation to perform analysis on from the function call:
	  obsname=argv[i];
	  std::cout << "Observation under analysis is: " << obsname << std::endl;	  
        }
    }

  if (userIPH) {
    if (simulate_IPH) {
      std::cout << "Fixing background IPH multiplier at " << IPHb << " .\n";
      std::cout << std::endl;
    } else {
      std::cout << "Fixing background IPH at " << IPHb << " kR.\n";
      std::cout << std::endl;
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


  corona_simulator sim;
  std::cout << "Proceeding with interpolated computation.\n";
  sim.obs_import(obsname,simulate_IPH);//load the observation

  //here's an initial parameter guess:
  VecDoub parms(4);///mmm delicious parms
  parms[0]=30000.0;parms[1]=350.0;parms[2]=0.45;parms[3]=1.0;

  //this sets it up:
  Fitmrq testfit(sim.obs_vec,sim.I_obs,sim.DI_obs,parms,sim);

  //if the user specified values, fix this to the user value:
  if (userdens) {
    testfit.hold(0,nexo);
  }
  if (usertemp) {
    testfit.hold(1,Texo);
  }
  if (userIPH) {
    testfit.hold(2,IPHb);
  }
  if (usercal) {
    testfit.hold(3,cal);
  }
  
  //now do the fit, which updates all the contained objects:
  testfit.fit();

  //get the escape flux for these parameters
  double lambdac = G*mMars*mH/(kB*testfit.a[1]*rexo);
  double vth = std::sqrt(2*kB*parms[1]/mH);
  double escflux = testfit.a[0]*vth/(2*std::sqrt(pi))
                           *(1+lambdac)*std::exp(-lambdac);

  
  //print out the results:
  std::cout << "\nHere's the result of the fitting routine:\n"
	    << "Chi-squared is: " << testfit.chisq
	    << " on " << sim.nobs-parms.size() << " DOF."<< std::endl;

  std::cout << "\nBest fit model parameters are: \n";
  for (int i = 0; i < testfit.a.size(); i++)
    std::cout << testfit.a[i] << ",\t";
  std::cout << "\b\b\n";

  std::cout << "\nComputed escape flux is: \n";
  std::cout << escflux << std::endl;
    
  std::cout << "\nFormal fit parameter covariance matrix is as follows: \n";
  for (int i = 0; i < testfit.a.size(); i++) {
    std::cout << " [ ";
    for (int j = 0; j < testfit.a.size(); j++) {
      std::cout.width(10);
      std::cout << testfit.covar[i][j] << ",\t";
    }
    std::cout << "\b\b ]" << std::endl;
  }

  if (printsim) { //print out the simulated values at the best fit
    double I_calc[sim.nobs];
    
    for (int i = 0; i < sim.nobs; i++) {
      I_calc[i] = testfit.a[3]*sim.interp_iobs(sim.obs_vec[i],
					       testfit.a[0],
					       testfit.a[1],
					       testfit.a[2]);
    }

    //print out the results:
    std::cout << "\nHere's the best fit:\n"
	      << "\nI_calc = [ ";
    for (int i = 0; i < sim.nobs-1; i++)
      std::cout << I_calc[i] << ", ";
    std::cout << I_calc[sim.nobs-1] << " ]\n";

    if (simulate_IPH) {
      std::cout << "\nHere's the computed and scaled IPH vector:\n"
		<< "\nIPH = [ ";
      for (int i = 0; i < sim.nobs-1; i++)
	std::cout << testfit.a[3]*testfit.a[2]*sim.IPHb_model[i] << ", ";
      std::cout << testfit.a[3]*testfit.a[2]*sim.IPHb_model[sim.nobs-1] << " ]\n";
    }
    
  }

  std::cout << std::endl << "\n\nYou have reached the end of the program!\nGoodbye.\n\n";

  return 0; 
}

