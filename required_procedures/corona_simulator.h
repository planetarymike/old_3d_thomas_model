//corona_simulator.h -- generate modeled spacecraft intensities for a given observation

#ifndef __OBS_SIM_H
#define __OBS_SIM_H

#include <iostream> // for file output and dialog
#include <fstream>  // ""
#include <cmath>    // for cos and sin
#include <cstdlib>
#include "physical.h" // physical definitions, nH, nCO2, KK, etc.
#include "quadrature.h" // to compute line integrals
#include "definitions.h" // basic parameter definitions
#include "nr3.h" // type VecDoub
#include "los.h"       //line of sight integration routines
#include "rad_trans.h" // two-point column density, holstein t and g functions
#include "iph_sim.h" 
#include "generate_S.h"

struct corona_simulator {

  bool obsinit;
  string obsname;
  string obsfname;
  int nobs;
  double aMars; // Sun-Mars distance at time of obs, AU
  double Fsun_earth; // Solar flux at line center at Earth
  double Fsun_mars; // Solar flux at line center at Mars, possibly calculated from above
  char datetime[1024];

  VecDoub I_obs, DI_obs;  // observed intensity (kR) and its uncertainty
  VecDoub alttan, SZAtan; // radius and SZA of tangent point
  VecDoub ra, dec, scvel; //ra and dec of los, and S/C SSB vel projected onto this direction
  MatDoub pos, dir; // position and LOS direction of spacecraft.
  IPHsim IPH; //IPH simulator

  string outfilename;
  VecDoub rpts, tpts, ppts;

  double lineintcoef;//intensity line integral coefficient ~"g-value"

  Sobj current_S; // current Sobj
  atmointerp current_atmointerp; //array of atmointerps for each nH, T
  bool Sinit; //also reflects status of atmointerp b/c these are loaded together

  LOS_integrator LOS; //object to do line integration given S and tabulated HolT values

  VecDoub I_calc; //storage for calculated radiances
  
  //constructor
  corona_simulator() : LOS(HolTfilename) {
    //nothing is initialized at first
    obsinit  = 0;
    Sinit    = 0;
  };
  
  void obs_import(string obsname) {
    //loads data corresponding to an observation
    fstream obsfile;
    obsfile.open(obsname.c_str());//this assumes that the code is executing 
    //from the directory containing this file, or the filename is absolute
    
    if (obsfile.is_open()) {
      //now read in the data:
      //read the preliminary data in
      obsfile >> nobs; // number of observations in this file
      /* std::cout << "nobs = " << nobs << std::endl; */
      obsfile >> aMars; // Sun-Mars distance at time of obs, AU
      /* std::cout << "aMars = " << aMars << std::endl; */
      obsfile >> Fsun_earth; // photons/s/m2/Angstrom Solar flux at line center at Earth
      Fsun_earth /= 1e4; // photons/s/cm2/Angstrom
      Fsun_earth *= (1215.)*(1215e-8)/(3e10/*cm/s*/);// photons/s/cm2/Hz
      Fsun_mars = Fsun_earth/aMars/aMars;
      /* std::cout << "Fsun_mars = " << Fsun_mars << std::endl; */
      obsfile.getline(datetime,1024); // capture time of file (irrelevant for this code)

      //prepare vector parameters for read-in
      I_obs.resize(nobs);
      DI_obs.resize(nobs);
      alttan.resize(nobs);
      SZAtan.resize(nobs);
      pos.resize(nobs,3);
      dir.resize(nobs,3);
      ra.resize(nobs);
      dec.resize(nobs);
      scvel.resize(nobs);
     
      /* X axis points toward Sun!! */
      //loop over the number of observations, reading each in to the appropriate variable
      for (int row = 0; row < nobs; row++) {
	//get data from file
	obsfile >> row;
	obsfile >> I_obs[row];
	obsfile >> DI_obs[row];
	obsfile >> pos[row][0];// km
	pos[row][0] *= 1e5;// cm
	obsfile >> pos[row][1];
	pos[row][1] *= 1e5;// cm
	obsfile >> pos[row][2];
	pos[row][2] *= 1e5;// cm
	obsfile >> dir[row][0];
	obsfile >> dir[row][1];
	obsfile >> dir[row][2];
	obsfile >> alttan[row];
	obsfile >> SZAtan[row];
	obsfile >> ra[row];
	obsfile >> dec[row];
	obsfile >> scvel[row];
			 
	// std::cout << "On row " << row << ":\n";
	// std::cout << "I_obs = " << I_obs[row] << " +- " << DI_obs[row] << " kR,\n";
	// std::cout << "pos = { " << pos[row][0] << " , "
	// 	 << pos[row][1] << " , " << pos[row][2] << " },\n";
	// std::cout << "r_pos = "
	// 	 << sqrt(pos[row][0]*pos[row][0]
	// 		 + pos[row][1]*pos[row][1]
	// 		 + pos[row][2]*pos[row][2]) << " cm\n";
	// std::cout << "dir = { " << dir[row][0]
	// 	 << " , " << dir[row][1]
	// 	 << " , " << dir[row][2] << " },\n";
	// std::cout << "alttan = " << alttan[row] << " , SZAtan = " << SZAtan[row] << ".\n";
	// std::cout << "ra = " << ra[row] << " , dec = " << dec[row] << ".\n";
	// std::cout << "scvel = " << scvel[row] << ".\n";
	// std::cin.get();
      }
      obsfile.close();

      //set up IPH simulator
      IPH = IPHsim(nobs, ra, dec, scvel, pos, dir);

      obsinit = 1;//file successfully loaded.
      
    } else {//obsfile not found!
      std::cout << "Observation file " << obsfile << " not found!\n" 
		<< "Check that file path is correct.\n";
      
    }
  }

  void get_S(double nH, double T,bool forcesim=FALSE) {
    /* tries to load S from file; if file not found, calculates and
       produces file for future use */
    
    if (!forcesim) { 
	//what would the filename look like?
	string Sfname=Sfilename(nH,T);
	ifstream Sfile(Sfname.c_str());
	if (Sfile.good()) {
	  Sfile.close();
	  current_S=Sobj(Sfname);
	  //if S file exists, atmointerp does too.
	  current_atmointerp = atmointerp(current_S.nH,current_S.T,nphyspts);
	}
	Sfile.close();
    } else {
      //Sfile does not exist, we must create it.
      generate_S(nH, T, current_S, current_atmointerp);
    }
    Sinit=1;
  }

  double simulate_iobs(int iobs, Sobj &thisS, atmointerp &thisatmointerp, double IPHb=0.0) {
    if (!Sinit||!obsinit)
      throw("Load observation and source function before calling simulate()!");

    //simulate intensity for a single observation using the current
    //Sobj and the specicified IPH background.
    /* std::cout << "Simulating coordinate iobs = " << iobs  */
    /* 	      << ", S_nH =" << current_S.nH << ", S_T = " << current_S.T << ".\n"; */

    double nH = current_S.nH;
    double T = current_S.T;

    //from the input data, compute lineint = Fsun * doppler width (Hz) * sqrt(pi)
    lineintcoef = Fsun_mars*sHtot/(4*pi);
    //^^^^ photons/s/cm2/Hz * Hz = photons/cm2/s.
    /* std::cout << "doppler width* sqrt(pi) = " << sqrt(2*pi*kB*T/mH) << std::endl; */
    /* std::cout << "lineintcoef = " << lineintcoef << std::endl; */

    
    double tIcalc = 0.0;
    VecDoub posvec(3,pos[iobs]);
    VecDoub dirvec(3,dir[iobs]);
    tIcalc = LOS.integrate(thisS,
			   thisatmointerp,
			   posvec,
			   dirvec,
			   lineintcoef);// photons/cm2/s
    //convert to rayleighs
    tIcalc /= 1e6;// megaphoton/cm2/s
    tIcalc *= 4*pi/1e3; // now we're in kR; see C&H pg. 280-282
    
    //now add in the IPH
    //    tIcalc += IPH.sim(IPHb,nH,T,iobs);
    
    return tIcalc;
      
    //	  std::cout << "tIcalc = " << tIcalc << "\n";
    //	  std::cout << "I_obs = " << thisobs.I_obs[iobs]
    //		    << " +- " << thisobs.DI_obs[iobs] << " kR,\n";
    //	  std::cin.get();
  }

  void calc_I(double IPHb=0.0) {
    //simulates for all iobs and stores the result in I_calc
    if (!Sinit||!obsinit)
      throw("Load observation and source function before calling simulate()!");
    I_calc.resize(nobs);
    for (int iobs=0; iobs<nobs; iobs++)
      I_calc[iobs]=simulate_iobs(iobs,current_S,current_atmointerp,IPHb);
  }
  
};

#endif
