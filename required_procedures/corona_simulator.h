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
  VecDoub IPHb; //storage for calculated IPH brightnesses
  
  VecDoub rpts, tpts, ppts;

  double lineintcoef;//intensity line integral coefficient ~"g-value"

  //for a single number density and temperature, a single S, atmointerp, Icalc
  Sobj current_S; // current Sobj
  atmointerp current_atmointerp; //current atmointerp
  bool current_Sinit; //also reflects status of atmointerp b/c these are loaded together
  VecDoub current_I_calc; //storage for calculated radiances
  VecDoub current_IPHtrans; //storage for calculated IPH transmission factors

  //for interpolation on a grid of S in number density and temperature:
  //  int nnH, nT;//number of grid points -- defined in definitions.h
  VecDoub nH_vec, T_vec;//location of grid points
  Linear_interp nH_terp, T_terp;//interpolation objects

  Sobj** S_grid; // grid of S at each nH, T
  atmointerp** atmointerp_grid; //grid of atmospheres at each nH, T
  bool** Sinit_grid; //grid reflecting initialization state of each atmosphere
  VecDoub** I_calc_grid; //storage for calculated radiances
  VecDoub** IPHtrans_grid; //storage for calculated radiances
  
  LOS_integrator LOS; //object to do line integration given S and tabulated HolT values
  
  //constructor
  corona_simulator() : LOS(HolTfilename) {
    //nothing is initialized at first
    obsinit  = 0;
    current_Sinit = 0;
    
    nH_vec.resize(nnH);
    T_vec.resize(nT);
    for (int inH=0;inH<nnH;inH++) {
      if (nHgridlog) {
	nH_vec[inH]=exp(log(nHi)+(log(nHf)-log(nHi))/(nnH-1)*inH);
      } else {
	nH_vec[inH]=nHi+(nHf-nHi)/(nnH-1)*inH;
      }
    }
    nH_terp = Linear_interp(nH_vec,nH_vec); 
    for (int iT=0;iT<nT;iT++) {
      if (Tgridlog) {
	T_vec[iT]=exp(log(Ti)+(log(Tf)-log(Ti))/(nT-1)*iT);
      } else {
	T_vec[iT]=Ti+(Tf-Ti)/(nT-1)*iT;
      }
    }
    T_terp = Linear_interp(T_vec,T_vec);      
    
    //set up the simulator grids and initialize 
    S_grid = new Sobj*[nnH];
    atmointerp_grid = new atmointerp*[nnH];
    Sinit_grid = new bool*[nnH];
    I_calc_grid = new VecDoub*[nnH];
    for (int inH=0;inH<nnH;inH++) {
      S_grid[inH] = new Sobj[nT];
      atmointerp_grid[inH] = new atmointerp[nT];
      Sinit_grid[inH] = new bool[nT];
      I_calc_grid[inH] = new VecDoub[nT];
      for (int iT=0;iT<nT;iT++)
	Sinit_grid[inH][iT]=FALSE;
    }
  };

  ~corona_simulator() {
    for (int inH = 0; inH < nnH; inH++) {
      delete [] Sinit_grid[inH];
      delete [] S_grid[inH];
      delete [] atmointerp_grid[inH];
      delete [] I_calc_grid[inH];
    }
    delete [] Sinit_grid;
    delete [] S_grid;
    delete [] atmointerp_grid;
    delete [] I_calc_grid;
  }
  
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
      throw("obsfile not found!");
    }
  }

  void get_S(const double nH, const double T,
	     Sobj& thisS,atmointerp& thisatmointerp,bool& thisSinit,
	     bool forcesim=FALSE) {
    /* tries to load S from file; if file not found, calculates and
       produces file for future use */

    if (!forcesim) { 
	//what would the filename look like?
	string Sfname=Sfilename(nH,T);
	ifstream Sfile(Sfname.c_str());
	if (Sfile.good()) {
	  Sfile.close();
	  thisS=Sobj(Sfname);
	  //if S file exists, atmointerp does too.
	  thisatmointerp = atmointerp(thisS.nH,thisS.T,nphyspts);
	  thisSinit=1;
	}
	Sfile.close();
    }
    if (!thisSinit) {
      //Sfile does not exist, we must create it.
      generate_S(nH, T, thisS, thisatmointerp);
      thisSinit=1;
    }

    return;
  }
  void get_S(const double nH, const double T, bool forcesim=FALSE) {
    get_S(nH, T, current_S, current_atmointerp, current_Sinit, forcesim);
  }
  
  double simulate_iobs(int iobs, 
		       Sobj &thisS, 
		       atmointerp &thisatmointerp,
		       bool thisSinit,
		       double IPHb=0.0) {
    if (!thisSinit||!obsinit) {
      std::cout << "thisSinit = " << thisSinit << std::endl;
      std::cout << "obsinit = " << obsinit << std::endl;
      
      throw("Load observation and source function before calling simulate()!");
    }
    //simulate intensity for a single observation using the current
    //Sobj and the specicified IPH background.
    /* std::cout << "Simulating coordinate iobs = " << iobs  */
    /* 	      << ", S_nH =" << current_S.nH << ", S_T = " << current_S.T << ".\n"; */

    double nH = thisS.nH;
    double T = thisS.T;

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
  void simulate_iobs(int iobs, double IPHb) {
    simulate_iobs(iobs, current_S, current_atmointerp, current_Sinit, IPHb);
  }

  void calc_I(Sobj &thisS,
	      atmointerp &thisatmointerp,
	      bool &thisSinit,
	      VecDoub &thisI_calc,
	      double IPHb=0.0) {
    //simulates for all iobs and stores the result in I_calc
    if ((!thisSinit)||(!obsinit)) {
      //error checking implemented for use of current_S only
      std::cout << "thisSinit = " << thisSinit << std::endl;
      std::cout << "obsinit = " << obsinit << std::endl;
      
      throw("Load observation and source function before calling simulate()!");
    }

    thisI_calc.resize(nobs);
    for (int iobs=0; iobs<nobs; iobs++)
      thisI_calc[iobs]=simulate_iobs(iobs,thisS,thisatmointerp,thisSinit,IPHb);
  }
  void calc_I(double IPHb=0.0) {
    std::cout << "calc_I current_Sinit = " << current_Sinit << "\n";
    calc_I(current_S, current_atmointerp, current_Sinit, current_I_calc, IPHb);
  }
  
  //check that intensities exist for interpolation, fill them in if they don't
  void grid_init(int inH, int iT,double IPHb,bool forcesim=FALSE)
  {
    for (int iinH=inH;iinH<=inH+1;iinH++) {
      for (int iiT=iT;iiT<=iT+1;iiT++) {
	if (!Sinit_grid[iinH][iiT]) {
	  //load or simulate the source function
	  get_S(nH_vec[iinH],
		T_vec[iiT],
		S_grid[iinH][iiT],
		atmointerp_grid[iinH][iiT],
		Sinit_grid[iinH][iiT],
		forcesim);
	  //calculate the intensities
	  calc_I(S_grid[iinH][iiT],
		 atmointerp_grid[iinH][iiT],
		 Sinit_grid[iinH][iiT],
		 I_calc_grid[iinH][iiT],
		 IPHb);
	}
      }
    }
  }

  
  double interp_iobs(int iobs, double nHp, double Tp, double IPHb=0.0,bool forcesim=FALSE) {
    double iIcalc, iIPHtrans, anH, aT;
    //find the grid square:
    /* std::cout << "nH = " << nHp << std::endl; */
    int inH = (nH_terp).cor ? (nH_terp).hunt(nHp) : (nH_terp).locate(nHp);
    // std::cout << "inH = " << inH << std::endl;
    // std::cout << "nH[" << inH << "] = " << (*nH_terp).xx[inH] << std::endl;
    // std::cout << "nH[" << inH + 1 << "] = " << (*nH_terp).xx[inH + 1] << std::endl;
    // std::cout << "T = " << Tp << std::endl;
    int iT = (T_terp).cor ? (T_terp).hunt(Tp) : (T_terp).locate(Tp);
    // std::cout << "iT = " << iT << std::endl;
    // std::cout << "T[" << iT << "] = " << (*T_terp).xx[iT] << std::endl;
    // std::cout << "T[" << iT + 1 << "] = " << (*T_terp).xx[iT + 1] << std::endl;
    
    // std::cout << "inH = " << inH << std::endl;
    // std::cout << "iT = " << iT << std::endl;    
    
    grid_init(inH,iT,IPHb,forcesim);//simulate or load intensities
    
    //interpolate:
    anH = (nHp-(nH_terp).xx[inH])/((nH_terp).xx[inH+1]-(nH_terp).xx[inH]);
    aT  = ( Tp-( T_terp).xx[ iT])/(( T_terp).xx[ iT+1]-( T_terp).xx[ iT]);

    /* std::cout << "aobs = " << aobs << std::endl; */
    /* std::cout << "anH = " << anH << std::endl; */
    /* std::cout << "aT = " << aT << std::endl;     */
    
    iIcalc = (1.-anH)*(1.-aT)*I_calc_grid[inH  ][iT  ][iobs]
            +    anH *(1.-aT)*I_calc_grid[inH+1][iT  ][iobs]
            +(1.-anH)*    aT *I_calc_grid[inH  ][iT+1][iobs]
            +    anH*     aT *I_calc_grid[inH+1][iT+1][iobs];

    //now add in the IPH for points with tangent altitudes above 200km
    //if iph is simulated using the Quemerais code, it is multiplied
    //in here. Otherwise, this factor was set to 1 in the initial load
    //of the obs_data and IPHb represents the brightness of the IPH
    //and not a multiplicative scaling factor.
    // iIPHtrans =  (1.-anH)*(1.-aT)*corona_simulator_grid[inH  ][iT  ].IPHtrans[iobs]
    //            +     anH *(1.-aT)*corona_simulator_grid[inH+1][iT  ].IPHtrans[iobs]
    //            + (1.-anH)*    aT *corona_simulator_grid[inH  ][iT+1].IPHtrans[iobs]
    //            +     anH *    aT *corona_simulator_grid[inH+1][iT+1].IPHtrans[iobs]
    //    iIcalc += IPHb*iIPHtrans;
    
    return iIcalc;
  }

  void interp_I(double nHp, double Tp, 
		VecDoub &thisI_calc,double IPHb=0.0, bool forcesim=FALSE) {

    thisI_calc.resize(nobs);
    for (int iobs=0; iobs<nobs; iobs++)
      thisI_calc[iobs]=interp_iobs(iobs, nHp, Tp, IPHb, forcesim);
  }


  
};

#endif
