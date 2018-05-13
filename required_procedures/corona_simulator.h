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
#include "iph_model_interface.h"

struct corona_simulator {

  bool obsinit;
  bool forcesim;
  string obsname;
  string obsfname;
  int nobs;
  double aMars; // Sun-Mars distance at time of obs, AU
  double Fsun_earth; // Solar flux at line center at Earth
  double Fsun_quemerais; //solar flux for quemerais code
  double Fsun_mars; // Solar flux at line center at Mars, possibly calculated from above
  VecDoub marspos;
  char datetime[1024];

  VecDoub obs_vec;
  VecDoub I_obs, DI_obs;  // observed intensity (kR) and its uncertainty
  VecDoub alttan, SZAtan; // radius and SZA of tangent point
  VecDoub ra, dec, scvel; //ra and dec of los, and S/C SSB vel projected onto this direction
  MatDoub pos, dir; // position and LOS direction of spacecraft.
  IPHsim IPH; //IPH simulator
  VecDoub IPHb_model; //storage for calculated IPH brightnesses
  
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
  corona_simulator(bool forcesimm=FALSE) : LOS(HolTfilename), forcesim(forcesimm) {
    //nothing is initialized at first
    obsinit  = 0;
    current_Sinit = 0;
    
    nH_vec.resize(nnH);
    if (nHgridextend) {
        //hack to extend grid to larger densities
        for (int inH=0;inH<nnH;inH++) {
            if (inH<nnH-n_nHgridextend) {
                if (nHgridlog) {
                    nH_vec[inH]=exp(log(nHi)+(log(nHlogf)-log(nHi))/(nnH-n_nHgridextend-1)*inH);
                } else {
                    nH_vec[inH]=nHi+(nHf-nHi)/(nnH-1)*inH;
                }
                //	  std::cout << "nH_vec[" << inH << "] = " << nH_vec[inH] << std::endl;
            } else {
                nH_vec[inH]=nHlogf*(inH-(nnH-n_nHgridextend)+2);
                //	  std::cout << "nH_vec[" << inH << "] = " << nH_vec[inH] << std::endl;
            }
        }
        //      std::cin.get();
    } else {
        nH_vec.resize(nnH);
        for (int inH=0;inH<nnH;inH++) {
            if (nHgridlog) {
                nH_vec[inH]=exp(log(nHi)+(log(nHf)-log(nHi))/(nnH-1)*inH);
            } else {
                nH_vec[inH]=nHi+(nHf-nHi)/(nnH-1)*inH;
            }
        }
    }
    nH_terp = Linear_interp(nH_vec,nH_vec);
      
    T_vec.resize(nT);
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
    IPHtrans_grid = new VecDoub*[nnH];
    for (int inH=0;inH<nnH;inH++) {
      S_grid[inH] = new Sobj[nT];
      atmointerp_grid[inH] = new atmointerp[nT];
      Sinit_grid[inH] = new bool[nT];
      I_calc_grid[inH] = new VecDoub[nT];
      IPHtrans_grid[inH] = new VecDoub[nT];
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
      delete [] IPHtrans_grid[inH];
    }
    delete [] Sinit_grid;
    delete [] S_grid;
    delete [] atmointerp_grid;
    delete [] I_calc_grid;
    delete [] IPHtrans_grid;
  }
  
  void obs_import(string obsname,bool simulate_iph=TRUE) {
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
      Fsun_quemerais=Fsun_earth;
      Fsun_earth *= (1215.)*(1215e-8)/(3e10/*cm/s*/);// photons/s/cm2/Hz
      Fsun_mars = Fsun_earth/aMars/aMars;
      /* std::cout << "Fsun_mars = " << Fsun_mars << std::endl; */
      marspos.resize(3);
      obsfile >> marspos[0]; //ecliptic coordinates of Mars at the time of the observation
      obsfile >> marspos[1]; //      these coordinates are used to determine the model 
      obsfile >> marspos[2]; //      IPH brightness that seeds the background subtraction.
      std::cout << "mars_pos = [" 
		<< marspos[0] << ", " 
		<< marspos[1] << ", "
		<< marspos[2] << "]"<< std::endl;
      obsfile.getline(datetime,1024); // capture time of file (irrelevant for this code)

      //prepare vector parameters for read-in
      obs_vec.resize(nobs);
      I_obs.resize(nobs);
      DI_obs.resize(nobs);
      alttan.resize(nobs);
      SZAtan.resize(nobs);
      pos.resize(nobs,3);
      dir.resize(nobs,3);
      ra.resize(nobs);
      dec.resize(nobs);
      scvel.resize(nobs);
      IPHb_model.resize(nobs);
      
      /* X axis points toward Sun!! */
      //loop over the number of observations, reading each in to the appropriate variable
      for (int row = 0; row < nobs; row++) {
	//get data from file
	obsfile >> row;
	obs_vec[row] = row;
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

	if(simulate_iph) {
	  //get the iph brightness from the Quemerais model
	  /* std::cout << "Getting brightness:\n" */
	  /* 	    << "Fsun_earth = " << Fsun_quemerais << " photons/s/cm2/Angstrom\n" */
	  /* 	    << "ra  [" << row << "] = " << ra[row] << "\n" */
	  /* 	    << "dec [" << row << "] = " << dec[row] << "\n"; */
	  IPHb_model[row] = quemerais_iph_model(Fsun_quemerais, marspos, ra[row], dec[row]);
	  /* std::cout << " Brightness = " << IPHb_model[row] << std::endl; */
	} else {
	  //assume IPH is uniform across all observations in this file
	  IPHb_model[row] = 1.0;
	}
			 
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
      toss("obsfile not found!");
    }
  }

  void get_S(const double nH, const double T,
	     Sobj& thisS,atmointerp& thisatmointerp,bool& thisSinit) {
    /* tries to load S from file; if file not found, calculates and
       produces file for future use */

    if (!forcesim) { 
	//what would the filename look like?
	string Sfname=Sfilename(nH,T);
        std::cout << "trying to load from " << Sfname << std::endl;
	ifstream Sfile(Sfname.c_str());
	if (Sfile.good()) { //file named like this exists?
	  Sfile.seekg(0, ios::end);
	  if (Sfile.tellg() != 0) { //file length isn't zero?
	    Sfile.close();
	    thisS=Sobj(Sfname);
	    //if S file exists, atmointerp does too.
	    thisatmointerp = atmointerp(thisS.nH,thisS.T,nphyspts);
	    thisSinit=1;
	  }
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
  void get_S(const double nH, const double T) {
    get_S(nH, T, current_S, current_atmointerp, current_Sinit);
  }
  
  double simulate_IPH_transmission(int iobs, double nH, double T) {
    //simulate IPH transmission for a single point in altitude, density, temperature
    /* std::cout << "Simulating coordinate iobs = " << iobs  */
    /* 	      << ", inH =" << inH << ", iT = " << iT << ".\n"; */

    double ttransfrac;
    if (alttan[iobs]>120) {
      ttransfrac = IPH.sim(1.0,nH,T,iobs);
    } else {
      ttransfrac = 0.0;
    }
    return ttransfrac;
  }


  double simulate_iobs(int iobs, 
		       Sobj &thisS, 
		       atmointerp &thisatmointerp,
		       bool thisSinit,
		       double IPHb=0.0) {
    if (!thisSinit||!obsinit) {
      std::cout << "thisSinit = " << thisSinit << std::endl;
      std::cout << "obsinit = " << obsinit << std::endl;
      
      toss("Load observation and source function before calling simulate()!");
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
    tIcalc += IPHb*IPHb_model[iobs]*simulate_IPH_transmission(iobs, nH, T);

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
      
      toss("Load observation and source function before calling simulate()!");
    }

    thisI_calc.resize(nobs);
    for (int iobs=0; iobs<nobs; iobs++)
      thisI_calc[iobs]=simulate_iobs(iobs,thisS,thisatmointerp,thisSinit,IPHb);
  }
  void calc_I(double IPHb=0.0) {
    std::cout << "calc_I current_Sinit = " << current_Sinit << "\n";
    calc_I(current_S, current_atmointerp, current_Sinit, current_I_calc, IPHb);
  }

  void calc_IPHb_transmission(double nH, double T, VecDoub &thisIPHtrans){
    thisIPHtrans.resize(nobs);
    for (int iobs=0;iobs<nobs;iobs++)
      thisIPHtrans[iobs]=IPHb_model[iobs]*simulate_IPH_transmission(iobs, nH, T);
  }
  
  //check that intensities exist for interpolation, fill them in if they don't
  void grid_init(int inH, int iT)
  {
    for (int iinH=inH;iinH<=inH+1;iinH++) {
      for (int iiT=iT;iiT<=iT+1;iiT++) {
	if (!Sinit_grid[iinH][iiT]) {
	  //load or simulate the source function
	  get_S(nH_vec[iinH],
		T_vec[iiT],
		S_grid[iinH][iiT],
		atmointerp_grid[iinH][iiT],
		Sinit_grid[iinH][iiT]);
	  //calculate the intensities
	  calc_I(S_grid[iinH][iiT],
		 atmointerp_grid[iinH][iiT],
		 Sinit_grid[iinH][iiT],
		 I_calc_grid[iinH][iiT],
		 0.0);//IPHb=0 here because it is added in later by interp_iobs
	  //calculate the IPH reduction factor
	  calc_IPHb_transmission(nH_vec[iinH],
				 T_vec[iiT],
				 IPHtrans_grid[iinH][iiT]);
	}
      }
    }
  }

  
  double interp_iobs(int iobs, double nHp, double Tp, double IPHb=0.0) {
    if (nHp<nHi||nHp>nHf||Tp<Ti||Tp>Tf) {
      std::cout << "nH or T outside of parameter grid!\n"
		<< " nH = " << nHp << std::endl
		<< "  T = " <<  Tp << std::endl;
      return -1;
    } else {
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
      
      grid_init(inH,iT);//simulate or load intensities
      
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
      iIPHtrans =  (1.-anH)*(1.-aT)*IPHtrans_grid[inH  ][iT  ][iobs]
                 +     anH *(1.-aT)*IPHtrans_grid[inH+1][iT  ][iobs]
                 + (1.-anH)*    aT *IPHtrans_grid[inH  ][iT+1][iobs]
                 +     anH *    aT *IPHtrans_grid[inH+1][iT+1][iobs];
      iIcalc += IPHb*iIPHtrans;
      
      return iIcalc;
    }
  }

  void interp_I(double nHp, double Tp, 
		VecDoub &thisI_calc,double IPHb=0.0) {

    thisI_calc.resize(nobs);
    for (int iobs=0; iobs<nobs; iobs++)
      thisI_calc[iobs]=interp_iobs(iobs, nHp, Tp, IPHb);
  }

  double IPHtrans_interp(int iobs, double nHp, double Tp) {
    int inH, iT;
    double iIPHtrans, anH, aT;
    //find the grid square:
    inH = (nH_terp).cor ? (nH_terp).hunt(nHp) : (nH_terp).locate(nHp);
    // std::cout << "inH = " << inH << std::endl;
    // std::cout << "nH[" << inH << "] = " << (*nH_terp).xx[inH] << std::endl;
    // std::cout << "nH[" << inH + 1 << "] = " << (*nH_terp).xx[inH + 1] << std::endl;
    // std::cout << "T = " << Tp << std::endl;
    iT = (T_terp).cor ? (T_terp).hunt(Tp) : (T_terp).locate(Tp);
    // std::cout << "iT = " << iT << std::endl;
    // std::cout << "T[" << iT << "] = " << (*T_terp).xx[iT] << std::endl;
    // std::cout << "T[" << iT + 1 << "] = " << (*T_terp).xx[iT + 1] << std::endl;
    
    // std::cout << "inH = " << inH << std::endl;
    // std::cout << "iT = " << iT << std::endl;    

    grid_init(inH,iT);//simulate or load intensities
    
    //interpolate:
    anH = (nHp-(nH_terp).xx[inH])/((nH_terp).xx[inH+1]-(nH_terp).xx[inH]);
    aT = (Tp-(T_terp).xx[iT])/((T_terp).xx[iT+1]-(T_terp).xx[iT]);

    /* std::cout << "anH = " << anH << std::endl; */
    /* std::cout << "aT = " << aT << std::endl;     */
    
    //now add in the IPH for points with tangent altitudes above 200km
    //if iph is simulated using the Quemerais code, it is multiplied
    //in here. Otherwise, this factor was set to 1 in the initial load
    //of the obs_data and IPHb represents the brightness of the IPH
    //and not a multiplicative scaling factor.
    iIPHtrans =  (1.-anH)*(1.-aT)*IPHtrans_grid[inH  ][iT  ][iobs]
               +     anH *(1.-aT)*IPHtrans_grid[inH+1][iT  ][iobs]
               + (1.-anH)*    aT *IPHtrans_grid[inH  ][iT+1][iobs]
               +     anH *    aT *IPHtrans_grid[inH+1][iT+1][iobs];
    
    return iIPHtrans;
  }

  void derivs(int iobs, double nHp, double Tp, double IPHb, VecDoub_O &dIcalc) {
    //find the grid square:
    int inH = (nH_terp).cor ? (nH_terp).hunt(nHp) : (nH_terp).locate(nHp);
    int iT = (T_terp).cor ? (T_terp).hunt(Tp) : (T_terp).locate(Tp);

    grid_init(inH,iT);//simulate or load intensities
    
    double nH1=(nH_terp).xx[inH];
    double nH2=(nH_terp).xx[inH+1];
    double T1=(T_terp).xx[iT];
    double T2=(T_terp).xx[iT+1];

    double fnH=(nHp-nH1)/(nH2-nH1);
    double fT=(Tp-T1)/(T2-T1);

    double I11=I_calc_grid[inH  ][iT  ][iobs];
    double I12=I_calc_grid[inH  ][iT+1][iobs];
    double I21=I_calc_grid[inH+1][iT  ][iobs];
    double I22=I_calc_grid[inH+1][iT+1][iobs];

    double trans11=IPHtrans_grid[inH  ][iT  ][iobs];
    double trans12=IPHtrans_grid[inH  ][iT+1][iobs];
    double trans21=IPHtrans_grid[inH+1][iT  ][iobs];
    double trans22=IPHtrans_grid[inH+1][iT+1][iobs];

    //we must bilinearly interpolate the derivatives
    //intensity
    double dIdnH, dIdT;
    dIdnH  = (1.-fT)*(I21-I11);
    dIdnH +=     fT *(I22-I12);
    dIdnH /= (nH2-nH1);

    dIdT  = (1.-fnH)*(I12-I11);
    dIdT +=     fnH *(I22-I21);
    dIdT /= (T2-T1);

    //transmission
    double dtransdnH, dtransdT;
    dtransdnH  = (1.-fT)*(trans21-trans11);
    dtransdnH +=     fT *(trans22-trans12);
    dtransdnH /= (nH2-nH1);

    dtransdT  = (1.-fnH)*(trans12-trans11);
    dtransdT +=     fnH *(trans22-trans21);
    dtransdT /= (T2-T1);

    dIcalc[0] = dIdnH + IPHb*dtransdnH;
    dIcalc[1] = dIdT  + IPHb*dtransdT;
  }

  //this is used for the chi-square minimzation by L-M
  void operator()(const double obsp, VecDoub_I &a, double &iIcalc, VecDoub_O &dIcalcda)
  {
    double nHp = a[0];
    double Tp = a[1];
    double IPHb = a[2];
    double cal = a[3];
    
    double tempcalc = interp_iobs(obsp, nHp, Tp, IPHb);
    iIcalc = tempcalc*cal;

    VecDoub dItemp(2);
    derivs(obsp, nHp, Tp, IPHb, dItemp);
    dIcalcda[0]=dItemp[0];
    dIcalcda[1]=dItemp[1];
    //the derivative wrt the IPH parameter is computed as the product
    //of the transmission and the iphb_calc parameter. If IPHb is a
    //multiplicative factor on the simulated brightness, this
    //multplies in the brightness simulated. If IPHb is a scalar
    //background, iphb_calc is set to 1 in the initial load of the
    //obs_data.
    dIcalcda[2]=IPHtrans_interp(obsp, nHp, Tp);
    dIcalcda[3]=tempcalc;
  }



  
};

#endif
