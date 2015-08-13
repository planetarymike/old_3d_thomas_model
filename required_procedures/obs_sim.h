//obs_sim.h -- generate modeled spacecraft intensities for a given observation

#ifndef __OBS_SIM_H
#define __OBS_SIM_H

#include <iostream> // for file output and dialog
#include <fstream>  // ""
#include <cmath>    // for cos and sin
#include <cstdlib>
#include "physical.h" // physical definitions, nH, nCO2, KK, etc.
#include "quadrature.h" // to compute line integrals
#include "definitions.h" // basic parameter definitions
#include "pipeline_los_defs.h" // definitions for LOS integration
#include "nr3.h" // type VecDoub
//#include "nHinterp.h" // for interpolated hydrogen function
//#include "interpgen.h" // runtime hydrogen interpolation
#include "los.h"       //line of sight integration routines
#include "obs_data.h"
#include "rad_trans.h" // two-point column density, holstein t and g functions
#include "iph_sim.h"

struct obs_simulator {
  obs_data thisobs;
  string outfilename;
  VecDoub rpts, tpts, ppts;
  string src_fn_folder;
  string ptsfile;
  string solnTptsfile;
  int nsolnH;
  VecDoub solnH;
  int nsolT;
  VecDoub solT;
  double lineintcoef;//intensity line integral coefficient "g-value"
  Sobj ***Sobjs; // array of Sobjs for each nH, T
  atmointerp ***atmointerps; //array of atmointerps for each nH, T
  bool **filesinit; //bool to see if Sobj, atmointerps have been initialized for this nH, T 
  LOS_integrator *LOS; //object to do line integration using tabulated HolT values
  IPHsim *IPH;

  
  //constructor. takes an obs_data
  obs_simulator(obs_data &thisobss) : thisobs(thisobss) {

    IPH = new IPHsim(thisobs);
    
    src_fn_folder = SRCFNSLOC;//SRCFNSLOC is passed by a macro at compile-time
    std::cout << "Looking for Model Source Functions in: " << src_fn_folder << std::endl;
    ptsfile = src_fn_folder + "pts74x19x37.dat";
    solnTptsfile = src_fn_folder + "nTpts.dat";

    //set up the line of sight integrator
    LOS = new LOS_integrator(src_fn_folder+"HolTinterp.dat");
    
    //echo observation name to terminal
    std::cout << "Observation # = " << thisobs.obsname << std::endl;

    //describe input and destination file to user
    std::cout << "Input file = " << thisobs.fname << std::endl;
    //output file
    outfilename = thisobs.obsname + "_sim.dat";
    std::cout << "Output file = " << outfilename << std::endl;

    //pause
    //  std::cin.get();
  
    //first we need to import the radial and theta points:
    getfilepts(rpts,tpts,ppts,ptsfile);//ptsfile is defined in pipeline_los_defs.h

    ifstream solnTpts;
    solnTpts.open(solnTptsfile.c_str());
    solnTpts >> nsolnH;
    //        std::cout << "nsolnH = " << nsolnH << "\n";
    //        std::cin.get();
    solnH.resize(nsolnH);
    for (int i = 0; i < nsolnH; i++){
      solnTpts >> solnH[i];
      //          std::cout << "solnH[" << i << "] = " << solnH[i] << "\n";
      //          std::cin.get();
    }

    solnTpts >> nsolT;
    //    std::cout << "nsolT = " << nsolT << "\n";
    solT.resize(nsolT);
    for (int i = 0; i < nsolT; i++) {
      solnTpts >> solT[i];
      //    std::cout << "solT[" << i << "] = " << solT[i] << "\n";
      //    std::cin.get();
    }

    //set up the Sobjs and atmointerps
    Sobjs = new Sobj**[nsolnH];
    atmointerps = new atmointerp**[nsolnH];
    filesinit = new bool*[nsolnH];
    for (int i = 0; i < nsolnH; i++) {
      Sobjs[i] = new Sobj*[nsolT];
      atmointerps[i] = new atmointerp*[nsolT];
      filesinit[i] = new bool[nsolT];
      for (int j = 0; j < nsolT; j++) {
	filesinit[i][j]=false;
      }
    }

  };

  ~obs_simulator() {
    //free the Sobjs and atmointerps
    //    std::cout << "obs_simulator destructor is being called.\n";
    for (int i = 0; i < nsolnH; i++) {
      for (int j = 0; j < nsolT; j++) {
	//free only the Sobjs and atmointerps loaded in from file
	if (filesinit[i][j]) {
	  //	  std::cout << "i, j = " << i << ", " << j << std::endl; 
  	  delete Sobjs[i][j];
  	  delete atmointerps[i][j];
	}
      }
      delete[] Sobjs[i];
      delete[] atmointerps[i];
      delete[] filesinit[i];
    }
    delete[] Sobjs;
    delete[] atmointerps;
    delete[] filesinit;

    //free the LOS integrator
    delete LOS;
  }

  void writecoords() {
    //writes the coordinates used in the simulation to file
    //uncomment some code to have the file checked before rewriting.


    //altitudes
    /* ifstream outfile; */
    string coordfilename;
    coordfilename = thisobs.obsname + "_sim_coords.dat";
    /* outfile.open(coordfilename.c_str()); */

    /* if(outfile.good()) {//outfile already exists, do not rewrite */
    /*   std::cout << coordfilename << " already exists, will not rewrite.\n"; */
    /*   //      std::cin.get(); */
    /*   outfile.close(); */
    /* } */
    /* else { //outfile does not exist, write it */
    /* outfile.close(); */
    ofstream outfile(coordfilename.c_str());

    //output altitudes
    outfile << "//alt (km); each block is at one altitude\n";
    outfile << thisobs.alttan.size() << std::endl;
    for (int i = 0; i < thisobs.alttan.size(); i++) {
      outfile.width(20);
      outfile << thisobs.alttan[i];
    }
    outfile << std::endl;

    //output densities
    outfile << "//nH (cm^-3); each row is at one density\n";
    outfile << nsolnH << std::endl;
    for (int i = 0; i < nsolnH; i++) {
      outfile.width(20);
      outfile << solnH[i];
    }
    outfile << std::endl;

    //output densities
    outfile << "//Temp (K); each column is at one temperature\n";
    outfile << nsolT << std::endl;
    for (int i = 0; i < nsolT; i++) {
      outfile.width(20);
      outfile << solT[i];
    }
    outfile << std::endl;

    outfile.close();
    /* } */
  }

  double simulate_gridpoint(int iobs, int inH, int iT, double IPHb=0.0) {
    //simulate intensity observed for a single point in altitude, density, temperature
    /* std::cout << "Simulating coordinate iobs = " << iobs  */
    /* 	      << ", inH =" << inH << ", iT = " << iT << ".\n"; */

    double nH = solnH[inH];
    double T = solT[iT];

    //from the input data, compute lineint = Fsun * doppler width (Hz) * sqrt(pi)
    lineintcoef = thisobs.Fsun_mars*(1/121.6e-7)*sqrt(2*pi*kB*T/mH); 
    //^^^^ photons/s/cm2/Hz * Hz = photons/cm2/s.

    if (!(filesinit[inH][iT])) {
      //if files have not been loaded, load them
      char Sfile[200];
      sprintf(Sfile, "%sS_nH%d_T%d.dat", src_fn_folder.c_str(),(int) nH,(int) T);
      //	  std::cout << "Sfile = " << Sfile << "\n";
      Sobjs[inH][iT] = new Sobj(rpts,tpts,ppts,Sfile);
      char atmofile[200];
      sprintf(atmofile, "%satm_nH%d_T%d.dat", src_fn_folder.c_str(),(int) nH,(int) T);
      atmointerps[inH][iT] = new atmointerp(nH, T, nphyspts, atmofile);
      filesinit[inH][iT]=true;
    }
          
    double tIcalc = 0.0;
    VecDoub posvec(3,thisobs.pos[iobs]);
    VecDoub dirvec(3,thisobs.dir[iobs]);
    tIcalc = (*LOS).integrate(*(Sobjs[inH][iT]),
			   *(atmointerps[inH][iT]),
			   posvec,
			   dirvec,
			   lineintcoef);// photons/cm2/s
    //convert to rayleighs
    tIcalc /= 1e6;// megaphoton/cm2/s
    tIcalc *= 4*pi/1e3; // now we're in kR; see C&H pg. 280-282

    //now add in the IPH
    tIcalc += (*IPH).sim(IPHb,nH,T,iobs);
    
    return tIcalc;
      
    //	  std::cout << "tIcalc = " << tIcalc << "\n";
    //	  std::cout << "I_obs = " << thisobs.I_obs[iobs]
    //		    << " +- " << thisobs.DI_obs[iobs] << " kR,\n";
    //	  std::cin.get();
  }

  
  void simulate_allcords(double IPHb=0.0) {
    //simulates for all grid points and writes simulation output to file

    //check to see whether simulation output file exists before re-simulating
    ifstream checkoutfile(outfilename.c_str());

    if(checkoutfile.good()) { //file exists, do not rewrite
      std::cout << outfilename << " already exists, will not rewrite.\n";
      //      std::cin.get();
      checkoutfile.close();
    } else { //file does not exist, simulate away!
      std::cout << "Simulation file does not exist, beginning simulation.\n";
      checkoutfile.close();

      //all of the data has been read in already to the obs_data object thisobs.
      //all we need to do is call the appropriate member object
      
      //allocate the storage array for Icalc
      double ***I_calc;
      I_calc = new double**[thisobs.ndat];
      for (int iobs = 0; iobs < thisobs.ndat; iobs++) {
	I_calc[iobs] = new double*[nsolnH];
	for (int inH = 0; inH < nsolnH; inH++) {
	  I_calc[iobs][inH] = new double[nsolT];
	}
      }

      std::cout << "Beginning simulation for each altitude, density, and temperature." << std::endl;
    
      /* X axis points toward Sun!! */
      //this parallelizes the loop:
      //#pragma omp parallel for
      for (int iobs = 0; iobs < thisobs.ndat; iobs++) {
	for (int inH = 0; inH < nsolnH; inH++) {
	  for (int iT = 0; iT < nsolT; iT++) {
	    I_calc[iobs][inH][iT] = simulate_gridpoint(iobs,inH,iT,IPHb);
	  }
	}
	std::cout << "Done with altitude #" << iobs+1 << " of " << thisobs.ndat << std::endl;
      }

      //write the data to file
      ofstream outfile;
      //set up output file to do some chi-sq
      outfile.open(outfilename.c_str());
      outfile.precision(10);
      outfile.setf(ios_base::scientific);

      for (int iobs = 0; iobs < thisobs.ndat; iobs++) {
	outfile << (iobs + 1) << std::endl;
	for (int inH = 0; inH < nsolnH; inH++) {
	  for (int iT = 0; iT < nsolT; iT++) {
	    outfile.width(20);
	    outfile << I_calc[iobs][inH][iT];
	  }
	  outfile << std::endl;
	}
      }
      outfile.close();

      //free the Icalc array
      for (int iobs = 0; iobs < thisobs.ndat; iobs++) {
	for (int inH = 0; inH < nsolnH; inH++) {
	  delete [] I_calc[iobs][inH];
	}
	delete [] I_calc[iobs];
      }
      delete [] I_calc;
    }
  }

};

#endif
