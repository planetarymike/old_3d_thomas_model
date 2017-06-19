//definitions.h -- model parameters that are global in nature; physical constants, etc.

#ifndef __DEF_HMARS_H
#define __DEF_HMARS_H

#include <cmath>

//basic model parameters
const double rMars = 3395e5;// cm, radius of Mars
const double mMars = .1076 * 5.98e27;// gm, Mass of Mars (from fraction of Earth mass)
const double nr0 = 2.6e13; // cm^-3, number density of CO2 atmosphere at surface
const double rexo = 3595e5;// cm, exobase altitude

const double rminatm = rMars + 80e5; // cm, minimum altitude in model atmosphere
const double rmax = 50000e5+rMars;

//number of interpolation points in the hydrogen atmosphere
const int nphyspts = 10000;//> 366!! or segfault!

//physical  constants
const double G = 6.67e-8;// dyn cm^2 gm^-2, Cavendish constant
const double kB = 1.38e-16;// erg K^-1, Boltzmann constant
const double clight = 3e10;// cm s^-1, speed of light
const double mH = 1.673e-24;// gm, mass of hydrogen atom
const double mCO2 = 44 * mH;// gm, mass of CO2 molecule

//mathematical constants
const double pi = 3.1415926535898;

//radiative transfer parameters
const double sHtot = 2.647e-2 * 0.416;// cm^2 Hz,
                                      // total cross section of Ly alpha pi*e^2/(m c) * f

const double sCO2 = 6.3e-20; //cm^2 CO2 cross section at Ly alpha

//parallel processing definitions:
const int master = 0; // rank of the master process; needed in both main() and interpgen.h

//basic boolean definitions:
const bool TRUE = 1;
const bool FALSE = 0;

//______________________________________________________
//------------------CONTROL PARAMETERS------------------
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
//smallest resolvable length
const double dsfrac = 0.05;//fraction of the smallest box dimension

//differential optical depth
const double dtau = 0.005;

//maximum tau for which HolG and HolT have been precomputed:
const double taumax = 100000.;

//max and min path length
const double dsmin = 0.1e5; // 0.1km
const double dsmax = 100e5; // 100km
  
//maximum number of iterations
const int maxit = 2000000;

//define the number of boundaries in each dimension

//select the method for distributing the radial points in
//altitude. Options are: taupts, logpts, loglinpts.
const string rmethod = "loglinpts";
const int nrpts = 160;//includes boundaries at rmin and rmax
                      //strongest dependence of computed intensities is here

//select the method for distributing the theta points in
//angle. Options are: uniform, doubleterm (more points near the terminator)
const string tmethod = "uniform";//"doubleterm";
const int nthetapts = 17;//includes boundaries at 0 and 180
                         //if using method doubleterm, must be of the form
                         //4N+5 for N=0,1,2,
                         //ie 5, 9, 13, 17, 21, 25, 29, 33, 37, ...

const int nphipts = 2;//includes boundaries at 0 and 360 so 2 implies
                      //no phi boundaries (azimuthially symmetric)

//select the method for distributing the rays in theta angle.
//options are "gaussian" and "uniform"
//phi rays are always uniform
const string raymethod = "gaussian";
//number of rays
const int ntrays = 5;//weak dependence on number of rays
const int nprays = 10;

//compute CO2 absoprtion self-consistently in the inner loop?
const bool inloop=TRUE;//bad things happen if FALSE

//use chaufray hydrogen density?
const bool usechauH = FALSE;

//grid definitions (for interpolating source functions)
//grid is populated on demand as source functions are requested.
const bool nHgridextend=TRUE;//this is a hack to extend the grid at
			     //high densities because some fits
			     //wandered off there.
const int  n_nHgridextend=19;

const int        nnH = (nHgridextend) ? 100+n_nHgridextend : 100;
const double     nHi = 10000;
const double     nHlogf = 1000000;//getting really hacky here to
				  //maintain backwards compatibility
const double     nHf = (nHgridextend) ? 1e6*(n_nHgridextend+1) : nHlogf;
const bool nHgridlog = TRUE;

const int         nT=191;
const double      Ti=100;
const double      Tf=2000;
const bool  Tgridlog=FALSE;

  
//define the interpolation data locations
const string tabdataloc="./tabulated_data/";
const string HolTfilename=tabdataloc+"HolTinterp.dat";
const string HolGfilename=tabdataloc+"HolGinterp.dat";
const string eff_filename=tabdataloc+"eff_interp.dat";

//  Holinterp HolG_lookup("./tabulated_data/HolGinterp.dat");

//define the LOS profile filename
const string losproffname=tabdataloc+"H_LOS_prof.dat";//SRCFNSLOC is
						      //passed by a
						      //macro at
						      //compile-time

//define the location to search for and store generated source functions
const string srcfnsloc="./source_functions/";
  
//______________________________________________________
//----------------END CONTROL PARAMETERS----------------
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 




#endif
