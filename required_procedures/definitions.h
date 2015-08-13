//definitions.h -- model parameters that are global in nature; physical constants, etc.

#ifndef __DEF_HMARS_H
#define __DEF_HMARS_H

#include <cmath>

//basic model parameters
const double rMars = 3395e5;// cm, radius of Mars
const double mMars = .1076 * 5.98e27;// gm, Mass of Mars (from fraction of Earth mass)
const double nr0 = 2.6e13; // cm^-3, number density of CO2 atmosphere at surface
const double rexo = 3595e5;// cm, exobase altitude

//const double Tinf = 240.0;// K, temperature of the atmosphere at infinity; superseded by command line input
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
//const double Hn0test = 2.5e5; // fiducial value for hydrogen density; superseded by command line inputs
const double sHtot = 2.647e-2 * 0.416;// cm^2 Hz,
                                      // total cross section of Ly alpha pi*e^2/(m c) * f

const double sCO2 = 6.3e-20; //cm^2 CO2 cross section at Ly alpha
//const double sH = 5.96e-12/sqrt(Tinf); //effective H cross section at Ly alpha; function of
                                         //temperature throughout atmosphere; now defined in physical.h

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
const double dtau = 0.01;

//max and min path length
const double dsmin = 0.1e5; // 0.1km
const double dsmax = 100e5; // 100km
  
//maximum number of iterations
const int maxit = 100000;

//define the number of boundaries in each dimension
const int nrpts = 74;//includes boundaries at rmin and rmax
const int nthetapts = 19;//includes boundaries at 0 and 180
const int nphipts = 2;//includes boundaries at 0 and 360
  
//number of rays
const int ntrays = 6;
const int nprays = 12;

//use chaufray hydrogen density?
const bool usechauH = FALSE;

//select the method for distributing the radial points in
//altitude. Options are: taupts, chaupts, logpts, loglinpts.
const string rmethod = "taupts";

//select the method for distributing the rays in angle.
//default is gauss-jordan quadrature in theta, uniform in phi
const string raymethod = "gaussian";
  
//define the HolG and HolT interpolation function locations
//SRCFNSLOC is passed by a macro at compile time
const string srcfnsloc=SRCFNSLOC;
const string HolTfilename=srcfnsloc+"HolTinterp.dat";
//  Holinterp HolG_lookup("./tabulated_data/HolGinterp.dat");

//define the LOS profile filename
const string losproffname=srcfnsloc+"H_LOS_prof.dat";//SRCFNSLOC is passed by a macro at compile-time
  

  
//some radiative transfer and atmospheric parameters are defined in definitions.h

//______________________________________________________
//----------------END CONTROL PARAMETERS----------------
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 




#endif
