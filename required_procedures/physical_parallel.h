//physical.h -- the physical atmosphere model, based on Krasnopolski 2002

#ifndef __PHYSICAL_H_
#define __PHYSICAL_H_

#include <cmath>
#include "quadrature.h"
#include "interp_1d.h"
#include "definitions.h"
#include <iostream>
#include <fstream>
#include <string>
#include "nr3.h"
#include <gsl/gsl_sf_gamma.h>
//#include <mpi.h> // parallelization//must be included in top line of main file

using std::exp;

double Temp(const double &r, const double &Texo)
// computes the analytic thermospheric temperature as given by Krasnopolsky (2002)
{
  const double rdiff = r - (90e5 + rMars);
  return Texo - (Texo - 125.0)*exp(-rdiff*rdiff*1e-10/(11.4*Texo));
}

double Tprime(const double &r, const double &T,const double &Texo)
// computes the derivative of the analytic temperature profile
{
  const double rdiff = r - (90e5 + rMars);
  return ( Texo - T ) * ( 2*rdiff*1e-10 / (11.4*Texo) );
}

double oneoverHn(const double &r, const double &Tt, const double &Ttprime)
// computes the inverse CO2 scale height, used below the exobase for the bulk atmosphere
{
  const double f = G*mMars*mCO2/(kB*r*r);
  return (f + Ttprime)/Tt;
}

double oneoverHncomp(const double &r,const double &Texo)
// computes the inverse CO2 scale height, used below the exobase for
// the bulk atmosphere. Computes the temperature and its derivative
// instead of having these passed.
{
  const double f = G*mMars*mCO2/(kB*r*r);
  const double Tt = Temp(r,Texo);
  return f/Tt + Tprime(r,Tt,Texo)/Tt;
}


double oneoverHH(const double &r, const double &Tt, const double &Ttprime)
// computes the inverse hydrogen scale height, used above the exobase for hydrogen
{
  const double f = G*mMars*mH/(kB*r*r);
  return (f + Ttprime)/Tt;
}

struct Hninv {
  double Texo;
  Hninv(const double &Texoo) : Texo(Texoo) {}
  double operator()(double r) {
    return oneoverHncomp(r,Texo);
  }
};

double nCO2(const double &r,const double &Texo)
//computes CO2 number density by integrating scale height
{
  Hninv invHn(Texo);
  return nr0 * exp(-qsimp(invHn, rminatm, r));
}

double DH(const double &r, const double &T, const double &dCO2)
// diffusion coefficient of hydrogen
{
  const double DH0 = 8.4e17;// cm^2 s^-1
  const double s = 0.6;
  return std::pow(T,s) * DH0/dCO2;
}

inline double KK(const double &dCO2,const double &Texo) { return 1.2e12 * std::sqrt(Texo/dCO2); } // eddy diffusion coefficient

double fone(const double &r,const double &Texo)
// first integration function
{
  const double Tt = Temp(r,Texo);
  const double Ttprime = Tprime(r, Tt, Texo);
  const double dCO2 = nCO2(r,Texo);

  const double tDH = DH(r, Tt, dCO2);
  const double tKK = KK(dCO2,Texo);
  
  const double aT = -0.25; //thermal coefficient
  
  const double toneoverHH = oneoverHH(r, Tt, Ttprime);
  const double toneoverHn = oneoverHn(r, Tt, Ttprime);

  const double Trat = Ttprime/Tt;

  return ( 1.0/(tDH + tKK) ) * ( tDH*( toneoverHH + (1.0 + aT)*Trat ) + tKK*( toneoverHn + Trat ) );
}

struct foneint {
  double Texo;
  foneint(const double &Texoo) : Texo(Texoo) {}
  double operator()(double r) {
    return fone(r,Texo);
  }
};
  

double f2(const double &r,const double &Texo)
// second integration function
{
  const double Tt = Temp(r,Texo);
  const double Ttprime = Tprime(r, Tt, Texo);
  const double dCO2 = nCO2(r,Texo);
  
  return rexo*rexo/(r*r*(DH(r, Tt, dCO2) + KK(dCO2,Texo)));
}

// struct f2int {
//   double Texo;
//   f2int(double &Texoo) : Texo(Texoo) {}
//   double operator()(double r) {
//     return f2(r,Texo);
//   }
// };

double f3(const double &r,const double &Texo)
// third integration function
{
  foneint fone_int(Texo);
  return qsimp(fone_int, rexo, r);
}

// struct f3int {
//   double Texo;
//   f3int(double &Texoo) : Texo(Texoo) {}
//   double operator()(double r) {
//     return f3(r,Texo);
//   }
// };

double f4(const double &r,const double &Texo)
//fourth integration function
{
  return f2(r,Texo)*exp(f3(r,Texo));
}

struct f4int {
  double Texo;
  f4int(const double &Texoo) : Texo(Texoo) {}
  double operator()(double r) {
    return f4(r,Texo);
  }
};

double f5(const double &r,const double &Texo)
//fifth integration function
{
  f4int f4_int(Texo);
  return qsimp(f4_int, rexo, r);
}

double nHltexo(const double &r, const double &nHexo, const double &Texo)
// computes hydrogen number density as a function of altitude by
// solving the diffusion equation
{
  const double lambdac = G*mMars*mH/(kB*Texo*rexo);//chamberlain lambda
  const double uTe = 0.5 * sqrt( 2.0*kB*Texo / (mH*pi) ) * (1.0 + lambdac) * exp(-lambdac); // effusion velocity

  return nHexo*exp(-f3(r,Texo))*(1.0 - uTe*f5(r,Texo));
}

//inline double nHltexo(const double &r) { return nHltexo(r, 1.0, Texo); } // provide default values

double nHgtexo(const double &r, const double &nHexo, const double &Texo)
// computes hhydrogen number density as a function of altitude above
// the exobase, assuming a chamberlain exosphere w/o satellite particles
{
  const double lambda = G*mMars*mH/(kB*Texo*r);//chamberlain lambda
  const double lambdac = G*mMars*mH/(kB*Texo*rexo);//chamb lambda @ rexo
  const double psione = lambda*lambda/(lambda+lambdac);
  
  const double g1 = gsl_sf_gamma_inc_P(1.5, lambda);
  const double g1c = gsl_sf_gamma_inc_P(1.5, lambdac);
  const double g2 = gsl_sf_gamma_inc_P(1.5, lambda - psione);
  
  // calculate the fraction of the exobase density at this altitude:
  // this looks different than the formula in C&H b/c the incomplete
  // gamma function in the gsl is normalized; to get gamma(3/2,x) from
  // C&H we have to multiply g1 by Gamma(3/2). This enables us to pull
  // the complete Gamma function outside the parens.
  double frac = ( 1.0 + g1 - sqrt( 1.0 - lambda*lambda/(lambdac*lambdac) ) * exp(-psione) * ( 1.0 + g2 ) );  
  // normalize to the fraction at the exobase:
  double norm = (1.0 + g1c);
  frac /= norm;

  // multiply by the exobase density and return
  return nHexo*frac*exp(lambda-lambdac);
}

//inline double nHgtexo(const double &r) { return nHgtexo(r, 1.0, Texo); }// provide default values

double nH(const double &r, const double &nHexo, const double &Texo)
// combine all the hydrogen regions into one function
{
  const double km120 = 120e5+rMars;
  
  if (r >= rexo)
    return nHgtexo(r, nHexo, Texo);
  else if (r > km120)
    return nHltexo(r, nHexo, Texo);
  else
    return nHltexo(km120, nHexo, Texo);
}

//inline double nH(const double &r) { return nH(r, 1.0, Texo); }// provide default values
//inline double nHonepar(const double &r) { return nH(r, 1.0, Texo); }// provide an integrable function


///replacing the above expensive iterated integrals with a table-based approach


struct tabular_atmo {
  //structure to perform 
  double nexo, Texo;
  int n_interp_pts;
  int n_quad_pts;
  int quad_to_interp_ratio;
  double *nHlist, *nCO2list, *f2list, *f1list, *D, *K, *Hinv, *H0inv, *T, *Tp, *alt, *logalt;
  double f3, f4, f5;
  double dalt;
  VecDoub logaltvec, lognHvec, lognCO2vec;
  Linear_interp *lognHinterp;
  Linear_interp *lognCO2interp;
  
  //constructor performs all of the tabulation and integration
tabular_atmo(double nexoo, double Texoo, int n_interp_ptss=10000, string filename="") : 
  nexo(nexoo), Texo(Texoo), n_interp_pts(n_interp_ptss) {

    quad_to_interp_ratio=100;
    n_quad_pts=quad_to_interp_ratio*n_interp_pts;

    nHlist = new double[n_quad_pts];
    nCO2list = new double[n_quad_pts];
    f2list = new double[n_quad_pts];
    f1list = new double[n_quad_pts];
    D = new double[n_quad_pts];
    K = new double[n_quad_pts];
    Hinv = new double[n_quad_pts];
    H0inv = new double[n_quad_pts];
    T = new double[n_quad_pts];
    Tp = new double[n_quad_pts];

    logalt = new double[n_quad_pts];
    alt = new double[n_quad_pts];
    const double altmin = (rminatm-rMars)/1e5; // minimum and maximum altitude
    const double altmax = 2*(rmax-rMars)/1e5;
    const double logaltmin = log(altmin);
    const double logaltmax = log(altmax);
    const double dlogalt = (logaltmax - logaltmin)/((double) n_quad_pts-1);
    double dalt = 0.0; 
    double lastalt=rminatm;
    int altexopos = 0;		/* counter for largest altitude less than exobase */

    /* variables for integration in the loop */
    double H0invint = 0.0;
    double lastH0inv = 0.0;

    /* Thermal diffusion coefficient */
    double aT=-0.25;


    for (int i = 0; i<n_quad_pts; i++) {
      logalt[i] = logaltmin + i*dlogalt;
      alt[i] = rMars+1e5*exp(
			     logalt[i]
			     );
      dalt = alt[i]-lastalt;
      lastalt = alt[i];

      T[i] = Temp(alt[i], Texo);
      Tp[i] = Tprime(alt[i], T[i], Texo);

      H0inv[i] = oneoverHn(alt[i], T[i], Tp[i]);
      H0invint += 0.5*(lastH0inv+H0inv[i])*dalt; /* Integration by
						    trapeziodal rule
						    at each step */
      lastH0inv = H0inv[i];

      nCO2list[i] = nr0 * exp(-H0invint);
      
      D[i] = DH(alt[i], T[i], nCO2list[i]);
      K[i] = KK(nCO2list[i], Texo);

      Hinv[i] = oneoverHH(alt[i], T[i], Tp[i]);    

      /* now, split the computation up for hydrogen; Chamberlain above
	 exobase, diffusion below. */
      if (alt[i] < rexo) {
	altexopos=i;

	/* do what we can on the way up */
	f1list[i] = (1.0/(D[i] + K[i]))
	  *(D[i]*(Hinv[i] + (1.0 + aT)*(Tp[i]/T[i]))
	    +K[i]*(H0inv[i] + Tp[i]/T[i]));

	f2list[i] = rexo*rexo/(alt[i]*alt[i]*(D[i] + K[i]));

      } else { /* alt >= exokm */
	nHlist[i] = nHgtexo(alt[i], nexo, Texo);
      }
    }

    /* Now loop back down from exobase to integrate H in the diffusion
       regime. */

    /* First get the values at rexo and integrate down to the first
       point on the grid, to begin the computation */
    double f1_exo=fone(rexo, Texo);
    double f2_exo=f2(rexo, Texo);
    double dtop=alt[altexopos]-rexo;

    double lastf1=f1_exo;
    f3=0.5*(lastf1+f1list[altexopos])*dtop;
    lastf1=f1list[altexopos];

    double lastf4 = f2_exo;
    f4=f2list[altexopos]*exp(f3);
    f5=0.5*(lastf4+f4)*dtop;
    lastf4=f4;

    const double lambdac = G*mMars*mH/(kB*Texo*rexo);//chamberlain lambda
    const double uTe = 0.5*(sqrt( 2.0*kB*Texo / (mH*pi) ) 
			    * (1.0 + lambdac) * exp(-lambdac)); // effusion velocity

    nHlist[altexopos] = nexo*exp(-f3)*(1.0 - uTe*f5);
    lastalt=alt[altexopos];
    
    int alt120pos=0;
    double km120=rMars+120e5;

    /* Now loop down to the next breakpoint */
    for (int i = altexopos-1; i >= 0; i--) {
      if (alt[i] < km120)	break;
      alt120pos=i;//count down to the largest altitude greater than 120km

      dalt=alt[i]-lastalt;
      lastalt=alt[i];

      f3+=0.5*(lastf1+f1list[i])*dalt;
      lastf1=f1list[i];

      f4=f2list[i]*exp(f3);
      f5+=0.5*(lastf4+f4)*dalt;
      lastf4=f4;
      
      nHlist[i] = nexo*exp(-f3)*(1.0 - uTe*f5);
    }

    /* fill the gap between the last position and 120*/
    double dgap=km120-alt[alt120pos];

    f3+=0.5*(f1list[alt120pos]+fone(km120,Texo))*dgap;
    f4=f2(km120,Texo)*exp(f3);
    f5+=0.5*(lastf4+f4)*dgap;

    double nH120 = nexo*exp(-f3)*(1.0 - uTe*f5);

    for (int i = alt120pos-1; i >= 0; i--) nHlist[i] = nH120;


    //set up interpolation vectors
    logaltvec.resize(n_interp_pts);
    lognHvec.resize(n_interp_pts);
    lognCO2vec.resize(n_interp_pts);
    for (int i = 0; i < n_interp_pts; i++) {
      logaltvec[i] = logalt[i*quad_to_interp_ratio];
      lognHvec[i] = log(nHlist[i*quad_to_interp_ratio]);
      lognCO2vec[i] = log(nCO2list[i*quad_to_interp_ratio]);
    }
    //proceed with interpolation
    lognHinterp = new Linear_interp(logaltvec, lognHvec); // interpolate
    lognCO2interp = new Linear_interp(logaltvec, lognCO2vec); // interpolate

    int rank=MPI::COMM_WORLD.Get_rank();
    int size=MPI::COMM_WORLD.Get_size();
    if (filename.length() > 0) { /* if a filename has been provided,
				    save to file */
      if (rank==master) {
	/* save to file */
	ofstream file;
	file.open(filename.c_str());
	int iw = ((int) log(n_interp_pts)/log(10))+3;
	int w = 15;
	file.width(iw);
	file << "n:";
	file.width(w);
	file << "logalt:";
	file.width(w);
	file << "alt (km)";
	file.width(w);
	file << "lognH:";
	file.width(w);
	file << "nH (cm^-3)";
	file.width(w);
	file << "lognCO2:";
	file.width(w);
	file << "nCO2 (cm^-3)"<< std::endl;
	for (int i = 0; i < n_interp_pts; i++) {
	  file.width(iw);
	  file << i;
	  file.width(w);
	  file << logaltvec[i];
	  file.width(w);
	  file << exp(logaltvec[i]);
	  file.width(w);
	  file << lognHvec[i];
	  file.width(w);
	  file << exp(lognHvec[i]);
	  file.width(w);
	  file << lognCO2vec[i];
	  file.width(w);
	  file << exp(lognCO2vec[i]) << std::endl;
	}
	file.close();
      }
    }
  }

  /* Member functions to retrieve H and CO2 densities */
  double nH(const double &r) const {
    double km = (r - rMars)/1e5;
    return exp((*lognHinterp).interp(log(km)));
  }
  double nCO2(const double &r) const {
    double km = (r - rMars)/1e5;
    return exp((*lognCO2interp).interp(log(km)));
  }
  
  ~tabular_atmo() {
    delete[] nHlist;
    delete[] nCO2list;
    delete[] f2list;
    delete[] f1list;
    delete[] D;
    delete[] K;
    delete[] Hinv;
    delete[] H0inv;
    delete[] T;
    delete[] Tp;
    delete[] logalt;
    delete[] alt;

    delete lognHinterp;
    delete lognCO2interp;
  }
};













#endif

