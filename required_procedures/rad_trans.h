//rad_trans.h -- routines for radiative transfer


#ifndef __RAD_TRANS_H
#define __RAD_TRANS_H

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include "quadrature.h"
#include "definitions.h"
#include "gauss_wgts.h" // obtain gaussian quadrature points
#include "physical.h" // physical definitions, nH, nCO2, KK, etc.
#include "roots.h"
#include "interpgen.h" // hydrogen interpolation at runtime

using std::log;
using std::exp;
using std::pow;


inline double sH(const double &T) { 
  /* std::cout << "sHcentercoef = " << sHcentercoef << std::endl; */
  /* std::cin.get(); */
  return sHcentercoef/sqrt(T); 
} //effective H cross section at Ly alpha line center

inline double dtau_H(const double &r, atmointerp &thisatmointerp) {
  double temp = sH(Temp(r,thisatmointerp.Texo))*thisatmointerp.nH(r);
  return temp;
}

inline double dtau_H(const double &r,const double &t, atmointerp &thisatmointerp) {
  double temp = sH(Temp(r,thisatmointerp.Texo))*thisatmointerp.nH(r,t);
  return temp;
}

inline double dtau_H(const double &r,
		     const double &t,
		     const double &p, 
		     atmointerp &thisatmointerp) {
  return sH(Temp(r,thisatmointerp.Texo))*thisatmointerp.nH(r,t,p);
}


struct dtau_H_int {
  atmointerp thisatmointerp;
  dtau_H_int(atmointerp thisatmointerpp) : thisatmointerp(thisatmointerpp) {}
  double operator()(double r) {
    return dtau_H(r,thisatmointerp);
  }
  double operator()(VecDoub r0) {
    double r=0.0;
    for (int i=0; i<r0.size();i++) {
      r+=r0[i]*r0[i];
    }
    r=sqrt(r);
    return dtau_H(r,thisatmointerp);
  }
};

inline double dtau_CO2(const double &r, atmointerp &thisatmointerp) {
  return sCO2*thisatmointerp.nCO2(r);
}

inline double dtau_CO2(const double &r, const double &t, atmointerp &thisatmointerp) {
  return sCO2*thisatmointerp.nCO2(r,t);
}

inline double dtau_CO2(const double &r, 
		       const double &t, 
		       const double &p,
		       atmointerp &thisatmointerp) {
  return sCO2*thisatmointerp.nCO2(r,t,p);
}


/* inline double dtau_CO2(const double &r, const double &Texo) { */
/*   return sCO2*nCO2(r,Texo); */
/* } */

struct dtau_CO2_int {
  atmointerp thisatmointerp;
  dtau_CO2_int(atmointerp thisatmointerpp) : thisatmointerp(thisatmointerpp) {}
  double operator()(double r) {
    return dtau_CO2(r,thisatmointerp);
  }
  double operator()(VecDoub r0) {
    double r=0.0;
    for (int i=0; i<r0.size();i++) {
      r+=r0[i]*r0[i];
    }
    r=sqrt(r);
    return dtau_CO2(r,thisatmointerp);
  }
};

/* struct dtau_CO2_int { */
/*   double Texo; */
/*   dtau_CO2_int(const double &Texoo) : Texo(Texoo) {} */
/*   double operator()(double r) { */
/*     return dtau_CO2(r,Texo); */
/*   } */
/* }; */

inline double dtau_tot(const double &r, atmointerp &thisatmointerp) {
  return dtau_H(r,thisatmointerp)+dtau_CO2(r,thisatmointerp);
}

struct dtau_tot_int {
  atmointerp thisatmointerp;
  dtau_tot_int(atmointerp thisatmointerpp) : thisatmointerp(thisatmointerpp) {}
  double operator()(double r) {
    return dtau_tot(r,thisatmointerp);
  }
};


struct HolGaux
//auxiliary functor for computing the Holstein T function
{
  double tau;
  HolGaux(const double tauu) : tau(tauu) {};
  double operator()(double x) {
    return exp(-2.0*x*x-tau*exp(-x*x));
  }
};

double HolG(const double tau)
{
  //  std::cout << tau << std::endl;
  HolGaux Gt(tau);
  return qsimpinf(Gt)/sqrt(pi);
}

struct HolTaux
//auxiliary functor for computing the Holstein T function
{
  double tau;
  HolTaux(const double tauu) : tau(tauu) {};
  double operator()(double x) {
    return exp(-x*x-tau*exp(-x*x));
  }
};

double HolT(const double tau)
{
  HolTaux Tt(tau);
  return qsimpinf(Tt)/sqrt(pi);
}

/* template <class T> */
/* double col_int_rsym(T &func, VecDoub &r1, VecDoub &r2) */
/*  // return column density of the passed function between the two */
/*  // points r1 and r2 */
/* { */
/*   //extract cartesian coordinates */
/*   double x1 = r1[0], y1 = r1[1], z1 = r1[2]; */
/*   double x2 = r2[0], y2 = r2[1], z2 = r2[2]; */
  
/*   if (r1[0] == r2[0] && r1[1] == r2[1] && r1[2] == r2[2]) {  */
/*     // column density is zero if the points are the same */
/*     return 0.0; */
/*   } else { */
/*     // compute differences */
/*     double dx = x2 - x1; */
/*     double dy = y2 - y1; */
/*     double dz = z2 - z1; */

/*     // here's the line parameter of the minimum point: */
/*     double smin = -( x1*dx + y1*dy +z1*dz )/( dx*dx + dy*dy + dz*dz ); */

/*     //compute the radius at the minimum point */
/*     double lineminr = (x1+smin*dx)*(x1+smin*dx) +  */
/*                       (y1+smin*dy)*(y1+smin*dy) +  */
/* 		      (z1+smin*dz)*(z1+smin*dz); */
/*     lineminr = sqrt(lineminr); */
    
/*     if ( smin > 0.0 && smin < 1.0 && lineminr < rminatm ) { */
/*       // if the minimum altitude is less than the specified minimum, */
/*       // the column density is infinite */
/*       return 1e99; */
/*     } else { // actually compute the column density */
/*       //std::cout << "r1 = { " << x1 << " , " << y1 << " , " << z1 << " }" << std::endl; */
/*       //std::cout << "r2 = { " << x2 << " , " << y2 << " , " << z2 << " }" << std::endl; */
/*       return lineqsimp_rsym(func, r1, r2, 1.0e-4); */
/*     } */
/*   } */
/* } */

/* template <class T> */
/* double col_int(T &func, VecDoub &r1, VecDoub &r2) */
/*  // return column density of the passed function between the two */
/*  // points r1 and r2 */
/* { */
/*   //extract cartesian coordinates */
/*   double x1 = r1[0], y1 = r1[1], z1 = r1[2]; */
/*   double x2 = r2[0], y2 = r2[1], z2 = r2[2]; */
  
/*   if (r1[0] == r2[0] && r1[1] == r2[1] && r1[2] == r2[2]) {  */
/*     // column density is zero if the points are the same */
/*     return 0.0; */
/*   } else { */
/*     // compute differences */
/*     double dx = x2 - x1; */
/*     double dy = y2 - y1; */
/*     double dz = z2 - z1; */

/*     // here's the line parameter of the minimum point: */
/*     double smin = -( x1*dx + y1*dy +z1*dz )/( dx*dx + dy*dy + dz*dz ); */

/*     //compute the radius at the minimum point */
/*     double lineminr = (x1+smin*dx)*(x1+smin*dx) +  */
/*                       (y1+smin*dy)*(y1+smin*dy) +  */
/* 		      (z1+smin*dz)*(z1+smin*dz); */
/*     lineminr = sqrt(lineminr); */
    
/*     if ( smin > 0.0 && smin < 1.0 && lineminr < rminatm ) { */
/*       // if the minimum altitude is less than the specified minimum, */
/*       // the column density is infinite */
/*       return 1e99; */
/*     } else { // actually compute the column density */
/*       //std::cout << "r1 = { " << x1 << " , " << y1 << " , " << z1 << " }" << std::endl; */
/*       //std::cout << "r2 = { " << x2 << " , " << y2 << " , " << z2 << " }" << std::endl; */
/*       return lineqsimp(func, r1, r2, 1.0e-4); */
/*     } */
/*   } */
/* } */



////three argument function useful for interpolation
//void getquadpts(VecDoub &rpts, VecDoub &tpts,const Hinterp &nHinterp) {
//  //create dummies for everything else
//  int nr = rpts.size(), nt = tpts.size();
//  VecDoub rwts(nr), twts(nt), ppts(nt), pwts(nt);
//  bool printout = 0;
//  getquadpts(rpts,rwts,tpts,twts,ppts,pwts,printout,nHinterp);
//  return;
//}
//
//
//void getquadpts_thetain(VecDoub &rpts, VecDoub &rwts, 
//			VecDoub &thetapts, VecDoub &thetawts,
//			VecDoub &phipts, VecDoub &phiwts, bool &printout, const Hinterp &nHinterp)
////gets quadrature points and weights and stores them in the passed
////vector objects. differs from the above because GL quadrature is used in theta.
//{
//
//  // get radial points that are uniformly spaced in optical depth:
//  taufracpts(rpts,rwts,printout,nHinterp);
//
//  //for theta, simply Gauss-Legendre quadrature
//  gauleg(0.0,pi,thetapts,thetawts);
//
//  if (printout) {
//    double tot = 0.0;
//    for (int i = 0; i < thetapts.size(); i++) {
//      std::cout << "i = " << i << std::endl;
//      std::cout << "thetapts[" << i << "] = " << thetapts[i] << std::endl;
//      std::cout << "thetawts[" << i << "] = " << thetawts[i] << std::endl;
//      tot += thetawts[i];
//    }
//    std::cout << "tot = " << tot << std::endl;
//  }
//  
//  //for phi, trapezoidal quadrature is sufficient since the problem is
//  //symmetric in this coordinate
//  //  tot = 0.0;
//  for (int i = 0; i < phipts.size(); i++) {
//    phipts[i] = i*2.0*pi/phipts.size();
//    phiwts[i] = 2.0*pi/phipts.size();
//    //    tot += phiwts[i];
//  }
//  //  std::cout << "tot = " << tot << std::endl;
//
//  //print points to file
//  printpts(rpts,rwts,thetapts,thetawts,phipts,phiwts,"pts_thetain");
//  return ;
//  
//}
//
////two argument function useful for interpolation
//void getquadpts_thetain(VecDoub &rpts, VecDoub &tpts, const Hinterp &nHinterp) {
//  //create dummies for everything else
//  int nr = rpts.size(), nt = tpts.size();
//  VecDoub rwts(nr), twts(nt), ppts(nt), pwts(nt);
//  bool printout = 0;
//  getquadpts_thetain(rpts,rwts,tpts,twts,ppts,pwts,printout,nHinterp);
//  return;
//}


void getraypts(VecDoub &thetapts, VecDoub &thetawts,
	       VecDoub &phipts, VecDoub &phiwts,
	       const bool &printout, const string & raymethod)
/* gets ray angular points and weights, using simple trapezoidal
   quadrature in phi and G-L quadrature in theta.*/
{

 //for theta, simply Gauss-Legendre quadrature
  if (raymethod=="gaussian") {
    gauleg(0.0,pi,thetapts,thetawts);
  } else if (raymethod=="uniform") {
    for (int i = 0; i < thetapts.size(); i++) {
      thetapts[i] = i*pi/(thetapts.size()-1);
      thetawts[i] = pi/(thetapts.size()-1);
    }
    thetawts[0]/=2;
    thetawts[thetapts.size()-1]/=2;
  } else {
  std::cout << "raymothod \"" << raymethod << "\" not supported!\n"; 
  throw("raymethod not supported.");
 }

 if (printout) {
   double tot = 0.0;
   for (int i = 0; i < thetapts.size(); i++) {
     std::cout << "i = " << i << std::endl;
     std::cout << "thetapts[" << i << "] = " << thetapts[i] << std::endl;
     std::cout << "thetawts[" << i << "] = " << thetawts[i] << std::endl;
     tot += thetawts[i];
   }
   std::cout << "tot = " << tot << std::endl;
   //   std::cin.get();
 }
 
 //for phi, trapezoidal quadrature
 for (int i = 0; i < phipts.size(); i++) {
   phipts[i] = i*2.0*pi/phipts.size();
   phiwts[i] = 2.0*pi/phipts.size();
 }
 
 if (printout) {
   double tot = 0.0;
   for (int i = 0; i < phipts.size(); i++) {
     std::cout << "i = " << i << std::endl;
     std::cout << "phipts[" << i << "] = " << phipts[i] << std::endl;
     std::cout << "phiwts[" << i << "] = " << phiwts[i] << std::endl;
     tot += phiwts[i];
   }
   std::cout << "tot = " << tot << std::endl;
 }


 //  std::cout << "tot = " << tot << std::endl;

 //print points to file
 // printpts(thetapts,thetawts,phipts,phiwts,"ray_pts");
 return ;
 
}



void logpts(VecDoub &rpts, VecDoub &rwts, const bool printout, atmointerp &thisatmointerp)
// gets radial points that are evenly spaced in log radius,
// and prints the spacing of the points
{
  //how many rpts are there?
  int nrpts = rpts.size();

  //log spacing throughout atmosphere
  double logmax = log(rmax-rMars);
  double logmin = log(rminatm-rMars);
  double logspace = (logmax-logmin)/((double) nrpts-1);
  
  //now get the rpoints:
  for (int i = 0; i < nrpts; i++) {
    rpts[i]=exp(logmax-(i)*logspace)+rMars;
    std::cout << "rpts["<<i<<"] = "<<(rpts[i]-rMars)/1e5<<"\n";
  }
  
  //if we've got the points, we can compute the coefficients, assuming
  //simple trapezoidal quadrature:
  rwts[0] = (rpts[0]-rpts[1])/2.0;
  rwts[nrpts-1] = (rpts[nrpts-2]-rpts[nrpts-1])/2.0;
  for (int i = 1; i < nrpts-1; i++) 
    rwts[i] = (rpts[i-1]-rpts[i+1])/2.0;
  for (int i = 0; i < nrpts; i++)
    rwts[i] /= (rmax-rminatm);

  if (printout) {

    //just for funsies let's see what the optical depth looks like at each point:
    dtau_tot_int dtau(thisatmointerp);
    double tautot = qsimp(dtau,rminatm,2.0*rmax,1.0e-4);
    if (printout) std::cout << "total tau through atmosphere is "
			    << tautot << std::endl;
    double tau[nrpts];
    tau[0] = 0.0;
    for (int i=1; i < nrpts; i++) {
      tau[i] = qsimp(dtau,rpts[i],2.0*rmax,1.0e-4);
      //   std::cout << "tau[" << i << "] = " << tau[i] << std::endl;
    }

    //print the altitudes being used:
    std::cout << "These are the altitudes in use for this run:" << std::endl;
    std::cout.width(15);
    std::cout<<"alt(km)";
    std::cout.width(15);
    std::cout<<"nH cm^-3";
    std::cout.width(15);
    std::cout<<"tau";
    std::cout.width(15);
    std::cout<<"T (K)";
    std::cout.width(15);
    std::cout<<"nC02 cm^-3";
    std::cout<<std::endl;
    for (int i = 0; i < nrpts; i++) {
      std::cout.width(15);
      std::cout << (rpts[i]-rMars)/1e5;
      std::cout.width(15);
      std::cout << thisatmointerp.nH(rpts[i]);
      std::cout.width(15);
      std::cout << tau[i];
      std::cout.width(15);
      std::cout << Temp(rpts[i],thisatmointerp.Texo);
      std::cout.width(15);
      std::cout << thisatmointerp.nCO2(rpts[i]);
      std::cout << std::endl;
    }
    //    std::cin.get();
  }
}

void loglinpts(VecDoub &rpts, VecDoub &rwts, const bool printout,atmointerp &thisatmointerp)
/* gets radial points that are split, with half linearly spaced below
   the exobase and half logarithmically spaced above. */
{
  //how many rpts are there?
  int nrpts = rpts.size();

  //need to put 1/2 of radial points below rexo;
  int nbelowrexo=nrpts/2;

  //log spacing above rexo
  double logmax = log(rmax-rMars);
  double logmin = log(rexo-rMars);
  double logspace = (logmax-logmin)/((double) nrpts-nbelowrexo);
  
  //now get the rpoints:
  for (int i = 0; i < nrpts-nbelowrexo; i++) {
    rpts[i]=exp(logmax-(i)*logspace)+rMars;
    std::cout << "rpts["<<i<<"] = "<<(rpts[i]-rMars)/1e5<<"\n";
  }

  //space linearly below the exobase
  double linspace = (rexo-rminatm)/((double) nbelowrexo-1);
  for (int i = nrpts-nbelowrexo; i<nrpts; i++) {
    rpts[i]=rexo-(i-(nrpts-nbelowrexo))*linspace;
    std::cout << "rpts["<<i<<"] = "<<(rpts[i]-rMars)/1e5<<"\n";
  }

  //if we've got the points, we can compute the coefficients, assuming
  //simple trapezoidal quadrature:
  rwts[0] = (rpts[0]-rpts[1])/2.0;
  rwts[nrpts-1] = (rpts[nrpts-2]-rpts[nrpts-1])/2.0;
  for (int i = 1; i < nrpts-1; i++) 
    rwts[i] = (rpts[i-1]-rpts[i+1])/2.0;
  for (int i = 0; i < nrpts; i++)
    rwts[i] /= (rmax-rminatm);

  if (printout) {
    //just for funsies let's see what the optical depth looks like at each point:
    dtau_tot_int dtau(thisatmointerp);
    double tautot = qsimp(dtau,rminatm,2.0*rmax,1.0e-4);
    if (printout) std::cout << "total tau through atmosphere is "
			    << tautot << std::endl;
    double tau[nrpts];
    tau[0] = 0.0;
    for (int i=1; i < nrpts; i++) {
      tau[i] = qsimp(dtau,rpts[i],2.0*rmax,1.0e-4);
      //   std::cout << "tau[" << i << "] = " << tau[i] << std::endl;
    }

    //print the altitudes being used:
    std::cout << "These are the altitudes in use for this run:" << std::endl;
    std::cout.width(15);
    std::cout<<"alt(km)";
    std::cout.width(15);
    std::cout<<"nH cm^-3";
    std::cout.width(15);
    std::cout<<"tau";
    std::cout.width(15);
    std::cout<<"T (K)";
    std::cout.width(15);
    std::cout<<"nC02 cm^-3";
    std::cout<<std::endl;
    for (int i = 0; i < nrpts; i++) {
      std::cout.width(15);
      std::cout << (rpts[i]-rMars)/1e5;
      std::cout.width(15);
      std::cout << thisatmointerp.nH(rpts[i]);
      std::cout.width(15);
      std::cout << tau[i];
      std::cout.width(15);
      std::cout << Temp(rpts[i],thisatmointerp.Texo);
      std::cout.width(15);
      std::cout << thisatmointerp.nCO2(rpts[i]);
      std::cout << std::endl;
    }
    //    std::cin.get();
  }
}

struct auxtau {
  //auxiliary function to compute the local importance of a point for
  //the radiative transfer, T(tau_H)*exp(-tau_CO2)
  dtau_H_int * dtauH;
  dtau_CO2_int * dtauCO2;
  auxtau(dtau_H_int * dtauHH, dtau_CO2_int * dtauCO22) : dtauH(dtauHH), dtauCO2(dtauCO22) {};
  double operator()(double & r) {
    return 1.0-HolT((*dtauH)(r));
  }
};


struct taufind
// auxiliary function for finding radial gridpoints
{
  double offset;
  auxtau * myauxtau;
  dtau_CO2_int * dtauCO2;
  taufind(double & offsett, 
	  auxtau * myauxtauu, 
	  dtau_CO2_int * dtauCO22) : offset(offsett), 
				     myauxtau(myauxtauu), 
				     dtauCO2(dtauCO22) {};
  double operator()(double & r) {
    return qsimp((*myauxtau),rminatm,r,1.0e-4)*(1-exp(-qsimp((*dtauCO2),rminatm,r,1.0e-4))) - offset;
  }
};

void taufracpts(VecDoub &rpts, VecDoub &rwts, const bool printout,atmointerp &thisatmointerp)
// gets radial points that are evenly spaced in RT importance,
// and prints some properties of the points
{
  //how many rpts are there?
  int nrpts = rpts.size();

  //exobase temperature is stored in interpolation object:
  const double Texo = thisatmointerp.Texo;

  //get total optical depth from infinity to rmin
  dtau_H_int dtauH(thisatmointerp);
  dtau_CO2_int dtauCO2(thisatmointerp);
  auxtau myauxtau(&dtauH,&dtauCO2);
  double tau = qsimp(myauxtau,rminatm,rmax,1.0e-4)*(1-exp(-qsimp(dtauCO2,rminatm,rmax,1.0e-4)));
  if (printout) std::cout << "total Int(T(tau_H)*exp(-tau_CO2)) through atmosphere is " 
			  << tau << std::endl;
  double taustep = tau/(nrpts-1);
  if (printout) std::cout << "spacing between radial points is " << taustep << std::endl;
  
  //set up endpoints
  rpts[0] = rmax;
  rpts[nrpts-1] = rminatm;

  //now get the rpoints moving inward:
  double ttau = tau;
  for (int i = 1; i <nrpts-1; i++) {
    std::cout << "finding point " << i << std::endl;
    ttau -= taustep;
    taufind thistau(ttau, &myauxtau, &dtauCO2);
    rpts[i] = zbrent(thistau,rminatm,rpts[i-1],1e-4);
    std::cout << "rpts[" << i << "] = " << (rpts[i]-rMars)/1e5 << std::endl;
  }

  //if we've got the points, we can compute the coefficients, assuming
  //simple trapezoidal quadrature:
  rwts[0] = (rpts[0]-rpts[1])/2.0;
  rwts[nrpts-1] = (rpts[nrpts-2]-rpts[nrpts-1])/2.0;
  for (int i = 1; i < nrpts-1; i++) 
    rwts[i] = (rpts[i-1]-rpts[i+1])/2.0;
  for (int i = 0; i < nrpts; i++)
    rwts[i] /= (rmax-rminatm);
  if (printout) {
    //print the altitudes being used:
    std::cout << "These are the altitudes in use for this run:" << std::endl;
    std::cout.width(15);
    std::cout<<"alt(km)";
    std::cout.width(15);
    std::cout<<"nH cm^-3";
    std::cout.width(15);
    std::cout<<"tau";
    std::cout.width(15);
    std::cout<<"T (K)";
    std::cout.width(15);
    std::cout<<"nC02 cm^-3";
    std::cout<<std::endl;
    ttau = 0;
    for (int i = 0; i < nrpts; i++) {
      std::cout.width(15);
      std::cout << (rpts[i]-rMars)/1e5;
      std::cout.width(15);
      std::cout << thisatmointerp.nH(rpts[i]);
      std::cout.width(15);
      std::cout << ttau;
      ttau += taustep;
      std::cout.width(15);
      std::cout << Temp(rpts[i],Texo);
      std::cout.width(15);
      std::cout << thisatmointerp.nCO2(rpts[i]);
      std::cout << std::endl;
    }
    //    std::cin.get();
  }
}

void doubletermtpts(VecDoub &thetapts, 
		    VecDoub &thetawts, 
		    const bool printout)
// gets theta points that are spaced twice as finely in the ~60
// degrees near the terminator as they are in the remainder of the
// space
{
  double termfrac=pi/6.;//fraction of SZA considered "near" the
                             //terminator

  int ntpts=thetapts.size();
  std::cout << "ntpts = " << ntpts << std::endl;
  if (ntpts<5||(ntpts-5)%4!=0) { 
    throw("for method doubletermpts, nthetapoints must be of the form 4N+5 for N>=0.");
  }
  thetapts[0]=0;
  thetapts[(ntpts-1)/2]=pi/2.;
  thetapts[ntpts-1]=pi;
  int nmidpts=(ntpts-5)/4;
  thetapts[nmidpts+1]=pi/2-termfrac/2;
  thetapts[ntpts-1-(nmidpts+1)]=pi/2+termfrac/2;
  for (int i = 1; i <= nmidpts; i++) {
    thetapts[i]=thetapts[nmidpts+1]/(nmidpts+1)*i;
    thetapts[ntpts-1-nmidpts-1+i]=pi/2+termfrac/2+thetapts[i];
    thetapts[nmidpts+1+i]=pi/2-termfrac/2+termfrac/2/(nmidpts+1)*i;
    thetapts[ntpts-1-2*(nmidpts+1)+i]=pi/2+termfrac/2/(nmidpts+1)*i;
  } 

  //now define the weights
  thetawts[0]=thetapts[1]/2;
  thetawts[ntpts-1]=thetapts[1]/2;
  for (int i = 1; i < ntpts-1; i++) {
    thetawts[i]=(thetapts[i+1]-thetapts[i-1])/2.;
  }

  if (printout) {
    double tot = 0.0;
    for (int i = 0; i < thetapts.size(); i++) {
      std::cout << "i = " << i << std::endl;
      std::cout << "thetapts[" << i << "] = " << thetapts[i] 
		<< " ("<< 180./pi*thetapts[i] << ")"<< std::endl;
      std::cout << "thetawts[" << i << "] = " << thetawts[i] << std::endl;
      tot += thetawts[i];
    }
    std::cout << "tot = " << tot << std::endl;
  }
  //  std::cin.get();

}



void getquadpts(VecDoub &rpts, VecDoub &rwts, 
		VecDoub &thetapts, VecDoub &thetawts,
		VecDoub &phipts, VecDoub &phiwts, const bool printout,
		atmointerp &thisatmointerp, 
		const string rmethod="loglinpts",const string tmethod="uniform")
//gets quadrature points and weights and stores them in the passed
//vector objects
{

  if (rmethod=="taupts") {
    // get radial points that are uniformly spaced in optical depth:
    taufracpts(rpts,rwts,printout,thisatmointerp);      
  } else if (rmethod=="logpts") {
    // use points with equal logarithmic spacing:
    logpts(rpts,rwts,printout,thisatmointerp);
  } else if (rmethod=="loglinpts") {
    // use points with split linear/logarithmic spacing:
    loglinpts(rpts,rwts,printout,thisatmointerp);
  } else {
    string msg="No method "+rmethod+" in getquadpts!";
    throw(msg.c_str());
  }
  
  if (tmethod=="doubleterm") {
    doubletermtpts(thetapts,thetawts,printout);
  } else if (tmethod=="uniform") {
    //spacing in SZA is uniform (trapezoidal quadrature):
    for (int i = 0; i < thetapts.size(); i++) {
      thetapts[i] = i*pi/(thetapts.size()-1);
      thetawts[i] = pi/(thetapts.size()-1);
    }
    thetawts[0] /= 2.0;
    thetawts[thetapts.size()-1] /= 2.0;
  } else {
    string msg="No method "+tmethod+" in getquadpts!";
    throw(msg.c_str());
  }

  //  double tot = 0.0;
  //  for (int i = 0; i < thetapts.size(); i++) {
  //    std::cout << "i = " << i << std::endl;
  //    std::cout << "thetapts[" << i << "] = " << thetapts[i] << std::endl;
  //    std::cout << "thetawts[" << i << "] = " << thetawts[i] << std::endl;
  //    tot += thetawts[i];
  //  }
  //  std::cout << "tot = " << tot << std::endl;
  
  
  //for phi, trapezoidal quadrature is sufficient since the problem is
  //symmetric in this coordinate
  //  tot = 0.0;
  for (int i = 0; i < phipts.size(); i++) {
    phipts[i] = i*2.0*pi/(phipts.size()-1);
    phiwts[i] = 2.0*pi/(phipts.size()-1);
    //    tot += phiwts[i];
  }
  //  std::cout << "tot = " << tot << std::endl;

  //print points to file
  //  printpts(rpts,rwts,thetapts,thetawts,phipts,phiwts,"pts");

  return ;
  
}

  
#endif
