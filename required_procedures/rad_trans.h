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


inline double sH(const double &T) { return 5.96e-12/sqrt(T); } //effective H cross section at Ly alpha line center

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
 }
 
 //for phi, trapezoidal quadrature
 //  tot = 0.0;
 for (int i = 0; i < phipts.size(); i++) {
   phipts[i] = i*2.0*pi/phipts.size();
   phiwts[i] = 2.0*pi/phipts.size();
   //    tot += phiwts[i];
 }
 //  std::cout << "tot = " << tot << std::endl;

 //print points to file
 // printpts(thetapts,thetawts,phipts,phiwts,"ray_pts");
 return ;
 
}


struct taufind
// auxiliary function for finding radial gridpoints
{
  double offset;
  dtau_tot_int * dtau;
  taufind(double & offsett, dtau_tot_int * dtauu) : offset(offsett), dtau(dtauu) {};
  double operator()(double & x) {
    return qsimp((*dtau),rminatm,x,1.0e-4) - offset;
  }
};

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


void taufracpts(VecDoub &rpts, VecDoub &rwts, const bool printout,atmointerp &thisatmointerp)
// gets radial points that are evenly spaced in optical depth,
// and prints the optical depth spacing of the points
{
  //how many rpts are there?
  int nrpts = rpts.size();

  //exobase temperature is stored in interpolation object:
  const double Texo = thisatmointerp.Texo;

  //get total optical depth from infinity to rmin
  dtau_tot_int dtau(thisatmointerp);
  double tau = qsimp(dtau,rminatm,rmax,1.0e-4);
  if (printout) std::cout << "total tau through atmosphere is " << tau << std::endl;
  double taustep = tau/(nrpts-1);
  if (printout) std::cout << "spacing in tau between radial points is " << taustep << std::endl;
  
  //set up endpoints
  rpts[0] = rmax;
  rpts[nrpts-1] = rminatm;

  //now get the rpoints moving inward:
  double ttau = tau;
  for (int i = 1; i < nrpts - 1; i++) {
    //    std::cout << "finding point " << i << std::endl;
    ttau -= taustep;
    taufind thistau(ttau, &dtau);
    rpts[i] = zbrent(thistau,rminatm,rpts[i-1],1e-4);
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
  }
}

void chaurpts(VecDoub &rpts, VecDoub &rwts, const bool printout,atmointerp &thisatmointerp)
// gets radial points that are evenly spaced in optical depth,
// and prints the optical depth spacing of the points
{
  std::cout << "Using chaufray's tabulated values for radial points" << std::endl;

  int nrpts = rpts.size();
  
  // Chaufray provides 74 r points. By interpolation (see
  // ./chaufrays_code/get_chau_r_pts.nb), I have created files with 2x
  // and 4x this number of points. Print an error if the number of
  // points is not one of these three numbers.
  string fname; 
  switch (nrpts)
    {
    case 74:
      fname = "chau_r_pts.dat";
      break;
    case 147:
      fname = "chau_r_pts_x2.dat";
      break;
    case 293:
      fname = "chau_r_pts_x4.dat";
    default:
      throw("Using chaurpts doesn't work unless you use 74, 147, or 293 points in r!");      
    }

  //set up the file input to read the points in:
  ifstream chaufile;
  chaufile.open(fname.c_str());
      
  //now get the rpoints from the file
  for (int i = 0; i < nrpts; i++) {
    chaufile >> rpts[i];
    //    std::cout << "rpts[" << i << "] = " << rpts[i] << std::endl;
    //convert from km in alt to cm in radius:
    rpts[i] = rpts[i]*1e5 + rMars;
    //    std::cout << "rpts[" << i << "] = " << rpts[i] << std::endl;
  }

  //if we've got the points, we can compute the coefficients, assuming
  //simple trapezoidal quadrature:
  rwts[0] = (rpts[0]-rpts[1])/2.0;
  rwts[nrpts-1] = (rpts[nrpts-2]-rpts[nrpts-1])/2.0;
  for (int i = 1; i < nrpts-1; i++) 
    rwts[i] = (rpts[i-1]-rpts[i+1])/2.0;
  for (int i = 0; i < nrpts; i++)
    rwts[i] /= (rmax-rminatm);

  //exobase temperature is stored in interpolation object:
  const double Texo = thisatmointerp.Texo;

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
  

  if (printout) {
    //print the altitudes being used:
    std::cout << "These are the altitudes in use for this run:" << std::endl;
    std::cout.width(15);
    std::cout<<"radius(cm)";
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
      std::cout << rpts[i];
      std::cout.width(15);
      std::cout << (rpts[i]-rMars)/1e5;
      std::cout.width(15);
      std::cout << thisatmointerp.nH(rpts[i]);
      std::cout.width(15);
      std::cout << tau[i];
      std::cout.width(15);
      std::cout << Temp(rpts[i],Texo);
      std::cout.width(15);
      std::cout << thisatmointerp.nCO2(rpts[i]);
      std::cout << std::endl;
    }
  }
}

// void getfname(ofstream &file, const char *tag,const int nrpts,const int ntpts, const int nppts)
// {
//   //file extension
//   char ext[50];
//   sprintf(ext,"%dx%dx%d.dat",nrpts,ntpts,nppts);

//   //get filename, open file
//   //std::cout << "Input a filename for the "<< tag <<" matrix (or Enter to use "<< tag << ext << "): ";
//   char name[50];
//   //std::cin.getline(name,50);
//   //if (!strcmp(name,""))
//     sprintf(name,"%s%s",tag,ext);
//   file.open(name);
//   file.precision(10);
//   file.setf(ios_base::scientific);
// }

// void printpts(const VecDoub &rpts, const VecDoub &rwts, 
// 	      const VecDoub &tpts, const VecDoub &twts,
// 	      const VecDoub &ppts, const VecDoub &pwts, const char* tag)
// {
//   int nrpts = rpts.size();
//   int ntpts = tpts.size();
//   int nppts = ppts.size();
//   ofstream ptsfile;
//   getfname(ptsfile,tag,nrpts,ntpts,nppts);

//   std::cout << "Points used for RT in atm are being saved to: "
// 	    << tag << nrpts << "x" << ntpts << "x" << nppts << ".dat" << std::endl;
  
//   ptsfile << "This file contains the radial and angular points and weights "
// 	  << "used for the calculation specified in the filename.\n\n";

//   ptsfile << "Radial points and weights: (in km)\n\n";
//   for (int i = 0; i < nrpts; i++) {
//     ptsfile.width(20);
//     ptsfile << (rpts[i]-rMars)/1e5;
//     ptsfile.width(20);
//     ptsfile << rwts[i];
//     ptsfile << std::endl;
//   }

//   ptsfile << "SZA points and weights: (in radians)\n\n";
//   for (int i = 0; i < ntpts; i++) {
//     ptsfile.width(20);
//     ptsfile << tpts[i];
//     ptsfile.width(20);
//     ptsfile << twts[i];
//     ptsfile << std::endl;
//   }

//   ptsfile << "Phi points and weights: (in radians)\n\n";
//   for (int i = 0; i < nppts; i++) {
//     ptsfile.width(20);
//     ptsfile << ppts[i];
//     ptsfile.width(20);
//     ptsfile << pwts[i];
//     ptsfile << std::endl;
//   }

//   ptsfile.close();
//   return ;
// }

void getquadpts(VecDoub &rpts, VecDoub &rwts, 
		VecDoub &thetapts, VecDoub &thetawts,
		VecDoub &phipts, VecDoub &phiwts, const bool printout,
		atmointerp &thisatmointerp, const string rmethod="logpts")
//gets quadrature points and weights and stores them in the passed
//vector objects
{

  if (rmethod=="taupts") {
    // get radial points that are uniformly spaced in optical depth:
    taufracpts(rpts,rwts,printout,thisatmointerp);      
  } else if (rmethod=="chaupts") {
    // use chaufray's tabulated points
    chaurpts(rpts,rwts,printout,thisatmointerp);
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
  
  //spacing in SZA is uniform (trapezoidal quadrature):
  for (int i = 0; i < thetapts.size(); i++) {
    thetapts[i] = i*pi/(thetapts.size()-1);
    thetawts[i] = pi/(thetapts.size()-1);
  }
  thetawts[0] /= 2.0;
  thetawts[thetapts.size()-1] /= 2.0;

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
