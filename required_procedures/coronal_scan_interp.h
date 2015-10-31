//coronal_scan_interp.h -- routines to compute coronal scan
//                         intensities from an interpolated grid of
//                         source functions, computed on demand.

#ifndef __CORONAL_SCAN_INTERP_H

#include <iostream>
#include <fstream>
#include "nr3.h"
#include "definitions.h" //grid definitions
#include <cmath>
#include "corona_simulator.h" //simulation routines for synthetic atmospheres

using std::log;
using std::exp;

struct coronal_scan_interp
/* this object provides utilities for fitting observations. It
   interpolates intensities from generated source functions on a
   specified grid in temperature and density. Interpolation is linear
   in each dimension. This object will also return the derivative in
   each dimension. */
{
  int nobs, nnH, nT;
  VecDoub nH_vec, T_vec;
  Linear_interp nH_terp;
  Linear_interp T_terp;

  string obsname;
  bool forcesim;
  double *iIcalc;
 
  corona_simulator **corona_simulator_grid;

  //constructor takes the name of the observation to be loaded as an argument
  coronal_scan_interp(string obsnamee,bool forcesimm=FALSE): 
    obsname(obsnamee), forcesim(forcesimm) {
    //now set up the number density and temperature grid, based on
    //numbers in definitions.h
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

    //set up the coronal simulator grid
    corona_simulator_grid = new corona_simulator*[nnH];
    for (int inH=0;inH<nnH;inH++)
      corona_simulator_grid[nnH] = new corona_simulator[nT];
  }


  ~coronal_scan_interp() {
    for (int inH = 0; inH < nnH; inH++)
      delete [] corona_simulator_grid[nnH];
    delete [] corona_simulator_grid;
  }
  
  //check that intensities exist for interpolation, fill them in if they don't
  void check_intensities(int inH, int iT)
  {
    for (int iinH=inH;iinH<=inH+1;iinH++) {
      for (int iiT=iT;iiT<=iT+1;iiT++) {
	if (!corona_simulator_grid[iinH][iiT].obsinit) {
	  //load the observation
	  corona_simulator_grid[iinH][iiT].obs_import(obsname);
	  //load or simulate the source function
	  corona_simulator_grid[iinH][iiT].get_S(nH_vec[iinH],T_vec[iiT],forcesim);
	  //calculate the intensities
	  corona_simulator_grid[iinH][iiT].calc_I();
	  //initialize or reset the number of observation points
	  nobs=corona_simulator_grid[iinH][iiT].nobs;
	}
      }
    }
  }
  
  double interp(int iobs, double nHp, double Tp, double IPHb=0.0) {
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

    check_intensities(inH,iT);//check that these intensities have been simulated
    
    //interpolate:
    anH = (nHp-(nH_terp).xx[inH])/((nH_terp).xx[inH+1]-(nH_terp).xx[inH]);
    aT  = ( Tp-( T_terp).xx[ iT])/(( T_terp).xx[ iT+1]-( T_terp).xx[ iT]);

    /* std::cout << "aobs = " << aobs << std::endl; */
    /* std::cout << "anH = " << anH << std::endl; */
    /* std::cout << "aT = " << aT << std::endl;     */
    
    iIcalc = (1.-anH)*(1.-aT)*corona_simulator_grid[inH  ][iT  ].I_calc[iobs]
            +    anH *(1.-aT)*corona_simulator_grid[inH+1][iT  ].I_calc[iobs]
            +(1.-anH)*    aT *corona_simulator_grid[inH  ][iT+1].I_calc[iobs]
            +    anH*     aT *corona_simulator_grid[inH+1][iT+1].I_calc[iobs];

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

  // double IPHtrans_interp(int iobs, double nHp, double Tp) {
  //   double iIPHtrans, anH, aT;
  //   //find the grid square:
  //   /* std::cout << "nH = " << nHp << std::endl; */
  //   int inH = (nH_terp).cor ? (nH_terp).hunt(nHp) : (nH_terp).locate(nHp);
  //   // std::cout << "inH = " << inH << std::endl;
  //   // std::cout << "nH[" << inH << "] = " << (*nH_terp).xx[inH] << std::endl;
  //   // std::cout << "nH[" << inH + 1 << "] = " << (*nH_terp).xx[inH + 1] << std::endl;
  //   // std::cout << "T = " << Tp << std::endl;
  //   int iT = (T_terp).cor ? (T_terp).hunt(Tp) : (T_terp).locate(Tp);
  //   // std::cout << "iT = " << iT << std::endl;
  //   // std::cout << "T[" << iT << "] = " << (*T_terp).xx[iT] << std::endl;
  //   // std::cout << "T[" << iT + 1 << "] = " << (*T_terp).xx[iT + 1] << std::endl;
    
  //   // std::cout << "inH = " << inH << std::endl;
  //   // std::cout << "iT = " << iT << std::endl;    

  //   check_intensities(inH,iT);//check that these intensities have been simulated
    
  //   //interpolate:
  //   anH = (nHp-(nH_terp).xx[inH])/((nH_terp).xx[inH+1]-(nH_terp).xx[inH]);
  //   aT = (Tp-(T_terp).xx[iT])/((T_terp).xx[iT+1]-(T_terp).xx[iT]);

  //   /* std::cout << "anH = " << anH << std::endl; */
  //   /* std::cout << "aT = " << aT << std::endl;     */
    
  //   //now add in the IPH for points with tangent altitudes above 200km
  //   //if iph is simulated using the Quemerais code, it is multiplied
  //   //in here. Otherwise, this factor was set to 1 in the initial load
  //   //of the obs_data and IPHb represents the brightness of the IPH
  //   //and not a multiplicative scaling factor.
  //   iIPHtrans =  (1.-anH)*(1.-aT)*corona_simulator_grid[inH  ][iT  ].IPHtrans[iobs]
  //              +     anH *(1.-aT)*corona_simulator_grid[inH+1][iT  ].IPHtrans[iobs]
  //              + (1.-anH)*    aT *corona_simulator_grid[inH  ][iT+1].IPHtrans[iobs]
  //              +     anH *    aT *corona_simulator_grid[inH+1][iT+1].IPHtrans[iobs]
    
  //   return iIPHtrans;
  // }

  void derivs(int iobs, double nHp, double Tp, double IPHb, VecDoub_O &dIcalc) {
    //find the grid square:
    int inH = (nH_terp).cor ? (nH_terp).hunt(nHp) : (nH_terp).locate(nHp);
    int iT = (T_terp).cor ? (T_terp).hunt(Tp) : (T_terp).locate(Tp);

    check_intensities(inH,iT);//check that these intensities have been simulated
    
    double nH1=(nH_terp).xx[inH];
    double nH2=(nH_terp).xx[inH+1];
    double T1=(T_terp).xx[iT];
    double T2=(T_terp).xx[iT+1];

    double fnH=(nHp-nH1)/(nH2-nH1);
    double fT=(Tp-T1)/(T2-T1);

    double I11=corona_simulator_grid[inH  ][iT  ].I_calc[iobs];
    double I12=corona_simulator_grid[inH  ][iT+1].I_calc[iobs];
    double I21=corona_simulator_grid[inH+1][iT  ].I_calc[iobs];
    double I22=corona_simulator_grid[inH+1][iT+1].I_calc[iobs];

    // double trans11=corona_simulator_grid[inH  ][iT  ].IPHtrans[iobs];
    // double trans12=corona_simulator_grid[inH  ][iT+1].IPHtrans[iobs];
    // double trans21=corona_simulator_grid[inH+1][iT  ].IPHtrans[iobs];
    // double trans22=corona_simulator_grid[inH+1][iT+1].IPHtrans[iobs];

    //we must bilinearly interpolate the derivatives
    //intensity
    double dIdnH, dIdT;

    dIdnH =  (1.-fT)*(I21-I11);
    dIdnH +=     fT *(I22-I12);
    dIdnH /= (nH2-nH1);
    dIcalc[0]=dIdnH;
    
    dIdT =  (1.-fnH)*(I12-I11);
    dIdT +=     fnH *(I22-I21);
    dIdT /= (T2-T1);
    dIcalc[1]=dIdT;

    //transmission
    // double dtransdnH, dtransdT;
    // dtransdnH =  (1.-fT)*(trans21-trans11);
    // dtransdnH +=     fT *(trans22-trans12);
    // dtransdnH /= (nH2-nH1);
    // dIcalc[0]+=IPHb*dtransdnH;

    // dtransdT =  (1.-fnH)*(trans12-trans11);
    // dtransdT +=     fnH *(trans22-trans21);
    // dtransdT /= (T2-T1);
    // dIcalc[1]+=IPHb*dtransdT;
  }

  //this is used for the chi-square minimzation by L-M
  void operator()(const int iobs, VecDoub_I &a, double &iIcalc, VecDoub_O &dIcalcda)
  {
    double nHp = a[0];
    double Tp = a[1];
    double IPHb = a[2];
    double cal = a[3];
    
    double tempcalc = interp(iobs, nHp, Tp, IPHb);
    iIcalc = tempcalc*cal;

    VecDoub dItemp(2);
    derivs(iobs, nHp, Tp, IPHb, dItemp);
    dIcalcda[0]=cal*dItemp[0];
    dIcalcda[1]=cal*dItemp[1];
    //the derivative wrt the IPH parameter is computed as the product
    //of the transmission and the iphb_calc parameter. If IPHb is a
    //multiplicative factor on the simulated brightness, this
    //multplies in the brightness simulated. If IPHb is a scalar
    //background, iphb_calc is set to 1 in the initial load of the
    //obs_data.
    //    dIcalcda[2]=cal*IPHtrans_interp(iobs, nHp, Tp);
    dIcalcda[3]=tempcalc;
  }
  
};

#endif
