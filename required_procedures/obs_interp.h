/* obs_interp.h -- structure to read in and interpolate over simulated
   observation data for a given observation. Simulations of the source
   function are integrated on demand as required. */

#ifndef __OBS_INTERP_H
#define __OBS_INTERP_H

#include "nr3.h"
#include "interp_1d.h"
#include <iostream> // for file output and dialog
#include <fstream>  // ""
#include "obs_sim.h"

struct obs_interp
/* this object provides utilities for fitting observations. It
   interpolates in temperature and density, and also interpolates in
   observation number (though this is mostly as a
   convenience). Interpolation is linear in each dimension. This
   object will also return the derivative in each dimension. */
{
  obs_simulator *obs_sim;
  string obsname;
  string fname;
  fstream file;
  int iobs, inH, iT;
  double ***Icalc;
  int nobs, nnH, nT;
  double *obs_list, *nH_list, *T_list;
  VecDoub obs_vec, nH_vec, T_vec;
  Linear_interp *obs_terp;
  Linear_interp *nH_terp;
  Linear_interp *T_terp;

  //constructor takes the name of the observation to be loaded as an argument
  obs_interp(obs_simulator *obs_simm): obs_sim(obs_simm) {
    //first, load up all the coordinates information from the coordinates file
    obsname = (*obs_sim).thisobs.obsname;
    fname = obsname + "_sim_coords.dat";
    std::cout << "Reading computed intensity geometrical data from " << fname << std::endl;
    file.open(fname.c_str());//this assumes that the code is executing 
    //from the directory containing this file
    //now read in the observation number, densities, and temperatures
    //OBS NUMBER
    string dum; //dummy string to skip comment lines
    std::getline(file,dum);
    //     std::cout << dum;
    file >> nobs;
    //     std::cout << "nobs = " << nobs << std::endl;
    obs_list = new double[nobs];
    for (iobs = 0; iobs < nobs; iobs++) {
      obs_list[iobs]=iobs;
      file >> dum; // need to stuff file input somewhere
      //      std::cout << "obs_list[" << iobs << "] = " << obs_list[iobs] << std::endl;
      /* std::cin.get(); */
    }
    //DENSITIES
    std::getline(file,dum);//clear the carraige return
    std::getline(file,dum);
    //     std::cout << dum;
    file >> nnH;
    //     std::cout << "nnH = " << nnH << std::endl;
    nH_list = new double[nnH];
    for (inH = 0; inH < nnH; inH++) {
      file >> nH_list[inH];
      /* std::cout << "nH_list[" << inH << "] = " << nH_list[inH] << std::endl; */
      /* std::cin.get(); */
    }
    //TEMPERATURES
    std::getline(file,dum);//clear the carraige return
    std::getline(file,dum);
    //
    //     std::cout << dum;
    file >> nT;
    //     std::cout << "nT = " << nT << std::endl;
    T_list = new double[nT];
    for (iT = 0; iT < nT; iT++) {
      file >> T_list[iT];
      /* std::cout << "T_list[" << iT << "] = " << T_list[iT] << std::endl; */
      /* std::cin.get(); */
    }
    file.close();
     
    //now initialize the intensities:
    Icalc = new double**[nobs];
    for (iobs = 0; iobs < nobs; iobs++) {
      Icalc[iobs] = new double*[nnH];
      for (inH = 0; inH < nnH; inH++) {
	Icalc[iobs][inH] = new double[nT];
	for (iT = 0; iT < nT; iT++) {
	  Icalc[iobs][inH][iT] = -1;
	}
      }
    }

    //load the coordinates into interpolation objects for coord finding
    //OBS NO
    obs_vec.resize(nobs);
    for (int i = 0; i < nobs; i++) obs_vec[i]=obs_list[i];
    obs_terp = new Linear_interp(obs_vec,obs_vec);
    // double obsp;
    // std::cout << "Enter an obsitude: ";
    // while (std::cin >> obsp) {
    //   iobs = (*obs_terp).cor ? (*obs_terp).hunt(obsp) : (*obs_terp).locate(obsp);
    //   std::cout << "Obs grid is " << iobs << std::endl;
    //   std::cout << "Enter an obs number: ";
    // }
    //DENSITY
    nH_vec.resize(nnH);
    for (int i = 0; i < nnH; i++) nH_vec[i]=nH_list[i];
    nH_terp = new Linear_interp(nH_vec,nH_vec); 
    // double nHp;
    // std::cout << "Enter a density: ";
    // while (std::cin >> nHp) {
    //   inH = (*nH_terp).cor ? (*nH_terp).hunt(nHp) : (*nH_terp).locate(nHp);
    //   std::cout << "nH grid is " << inH << std::endl;
    //   std::cout << "Enter a density: ";
    // }
    //TEMPERATURE
    T_vec.resize(nT);
    for (int i = 0; i < nT; i++) T_vec[i]=T_list[i];
    T_terp = new Linear_interp(T_vec,T_vec); 
    //     std::cout << "T[5] = " << (*T_terp).xx[5] << std::endl;
    // double Tp;
    // std::cout << "Enter a temperature: ";
    // while (std::cin >> Tp) {
    //   iT = (*T_terp).cor ? (*T_terp).hunt(Tp) : (*T_terp).locate(Tp);
    //   std::cout << "T grid is " << iT << std::endl;
    //   std::cout << "Enter a temperature: ";
    // }
  }


  ~obs_interp() {
    //destructor frees memory allocated to hold interpolation objects
    delete [] obs_list;
    delete [] nH_list;
    delete [] T_list;

    for (int iobs = 0; iobs < nobs; iobs++) {
      for (int inH = 0; inH < nnH; inH++) {
	delete [] Icalc[iobs][inH];
      }
      delete [] Icalc[iobs];
    }
    delete [] Icalc;

    delete obs_terp;
    delete nH_terp;
    delete T_terp;
  }

  //check that intensities exist for interpolation, fill them in if they don't
  void check_intensities(int iobs, int inH, int iT, double IPHb)
  {
     for (int iiobs = iobs; iiobs <= iobs+1; iiobs++) {
       for (int iinH = inH; iinH <= inH+1; iinH++) {
	 for (int iiT = iT; iiT <= iT+1; iiT++) {
	   if(Icalc[iiobs][iinH][iiT] == -1)//intensity was initialized but not computed
	     {
	       Icalc[iiobs][iinH][iiT] = (*obs_sim).simulate_gridpoint(iiobs,iinH,iiT,IPHb);
	     };
	 }
       }
     }
  }

  double interp(double obsp, double nHp, double Tp, double IPHb=0.0) {
    double iIcalc, aobs, anH, aT;
    //find the grid square:
    /* std::cout << "obs = " << obsp << std::endl; */
    iobs = (*obs_terp).cor ? (*obs_terp).hunt(obsp) : (*obs_terp).locate(obsp);
    /* std::cout << "iobs = " << iobs << std::endl; */
    /* std::cout << "obs[" << iobs << "] = " << (*obs_terp).xx[iobs] << std::endl; */
    /* std::cout << "obs[" << iobs + 1 << "] = " << (*obs_terp).xx[iobs + 1] << std::endl; */
    /* std::cout << "nH = " << nHp << std::endl; */
    inH = (*nH_terp).cor ? (*nH_terp).hunt(nHp) : (*nH_terp).locate(nHp);
    // std::cout << "inH = " << inH << std::endl;
    // std::cout << "nH[" << inH << "] = " << (*nH_terp).xx[inH] << std::endl;
    // std::cout << "nH[" << inH + 1 << "] = " << (*nH_terp).xx[inH + 1] << std::endl;
    // std::cout << "T = " << Tp << std::endl;
    iT = (*T_terp).cor ? (*T_terp).hunt(Tp) : (*T_terp).locate(Tp);
    // std::cout << "iT = " << iT << std::endl;
    // std::cout << "T[" << iT << "] = " << (*T_terp).xx[iT] << std::endl;
    // std::cout << "T[" << iT + 1 << "] = " << (*T_terp).xx[iT + 1] << std::endl;
    
    // std::cout << "iobs = " << iobs << std::endl;
    // std::cout << "inH = " << inH << std::endl;
    // std::cout << "iT = " << iT << std::endl;    

    check_intensities(iobs,inH,iT,IPHb);//check that these intensities have been simulated
    
    //interpolate:
    aobs = (obsp-(*obs_terp).xx[iobs])/((*obs_terp).xx[iobs+1]-(*obs_terp).xx[iobs]);
    anH = (nHp-(*nH_terp).xx[inH])/((*nH_terp).xx[inH+1]-(*nH_terp).xx[inH]);
    aT = (Tp-(*T_terp).xx[iT])/((*T_terp).xx[iT+1]-(*T_terp).xx[iT]);

    /* std::cout << "aobs = " << aobs << std::endl; */
    /* std::cout << "anH = " << anH << std::endl; */
    /* std::cout << "aT = " << aT << std::endl;     */
    
    iIcalc = (1.-aobs)*(1.-anH)*(1.-aT)*Icalc[iobs][inH][iT]
      + aobs*(1.-anH)*(1.-aT)*Icalc[iobs+1][inH][iT]
      + (1.-aobs)*anH*(1.-aT)*Icalc[iobs][inH+1][iT]
      + (1.-aobs)*(1.-anH)*aT*Icalc[iobs][inH][iT+1]
      + aobs*anH*(1.-aT)*Icalc[iobs+1][inH+1][iT]
      + (1.-aobs)*anH*aT*Icalc[iobs][inH+1][iT+1]
      + aobs*(1.-anH)*aT*Icalc[iobs+1][inH][iT+1]
      + aobs*anH*aT*Icalc[iobs+1][inH+1][iT+1];

    return iIcalc;
  }

  void derivs(double obsp, double nHp, double Tp, double IPHb, VecDoub_O &dIcalc) {
    //find the grid square:
    iobs = (*obs_terp).cor ? (*obs_terp).hunt(obsp) : (*obs_terp).locate(obsp);
    inH = (*nH_terp).cor ? (*nH_terp).hunt(nHp) : (*nH_terp).locate(nHp);
    iT = (*T_terp).cor ? (*T_terp).hunt(Tp) : (*T_terp).locate(Tp);

    check_intensities(iobs,inH,iT,IPHb);//check that these intensities have been simulated
    
    double obs1=(*obs_terp).xx[iobs];
    double obs2=(*obs_terp).xx[iobs+1];
    double nH1=(*nH_terp).xx[inH];
    double nH2=(*nH_terp).xx[inH+1];
    double T1=(*T_terp).xx[iT];
    double T2=(*T_terp).xx[iT+1];

    double fobs=(obsp-obs1)/(obs2-obs1);
    double fnH=(nHp-nH1)/(nH2-nH1);
    double fT=(Tp-T1)/(T2-T1);

    double I111=Icalc[iobs][inH][iT];
    double I112=Icalc[iobs][inH][iT+1];
    double I121=Icalc[iobs][inH+1][iT];
    double I122=Icalc[iobs][inH+1][iT+1];
    double I211=Icalc[iobs+1][inH][iT];
    double I212=Icalc[iobs+1][inH][iT+1];
    double I221=Icalc[iobs+1][inH+1][iT];
    double I222=Icalc[iobs+1][inH+1][iT+1];

    //we must bilinearly interpolate the derivatives
    double dIdobs, dIdnH, dIdT;
    dIdobs =  (1.-fnH)*(1.-fT)*(I211-I111);
    dIdobs += fnH*(1.-fT)*(I221-I121);
    dIdobs += (1.-fnH)*fT*(I212-I112);
    dIdobs += fnH*fT*(I222-I122);
    dIdobs /= (obs2-obs1);

    dIdnH =  (1.-fobs)*(1.-fT)*(I121-I111);
    dIdnH += fobs*(1.-fT)*(I221-I211);
    dIdnH += (1.-fobs)*fT*(I122-I112);
    dIdnH += fobs*fT*(I222-I212);
    dIdnH /= (nH2-nH1);

    dIdT =  (1.-fobs)*(1.-fnH)*(I112-I111);
    dIdT += fobs*(1.-fnH)*(I212-I211);
    dIdT += (1.-fobs)*fnH*(I122-I121);
    dIdT += fobs*fnH*(I222-I221);
    dIdT /= (T2-T1);

    dIcalc[0]=dIdobs;
    dIcalc[1]=dIdnH;
    dIcalc[2]=dIdT;
  }

  //this is used for the chi-square minimzation by L-M
  void operator()(const double obsp, VecDoub_I &a, double &iIcalc, VecDoub_O &dIcalcda)
  {
    //    std::cout << "T[5] = " << (*T_terp).xx[5] << std::endl;
    double nHp = a[0];
    double Tp = a[1];
    double IPHb = a[2];

    iIcalc = interp(obsp, nHp, Tp, IPHb);

    VecDoub dItemp(3);
    derivs(obsp, nHp, Tp, IPHb, dItemp);
    dIcalcda[0]=dItemp[1];
    dIcalcda[1]=dItemp[2];
    dIcalcda[2]=(*obs_sim).IPH->sim(1.0,nHp,Tp,obsp);
  }
  
};

#endif
