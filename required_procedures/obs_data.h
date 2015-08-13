//obs_data.h -- object to read observation data from file into memory

#ifndef __OBS_DATA_H
#define __OBS_DATA_H


#include <iostream>
#include <fstream>
#include "nr3.h"


struct obs_data {

  string obsname;
  string fname;
  int ndat;
  double aMars; // Sun-Mars distance at time of obs, AU
  double Fsun_earth; // Solar flux at line center at Earth
  double Fsun_mars; // Solar flux at line center at Mars, possibly calculated from above
  char datetime[1024];

  VecDoub I_obs, DI_obs;  // observed intensity (kR) and its uncertainty
  VecDoub alttan, SZAtan; // radius and SZA of tangent point
  VecDoub ra, dec, scvel; //ra and dec of los, and S/C SSB vel projected onto this direction
  MatDoub pos, dir; // position and LOS direction of spacecraft.

  //constructor takes the name of the observation as an argument and loads the data
  obs_data(string obsnamee): obsname(obsnamee) {
     fname = obsname + "_obs.dat";
     fstream obsfile;
     obsfile.open(fname.c_str());//this assumes that the code is executing 
                              //from the directory containing this file

     //now read in the data:
       //read the preliminary data in
     obsfile >> ndat; // number of observations in this file
         /* std::cout << "ndat = " << ndat << std::endl; */
     obsfile >> aMars; // Sun-Mars distance at time of obs, AU
         /* std::cout << "aMars = " << aMars << std::endl; */
     obsfile >> Fsun_earth; // photons/s/m2/Angstrom Solar flux at line center at Earth
     Fsun_earth /= 1e4; // photons/s/cm2/Angstrom
     Fsun_earth *= (1215.)*(1215e-8)/(3e10/*cm/s*/);// photons/s/cm2/Hz
     Fsun_mars = Fsun_earth/aMars/aMars;
         /* std::cout << "Fsun_mars = " << Fsun_mars << std::endl; */
     obsfile.getline(datetime,1024); // capture time of file (irrelevant for this code)

     //prepare vector parameters for read-in
     I_obs.resize(ndat);
     DI_obs.resize(ndat);
     alttan.resize(ndat);
     SZAtan.resize(ndat);
     pos.resize(ndat,3);
     dir.resize(ndat,3);
     ra.resize(ndat);
     dec.resize(ndat);
     scvel.resize(ndat);
     
     /* X axis points toward Sun!! */
     //loop over the number of observations, reading each in to the appropriate variable
     for (int row = 0; row < ndat; row++) {
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
			 
       /* std::cout << "On row " << row << ":\n"; */
       /* std::cout << "I_obs = " << I_obs[row] << " +- " << DI_obs[row] << " kR,\n"; */
       /* std::cout << "pos = { " << pos[row][0] << " , " */
       /* 	 << pos[row][1] << " , " << pos[row][2] << " },\n"; */
       /* std::cout << "r_pos = " */
       /* 	 << sqrt(pos[row][0]*pos[row][0] */
       /* 		 + pos[row][1]*pos[row][1] */
       /* 		 + pos[row][2]*pos[row][2]) << " cm\n"; */
       /* std::cout << "dir = { " << dir[row][0] */
       /* 	 << " , " << dir[row][1] */
       /* 	 << " , " << dir[row][2] << " },\n"; */
       /* std::cout << "alttan = " << alttan[row] << " , SZAtan = " << SZAtan[row] << ".\n"; */
       /* std::cout << "ra = " << ra[row] << " , dec = " << dec[row] << ".\n"; */
       /* std::cout << "scvel = " << scvel[row] << ".\n"; */
       /* /\* std::cin.get(); *\/ */
     }
     obsfile.close();
  }


};


#endif
