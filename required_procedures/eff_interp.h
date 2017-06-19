//eff_interp.h -- routines to convert between effusion velocity and temperature.

#ifndef __EFF_INTERP__
#define __EFF_INTERP__

#include <iostream> // for file output and dialog
#include <fstream>  // ""
#include <cmath>    // for exp and log
#include <cstdlib>
#include "definitions.h" // basic parameter definitions
#include "interp_1d.h" //interpolation

struct eff_interp
{
  /* Structure to load tabular Holstein functions */
  string efffilename;
  bool silent;
  int npts;
  VecDoub temp_vec, eff_vec;
  double effmin,effmax;
  Linear_interp effinterp;

  eff_interp(string efffilename=eff_filename, bool silent=FALSE) {
    if (!silent)
      std::cout << "Reading interpolated Hol values from " << efffilename << std::endl;
    //load the data from the file
    ifstream efffile;
    efffile.open(efffilename.c_str());
    //read number of points
    efffile >> npts;
    /* std::cout << "npts = " << npts << std::endl; */
    temp_vec.resize(npts);
    eff_vec.resize(npts);
    // clear the header
    char dumline[100];
    efffile.getline(dumline,100);
    efffile.getline(dumline,100);
    //    std::cout << dumline << std::endl;
    for (int i = 0; i < npts; i++) {
      efffile >> temp_vec[i];
      //      std::cout << "temp_vec[" << i << "] = " << temp_vec[i] << std::endl;
      efffile >> eff_vec[i];
      //      std::cout << "eff_vec[" << i << "] = " << eff_vec[i] << std::endl;
      //      std::cin.get();
    }
    efffile.close();
    effmin= eff_vec[0] > eff_vec[npts-1] ? eff_vec[npts-1] : eff_vec[0];
    effmax= eff_vec[0] > eff_vec[npts-1] ? eff_vec[0] : eff_vec[npts-1];
    effinterp = Linear_interp(eff_vec,temp_vec);
  }

  double interp(const double eff) {
    if (eff > effmax || eff < effmin)
      throw("eff out of range!");
    return effinterp.interp(eff);
  }
};

#endif
