/* interpgen.h -- generate interpolating function for hydrogen and CO2
   and return it. Important, because computing the physical atmosphere
   dominates runtime without interpolation. */

#ifndef __INTERPGEN__
#define __INTERPGEN__

#include <iostream> // for file output and dialog
#include <fstream>  // ""
#include <cmath>    // for exp and log
#include <cstdlib>
#include "physical.h" // physical definitions, nH, nCO2, KK, etc.
#include "definitions.h" // basic parameter definitions
#include "interp_1d.h" //interpolation
//#include <mpi.h> // parallelization//must be included in top line of main file

using std::log;
using std::exp;

/* function to interpolate across the atmosphere. Either takes a file
   for input or computes directly. It's kludegy, but this is required
   for backward compatibility with much older code. Writing/reading
   files is good because it allows for line integrations to use
   exactly the same densities as the source function generation.*/
struct atmointerp {
  double nexo, Texo;
  int npts;
  string fname;
  VecDoub logalt, lognH, lognCO2;
  Linear_interp lognHinterp;
  Linear_interp lognCO2interp;
  bool init;

  //default constructor
  atmointerp() {init=0;}

  atmointerp(double nexoo, 
	     double Texoo,
	     int nptss, 
	     string fnamee = "") : nexo(nexoo), 
				   Texo(Texoo), 
				   npts(nptss), 
				   fname(fnamee) 
  {
    if (fname.size() == 0) {
      char filename[100];
      sprintf(filename,"atm_nH%G_T%G_npts=%d.dat",nexo,Texo,npts);
      fname=srcfnsloc+filename;
    }
    logalt.resize(npts);
    lognH.resize(npts);
    lognCO2.resize(npts);
    //    std::cout << "npts = " << npts << "\n";
    //two ways to get interp points: from a file, or by direct calculation
    ifstream file;
    file.open(fname.c_str());
    if (file.good()) { // if the appropriate file exists, read from that
      std::cout << "Reading atmospheric density values from " << fname << std::endl;

      // clear the header
      char dumline[100];
      file.getline(dumline,100);
      //file structure looks like this:
      //n:     logalt:    alt (km)    lognH:   nH (cm^-3)    lognCO2:   nCO2 (cm^-3)
      //we only need logalt, lognH, and lognCO2
      double dum; // dummy variable for readin only
      for (int i = 0; i < npts; i++) {
	file >> dum;
	file >> logalt[i];
	// std::cout << "logalt[" << i << "] = " << logalt[i] << "\n";
	file >> dum;
	file >> lognH[i];
 	// std::cout << "lognH[" << i << "] = " << lognH[i] << "\n";
	file >> dum;
	file >> lognCO2[i];
	// std::cout << "lognCO2[" << i << "] = " << lognCO2[i] << "\n";
	file >> dum;
	// cin.get();
      }
      file.close();
      init=1;
    } else { // no physical interpolation file exists, calculate directly
      file.close();

      tabular_atmo atmotable(nexo,Texo,npts,fname);

      for (int i = 0; i < npts; i++) {
	logalt[i] = (atmotable.logaltvec)[i];
	//	std::cout << "logalt[" << i << "] = " << logalt[i] << " , ";
	lognH[i] = (atmotable.lognHvec)[i];
	//	std::cout << "lognH[" << i << "] = " << lognH[i] << ", ";
	lognCO2[i] = (atmotable.lognCO2vec)[i];
	//	std::cout << "lognCO2[" << i << "] = " << lognCO2[i] << "\n";
	//	std::cin.get();
      }
      //      getatmointerp(logalt,lognH,lognCO2,nexo,Texo,npts);
    }
    //    std::cout << "npts = " << npts << "\n";

    //proceed with interpolation
    lognHinterp = Linear_interp(logalt, lognH); // interpolate
    lognCO2interp = Linear_interp(logalt, lognCO2); // interpolate
    init=1;
  }
  double nH(const double &r) {
    if (!init)
      throw("atmointerp not initialized!");
    double km = (r - rMars)/1e5;
    return exp(lognHinterp.interp(log(km)));
  }
  double nCO2(const double &r) {
    if (!init)
      throw("atmointerp not initialized!");
    double km = (r - rMars)/1e5;
    /* if (km > 300) { */
    /*   return 0; */
    /* } else { */
      return exp(lognCO2interp.interp(log(km)));
    /* } */
  }
  double nH(const double &r, const double &t) {
    return nH(r);
  }
  double nCO2(const double &r, const double &t) {
    return nCO2(r);
  }
  double nH(const double &r, const double &t, const double &p) {
    return nH(r);
  }
  double nCO2(const double &r, const double &t, const double &p) {
    return nCO2(r);
  }

  double nH(const VecDoub &rpt) {
    double thisr=0.0;
    for (int i=0; i<rpt.size();i++)
      thisr+=rpt[i]*rpt[i];
    thisr=sqrt(thisr);
    return nH(thisr);
  }
  double nCO2(const VecDoub &rpt) {
    double thisr=0.0;
    for (int i=0; i<rpt.size();i++)
      thisr+=rpt[i]*rpt[i];
    thisr=sqrt(thisr);
    return nCO2(thisr);
  }
};

void makeCO2file(string altfname, string outfname, int npts, double Texo) {
  /* Function to generate CO2 densities at the altitudes specified by
     input altitude file and write them out to outfile.*/

  double *logalt;
  double *alt;
  double *lognCO2;
  double *nCO2vec;

  logalt = new double[npts];
  alt = new double[npts];
  lognCO2 = new double[npts];
  nCO2vec = new double[npts];

  /* Read in the existing values, and add CO2 to the mix */
  std::cout << "Reading atmospheric altitudes from " << altfname << std::endl;
  ifstream altfile;
  altfile.open(altfname.c_str());

  // clear the header
  char dumline[100];
  altfile.getline(dumline,100);

  double dum; // dummy variable for readin only
  for (int i = 0; i < npts; i++) {
    //load the altitudes and compute the CO2 densities
    //altfile structure looks like this:
    //    n:        logalt:       alt (km)         lognH:     nH (cm^-3)
    altfile >> dum;
    altfile >> logalt[i];
    altfile >> alt[i];
    //    std::cout << "alt[" << i << "] = " << alt[i] << " , ";
    altfile >> dum;
    altfile >> dum;

    nCO2vec[i]=nCO2(rMars+1e5*alt[i],Texo);
    //    std::cout << "nCO2vec[" << i << "] = " << nCO2vec[i] << std::endl;
    lognCO2[i]=log(nCO2vec[i]);
    //    std::cin.get();

  }
  altfile.close();

  
  //save the values to outfile
  ofstream outfile;
  outfile.open(outfname.c_str());
  int iw = ((int) log(npts)/log(10))+3;
  int w = 15;
  outfile.width(iw);
  outfile << "n:";
  outfile.width(w);
  outfile << "logalt:";
  outfile.width(w);
  outfile << "alt (km)";
  outfile.width(w);
  outfile << "lognCO2:";
  outfile.width(w);
  outfile << "nCO2 (cm^-3)"<< std::endl;
  for (int i = 0; i < npts; i++) {
    outfile.width(iw);
    outfile << i;
    outfile.width(w);
    outfile << logalt[i];
    outfile.width(w);
    outfile << alt[i];
    outfile.width(w);
    outfile << lognCO2[i];
    outfile.width(w);
    outfile << nCO2vec[i] << std::endl;
  }
  outfile.close();

  return;
}



void addCO2(string innHname, string innCO2name, string outfname, int npts) {
  /* Function to add CO2 density values to preexisting nH density
     files. Historical curiousity only, because the code
     auto-generates files with CO2 included now. Only previous
     versions did not. */
  
  double *logalt;
  double *alt;
  double *lognH;
  double *nH;
  double *lognCO2;
  double *nCO2;

  logalt = new double[npts];
  alt = new double[npts];
  lognH = new double[npts];
  nH = new double[npts];
  lognCO2 = new double[npts];
  nCO2 = new double[npts];

  /* Read in the existing values, and add CO2 to the mix */
  std::cout << "Reading atmospheric density values from " << innHname 
	    << " and " << innCO2name << " ." << std::endl;
  ifstream innHfile;
  innHfile.open(innHname.c_str());
  ifstream innCO2file;
  innCO2file.open(innCO2name.c_str());
  // clear the header
  char dumline[100];
  innHfile.getline(dumline,100);
  innCO2file.getline(dumline,100);

  double dum; // dummy variable for readin only
  for (int i = 0; i < npts; i++) {
    //first load the H data
    //innHfile structure looks like this:
    //    n:        logalt:       alt (km)         lognH:     nH (cm^-3)
    innHfile >> dum;
    innHfile >> logalt[i];
    //	std::cout << "logalt[" << i << "] = " << logalt[i] << " , ";
    innHfile >> alt[i];
    innHfile >> lognH[i];
    //	std::cout << "lognH[" << i << "] = " << lognH[i] << "\n";
    innHfile >> nH[i];

    //now the CO2 data
    //innCO2file looks like this:
    //    n:        logalt:       alt (km)         lognCO2:   nCO2 (cm^-3)
    innCO2file >> dum;
    innCO2file >> dum;
    if ((dum-logalt[i])/logalt[i] > 1e-4) { throw("Altitudes do not match!"); }
    //	std::cout << "logalt[" << i << "] = " << logalt[i] << " , ";
    innCO2file >> dum;
    innCO2file >> lognCO2[i];
    //	std::cout << "lognH[" << i << "] = " << lognH[i] << "\n";
    innCO2file >> nCO2[i];
  }
  innHfile.close();
  innCO2file.close();
  
  //save the values to outfile
  ofstream outfile;
  outfile.open(outfname.c_str());
  int iw = ((int) log(npts)/log(10))+3;
  int w = 15;
  outfile.width(iw);
  outfile << "n:";
  outfile.width(w);
  outfile << "logalt:";
  outfile.width(w);
  outfile << "alt (km)";
  outfile.width(w);
  outfile << "lognH:";
  outfile.width(w);
  outfile << "nH (cm^-3)";
  outfile.width(w);
  outfile << "lognCO2:";
  outfile.width(w);
  outfile << "nCO2 (cm^-3)"<< std::endl;
  for (int i = 0; i < npts; i++) {
    outfile.width(iw);
    outfile << i;
    outfile.width(w);
    outfile << logalt[i];
    outfile.width(w);
    outfile << alt[i];
    outfile.width(w);
    outfile << lognH[i];
    outfile.width(w);
    outfile << nH[i];
    outfile.width(w);
    outfile << lognCO2[i];
    outfile.width(w);
    outfile << nCO2[i] << std::endl;
  }
  outfile.close();

  return;
}

#endif
  

