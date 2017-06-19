#include <cmath>
#include <iostream>
#include <fstream>
#include "rad_trans.h"
#include <limits> //sizes of double
typedef std::numeric_limits< double > dbl;

int main() {
  // a program to generate intepolation points for the function HolT
  double tempbottom=100.0;
  double temptop=1600.0;
  double tempstep=1.0;
  int ntemps=(temptop-tempbottom)/tempstep+1;
  double lc, lcsqrt, u, eff;
  
  ofstream efffile;
  efffile.open("../../tabulated_data/eff_interp.dat");
  int w = 15;
  efffile.width(w);
  efffile << ntemps << std::endl;
  efffile.width(w);
  efffile << "temp [K]";
  efffile.width(w);
  efffile << "eff_v [cm/s]" << std::endl;
  for(double temp=tempbottom; temp <= temptop; temp+=tempstep) {
    efffile.width(dbl::digits10+10);
    efffile.precision(dbl::digits10);
    efffile << temp;
    efffile << "\t";
    efffile.width(dbl::digits10+10);
    efffile.precision(dbl::digits10);
    lc = G*mMars*mH/(kB*rexo*temp);
    lcsqrt=std::sqrt(lc);
    u = std::sqrt(2*kB*temp/mH);
    eff = u*exp(-lc)*(lc+1)/2.0/std::sqrt(pi);
    efffile << eff << std::endl;
  }
  efffile.close();


  return 0;
}
