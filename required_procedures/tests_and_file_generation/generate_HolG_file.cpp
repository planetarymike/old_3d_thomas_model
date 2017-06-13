#include <cmath>
#include <iostream>
#include <fstream>
#include "rad_trans.h"
#include <limits> //sizes of double
typedef std::numeric_limits< double > dbl;

int main() {
  // a program to generate intepolation points for the function HolT
  double taubottom=0.0;
  double tautop=100000.0;
  double taustep=0.05;
  int ntaus=(tautop-taustep)/taustep+1;

  ofstream HolGfile;
  HolGfile.open("../../tabulated_data/HolGinterp.dat");
  int w = 15;
  HolGfile.width(w);
  HolGfile << ntaus << std::endl;
  HolGfile.width(w);
  HolGfile << "tau";
  HolGfile.width(w);
  HolGfile << "HolG(tau)" << std::endl;
  for(double tau=taubottom; tau < tautop+taustep; tau+=taustep) {
    HolGfile.width(dbl::digits10+10);
    HolGfile.precision(dbl::digits10);
    HolGfile << tau;
    HolGfile << "\t";
    HolGfile.width(dbl::digits10+10);
    HolGfile.precision(dbl::digits10);
    HolGfile << HolG(tau) << std::endl;
  }
  HolGfile.close();


  return 0;
}
