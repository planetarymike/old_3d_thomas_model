#include <cmath>
#include <iostream>
#include <fstream>
#include "rad_trans.h"
#include <limits> //sizes of double
typedef std::numeric_limits< double > dbl;

int main() {
  // a program to generate intepolation points for the function HolT
  double taubottom=0.0;
  double tautop=10000.0;
  double taustep=0.01;
  int ntaus=(tautop-taustep)/taustep+1;

  ofstream HolTfile;
  HolTfile.open("../../tabulated_data/HolTinterp.dat");
  int w = 15;
  HolTfile.width(w);
  HolTfile << ntaus << std::endl;
  HolTfile.width(w);
  HolTfile << "tau";
  HolTfile.width(w);
  HolTfile << "HolT(tau)" << std::endl;
  for(double tau=taubottom; tau < tautop+taustep; tau+=taustep) {
    HolTfile.width(dbl::digits10+10);
    HolTfile.precision(dbl::digits10);
    HolTfile << tau;
    HolTfile << "\t";
    HolTfile.width(dbl::digits10+10);
    HolTfile.precision(dbl::digits10);
    HolTfile << HolT(tau) << std::endl;
  }
  HolTfile.close();


  return 0;
}
