#include <cmath>
#include <iostream>
#include <fstream>
#include "rad_trans.h"


int main() {
  // a program to generate intepolation points for the function HolT
  double taubottom=0.0;
  double tautop=100.0;
  double taustep=0.001;
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
    HolGfile.width(w);
    HolGfile << tau;
    HolGfile.width(w);
    HolGfile << HolG(tau) << std::endl;
  }
  HolGfile.close();


  return 0;
}
