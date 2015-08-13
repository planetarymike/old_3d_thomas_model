#include <mpi.h> // parallelization//must be included in top line of main file
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

  ofstream HolTfile;
  HolTfile.open("../Model_Source_Functions/HolTinterp.dat");
  int w = 15;
  HolTfile.width(w);
  HolTfile << ntaus << std::endl;
  HolTfile.width(w);
  HolTfile << "tau";
  HolTfile.width(w);
  HolTfile << "HolT(tau)" << std::endl;
  for(double tau=taubottom; tau < tautop+taustep; tau+=taustep) {
    HolTfile.width(w);
    HolTfile << tau;
    HolTfile.width(w);
    HolTfile << HolT(tau) << std::endl;
  }
  HolTfile.close();


  return 0;
}
