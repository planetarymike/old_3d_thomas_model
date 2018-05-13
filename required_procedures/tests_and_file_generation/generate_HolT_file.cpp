#include <cmath>
#include <iostream>
#include <fstream>
#include "rad_trans.h"
#include <limits> //sizes of double
typedef std::numeric_limits< double > dbl;

using std::log;
using std::exp;

int main() {
  // a program to generate intepolation points for the function HolT
  double taubottom=1e-4;
  double tautop=1e5;
  int ntaus=500;

  double lb=log(taubottom);
  double lt=log(tautop);
  double space=(lt-lb)/(ntaus-1.0);

  ofstream HolTfile;
  HolTfile.open("../../tabulated_data/HolTinterp.dat");
  int w = 15;
  HolTfile.width(w);
  HolTfile << ntaus << std::endl;
  HolTfile.width(w);
  HolTfile << "log(tau)";
  HolTfile.width(w);
  HolTfile << "log(HolT(tau))" << std::endl;
  for(double itau=0;itau<ntaus;itau++) {
    HolTfile.width(dbl::digits10+10);
    HolTfile.precision(dbl::digits10);
    double logtau=lb+space*itau;
    HolTfile << logtau;
    HolTfile << "\t";
    HolTfile.width(dbl::digits10+10);
    HolTfile.precision(dbl::digits10);
    HolTfile << log(HolT(exp(logtau))) << std::endl;
  }
  HolTfile.close();


  return 0;
}
