#include <cmath>
#include <iostream>
#include <fstream>
#include "rad_trans.h"
#include <limits> //sizes of double
typedef std::numeric_limits< double > dbl;

using std::log;
using std::exp;

int main() {
  // a program to generate intepolation points for the function HolG
  double taubottom=1e-4;
  double tautop=1e5;
  int ntaus=500;

  double lb=log(taubottom);
  double lt=log(tautop);
  double space=(lt-lb)/(ntaus-1.0);

  ofstream HolGfile;
  HolGfile.open("../../tabulated_data/HolGinterp.dat");
  int w = 15;
  HolGfile.width(w);
  HolGfile << ntaus << std::endl;
  HolGfile.width(w);
  HolGfile << "log(tau)";
  HolGfile.width(w);
  HolGfile << "log(HolG(tau))" << std::endl;
  for(double itau=0;itau<ntaus;itau++) {
    HolGfile.width(dbl::digits10+10);
    HolGfile.precision(dbl::digits10);
    double logtau=lb+space*itau;
    HolGfile << logtau;
    HolGfile << "\t";
    HolGfile.width(dbl::digits10+10);
    HolGfile.precision(dbl::digits10);
    HolGfile << log(HolG(exp(logtau))) << std::endl;
  }
  HolGfile.close();


  return 0;
}
