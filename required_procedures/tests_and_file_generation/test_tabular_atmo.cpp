// Program to test the speed and accuracy of the tabular atmosphere generation code


#include <iostream> // for file output and dialog
#include <cmath>    // for cos and sin
#include <cstdlib>
#include "nr3.h" // type VecDoub


#include "definitions.h" // basic parameter definitions
#include "physical.h" // physical definitions, nH, nCO2, KK, etc.
#include "quadrature.h" // to compute line integrals



int main(int argc, char* argv[]) {

  //get nexo and Texo from the command line
  double nexo = atof(argv[1]);
  std::cout << "nexo = " << nexo << std::endl;
  double Texo = atof(argv[2]);
  std::cout << "Texo = " << Texo << std::endl;

  tabular_atmo thisatmointerp(nexo, Texo, 1000000);

  double alt;
  std::cout << "Enter an altitude: ";
  while (std::cin >> alt) {
    std::cout << "nCO2(" << alt << ") = " << thisatmointerp.nCO2(rMars+alt*1e5) << std::endl;
    std::cout << "nH(" << alt << ") = " << thisatmointerp.nH(rMars+alt*1e5) << std::endl;
    std::cout << "Enter an altitude: ";
  }
  
  return 0;
  
}
