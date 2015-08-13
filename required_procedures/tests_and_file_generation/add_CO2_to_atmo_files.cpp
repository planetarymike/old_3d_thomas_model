#include <mpi.h> // parallelization//must be included in top line of main file
#include <iostream> // for file output and dialog
#include <fstream>  // ""
#include <cmath>    // for exp and log
#include <cstdlib>
#include "interpgen.h" //all routines to read and write CO2 and nH interpolation data files


int main() {
  // A program to build synthetic atmosphere files from seperate CO2
  // and H interpolation functions.

  char nHfile[100];
  char CO2file[100];
  char outfile[100];

  for (double nHexo=10000; nHexo<730000; nHexo+=20000) {
    for (double Texo=100; Texo<820; Texo+=20) {
      sprintf(nHfile,"../Model_Source_Functions/nH%G_T%G.dat",nHexo,Texo);
      sprintf(CO2file,"../Model_Source_Functions/nCO2_T%G.dat",Texo);
      sprintf(outfile,"../Model_Source_Functions/atm_nH%G_T%G.dat",nHexo,Texo);
    
      addCO2(nHfile,
	     CO2file,
	     outfile,
	     10000);
    }
  }


  return 0;
}


// int main() {
//   // This program generates CO2 files at the altitudes corresponding
//   // to those for which nH is tabulated.

//   char altfile[100];
//   char CO2file[100];

//   for (double Texo=100; Texo<820; Texo+=20) {
//     sprintf(altfile,"../Model_Source_Functions/nH10000_T%G.dat",Texo);
//     sprintf(CO2file,"../Model_Source_Functions/nCO2_T%G.dat",Texo);
    
//     makeCO2file(altfile,
// 		CO2file,
// 		10000,
// 		Texo);
//   }
//   return 0;
// }
