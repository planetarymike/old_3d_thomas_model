//getgausswts.cpp -- get and display weights for gaussian quadrature

#include <iostream>
#include <cmath>
#include <vector>
#include "gauss_wgts.h"

using std::vector;
using std::cout;
using std::endl;

int main()
{
  int NPTS = 30; // the number of abcissas and weights to return

  vector<double> legx (NPTS), legw (NPTS), lagx (NPTS), lagw (NPTS);
  
  gauleg(-1.0,1.0,legx,legw);
  gaulag(lagx,lagw,0.0);

  cout << "Gauss-Legendre Quadrature Weights" << endl;
  for(int i = 0; i < NPTS; i++)
    cout << "Abcissa " << i << ": " << legx[i] << ", Weight " << i << ": " << legw[i] << endl;

  cout << "Gauss-Laguerre Quadrature Weights" << endl;
  for(int i = 0; i < NPTS; i++)
  cout << "Abcissa " << i << ": " << lagx[i] << ", Weight " << i << ": " << lagw[i] << endl;
 
  return 0;
}
