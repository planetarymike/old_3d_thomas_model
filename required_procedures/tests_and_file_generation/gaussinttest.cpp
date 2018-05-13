//gaussinttest.cpp -- test the gauss-legendre and gauss-laguerre quadrature against analytic results

#include <iostream>
#include <cmath>
#include <vector>
#include "gauss_wgts.h"

using std::vector;
using std::cout;
using std::endl;

const double pi = 3.1415926535898;

//function prototypes
double f1(double x);
double f2(double x);


int main()
{
  int NPTS = 32; // the number of abcissas and weights to return

  //create the vectors to hold the weights and abcissae
  vector<double> legx (NPTS), legw (NPTS), lagx (NPTS), lagw (NPTS);
  
  gauleg(-1.0,1.0,legx,legw); // get x and w for gauss-legendre
  gaulag(lagx,lagw,2.0);      // same for gauss-laguerre

  double sum1 = 0;
  double sum2 = 0;
  for (int i = 0; i < NPTS; i++)
    {
      sum1 += legw[i]*f1(legx[i]);
      sum2 += lagw[i]*f2(lagx[i]);
    }

  cout << "With " << NPTS << " quadrature points:" << endl;
  cout << "I1 = " << sum1 << endl;
  cout << "I2 = " << sum2 << endl;
  
  return 0;
}


double f1(double x)
{
  return std::cos(x);
}

double f2(double x)
{
  double c = std::cos(x);
  return c*c;
}
