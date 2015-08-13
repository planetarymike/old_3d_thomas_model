//simpletrap.cpp -- basic integrations with the trapezoidal rule
#include <iostream>
#include <cmath>
#include "quadrature.h"

const double pi = 3.1415926535898;

struct Sine {
  double p, o, a;
  Sine(const double pp = 2*pi, const double oo = 0.0, const double aa = 1.0) 
    : p(pp), o(oo), a(aa) {}
  double operator()(const double x) {return a*sin(x*p/(2*pi)+o);
  }
};

int main()
{
  Sine asin;
  Trapzd<Sine> s(asin,0,pi);

  int m = 12;
  double val = 0.0;
  for(int j = 1; j <= m+1; j++)
    {
      val = s.next();
      std::cout << "After " << j << " iterations, val = " << val << std::endl;
    }

  //  int N = 21;
  //  int i = 0;
  //  double x = 0;
  //  for(i = 0; i < 21; i++)
  //    {
  //     std::cout << "sin(" << x << ") = " << asin(x) << std::endl;
  //     x += 2*pi/N;
  //    }

  return 0;
}
