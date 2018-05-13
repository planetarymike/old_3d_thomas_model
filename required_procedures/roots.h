//roots.h -- root finding from NR, Ch. 9

#include <vector>
#include <cmath>
#include "nr3.h"
#include <iostream>

template <class T>
double zbrent(T &func, const double x1, const double x2, const double tol)
// using Brent's method, return the root of a function or functor func
// known to lie between x1 and x2. The root is refined until its
// accuracy is tol.
{
  const int itmax = 100; // maximum number of iterations
  const double EPS=numeric_limits<Doub>::epsilon(); // machine precision
  double a=x1, b=x2, c=x2, d, e, fa=func(a), fb=func(b), fc, p, q, r, s, tol1, xm;
  if  ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb <0.0))
    toss("Root must be bracketed in zbrent");
  fc = fb;
  for (int iter = 0; iter < itmax; iter++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      // rename a, b, c, and adjust bounding interval d
      c = a;
      fc = fa;
      e = d = (b - a);
    }
    if (abs(fc) < abs(fb)) {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
    tol1 = 2.0*EPS*abs(b)+0.5*tol; // convergence check
    //    std::cout << "b = " << b << std::endl;
    xm = 0.5*(c-b);
    if (abs(xm) <= tol1 || fb == 0.0) return b;
    if (abs(e) >= tol1 && abs(fa) > abs(fb)) { // attempt inverse quadratic interpolation
      s = fb/fa;
      if (a == c) {
	p = 2.0*xm*s;
	q = 1.0-s;
      } else {
	q = fa/fc;
	r = fb/fc;
	p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
	q = (q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0) q = -q; // check whether interpolation returns in bounds value
      p = abs(p);
      double min1 = 3.0*xm*q - abs(tol1*q);
      double min2 = abs(e*q);
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
	e = d; // accept interpolation
	d = p/q;
      } else {
	d = xm;// interpolation failed; use bisection
	e = d;
      }
    } else { // bounds decreasing too slowly; use bisection
      d = xm;
      e = d;
    }
    a = b; // move last best guess to a;
    fa = fb;
    if (abs(d) > tol1) // evaluate new trial root
      b += d;
    else
      b += SIGN(tol1,xm);
    fb = func(b);
  }
  toss("Maximum number of iterations exceed in zbrent");
}
