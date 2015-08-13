//gauss_wgts.h -- header file to determine abcissas and weights for gauss-jordan quadrature

#include <vector>
#include <cmath>
#include <iostream>
#include <gsl/gsl_sf_gamma.h>

#ifndef __GAUSS_WGTS_H
#define __GAUSS_WGTS_H

using std::vector;
using std::abs;

void gauleg(const double x1, const double x2, VecDoub &x, VecDoub &w)
// Given a lower and upper limit, this returns arrays x and w (of length n),
// containing the abcissas and wieghts of the Gauss-Legendre n-point quadrature 
// formula (see NR pg 183)
{
  const double EPS = 1.0e-14; // EPS is the relative precision
  double z1, z, xm, xl, pp, p3, p2, p1;
  int n = x.size(); // the number of abcissas and weights
  int m = (n+1)/2; // number of roots to find since roots are symmetric
  xm = 0.5*(x2+x1);
  xl = 0.5*(x2-x1);


  for (int i = 0; i < m; i++) // loop to find roots
    {
      z = cos(3.141592654 * (i + 0.75) / (n + 0.5)); // initial guess for root
      // use Newton's method to find the root at the desired relative precision
      do {
	p1 = 1.0;
	p2 = 0.0;
	for (int j = 0; j < n; j++) // use the Legendre polynomial recurrence
	  {                         //   to evaluate the polynomial at z
	    p3 = p2;                //   uses recurrence relation
	    p2 = p1;                // (j+1)*P_(j+1) = (2j+1)*x*P_j - j*P_(j-1)
	    p1 = ((2*j+1)*z*p2 - j*p3)/(j+1);
	  }
	// p1 is now the legendre polynomial evaluated at z.
	// now get the derivative pp, also by a standard recurrence relation
	pp = n * (z*p1 - p2)/(z*z - 1.0);
	z1 = z;
	z = z1 - p1/pp; // Newton's method
      } while (abs(z-z1) > EPS);


      x[i] = xm - xl*z; // scale the root to the interval
      x[n-1-i] = xm + xl*z; // and put in its symmetric counterpart
      w[i] = 2.0*xl / ((1.0 - z*z)*pp*pp); // compute the weight (NR pg. 183)
      w[n-1-i] = w[i]; // and its symmetric counterpart
    }
}

void gaulag(VecDoub &x, VecDoub &w, const double alf)
// Given alf, the legendre polynomial parameter, this returns arrays x and w
// (of length n) containing the abcissas and weights of the n-point Gauss-
// Laguerre quadrature formula. Smallest abcissa is x[0], largest x[n-1].
// (see NR pg 184)
{
  const int MAXIT = 10;
  const double EPS = 1.0e-14; // EPS is relative precision
  int i, its, j;
  double ai, p1, p2, p3, pp, z, z1;
  int n = x.size();
  for (i = 0; i < n; i++) // loop to find roots
    {
      // make initial guesses
      if (i == 0) { // initial guess for smallest root
	z = (1.0 + alf) * (3.0 + 0.92*alf) / (1.0 + 2.4*n + 1.8 * alf);
      } else if (i == 1) { // initial guess for second root
	z += (15.0 + 6.25*alf) / (1.0 + 0.9*alf + 2.5*n);
      } else { // initial guess for other roots
	ai = i - 1;
	z += ((1.0 + 2.55*ai)/(1.9*ai) + 1.26*ai*alf/(1.0 + 3.5*ai)) 
	     * 
	     (z - x[i-2])/(1.0 + 0.3*alf);
      }
      for (its = 0; its < MAXIT; its++) // refine by Newton's method
	{
	  p1 = 1.0;
	  p2 = 0.0;
	  for (j = 0; j < n; j++) //loop up recurrence relation to evalute
	    {                     //laguerre polynomial at z
	      p3 = p2;
	      p2 = p1;
	      p1 = ((2*j + 1 + alf - z)*p2 - (j + alf)*p3) / (j+1);
	    }
	  // p1 is now the desired laguerre polynomial at z. now get the
	  // derivative by recurrence as well:
	  pp = (n*p1 - (n + alf)*p2) / z;
	  z1 = z;
	  z = z1 - p1/pp; // Newton's method
	  if (abs(z - z1) <= EPS) break;
	}
      if (its >= MAXIT) std::cout << "too many iterations in gaulag" << std::endl;
      x[i] = z; //store the root and weight
      w[i] = -exp(gsl_sf_lngamma(alf + n) - gsl_sf_lngamma(double(n)))/(pp*n*p2);
    }
}

#endif
