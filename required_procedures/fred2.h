// fred2.h -- routine for solving a one-dimensional fredholm equation of the seconds kind

// requires external functions to provide g(t) and lambda*Kern(ij)
// doub g(const double t) { ... }
// doub ak(const double t, const doub s) { ... }
//
// instantitiate like
// Fred2<double (double), double (double, double)> fred2(a,b,n,g,ak);
//
// call like
// double ans = fred2.fredin(x)

#include "nr3.h"
#include "ludcmp.h"

template <class G, class K>
struct Fred2
// solves a linear fredholm equation of the second kind
{
  const double a, b;
  const int n;
  G &g;
  K &ak;
  VecDoub t, f, w;
  Fred2(const double aa, const double bb, const int nn, G &gg, K &akk) :
    a(aa), b(bb), n(nn), g(gg), ak(akk), t(n), f(n), w(n)
    // quantities a and b are input as limits of integration. The
    // quantity n is the number of points to use in the Gaussian
    // quadrature. g and ak are user-supplied functions or functors
    // that return g(t) and lambda*K(t,s). The constructor computes
    // arrays t[0..n-1] and f[0..n-1] containing the abcissas t_i of
    // the quadrature and the solution f at those abcissas. Also
    // computed is the array w[0..n-1] of Gaussian weights for the
    // Nystrom interpolation routine fredin.
  {
    MatDoub omk(n,n);
    gauleg(a,b,t,w);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++)
	omk[i][j] = double(i == j) - ak(t[i],t[j])*w[j];
      f[i] = g(t[i]);
    }
    LUdcmp alu(omk);
    alu.solve(f,f);
  }

  double fredin(const double x)
  // Given arrays t[0..n-1] and w[0..n-1] containing the abcissas and
  // weights of the Gaussian quadrature, and given the solution array
  // f[0..n-1], this function returns the value of the function f at x
  // using the Nystrom interpolation formula.
  {
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += ak(x,t[i])*w[i]*f[i];
    return g(x) + sum;
  }
};
