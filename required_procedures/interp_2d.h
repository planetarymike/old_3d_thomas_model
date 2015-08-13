// interp_2d.h -- routines from NR for bilinear interpolation

#ifndef __INTERP_2D_H
#define __INTERP_2D_H

#include "nr3.h"
#include "interp_1d.h"

struct Bilin_interp
// object for bilinear interpolation on a matrix. Construct with a
// vector of x1 values, a vector of x2 values, and a matrix of
// tabulated function values yij. Then call interp for interpolated
// values.
{
  int m, n;
  const MatDoub &y;
  Linear_interp x1terp, x2terp;

  Bilin_interp(VecDoub_I &x1v, VecDoub_I &x2v, MatDoub_I &ym)
    : m(x1v.size()), n(x2v.size()), y(ym), x1terp(x1v,x1v), x2terp(x2v,x2v) {}
    // we need dummy 1 dim interp objects for their locate and hunt functions

  double interp(double x1p, double x2p) {
    int i, j;
    double yy, t, u;
    //find the grid square:
    i = x1terp.cor ? x1terp.hunt(x1p) : x1terp.locate(x1p);
    j = x2terp.cor ? x2terp.hunt(x2p) : x2terp.locate(x2p);

    //interpolate:
    t = (x1p-x1terp.xx[i])/(x1terp.xx[i+1]-x1terp.xx[i]);
    u = (x2p-x2terp.xx[i])/(x2terp.xx[i+1]-x2terp.xx[i]);
    yy = (1.-t)*(1.-u)*y[i][j] + t*(1.-u)*y[i+1][j] + (1.-t)*u*y[i][j+1] + t*u*y[i+1][j+1];
    return yy;
  }

};

#endif
