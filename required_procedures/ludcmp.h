//ludcmp.h -- routines to perform LU decomposition and solution of a matrix

#ifndef __LUDCMP_H
#define __LUDCMP_H

#include "nr3.h"

struct LUdcmp
// object for solving linear equations A.x = b using LU decomposition,
// and related functions.
{
  int n;
  MatDoub lu; // stores the decomposition
  VecInt indx; // stores the permutation (row switching)
  double d;  //used by det
  LUdcmp(MatDoub_I &a); // constructor takes matrix argument
  void solve(VecDoub_I &b, VecDoub_O &x); // solve for a single RHS
  void solve(MatDoub_I &b, MatDoub_O &x); // solve for multiple RHS
  void inverse(MatDoub_O &ainv); // calculate matrix inverse
  double det(); // return matrix determinant
  void mprove(VecDoub_I &b, VecDoub_IO &x); // iteratively improve
					    // solution
  MatDoub_I &aref; // used only by mprove
};

LUdcmp::LUdcmp(MatDoub_I &a) : n(a.nrows()), lu(a), indx(n), aref(a) {
  // Given a matrix a[0..n-1][0..n-1], this routine replaces it by the
  // LU decomposition of a rowwise permitation of itself. a is
  // input. On output, it is arranged as in NR Eqn 2.3.14;
  // indx[0..n-1] is an output vector that records the row
  // permutation; d is output as +/-1 depending on the signature of
  // the permutation. This routine is used with solve to solve linear
  // equations.
  const double TINY = 1.0e-40; // a small number
  int i, imax, j, k;
  double big, temp;
  VecDoub vv(n); // stores the scaling of each row
  d = 1.0; // no rows interchanged yet
  for (i = 0; i < n; i++){
    big = 0.0;
    for (j = 0; j < n; j++)
      if ( (temp = abs(lu[i][j])) > big ) big = temp;
    if (big == 0.0) throw("Singular Matrix in LUdcmp"); // No nonzero
							// largest
							// element
    vv[i] = 1.0/big; // save the scaling
  }
  for (k = 0; k < n; k++) {// outermost kij loop
    big = 0.0; // initialize the search for the largest pivot element
    imax = k;
    for (i = k; i < n; i++) {
      temp = vv[i]*abs(lu[i][k]);
      if (temp > big) { // is this pivot better?
	big = temp;
	imax = i;
      }
    }
    if (k != imax) { // do we need to interchange rows?
      for (j = 0; j < n; j++) { // yes we do
	temp = lu[imax][j];
	lu[imax][j] = lu[k][j];
	lu[k][j] = temp;
      }
      d = -d; // also change parity of d
      vv[imax] = vv[k]; // also interchange scale factor
    }
    indx[k] = imax;
    if (lu[k][k] == 0.0) lu[k][k] = TINY; //for some applications tiny
					  //is better than zero
    for (i = k+1; i < n; i++) {
      temp = lu[i][k] /= lu[k][k]; //divide by the pivot element
      for (j = k+1; j < n; j++)
	lu[i][j] -= temp*lu[k][j];
    }
  }
}

void LUdcmp::solve(VecDoub_I &b, VecDoub_O &x)
// solves the set of n linear equations A.x=b using the stored lu
// decomposition of A. b[0..n-1] is input as the right-hand side
// vector b, while x returns the solutions vector x; b and a may
// reference the same vector, in which case the solution overwrites
// the input. This routine takes into account the possibility that b
// will begin with many zero elements, so it is effeciant for use in
// matrix inversion.
{
  int i, ii = 0, ip, j;
  double sum;
  if (b.size() != n || x.size() != n)
    throw("LUdcmp::solve bad sizes");
  for (i = 0; i < n; i++) x[i] = b[i];
  for (i = 0; i < n; i++) {
    // When ii is set to a positive value, it will become the index of
    // the first nonvanishing element of b. We now do the forward
    // substitution, NR Eqn 2.3.6. The only new wrinkle is to
    // unscramble the permutation as we go.
    ip = indx[i];
    sum = x[ip];
    x[ip] = x[i];
    if (ii != 0)
      for (j = ii-1; j < i; j++) sum -= lu[i][j]*x[j];
    else if (sum != 0.0)
      ii = i+1;
    x[i] = sum;
  }
  for (i = n-1; i >= 0; i--) { // now backsubstitution, NR Eqn 2.3.7
    sum = x[i];
    for (j = i+1; j < n; j++) sum -= lu[i][j]*x[j];
    x[i]=sum/lu[i][i];
  }
}

void LUdcmp::solve(MatDoub_I &b, MatDoub_O &x)
// Solves m sets of n linear equations analagously to solve above. a
// and b may reference the same matrix, in which case the solution
// overwrites the input
{
  int i, j, m = b.ncols();
  if (b.nrows() != n || x.nrows() != n || b.ncols() != x.ncols())
    throw("LUdcmp::solve bad sizes");
  VecDoub xx(n);
  for (j = 0; j < m; j++) { // copy and solve each column in turn
    for (i = 0; i < n; i++) xx[i] = b[i][j];
    solve(xx,xx);
    for (i = 0; i < n; i++) x[i][j] = xx[i];
  }
}


void LUdcmp::inverse(MatDoub_O &ainv)
// using the stored LU decomposition, return the inverse matrix
{
  int i, j;
  ainv.resize(n,n);
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) ainv[i][j] = 0.0;
    ainv[i][i] = 1.0;
  }
  solve(ainv,ainv);
}

double LUdcmp::det()
// using the stored LU decomposition, return the matrix determinant
{
  double dd = d;
  for (int i = 0; i < n; i++) dd *= lu[i][i];
  return dd;
}

void LUdcmp::mprove(VecDoub_I &b, VecDoub_IO &x)
// improves a solution vector x[0..n-1] of the linear set of equations
// A.x = b. the vectors b[0..n-1] and x[0..n-1] are input. On output,
// x[0..n-1] is modified to an improved set of values
{
  int i, j;
  VecDoub r(n);
  for (i = 0; i < n; i++) { // calculate the RHS of NR Eqn 2.5.4
    long double sdp = -b[i];
    for (j = 0; j < n; j++)
      sdp += (long double) aref[i][j] * (long double) x[j];
    r[i] = sdp;
  }
  solve(r,r); // solve for the error term
  for (i = 0; i < n; i++) x[i] -= r[i];
}


#endif
