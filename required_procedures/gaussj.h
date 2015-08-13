//gaussj.h -- routines for solving a linear system via Gauss-Jordan elimination

#ifndef __GAUSSJ_H
#define __GAUSSJ_H

#include "nr3.h"
//#include <iostream>

void gaussj(MatDoub_IO &a, MatDoub_IO &b) 
//solves a linear system via gauss-jordan elimination (see p. 41-45
//in NR). The input matrix to invert is a[0..n-1][0..n-1].
//b[0..n-1][0..m-1] contains the RHS vectors to be solved for the
//input vectors. On output, a is replaced by its inverse and b is
//replaced by the solution vectors.
{
  int i, icol, irow, j, k, l, ll, n=a.nrows(), m=b.ncols();
  double big, dum, pivinv;
  VecInt indxc(n), indxr(n), ipiv(n); //pivoting bookkeepers
  for (j = 0; j < n; j++) ipiv[j] = 0;
  for (i = 0; i < n; i++) { //loop over columns to reduce
    big = 0.0;
    for (j = 0; j < n; j++) //search for a pivot element
      if (ipiv[j] != 1)
	for (k = 0; k < n; k++) {
	  if (ipiv[k] == 0) {
	    if (abs(a[j][k]) >= big) {
	      big = abs(a[j][k]);
	      irow = j;
	      icol = k;
	      //	      std::cout << "irow = " << irow << std::endl;
	      //	      std::cout << "icol = " << icol << std::endl;
	      //	      std::cin.get();
	    }
	  }
	}
    ++(ipiv[icol]);
    //we have found the largest element, so now interchange rows to
    //put it on the diagonal. also keep track of where it is moved:
    //indxc[i] is the column of the (i+1)th pivot element, and the
    //order in which the columns are reduced. indxr[i] is the row in
    //which the (i+1)th pivot was initially located. If indxr[i] !=
    //indxc[i], a column interchange is implied. The matrix will need
    //to be unscrambled after reduction.
    if (irow != icol) {
      for (l = 0; l < n; l++) SWAP(a[irow][l],a[icol][l]);
      for (l = 0; l < m; l++) SWAP(b[irow][l],b[icol][l]);
    }
    //we may now divide the pivot row by the pivot element, located at [irow][icol]
    indxr[i]=irow;
    indxc[i]=icol;
    if (a[icol][icol] == 0.0) throw("gaussj:: SINGULAR MATRIX");
    pivinv=1.0/a[icol][icol];
    a[icol][icol] = 1.0;
    for (l = 0; l < n; l++) a[icol][l] *= pivinv;
    for (l = 0; l < m; l++) b[icol][l] *= pivinv;
    for (ll = 0; ll < n; ll++) //now reduce the rows
      if(ll != icol) {         //except the pivot row
	dum = a[ll][icol];
	a[ll][icol] = 0.0;
	for (l = 0; l < n; l++) a[ll][l] -= a[icol][l]*dum;
	for (l = 0; l < m; l++) b[ll][l] -= b[icol][l]*dum;
      }
  }
  //and now we've reduced all the columns by finding n pivots. All we
  //need do now is unscramble the columns of the solution matrix:
  for (l = n-1; l >= 0; l--) 
    if (indxr[l] != indxc[l])
      for (k = 0; k < n; k++)
	SWAP(a[k][indxr[l]],a[k][indxc[l]]);
  //and done
}

void gaussj(MatDoub_IO &a)
//overloaded gaussj so that it can be called with no RHS b mat
{
  MatDoub b(a.nrows(),0);
  gaussj(a,b);
}





#endif
