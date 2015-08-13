//test.cpp -- testing qsimp on some integrals

#include <vector>
#include "quadrature.h"
#include <iostream>
#include <cmath>
#include "nr3.h"
#include "ludcmp.h"

// double ftest(double r)
// {
//   return r;
// }

int main()
{
  using std::cout; 
  using std::endl;

  const double mymat[9] = {1.0,10.0,3.0,1.0,-1.2,4.5,6.7,9.3,-10.0};
  
  MatDoub mat(3,3,mymat);

  cout << "Matrix: "<<endl;
  int row, col;
  cout.precision(5);
  for (row = 0; row < 3; row++){
    for (col = 0; col < 3; col++) {
      cout.width(10);
      cout << mat[row][col];
    }
    cout << endl;
  }
  
  LUdcmp sol(mat);

  cout << "Its LU decomposition: " << endl;
  for (row = 0; row < 3; row++){
    for (col = 0; col < 3; col++) {
      cout.width(10);
      cout << sol.lu[row][col];
    }
    cout << endl;
  }

  cout << "Its determinant: " << sol.det() << endl;

  VecDoub vsol(3);
  double myb[3] = {0.0,1.0,-1.0};
  VecDoub b(3,myb);
  cout << "Solution for rhs vector: {" << b[0];
  for (row = 1; row < 3; row++) cout << ", " << b[row];
  cout << "}" << endl;
  sol.solve(b,vsol);

  for (row = 0; row < 3; row++){
    cout.width(10);
    cout << vsol[row];
    cout << endl;
  }

  

  MatDoub inv(3,3);
  sol.inverse(inv);
  cout << "Its inverse: " << endl;
  for (row = 0; row < 3; row++){
    for (col = 0; col < 3; col++) {
      cout.width(10);
      cout << inv[row][col];
    }
    cout << endl;
  }


  
  //  const int dim = 3;

  //  double myvec1[dim] = {0.0,1.0,-1.0};
  //  double myvec2[dim] = {0.0,1.0,10.0};
  //  VecDoub r1(dim,myvec1), r2(dim,myvec2);

  //  cout << "The trap integral value is " << lineqtrap(ftest, r1, r2) << endl;
  //  cout << "The simp integral value is " << lineqsimp(ftest, r1, r2) << endl;
  //  //   cout << "The trans integral value is " << lineinttrans(ftest,r1,r2) << endl;


  // //  cout << "The integral value is " << qsimp(ftest, 0.0, 1.0e10) << endl;

  
}
