//test.cpp -- testing qsimp on some integrals

#include <vector>
#include "quadrature.h"
#include <iostream>
#include <cmath>
#include "nr3.h"
#include "ludcmp.h"

double ftest(VecDoub &r)
{
  // int dim=r.size();
  // double norm=0.0;
  // for (int i=0;i<dim;i++){ norm+=r[i]*r[i]; }
  // return sqrt(norm);
  return r[0]*cos(r[1])*r[2]*r[2];
}

int main()
{
  double r00[3]={0.0,0.0,0.0}, r11[3]={1.0,2.0,1.0};
  VecDoub r0(3,r00),r1(3,r11);

  cout << "The trap integral value is " << qlinetrap_2pt(ftest, r0, r1, 1.0e-8) << endl;
   //   cout << "The trans integral value is " << lineinttrans(ftest,r1,r2) << endl;


  //  cout << "The integral value is " << qsimp(ftest, 0.0, 1.0e10) << endl;

  
}

VecDoub line(const double &a) {
  // double r11[3]={-2.0,-1.0,0.0}, r22[3]={1.0,2.0,1.0};
  // VecDoub r1(3,r11),r2(3,r22);
  //  std::cout << "a = " << a << std::endl;

  double r[3];

  r[0]=cos(a);
  r[1]=sin(a);
  r[2]=3.0*a;
  
  // if (a<=0.0) {
  //   r[0]=a;
  //   //    std::cout<< "r[0] = " <<r[0] <<std::endl;
  //   r[1]=-1.0;
  //   r[2]=0.0;
  // } else if (a<=1.0) {
  //   r[0]=a;
  //   //    std::cout<< "r[0] = " <<r[0] <<std::endl;
  //   r[1]=a*a*a-1.0;
  //   r[2]=0.0;
  // } else if (a<=3.0) {
  //   r[0]=1.0;
  //   //    std::cout<< "r[0] = " <<r[0] <<std::endl;
  //   r[1]=a-1.0;
  //   r[2]=0.0;
  // }
  
  VecDoub rvec(3,r);
  // VecDoub rvec(3);

  // for (int i=0;i<3;i++) { rvec[i]=a*r2[i]+(1-a)*r1[i];}

  return rvec;

}


