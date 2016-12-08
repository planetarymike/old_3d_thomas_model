//bicubic.h -- 2D cubic interpolation on an irregular grid.

#ifndef __BICUBIC_H
#define __BICUBIC_H

#include <cmath>
#include "nr3.h"
#include "interp_1d.h"

double bicubic_iobs(const double x, const double y, 
		    Linear_interp& xgrid,
		    Linear_interp& ygrid,
		    VecDoub** F,
		    int iobs) {
      // find the grid square:
      int ix = xgrid.cor ? xgrid.hunt(x) : xgrid.locate(x);
      int nx = xgrid.n;
      int iy = ygrid.cor ? ygrid.hunt(y) : ygrid.locate(y);
      int ny = ygrid.n;

      //get the fractional coordinate values
      double ax, ay;
      ax = (x-xgrid.xx[ix])/(xgrid.xx[ix+1]-xgrid.xx[ix]);
      ay = (y-ygrid.xx[iy])/(ygrid.xx[iy+1]-ygrid.xx[iy]);      

      //to do bicubic interpolation we need the function values and
      //derivatives at the four corners of the box. 
      double f11, f12, f21, f22,
	f11x, f12x, f21x, f22x,
	f11y, f12y, f21y, f22y,
	f11xy, f12xy, f21xy, f22xy;
      f11=F[ix  ][iy  ][iobs];
      f12=F[ix  ][iy+1][iobs];
      f21=F[ix+1][iy  ][iobs];
      f22=F[ix+1][iy+1][iobs];
      //the tricky part is estimating the derivatives, in part because
      //we need to keep track of whether we're at the edge of the
      //parameter space

      //get the surrounding points for the purpose of calculating the
      //derivatives
      double f00, f01, f02, f03;
      double f10, f13;
      double f20, f23;
      double f30, f31, f32, f33;
      double f10x, f13x, f20x, f23x;
      double dx20, dx31;
      double dy20, dy31;

      //first the corner cases
      if (ix==0&&iy==0) {
	//get the values
	f00=F[ix  ][iy  ][iobs];
	f01=F[ix  ][iy  ][iobs];
	f02=F[ix  ][iy+1][iobs];
	f03=F[ix  ][iy+2][iobs];
	f10=F[ix  ][iy  ][iobs];
	f13=F[ix  ][iy+2][iobs];
	f20=F[ix+1][iy  ][iobs];
	f23=F[ix+1][iy+2][iobs];
	f30=F[ix+2][iy  ][iobs];
	f31=F[ix+2][iy  ][iobs];
	f32=F[ix+2][iy+1][iobs];
	f33=F[ix+2][iy+2][iobs];	

	//grid spacings
	dx20=xgrid.xx[ix+1]-xgrid.xx[ix  ];
	dx31=xgrid.xx[ix+2]-xgrid.xx[ix  ];
	dy20=ygrid.xx[iy+1]-ygrid.xx[iy  ];
	dy31=ygrid.xx[iy+2]-ygrid.xx[iy  ];
      } else if (ix==0&&iy==ny-2) {
	//get the values
	f00=F[ix  ][iy-1][iobs];
	f01=F[ix  ][iy  ][iobs];
	f02=F[ix  ][iy+1][iobs];
	f03=F[ix  ][iy+1][iobs];
	f10=F[ix  ][iy-1][iobs];
	f13=F[ix  ][iy+1][iobs];
	f20=F[ix+1][iy-1][iobs];
	f23=F[ix+1][iy+1][iobs];
	f30=F[ix+2][iy-1][iobs];
	f31=F[ix+2][iy  ][iobs];
	f32=F[ix+2][iy+1][iobs];
	f33=F[ix+2][iy+1][iobs];	

	//grid spacings
	dx20=xgrid.xx[ix+1]-xgrid.xx[ix  ];
	dx31=xgrid.xx[ix+2]-xgrid.xx[ix  ];
	dy20=ygrid.xx[iy+1]-ygrid.xx[iy-1];
	dy31=ygrid.xx[iy+1]-ygrid.xx[iy  ];
      } else if (ix==nx-2&&iy==0) {
	//get the values
	f00=F[ix-1][iy  ][iobs];
	f01=F[ix-1][iy  ][iobs];
	f02=F[ix-1][iy+1][iobs];
	f03=F[ix-1][iy+2][iobs];
	f10=F[ix  ][iy  ][iobs];
	f13=F[ix  ][iy+2][iobs];
	f20=F[ix+1][iy  ][iobs];
	f23=F[ix+1][iy+2][iobs];
	f30=F[ix+1][iy  ][iobs];
	f31=F[ix+1][iy  ][iobs];
	f32=F[ix+1][iy+1][iobs];
	f33=F[ix+1][iy+2][iobs];	

	//grid spacings
	dx20=xgrid.xx[ix+1]-xgrid.xx[ix-1];
	dx31=xgrid.xx[ix+1]-xgrid.xx[ix  ];
	dy20=ygrid.xx[iy+1]-ygrid.xx[iy  ];
	dy31=ygrid.xx[iy+2]-ygrid.xx[iy  ];
      } else if (ix==nx-2&&iy==ny-2) {
	//get the values
	f00=F[ix-1][iy-1][iobs];
	f01=F[ix-1][iy  ][iobs];
	f02=F[ix-1][iy+1][iobs];
	f03=F[ix-1][iy+1][iobs];
	f10=F[ix  ][iy-1][iobs];
	f13=F[ix  ][iy+1][iobs];
	f20=F[ix+1][iy-1][iobs];
	f23=F[ix+1][iy+1][iobs];
	f30=F[ix+1][iy-1][iobs];
	f31=F[ix+1][iy  ][iobs];
	f32=F[ix+1][iy+1][iobs];
	f33=F[ix+1][iy+1][iobs];	

	//grid spacings
	dx20=xgrid.xx[ix+1]-xgrid.xx[ix-1];
	dx31=xgrid.xx[ix+1]-xgrid.xx[ix  ];
	dy20=ygrid.xx[iy+1]-ygrid.xx[iy-1];
	dy31=ygrid.xx[iy+1]-ygrid.xx[iy  ];

	//then the edge cases
      } else if (ix==0) {     
	//get the values
	f00=F[ix  ][iy-1][iobs];
	f01=F[ix  ][iy  ][iobs];
	f02=F[ix  ][iy+1][iobs];
	f03=F[ix  ][iy+2][iobs];
	f10=F[ix  ][iy-1][iobs];
	f13=F[ix  ][iy+2][iobs];
	f20=F[ix+1][iy-1][iobs];
	f23=F[ix+1][iy+2][iobs];
	f30=F[ix+2][iy-1][iobs];
	f31=F[ix+2][iy  ][iobs];
	f32=F[ix+2][iy+1][iobs];
	f33=F[ix+2][iy+2][iobs];	

	//grid spacings
	dx20=xgrid.xx[ix+1]-xgrid.xx[ix  ];
	dx31=xgrid.xx[ix+2]-xgrid.xx[ix  ];
	dy20=ygrid.xx[iy+1]-ygrid.xx[iy-1];
	dy31=ygrid.xx[iy+2]-ygrid.xx[iy  ];
      } else if (ix==nx-2) {
	//get the values
	f00=F[ix-1][iy-1][iobs];
	f01=F[ix-1][iy  ][iobs];
	f02=F[ix-1][iy+1][iobs];
	f03=F[ix-1][iy+2][iobs];
	f10=F[ix  ][iy-1][iobs];
	f13=F[ix  ][iy+2][iobs];
	f20=F[ix+1][iy-1][iobs];
	f23=F[ix+1][iy+2][iobs];
	f30=F[ix+1][iy-1][iobs];
	f31=F[ix+1][iy  ][iobs];
	f32=F[ix+1][iy+1][iobs];
	f33=F[ix+1][iy+2][iobs];	

	//grid spacings
	dx20=xgrid.xx[ix+1]-xgrid.xx[ix-1];
	dx31=xgrid.xx[ix+1]-xgrid.xx[ix  ];
	dy20=ygrid.xx[iy+1]-ygrid.xx[iy-1];
	dy31=ygrid.xx[iy+2]-ygrid.xx[iy  ];
      } else if (iy==0) {
	//get the values
	f00=F[ix-1][iy  ][iobs];
	f01=F[ix-1][iy  ][iobs];
	f02=F[ix-1][iy+1][iobs];
	f03=F[ix-1][iy+2][iobs];
	f10=F[ix  ][iy  ][iobs];
	f13=F[ix  ][iy+2][iobs];
	f20=F[ix+1][iy  ][iobs];
	f23=F[ix+1][iy+2][iobs];
	f30=F[ix+2][iy  ][iobs];
	f31=F[ix+2][iy  ][iobs];
	f32=F[ix+2][iy+1][iobs];
	f33=F[ix+2][iy+2][iobs];	

	//grid spacings
	dx20=xgrid.xx[ix+1]-xgrid.xx[ix-1];
	dx31=xgrid.xx[ix+2]-xgrid.xx[ix  ];
	dy20=ygrid.xx[iy+1]-ygrid.xx[iy  ];
	dy31=ygrid.xx[iy+2]-ygrid.xx[iy  ];
      } else if (iy==ny-2) {
	//get the values
	f00=F[ix-1][iy-1][iobs];
	f01=F[ix-1][iy  ][iobs];
	f02=F[ix-1][iy+1][iobs];
	f03=F[ix-1][iy+1][iobs];
	f10=F[ix  ][iy-1][iobs];
	f13=F[ix  ][iy+1][iobs];
	f20=F[ix+1][iy-1][iobs];
	f23=F[ix+1][iy+1][iobs];
	f30=F[ix+2][iy-1][iobs];
	f31=F[ix+2][iy  ][iobs];
	f32=F[ix+2][iy+1][iobs];
	f33=F[ix+2][iy+1][iobs];	

	//grid spacings
	dx20=xgrid.xx[ix+1]-xgrid.xx[ix-1];
	dx31=xgrid.xx[ix+2]-xgrid.xx[ix  ];
	dy20=ygrid.xx[iy+1]-ygrid.xx[iy-1];
	dy31=ygrid.xx[iy+1]-ygrid.xx[iy  ];

	//then the interior
      } else if (ix>0&&ix<nx-2&&iy>0&&iy<ny-2) {      
	//get the values
	f00=F[ix-1][iy-1][iobs];
	f01=F[ix-1][iy  ][iobs];
	f02=F[ix-1][iy+1][iobs];
	f03=F[ix-1][iy+2][iobs];
	f10=F[ix  ][iy-1][iobs];
	f13=F[ix  ][iy+2][iobs];
	f20=F[ix+1][iy-1][iobs];
	f23=F[ix+1][iy+2][iobs];
	f30=F[ix+2][iy-1][iobs];
	f31=F[ix+2][iy  ][iobs];
	f32=F[ix+2][iy+1][iobs];
	f33=F[ix+2][iy+2][iobs];	

	//grid spacings
	dx20=xgrid.xx[ix+1]-xgrid.xx[ix-1];
	dx31=xgrid.xx[ix+2]-xgrid.xx[ix  ];
	dy20=ygrid.xx[iy+1]-ygrid.xx[iy-1];
	dy31=ygrid.xx[iy+2]-ygrid.xx[iy  ];
      } else {
	throw("can't locate x, y grid point in bicubic_iobs.");
      }

      //x derivatives
      f11x=(f21-f01)/dx20;
      f12x=(f22-f02)/dx20;
      f21x=(f31-f11)/dx31;
      f22x=(f32-f12)/dx31;

      //y derivatives
      f11y=(f12-f10)/dy20;
      f12y=(f13-f11)/dy31;
      f21y=(f22-f20)/dy20;
      f22y=(f23-f21)/dy31;

      //cross derivatives
      f10x=(f20-f00)/dx20;
      f13x=(f23-f03)/dx20;
      f20x=(f30-f10)/dx31;
      f23x=(f33-f13)/dx31;

      f11xy=(f12x-f10x)/dy20;
      f12xy=(f13x-f11x)/dy31;
      f21xy=(f22x-f20x)/dy20;
      f22xy=(f23x-f21x)/dy31;

      //now we can do the interpolation, which amounts to some matrix
      //operations to compute the coefficients and interpolated value
      double coefmat1[4][4] =
	{
	  { 1.0, 0.0, 0.0, 0.0},
	  { 0.0, 0.0, 1.0, 0.0},
	  {-3.0, 3.0,-2.0,-1.0},
	  { 2.0,-2.0, 1.0, 1.0}
	};
      double coefmat2[4][4] =
	{
	  { 1.0, 0.0,-3.0, 2.0},
	  { 0.0, 0.0, 3.0,-2.0},
	  { 0.0, 1.0,-2.0, 1.0},
	  { 0.0, 0.0,-1.0, 1.0}
	};
      double fmat[4][4] =
	{
	  {f11  , f12  , f11y , f12y },
	  {f21  , f22  , f21y , f22y },
	  {f11x , f12x , f11xy, f12xy},
	  {f21x , f22x , f21xy, f22xy}
	};

      double tempcoefmat[4][4];
      double coefmat[4][4];
      for (int ia=0;ia<4;ia++) {
	for (int ib=0;ib<4;ib++) {
	  tempcoefmat[ia][ib]=0.0;
	  coefmat[ia][ib]=0.0;
	}
      }

      for (int ia=0;ia<4;ia++) 
	for (int ib=0;ib<4;ib++)
	  for (int ic=0;ic<4;ic++)
	    tempcoefmat[ia][ib]+=fmat[ia][ic]*coefmat2[ic][ib];
      


      for (int ia=0;ia<4;ia++) 
	for (int ib=0;ib<4;ib++)
	  for (int ic=0;ic<4;ic++)
	    coefmat[ia][ib]+=coefmat1[ia][ic]*tempcoefmat[ic][ib];

      //now the interpolation coefficients are properly stored in
      //coefmat, so we can interpolate
      double retval=
	 coefmat[0][0]
	+coefmat[0][1]*ay
	+coefmat[0][2]*ay*ay
	+coefmat[0][3]*ay*ay*ay
	+ax*(
	     coefmat[1][0]
	     +coefmat[1][1]*ay
	     +coefmat[1][2]*ay*ay
	     +coefmat[1][3]*ay*ay*ay)
	+ax*ax*(
		coefmat[2][0]
		+coefmat[2][1]*ay
		+coefmat[2][2]*ay*ay
		+coefmat[2][3]*ay*ay*ay)
	+ax*ax*ax*(
		   coefmat[3][0]
		   +coefmat[3][1]*ay
		   +coefmat[3][2]*ay*ay
		   +coefmat[3][3]*ay*ay*ay);

return retval;
}


double bicubic_array(const double x, const double y, 
		     Linear_interp& xgrid,
		     Linear_interp& ygrid,
		     double** F) {
  VecDoub** myf;
  myf = new VecDoub*[xgrid.n];
  for (int i=0;i<xgrid.n;i++) {
    myf[i]=new VecDoub[ygrid.n];
    for (int j=0;j<ygrid.n;j++) {
      myf[i][j].resize(1);
      myf[i][j][0]=F[i][j];
    }
  }
  
  double retval=bicubic_iobs(x,y,xgrid,ygrid,myf,0);
  
  return retval;
}


#endif
