//iph_model_interface.cpp -- routine to call Quemerais fortran IPH model from C++

#include <cmath>

using std::cos;
using std::sin;


//fortran function prototype
extern "C" {
  void background_(float *lc,
		   float *x_pos,float *y_pos,float *z_pos,
		   float *x_look,float *y_look,float *z_look,
		   float *iph_b);
}

double quemerais_iph_model(double Fsun, VecDoub marspos, double ra, double dec) {
  //this routine computes the ecliptic look direction corresponding to
  //a given RA, dec, and feeds the appropriate parameters to the
  //Quemerais code, which returns a computed brightness.
  dec=pi/180*dec;
  ra=pi/180*ra;

  VecDoub j2000vec(3);
  j2000vec[0]=cos(dec)*cos(ra);
  j2000vec[1]=cos(dec)*sin(ra);
  j2000vec[2]=sin(dec);

  //ecliptic coordinates are related to J2000 coordinates via a
  //negative rotation about the x-axis equal to the obliquity of the
  //Earth.
  double eob=pi/180.*23.44;
  VecDoub eclipvec(3);
  eclipvec[0]=j2000vec[0];
  eclipvec[1]=j2000vec[1]*cos(-eob)-j2000vec[2]*sin(-eob);
  eclipvec[2]=j2000vec[2]*cos(-eob)+j2000vec[1]*sin(-eob);

  double iphb = 0.0;
  //fortran swap variables
  float lc_=Fsun;
  float x_pos_=marspos[0];
  float y_pos_=marspos[1];
  float z_pos_=marspos[2];
  float x_look_=eclipvec[0];
  float y_look_=eclipvec[1];
  float z_look_=eclipvec[2];
  float iphb_=iphb;

  /* std::cout << "float lc_=" <<lc_ << std::endl; */
  /* std::cout << "float x_pos_="<<x_pos_ << std::endl; */
  /* std::cout << "float y_pos_="<<y_pos_ << std::endl; */
  /* std::cout << "float z_pos_="<<z_pos_ << std::endl; */
  /* std::cout << "float x_look_="<<x_look_ << std::endl; */
  /* std::cout << "float y_look_="<<y_look_ << std::endl; */
  /* std::cout << "float z_look_="<<z_look_ << std::endl; */
  /* std::cout << "float iphb_="<<iphb_ << std::endl; */
  
  background_(&lc_,&x_pos_,&y_pos_,&z_pos_,&x_look_,&y_look_,&z_look_, &iphb_);

  /* std::cout << "float iphb_="<<iphb_ << std::endl; */

  iphb = (double) iphb_;

  return iphb/1000.; //iphb is in R, convert to kR
}

