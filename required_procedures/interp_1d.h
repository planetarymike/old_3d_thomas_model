//interp_1d.h -- base class for interpolation routines
#include <cmath>
#include "nr3.h"

#ifndef __Interp_H_
#define __Interp_H_

using std::pow;

struct Base_interp
// Abstract base class used by all interpolation routines in this
// chapter. Only interp is called by the user.
{
  int n, mmm, jsav, cor, dj;
  double *xx, *yy;
  Base_interp(VecDoub_I &x, VecDoub_I &y, int m)
  // constructor. set up for interpolating on a table of x's and y's
  // of length m. Normally called by derived class, not user.
    : n(x.size()), mmm(m), jsav(0), cor(0) {
    if (n==0) {
      xx = NULL;
      yy = NULL;
    } else {
      xx = new double[n];
      yy = new double[n];
    }
    for (int i=0;i<n;i++) {
      xx[i]=x[i];
      yy[i]=y[i];
    }
    
    dj = NR_MIN(1, (int) pow((double) n, 0.25) );
  }

  Base_interp(const Base_interp &B) {
    n=B.n;
    mmm=B.mmm;
    jsav=B.jsav;
    cor=B.cor;
    dj=B.dj;
    if (n==0) {
      xx = NULL;
      yy = NULL;
    } else {
      xx = new double[n];
      yy = new double[n];
    }
    for (int i=0; i<n; i++) {
      xx[i]=B.xx[i];
      yy[i]=B.yy[i];
    }
  }

  Base_interp&
  operator= (const Base_interp &B) {
    n=B.n;
    mmm=B.mmm;
    jsav=B.jsav;
    cor=B.cor;
    dj=B.dj;
    if (n==0) {
      xx = NULL;
      yy = NULL;
    } else {
      xx = new double[n];
      yy = new double[n];
    }
    for (int i=0; i<n; i++) {
      xx[i]=B.xx[i];
      yy[i]=B.yy[i];
    }
  }

  ~Base_interp() {
    if (n>0) {
      delete [] xx;
      delete [] yy;
    }
  }


  double interp(double x)
  // Given a value x, return an interpolated value, using data
  // pointed to by xx and yy.
  {
    int jlo = cor ? hunt(x) : locate (x);
    return rawinterp(jlo, x);
  }

  int index(double x)
  // Given a value x, return an interpolated value, using data
  // pointed to by xx and yy.
  {
    int jlo = cor ? hunt(x) : locate(x);
    return jlo;
  }

  int locate(const double x);
  int hunt(const double x);

  virtual double rawinterp(int jlo, double x) = 0;
  // derived classes provide this as the actual interpolation method used.

};

int Base_interp::locate(const double x)
// given a value x, return a value j such that x is (insofar as
// possible) centered in the subrange xx[j..j+mmm-1], where xx is the
// stored pointer. The values in xx must be monotonic, either
// increasing or decreasing. The returned value is not less than 0,
// nor greater than n-1).
{
  int ju, jm, jl;
  //  std::cout << "n = " << n << ", mmm = " << mmm << "\n";
  if (n < 2 || mmm < 2 || mmm > n) {
    //    std::cout << "n = " << n << ", mmm = " << mmm << std::endl;
    throw("locate size error");
  }

  bool ascnd = (xx[n-1] >= xx[0]); // true if ascending order of
				   // table, false otherwise
  jl = 0;//lower and upper size limits
  ju = n-1;
  while (ju-jl > 1) {  // until the appropriate range is found
    jm = (ju+jl) >> 1; // find the midpoint
    if (x >= xx[jm] == ascnd)
      jl = jm; // and replace either the lower limit
    else
      ju = jm; // or the upper limit, as appropriate
  }
  cor = abs(jl-jsav) > dj ? 0 : 1; // decide whether to use hunt() or
				   // locate() next time
  jsav = jl;
  return NR_MAX(0,NR_MIN(n-mmm,jl-((mmm-2)>>1)));
}


int Base_interp::hunt(const double x)
// given a value x, return a value j such that x is (insofar as
// possible) centered in the subrange xx[j..j+mmm-1], where xx is the
// stored pointer. The values in xx must be monotonic, either
// increasing or decreasing. The returned value is not less than 0,
// nor greater than n-1).
{
  int jl = jsav, jm, ju, inc = 1;
  //  std::cout << "n = " << n << ", mmm = " << mmm << std::endl;
  if (n < 2 || mmm < 2 || mmm > n) {
    std::cout << "n = " << n << ", mmm = " << mmm << std::endl;
    throw("hunt size error");
  }
  bool ascnd = (xx[n-1] >= xx[0]); // does table ascend?
  if (jl < 0 || jl > n-1) { // input guess not useful. go directly to
			    // bisection
    jl = 0;
    ju = n-1;
  } else {
    if (x >= xx[jl] == ascnd) { // hunt up:
      for(;;) {
	ju = jl + inc;
	if (ju >= n-1) { ju = n-1; break; } // off end of table
	else if (x < xx[ju] == ascnd) break; // found bracket
	else {
	  jl = ju;
	  inc += inc;
	}
      }
    } else { // hunt down:
      ju = jl;
      for (;;) {
	jl = jl - inc;
	if (jl <= 0) {jl = 0; break;} // off end of table
	else if (x >= xx[jl] == ascnd) break; // found bracket
	else {
	  ju = jl;
	  inc += inc;
	}
      }
    }
  }
  while (ju - jl > 1) { // hunt is done, so begin final bisection
    jm = (ju+jl) >> 1;
    if (x >= xx[jm] == ascnd)
      jl = jm;
    else
      ju = jm;
  }
  cor = abs(jl - jsav) > dj ? 0 : 1; // decide whether to use hunt or
				     // locate next time
  jsav = jl;
  return NR_MAX(0,NR_MIN(n-mmm,jl-((mmm-2)>>1)));
}

struct Linear_interp : Base_interp
// piecewise linear interpolation object. Construct with x and y
// vectors, then call interp for interpolated values
{
  Linear_interp(VecDoub_I &xv, VecDoub_I &yv) : Base_interp(xv,yv,2) {}
  Linear_interp() : Base_interp(VecDoub(0),VecDoub(0),2) {}
  
  double rawinterp(int j, double x) {
    if(xx[j] == xx[j+1]) return yy[j]; // table defective, but recover
    else return yy[j] + ((x-xx[j])/(xx[j+1]-xx[j]))*(yy[j+1]-yy[j]);
  }
};

struct Poly_interp : Base_interp
// Polynomial interpolation object. Construct with x and y vectors and
// the number M of points to be used locally (polynomial order plus
// one), then call interp for interpolated values.
{
  double dy;
  Poly_interp(VecDoub_I &xv, VecDoub_I &yv, int m)
    : Base_interp(xv, yv,m), dy(0.0) {}
  double rawinterp(int jl, double x);
};

double Poly_interp::rawinterp(int jl, double x)
// Given a value x, and using pointers to dara xx and yy, this routine
// returns an interpolated value y, and stores an error estimate
// dy. The returned value is obtained by mmm-point polynomial
// interpolation on the subrange xx[jl..jl+mmm-1]
{
  int i, m, ns = 0;
  double y, den, dif, dift, ho, hp, w;
  const double *xa = &xx[jl], *ya = &yy[jl];
  VecDoub c(mmm), d(mmm);
  dif = abs(x-xa[0]);
  for (i = 0; i < mmm; i++) { // find the index ns of the closest table
			     // entry
    if ((dift = abs(x-xa[i])) < dif) {
      ns = i;
      dif = dift;
    }
    c[i] = ya[i];
    d[i] = ya[i];
  }
  y = ya[ns--]; // initial approximataion to y
  for (m = 1; m < mmm; m++) { // for each column of the tableau
    for (i = 0; i < mmm-m; i++) { // we loop over the current C's and
				 // D's and update them
      ho = xa[i]-x;
      hp = xa[i+m]-x;
      w = c[i+1]-d[i];
      if ((den=ho-hp) == 0.0) throw("Poly_interp error");
      den = w/den;
      d[i] = hp*den;
      c[i] = ho*den;
    }
    y += (dy = (2*(ns+1) < (mmm-m) ? c[ns+1] : d[ns--]));
    // After each column in the tableau is completed, we secide which
    // correction, c or d, we want to add to our accumulating value od
    // y, i.e. which path to take throught the tableau --- forking up
    // or down. We do this in such a way as to take the most "straight
    // line" path through the tableau to its apex, updating ns
    // accordingly to keep track of where we are. This route keeps the
    // partial approximations centered (insofar as is possible) on the
    // target x. The last dy added is thus the error indication.
  }
  return y;
}

/* struct Poly_log_interp : Poly_interp { */
/*   // same thing as poly_interp, but in the log-log space */
/*   Poly_log_interp(VecDoub_I &xv, VecDoub_I &yv, int m) */
/*     : Poly_interp(xv, yv,m) {} */
  
/*   double rawinterp(int jl, double x) */
/*   // Given a value x, and using pointers to dara xx and yy, this routine */
/*   // returns an interpolated value y, and stores an error estimate */
/*   // dy. The returned value is obtained by mmm-point polynomial */
/*   // interpolation on the subrange xx[jl..jl+mmm-1] */
/*   { */
/*     int i, m, ns = 0; */
/*     double y, den, dif, dift, ho, hp, w; */
/*     const double *xa = &xx[jl], *ya = &yy[jl]; */
/*     VecDoub c(mmm), d(mmm); */
/*     dif = abs(x-xa[0]); */
/*     for (i = 0; i < mmm; i++) { // find the index ns of the closest table */
/*       // entry */
/*       if ((dift = abs(x=xa[i])) < dif) { */
/* 	ns = i; */
/* 	dif = dift; */
/*       } */
/*       c[i] = ya[i]; */
/*       d[i] = ya[i]; */
/*     } */
/*     y = ya[ns--]; // initial approximataion to y */
/*     for (m = 1; m < mmm; m++) { // for each column of the tableau */
/*       for (i = 0; i < mmm-m; i++) { // we loop over the current C's and */
/* 	// D's and update them */
/* 	ho = xa[i]-x; */
/* 	hp = xa[i+m]-x; */
/* 	w = c[i+1]-d[i]; */
/* 	if ((den=ho-hp) == 0.0) throw("Poly_interp error"); */
/* 	den = w/den; */
/* 	d[i] = hp*den; */
/* 	c[i] = ho*den; */
/*       } */
/*       y += (dy = (2*(ns+1) < (mmm-m) ? c[ns+1] : d[ns--])); */
/*       // After each column in the tableau is completed, we secide which */
/*       // correction, c or d, we want to add to our accumulating value od */
/*       // y, i.e. which path to take throught the tableau --- forking up */
/*       // or down. We do this in such a way as to take the most "straight */
/*       // line" path through the tableau to its apex, updating ns */
/*       // accordingly to keep track of where we are. This route keeps the */
/*       // partial approximations centered (insofar as is possible) on the */
/*       // target x. The last dy added is thus the error indication. */
/*     } */
/*     return std::exp(y); */

/*   } */
  
/*   double interp(double x) */
/*   // Given a value x, return an interpolated value, using data */
/*   // pointed to by xx and yy. */
/*   { */
/*     x = std::log(x); */
/*     int jlo = cor ? hunt(x) : locate (x); */
/*     return rawinterp(jlo, x); */
/*   } */
  
/* } */



#endif
