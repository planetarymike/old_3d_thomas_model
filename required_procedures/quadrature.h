//quadrature.h -- header file containing quadrature stuff from NR

#ifndef __QUADRATURE_NR_H
#define __QUADRATURE_NR_H

#include <vector>
#include <cmath>
#include "nr3.h"
#include <iostream>

//abstract base class for various quadrature algorithms
struct Quadrature{
  int n; //current level of refinement
  virtual double next(void) = 0; // returns the value of the integral at stage n
                                 // next() must be defined in the derived class
};

//Trapezoidal rule
template <class T>
struct Trapzd : Quadrature 
{
  double a, b, s; //limits of integration and current integral value
  T &func;
  //  Trapzd(...) {};
  Trapzd(T &funcc, const double aa, const double bb) : func(funcc), a(aa), b(bb) {n = 0;}
  double next() {
    double x, tnm, sum, del;
    int it, j;
    n++;
    if (n == 1) 
    	return (s = 0.5*(b-a)*(func(a)+func(b))); // initialize the integral value
    else 
      {
	for (it = 1, j = 1; j < n - 1; j++) it <<= 1; // i = number of sample points inside interval
	tnm = it;
	del = (b-a)/tnm; // width of each region
	x = a + 0.5 * del; // set x at the midpoint of the first interval
	for (sum = 0.0, j = 0; j < it; j++, x+=del) sum += func(x);
	s = 0.5*(s+(b-a)*sum/tnm);
	return s;
      }
  }
};


//midpoint rule for singularities on either end of the integral
template <class T>
struct Midpnt : Quadrature
{
  double a, b, s; // endpoints and current value of integral
  T &funk;
  Midpnt(T &funcc, const double aa, const double bb) :  // constructor takes endpoints 
    funk(funcc), a(aa), b(bb) {n=0;}                    // and function or functor to 
  double next()                                         // be integrated
  // the iterator. returns the nth stage of refinement; 
  // number of points is tripled with each call
  {
    int it, j;
    double x, tnm, sum, del, ddel;
    n++;
    if (n == 1) 
      return (s = (b-a)*func(0.5*(a+b)));
    else
      {
	for (it = 1, j = 1; j < n-1; j++) it *= 3;
	tnm = it;
	del = (b-a)/(3.0*tnm); // added points alternate in spacing between del and ddel
	ddel= del + del;
	x = a + 0.5*del;
	sum = 0.0;
	for (j = 0; j < it; j++)
	  {
	    sum += func(x);
	    x += ddel;
	    sum += func(x);
	    x += del;
	  }
	s = (s + (b-a)/tnm*sum)/3.0; // combine new sum with old to get refined value
	return s;
      }
  }
  virtual double func(const double x) {return funk(x);} 
  // identity mapping here generalizes to allow non-identity mappings 
  // for other kinds of mapped quadrature rules.
};

template <class T>
struct Midinf : Midpnt<T>
// exact replacement for Midpnt, except that now we evaluate the function 
// at points evenly spaced in 1/x instead of x. so that either the upper 
// limit is infinity or the lower is negative infinity
{
  double func(const double x){
    return Midpnt<T>::funk(1.0/x)/(x*x); // here's the change of variable
  }
  Midinf(T &funcc, const double aa, const double bb) :
    Midpnt<T>(funcc, aa, bb) {
    Midpnt<T>::a=1.0/bb; //set the limits of integration
    Midpnt<T>::b=1.0/bb;
  }
};

template <class T>
struct Midexp : Midpnt<T>
// Same thing as Midpnt, but used a mapped integration rule to integrate 
// to infinity with a falloff of type Exp(-x).
{
  double func(const double x){
    return Midpnt<T>::funk(-log(x))/x; // here's the change of variable
  }
  Midexp(T &funcc, const double aa, const double bb) :
    Midpnt<T>(funcc, aa, bb) {
    Midpnt<T>::a=0.0; //set the limits of integration
    Midpnt<T>::b=exp(-aa);
  }
};

template<class T>
struct DErule : Quadrature
//Implements the DE rule as explained on pp. 172-179 of NR
{
  double a, b, hmax, s;
  T &funk;
  DErule(T & funcc, const double aa = 0.0, const double bb = 1.0, const double hmaxx = 3.7) :
    funk(funcc), a(aa), b(bb), hmax(hmaxx) {n = 0;}
  //constructor. funcc is the function or functor that provides the
  //function to be integrated between limits aa and bb, also
  //input. The function operator in funcc takes two arguments, x and
  //delta, as described in NR. The range of integration in the
  //transformed variable t is (-hmaxx,hmaxx). Typical values of hmaxx
  //are 3.7 for logarithmic or milder singularities, and 4.3 for
  //square root singularities. See NR.

  double next()
  // simple trapezoidal integration in the transformed variable t.
  {
    double del, fact, q, sum, t, twoh;
    int it, j;
    n++;
    if (n == 1) {
      fact = 0.25;
      return s = hmax*2.0*(b-a)*fact*func(0.5*(b+a));
    } 
    else {
      for (it = 1, j = 1; j < n-1; j++) it <<= 1;
      twoh = hmax/it; //twice the spacing of the points to be added
      t = 0.5*twoh; // initialize t
      for (sum = 0.0, j = 0; j < it; j++)
	{
	  q = exp(-2.0*sinh(t));
	  del = (b-a)*q/(1.0 + q);
	  fact = q/((1.0+q)*(1.0+q))*cosh(t);
	  sum += fact*(func(a+del)+func(b-del));
	  t += twoh;
	}
      return  s = 0.5*s + (b-a)*twoh*sum;
    }
  }
  virtual double func(double x) {return funk(x);} 
  // identity mapping here generalizes to allow non-identity mappings 
  // for other kinds of mapped quadrature rules.
};


//performs the trapezoidal rule to a specified accuracy
template <class T>
double qtrap(T &func, const double a, const double b, const double EPS=1.0e-10)
// returns the integral of a function or functor from a to b. 
// EPS is the desired fractional accuracy and JMAX-1 the maximum 
// number of recursions allowed. Uses trapezoidal rule to integrate.
{
  const int JMAX = 20;
  double s, olds = 0.0;
  Trapzd<T> t(func, a, b); //define the integral object
  for (int j = 0; j < JMAX; j++)
    {
      s = t.next();
      if (j > 5) // avoid spurious early convergence
	if (abs(s-olds) < EPS*abs(olds) || (s == 0.0 && olds == 0.0)) return s;
      olds = s;
    }
  throw("Too many steps in routine qtrap");
};

//simpson's rule from trapezoidal rule:
template<class T>
double qsimp(T &func, const double a, const double b, const double EPS = 1.0e-10)
// same as qtrap, above, but uses simpson's rule as explained on NR pg. 165
{
  const int JMAX = 15;
  double s, st, ost = 0.0, os = 0.0;
  DErule<T> t(func, a, b);
  for (int j = 0; j < JMAX; j++)
    {
      st = t.next();
      s = (4.0*st - ost)/3.0;
      if (j > 5) // avoid spurious early convergence
	if (abs(s-os) < EPS*abs(os) || (s == 0.0 && os == 0.0)) return s;
      os = s;
      ost = st;      
    }
  throw("too many steps in routine qsimp");
}




//some semi-infinite ranges:

template<class T>
struct DEruleinf : Quadrature
//Implements the DE rule as explained on pp. 172-179 of NR
{
  double a, b, s; //limits of integration and current integral value
  T &funk;
 DEruleinf(T & funcc, const double aa = -4.3, const double bb = 4.3) :
  funk(funcc), a(aa), b(bb) {n = 0;}


  double next() {
    double x, tnm, sum, del;
    int it, j;
    n++;
    if (n == 1) 
    	return (s = 0.5*(b-a)*(func(a)+func(b))); // initialize the integral value
    else 
      {
	for (it = 1, j = 1; j < n - 1; j++) it <<= 1; // i = number of sample points inside interval
	tnm = it;
	del = (b-a)/tnm; // width of each region
	x = a + 0.5 * del; // set x at the midpoint of the first interval
	for (sum = 0.0, j = 0; j < it; j++, x+=del) sum += func(x);
	s = 0.5*(s+(b-a)*sum/tnm);
	return s;
      }
  }
  
  virtual double func(const double t) {
    double x = sinh(sinh(t));
    double dxdt = cosh(sinh(t))*cosh(t);
    return dxdt*funk(x);
  } 


  // identity mapping here generalizes to allow non-identity mappings 
  // for other kinds of mapped quadrature rules.
};

//simpson's rule from trapezoidal rule:
template<class T>
double qsimpinf(T &func, const double EPS = 1.0e-10)
// same as qtrap, above, but uses simpson's rule as explained on NR pg. 165
{
  const int JMAX = 15;
  double s, st, ost = 0.0, os = 0.0;
  Trapzd<T> t(func,-8.0,8.0);
  for (int j = 0; j < JMAX; j++)
    {
      st = t.next();
      s = (4.0*st - ost)/3.0;
      if (j > 5) // avoid spurious early convergence
	if (abs(s-os) < EPS*abs(os) || (s == 0.0 && os == 0.0)) return s;
      os = s;
      ost = st;      
    }
  throw("too many steps in routine qsimpinf");
}













//##############################
// here begins line integration:
//##############################

/* True line integrals require two things: a function f defined in
   N-dimensions, which is to be integrated along a line. We can
   implement this using all of the sophistication of the 1D
   integration techniques above, provided the user passes two
   functions to the intergrator: (1) the function pointer for a
   function to be integrated, and (2) a function defining the path of
   the line through the N-dimensional space. */
template<class Tfunc, class Tline>
struct linetrap : Quadrature
{
  double a, b, s; // endpoints and current value of integral
  Tfunc &funk;
  Tline &linefunc;
  linetrap(Tfunc & funcc,
	   Tline & linefuncc,
	   const double aa,
	   const double bb) : funk(funcc), linefunc(linefuncc), a(aa), b(bb) {n = 0;}
  double next() {
    double x0, x1, f0, f1, tnm, del;
    double rdist;
    VecDoub r0, r1;
    int it, j, i;
    n++;
    for (it = 1, j = 1; j < n; j++) it <<= 1; // it = number of sample points inside interval
    tnm = it;
    del = (b-a)/tnm; // width of each region
    //    std::cout<< "del = " << del << std::endl; 
    x0 = a, x1=a+del; //start at the very beginning
    s=0.0;
    r0=linefunc(x0);
    f0=funk(r0);
    for (j = 0; j < it; j++) {
      /* std::cout << "r0 = [" << r0[0] << ", " <<r0[1] << ", " << r0[2] <<"]" << std::endl; */
      //      std::cout<< "x1 = " << x1 <<std::endl;
      r1=linefunc(x1);
      /* std::cout << "r1 = [" << r1[0] << ", " <<r1[1] << ", " << r1[2] <<"]" << std::endl; */
      f1=funk(r1);
      /* std::cout << "f1 = " << f1 << std::endl; */
      rdist=0.0;
      for (i=0;i<r0.size();i++) { rdist+=(r1[i]-r0[i])*(r1[i]-r0[i]); }
      rdist=sqrt(rdist);
      /* std::cout << "rdist = " << rdist << std::endl; */
      s += rdist*0.5*(f0+f1);
      x0=x1;
      r0=r1;
      f0=f1;
      x1+=del;
      //      std::cin.get();
    }
    return s;
  }
  /* virtual double func(const double x) { */
  /*   return funk(linefunc(x)); */
  /* } */
};

//performs the trapezoidal rule to a specified accuracy
template <class Tfunc, class Tline>
double qlinetrap(Tfunc &func, Tline &linefunc, const double a, const double b, const double EPS=1.0e-10)
/* returns the integral of a function or functor on the line defined
   by linefunc from a to b. EPS is the desired fractional accuracy and
   JMAX-1 the maximum number of recursions allowed. Uses trapezoidal
   rule to integrate. */
{
  const int JMAX = 30;
  double s, olds = 0.0;
  linetrap<Tfunc,Tline> t(func, linefunc, a, b); //define the integral object
  for (int j = 0; j < JMAX; j++)
    {
      s = t.next();
      if (j > 5) // avoid spurious early convergence
	if (abs(s-olds) < EPS*abs(olds) || (s == 0.0 && olds == 0.0)) return s;
      olds = s;
      //      std::cout << "s = " << s << std::endl;
      
    }
  throw("Too many steps in routine qtrap");
};

template<class Tfunc>
struct linetrap_2pt : Quadrature
//integrates a function along an N-dimensional straight line.
{
  double a, b, s; // endpoints and current value of integral
  Tfunc &funk;
  VecDoub r0, r1;
  linetrap_2pt(Tfunc & funcc,
	       const VecDoub rr0,
	       const VecDoub rr1) : funk(funcc), r0(rr0), r1(rr1) {
    n = 0;
    a = 0;
    b = 0;
    for (int i = 0; i < r0.size(); i++)
      b+=(r1[i]-r0[i])*(r1[i]-r0[i]);
    b=sqrt(b);
  }
  VecDoub linefunc(double s) {
    double frac=s/b;
    VecDoub rs(r0.size());
    for (int i=0;i<r0.size();i++) {
      rs[i]=frac*(r1[i]-r0[i])+r0[i];
    }
    return rs;
  }
  double next() {
    double x0, x1, f0, f1, tnm, del;
    double rdist;
    VecDoub r0, r1;
    int it, j, i;
    n++;
    for (it = 1, j = 1; j < n; j++) it <<= 1; // it = number of sample points inside interval
    tnm = it;
    del = (b-a)/tnm; // width of each region
    //    std::cout<< "del = " << del << std::endl; 
    x0 = a, x1=a+del; //start at the very beginning
    s=0.0;
    r0=linefunc(x0);
    f0=funk(r0);
    for (j = 0; j < it; j++) {
      // std::cout << "r0 = [" << r0[0] << ", " <<r0[1] << ", " << r0[2] <<"]" << std::endl;
      // std::cout<< "x1 = " << x1 <<std::endl;
      r1=linefunc(x1);
      // std::cout << "r1 = [" << r1[0] << ", " <<r1[1] << ", " << r1[2] <<"]" << std::endl;
      f1=funk(r1);
      // std::cout << "f1 = " << f1 << std::endl;
      rdist=0.0;
      for (i=0;i<r0.size();i++) { rdist+=(r1[i]-r0[i])*(r1[i]-r0[i]); }
      rdist=sqrt(rdist);
      // std::cout << "rdist = " << rdist << std::endl;
      s += rdist*0.5*(f0+f1);
      if (s!=s) {
	std::cout << "r0 = [" << r0[0] << ", " <<r0[1] << ", " << r0[2] <<"]" << std::endl;
	std::cout << "x0 = " << x0 << std::endl;
	std::cout << "r0 = [" << r0[0] << ", " <<r0[1] << ", " << r0[2] <<"]" << std::endl;
	std::cout << "f0 = " << f0 << std::endl;
	std::cout << "x1 = " << x1 << std::endl;
	std::cout << "r1 = [" << r1[0] << ", " <<r1[1] << ", " << r1[2] <<"]" << std::endl;
	std::cout << "f1 = " << f1 << std::endl;
	std::cout << "rdist = " << rdist << std::endl;
      	std::cin.get();
      }
      x0=x1;
      r0=r1;
      f0=f1;
      x1+=del;
      //      std::cin.get();
    }
    return s;
  }
};

template <class Tfunc>
double qlinetrap_2pt(Tfunc &func, const VecDoub r0, VecDoub r1, const double EPS=1.0e-10)
/* returns the integral of a function or functor on the line defined
   by linefunc from a to b. EPS is the desired fractional accuracy and
   JMAX-1 the maximum number of recursions allowed. Uses trapezoidal
   rule to integrate. */
{
  const int JMAX = 30;
  double s, olds = 0.0;
  linetrap_2pt<Tfunc> t(func, r0, r1); //define the integral object
  for (int j = 0; j < JMAX; j++)
    {
      s = t.next();
      if (j > 5) // avoid spurious early convergence
	if (abs(s-olds) < EPS*abs(olds) || (s == 0.0 && olds == 0.0)) return s;
      olds = s;
      //      std::cout << "s = " << s << std::endl;
      
    }
  throw("Too many steps in routine qtrap");
};

 
/* //Trapezoidal rule */
/* template <class T> */
/* struct lineinttrap_rsym : Quadrature  */
/* { */
/*   double a, b, s; //limits of integration and current integral value */
/*   int dim; */
/*   VecDoub r1, r2; */
/*   T &funk; */
/*  lineinttrap_rsym(T &funcc, const VecDoub &r1p, const VecDoub &r2p) : funk(funcc) { */
/*     n = 0; */
/*     r1 = r1p; */
/*     r2 = r2p; */
/*     dim = r1p.size(); */
/*     a = 0.0; */
/*     b = 0.0; */
/*     for (int i = 0; i < dim; i++) */
/*       b += (r2[i]-r1[i])*(r2[i]-r1[i]); */
/*     b = sqrt(b); */
/*   } */
  
/*   double next() { */
/*     double x, tnm, sum, del; */
/*     int it, j; */
/*     n++; */
/*     if (n == 1)  */
/*     	return (s = 0.5*(b-a)*(func(a)+func(b))); // initialize the integral value */
/*     else  */
/*       { */
/* 	for (it = 1, j = 1; j < n - 1; j++) it <<= 1; // i = number of sample points inside interval */
/* 	tnm = it; */
/* 	del = (b-a)/tnm; // width of each region */
/* 	x = a + 0.5 * del; // set x at the midpoint of the first interval */
/* 	for (sum = 0.0, j = 0; j < it; j++, x+=del) sum += func(x); */
/* 	s = 0.5*(s+(b-a)*sum/tnm); */
/* 	return s; */
/*       } */
/*   } */

/*   double func(const double &x) */
/*   // returns the function value at a point which is a fraction x along */
/*   // the line between r1p and r2p: */
/*   { */
/*     double t1, t2, temp, r = 0.0; */
/*     for (int i = 0; i < dim; i++){ */
/*       t1 = r1[i]; */
/*       t2 = r2[i]; */
/*       temp = x/(b-a)*(t2-t1)+t1; */
/*       temp *= temp; */
/*       r += temp; */
/*     } */
/*     r = sqrt(r); */
    
/*     return funk(r); */
/*   } */
  
/* }; */

/* //performs the trapezoidal rule to a specified accuracy */
/* template <class T> */
/* double lineqtrap_rsym(T &func, const VecDoub r1, const VecDoub r2, const double EPS=1.0e-10) */
/* // returns the integral of a function or functor from a to b.  */
/* // EPS is the desired fractional accuracy and JMAX-1 the maximum  */
/* // number of recursions allowed. Uses trapezoidal rule to integrate. */
/* { */
/*   const int JMAX = 30; */
/*   double s, olds = 0.0; */
/*   lineinttrap_rsym<T> t(func, r1, r2); //define the integral object */
/*   for (int j = 0; j < JMAX; j++) */
/*     { */
/*       s = t.next(); */
/*       if (j > 5) // avoid spurious early convergence */
/* 	if (abs(s-olds) < EPS*abs(olds) || (s == 0.0 && olds == 0.0)) return s; */
/*       olds = s; */
/*     } */
/*   throw("Too many steps in routine qtrap"); */
/* }; */



/* template<class T> */
/* struct lineint_rsym : Quadrature */
/* // line integration of a function of radius alone between two points */
/* // with cartesian coordinatesp */
/* { */
/*   double a, b; */
/*   double hmax, s; */
/*   int dim; */
/*   VecDoub r1; */
/*   VecDoub r2; */
/*   T &funk; */
/*   lineint_rsym(T & funcc, const VecDoub &r1p, const VecDoub &r2p,  */
/* 	  const double hmaxx = 3.7) : */
/*     funk(funcc), hmax(hmaxx) { */
/*     n = 0; */
/*     r1 = r1p; */
/*     r2 = r2p; */
/*     dim = r1p.size(); */
/*     a = 0.0; */
/*     b = 0.0; */
/*     for (int i = 0; i < dim; i++) */
/*       b += (r2[i]-r1[i])*(r2[i]-r1[i]); */
/*     b = sqrt(b); */
/*   } */

/*   double next() */
/*   // simple trapezoidal integration in the transformed variable t. */
/*   { */
/*     double del, fact, q, sum, t, twoh; */
/*     int it, j; */
/*     n++; */
/*     if (n == 1) { */
/*       fact = 0.25; */
/*       return s = hmax*2.0*(b-a)*fact*func(0.5*(b+a)); */
/*     }  */
/*     else { */
/*       for (it = 1, j = 1; j < n-1; j++) it <<= 1; */
/*       twoh = hmax/it; //twice the spacing of the points to be added */
/*       t = 0.5*twoh;   // initialize t */
/*       for (sum = 0.0, j = 0; j < it; j++) */
/* 	{ */
/* 	  q = exp(-2.0*sinh(t)); */
/* 	  del = (b-a)*q/(1.0 + q); */
/* 	  fact = q/((1.0+q)*(1.0+q))*cosh(t); */
/* 	  sum += fact*(func(a+del)+func(b-del)); */
/* 	  t += twoh; */
/* 	} */
/*       return  s = 0.5*s + (b-a)*twoh*sum; */
/*     } */
/*   } */

/*   double func(const double &x) */
/*   // returns the function value at a point which is a fraction x along */
/*   // the line between r1p and r2p: */
/*   { */
/*     double t1, t2, temp, r = 0.0; */
/*     for (int i = 0; i < dim; i++){ */
/*       t1 = r1[i]; */
/*       t2 = r2[i]; */
/*       temp = x/(b-a)*(t2-t1)+t1; */
/*       temp *= temp; */
/*       r += temp; */
/*     } */
/*     r = sqrt(r); */
    
/*     return funk(r); */
/*   } */

/* }; */


/* //Transforms a line integral according to notes ``Line integrals of */
/* //radially symmetric functions'' into an integral in radius */
/* double _rmin,_rmax,_lsq,_cosxi;//swap variables */
/* double (*_funcpt)(double) = NULL; //swap function pointer */
/* double linefunc_g(const double &r); //function prototype */

/* double lineinttrans(double (*linefunc_f)(double), */
/* 		    const VecDoub &r1, const VecDoub &r2, */
/* 		    const double EPS = 1e-10) */
/* { */
/*   _funcpt = linefunc_f; // pass the function to the transformed integrand */
  
/*   double I; //current integral value */
/*   int dim = r1.size(); */

/*   double rr1 = 0.0, rr2 = 0.0, l = 0.0, lsq = 0.0, cosxi = 0.0; */
/*   // radii of initial and final points, length of line, lol squared, */
/*   // and r2*cos(xi), where xi is the angle between r1 and r2. */
  
/*   double srmin; // minimum radius along the line and path length at */
/* 		// which it occurs */

  
/*   // initialize radius-type parameters */
/*   for (int i = 0; i < dim; i++) { */
/*     rr1 += r1[i]*r1[i];// rr1 = |r1|^2 */
/*     rr2 += r2[i]*r2[i];// rr2 = |r2|^2 */
/*     cosxi += r1[i]*r2[i];// cosxi = r1.r2 */
/*   } */

/*   lsq = rr1 + rr2 - 2.0*cosxi;// lsq = |r1|^2 +| r2|^2 - 2*r1.r2 */
/*   _lsq = lsq; // pass lsq out */
/*   l = sqrt(lsq); */
/*   std::cout << "l = "<< l << std::endl; */
/*   srmin = (rr1 - cosxi)/l; // srmin = ( |r1|^2 - r1.r2 ) / l; */
/*   std::cout << "srmin = " << srmin << std::endl; */
  
/*   rr1 = sqrt(rr1);// rr1 = |r1| */
/*   std::cout << "rr1 = " << rr1 << std::endl; */
/*   rr2 = sqrt(rr2);// rr2 = |r2| */
/*   std::cout << "rr2 = " << rr2 << std::endl; */
/*   cosxi /= rr1*rr2;// cosxi = r1.r2/(|r1|*|r2|) */
/*   _cosxi = cosxi; // pass cosxi out */
/*   std::cout << "cosxi = " << cosxi << std::endl; */
  
/*   if (srmin <= 0) { // monotonically increasing */
/*     _rmin = rr1;//pass variables out */
/*     _rmax = rr2; */
/*     I = qsimp(linefunc_g,_rmin,_rmax,EPS); */
/*   } else if (srmin >= l) {// monotonically decreasing */
/*     _rmin = rr2; */
/*     _rmax = rr1; */
/*     I = qsimp(linefunc_g,_rmin,_rmax,EPS); */
/*   } else { // mixed; rerun this function with a split at the minimum */
/* 	   // radius. */
/*     VecDoub rmin(dim); // minimum radius point; calculate */
/*     for (int i = 0; i < dim; i++) rmin[i] = r1[i] + (r2[i] - r1[i])*srmin/l; */
/*     I = lineinttrans(linefunc_f, rmin, r1, EPS); */
/*     I += lineinttrans(linefunc_f,rmin,r2,EPS); */
/*   } */

/*   return I; */
/* } */

/* double linefunc_g(const double &r) */
/* // variable transformed function g for transformed line integration */
/* // (see notes ``Line integrals of radially symmetric functions'', pg. 5, */
/* // step 6). */
/* { */
/*   double rr1 = _rmin; */
/*   double rr2 = _rmax; */
/*   double lsq = _lsq; */
/*   double cxi = _cosxi; */

/*   return _funcpt(r)*1.0/sqrt( 1.0- */
/* 			  ( */
/* 			   1.0 - (rr2*cxi-rr1)*(rr2*cxi-rr1)/lsq */
/* 			   ) */
/* 			  *rr1*rr1/(r*r) */
/* 			  ); */


/* } */
    





/* // simpson line integrals */
/* template<class T> */
/* double lineqsimp(T &func, const VecDoub r1, const VecDoub r2, const double EPS = 1.0e-10) */
/* // same as qtrap, above, but uses simpson's rule as explained on NR pg. 165 */
/* { */
/*   const int JMAX = 30; */
/*   double s, st, ost = 0.0, os = 0.0; */
/*   lineint_rsym<T> t(func, r1, r2); */
/*   for (int j = 0; j < JMAX; j++) */
/*     { */
/*       st = t.next(); */
/*       s = (4.0*st - ost)/3.0; */
/*       if (j > 5) // avoid spurious early convergence */
/* 	if (abs(s-os) < EPS*abs(os) || (s == 0.0 && os == 0.0)) return s; */
/*       os = s; */
/*       ost = st;       */
/*     } */
/*   std::cout << "ERROR: too many steps in routine lineqsimp at "; */
/*   s = -1.0; */
/*   return s; */
/* } */


  
// VecDoub minpt(const VecDoub &r1, const VecDoub &r2, int dim)
// // returns the point between r1 and r2 inclusive which is closest to
// // the origin
// {
//   // figure out where along the line the minimum radius lies
//   double num = 0.0, denom = 0.0, t1, t2, crs;
//   int i;
//   for (i = 0; i < dim; i++) {
//     t1 = r1[i];
//     t2 = r2[i];
//     crs = t1*t2;
//     t1 *= t1;
//     t2 *= t2;
//     num += t1 - crs;
//     denom += t1 + t2 - 2.0*crs;
//     }
//   double minalpha = num / denom;
  
//   // return the point inside the range closest to the origin:
//   if (minalpha > 1.0)
//     return r2;
//   else if (minalpha < 0.0) 
//     return r1;
//   else {
//     VecDoub res dim;
//     for (i = 0; i < dim; i++){
//       t1 = r1[i];
//       t2 = r2[i];
//       res[i] = minalpha*(t2-t1)+t1;
//     }
//     return res;
//   }
// }

#endif
