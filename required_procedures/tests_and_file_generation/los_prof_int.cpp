//a program to compute line-of sight optical depths as a function of
//doppler shift through a chamberlain exosphere.

//This program generates H_LOS_prof.dat, required for correct
//computation of the IPH extinction.

#include "nr3.h"
#include "quadrature.h"
#include <cmath>
#include <iostream>

double pi=3.14159265358979;

//first, a function to define the limits of the integration in angle
void angles(double lc, double lb_cos_t, double v, double t, double cos_t, double a, VecDoub &anglevec) {
  //  double cos_t = cos(t);
  //  double lb_cos_t = lb*cos_t;
  double cx_sq = lb_cos_t*lb_cos_t/(lc+lb_cos_t);
  double vsq = v*v;
  if (a*a <= cx_sq && v*v <= cx_sq) {
    //circular section is entirely inside hyperbola, return the angles
    //that work when the circle is tangent to the hyperbola
    // double tan_t=tan(t);
    // double ret=-tan_t/sqrt(cx_sq/(a*a)-1);
    // anglevec[0]=acos(ret);
    // ret=tan_t/sqrt(cx_sq/(a*a)-1);
    // anglevec[1]=acos(ret);

    // or just return pi/2 for both angles
    double ret=pi/2;
    anglevec[0]=ret;
    anglevec[1]=ret;
    return;
  }
  else {
    //there is an intersection between sphere, plane, and hyperbola. Find it:
    double top1 = (vsq*(1+lb_cos_t/lc)-lb_cos_t*lb_cos_t/lc)*(1-lb_cos_t/lc);
    if (top1 < 0) {
      std::cout << "top1 less than zero!" << std::endl;
      std::cout << "lc = " << lc << std::endl;
      std::cout << "lb = " << lb_cos_t/cos_t << std::endl;
      std::cout << "v = " << v << std::endl;
      std::cout << "a = " << a << std::endl;
      std::cout << "t = " << t << std::endl;
      std::cout << "anglevec = [" << anglevec[0] << ", " << anglevec[1] << "]\n";
      std::cout << "top1 = " << top1 << std::endl;
      throw("top1 less than zero in angles()");
    }
    top1=sqrt(top1);
    double top2 = a*sin(t);
    double bottom = cos_t*sqrt(vsq-a*a);
    if (vsq < a*a) {
      std::cout << "vsq is less than a*a!" << std::endl;
      std::cout << "lc = " << lc << std::endl;
      std::cout << "lb = " << lb_cos_t/cos_t << std::endl;
      std::cout << "v = " << v << std::endl;
      std::cout << "a = " << a << std::endl;
      std::cout << "t = " << t << std::endl;
      std::cout << "anglevec = [" << anglevec[0] << ", " << anglevec[1] << "]\n";
      std::cout << "top1 = " << top1 << std::endl;
      std::cout << "top2 = " << top2 << std::endl;
      std::cout << "bottom = " << bottom << std::endl;
      throw("vsq less than a*a in angles()");
    }
    double ret=(top1-top2)/bottom;
    if (ret>1.0) ret=1.0;
    if (ret<-1.0) ret=-1.0;
    ret=acos(ret);
    anglevec[0]=ret;
    ret=(top1+top2)/bottom;
    if (ret>1.0) ret=1.0;
    if (ret<-1.0) ret=-1.0;
    ret=acos(ret);
    anglevec[1]=ret;
    return;
  }
}

double v_integrand(double lc, double lb_cos_t, double v, double t, double cos_t, double a) {
  VecDoub anglevec(2);
  angles(lc, lb_cos_t, v, t, cos_t, a, anglevec);
  if (v < abs(a)) {
    return 0.0;
  } else {
    double sqrt_lb_cos_t=sqrt(lb_cos_t);
    if (abs(a) >= sqrt_lb_cos_t) {
      return exp(-v*v)*v*2.0*(anglevec[0]);
    } else {
      return exp(-v*v)*v*2.0*(anglevec[0]+anglevec[1]);
    }
  }
}

using std::cout;
using std::cin;
using std::endl;

struct vint {
  double lc, lb, lb_cos_t, t, cos_t, a;
  vint(double lcc,double lbb,double tt, double aa) : lc(lcc), lb(lbb), t(tt), a(aa) {
    cos_t=cos(t);
    lb_cos_t=lb*cos_t;
  }
  double operator()(double v) {
    return v_integrand(lc,lb_cos_t,v,t,cos_t,a);
  }
};

double t_integrand(double lc, double lb, double t, double a) {
  vint innerint(lc,lb,t,a);
  double cos_t=cos(t);
  //  std::cout << "t = " << t << std::endl;  
  Midexp<vint> q(innerint,abs(a),1e99);
  return exp(lb*cos_t)/(cos_t*cos_t)*qromo(q,1e-6);
  //  return exp(lb*cos_t)/(cos_t*cos_t)*qsimpsemiinf(innerint,1e-4);
}

struct tint {
  double lc, lb, a;
  tint(double lcc,double lbb, double aa) : lc(lcc), lb(lbb), a(aa) {}
  double operator()(double t) {
    return t_integrand(lc,lb,t,a);
  }
};

double abs_a(double lc, double lb, double a, double t0) {
  tint outerint(lc,lb,a);
  return sqrt(lc*lc*lc)/lb*exp(-lc)*qsimp(outerint,t0,pi/2,1e-3);
}

//we need to multiply by the contants in front of the integral to
//return the optical depth. These factors are:
double s0=0.0110224;//cm2 Hz, lya natural line total cross section
double l0=121.6e-7;//cm, lya wavelength
double rc=(3396+200)*1e5;//cm, radius of exobase
double G=6.67e-8;//dyn cm2/gm2, Big G
double Mm=0.1076*5.98e27;//gm, mass of Mars
  
double taucoef=(s0*l0)*sqrt(rc*rc*rc/(2*pi*pi*pi*G*Mm));


//this function returns a profile of atmospheric optical depth with
//all factors accounted for except for H density at the
//exobase.
void tau_profile(double lc, double lb, double t0, VecDoub &avec, VecDoub &tauvec, double eps=1e-8,double astep=0.05,int maxsize=200) {
  double apos0=0.0;
  double taupos0=taucoef*abs_a(lc,lb,apos0,t0);
  // std::cout << "a = " << apos0 << std::endl;
  // std::cout << "tau = " << taupos0 << std::endl;
  if (taupos0 < eps) {
    // we are already below the absorption threshold
    avec.resize(1);
    avec[0]=apos0;
    tauvec.resize(1);
    tauvec[0]=taupos0;
    return;
  } else {
    double* apos = new double[maxsize];
    double* taupos = new double[maxsize];
    double* aneg = new double[maxsize];
    double* tauneg = new double[maxsize];
    int n_a=0;
    double tausize=taupos0;
    while (tausize > eps) {
      if (n_a == maxsize) {throw("too many steps in routine tau_profile; increase astep or increase maxsize");}
      apos[n_a]=(n_a+1)*astep;
      // std::cout << "a = " << apos[n_a] << std::endl;
      taupos[n_a]=taucoef*abs_a(lc,lb,apos[n_a],t0);
      // std::cout << "tau = " << taupos[n_a] << std::endl;
      aneg[n_a]=-(n_a+1)*astep;
      // std::cout << "a = " << aneg[n_a] << std::endl;
      tauneg[n_a]=taucoef*abs_a(lc,lb,aneg[n_a],t0);
      // std::cout << "tau = " << tauneg[n_a] << std::endl;
      tausize = taupos[n_a] > tauneg[n_a] ? taupos[n_a] : tauneg[n_a];
      n_a++;
    }
    avec.resize(2*n_a+1);
    tauvec.resize(2*n_a+1);
    avec[n_a]=apos0;
    tauvec[n_a]=taupos0;
    // std::cout << "n_a = " << n_a << std::endl;
    // std::cout << "avec[n_a] = " << avec[n_a] << std::endl;
    // std::cout << "tauvec[n_a] = " << tauvec[n_a] << std::endl;
    for (int i=0;i<n_a;i++) {
      // std::cout << "i = " << i << std::endl;
      // std::cout << "n_a-i-1 = " << n_a-i-1 << std::endl;
      // std::cout << "aneg[n_a-i] = " << aneg[n_a-i-1] << std::endl;
      // std::cout << "tauneg[n_a-i] = " << tauneg[n_a-i-1] << std::endl;
      avec[i]=aneg[n_a-i-1];
      tauvec[i]=tauneg[n_a-i-1];
      // std::cout << "avec[i] = " << avec[i] << std::endl;
      // std::cout << "tauvec[i] = " << tauvec[i] << std::endl;
      
      // std::cout << "n_a+i+1 = " << n_a+i+1 << std::endl;
      // std::cout << "apos[i+1] = " << apos[i] << std::endl;
      // std::cout << "taupos[i+1] = " << taupos[i] << std::endl;
      avec[n_a+i+1]=apos[i];
      tauvec[n_a+i+1]=taupos[i];
      // std::cout << "avec[n_a+i+1] = " << avec[n_a+i+1] << std::endl;
      // std::cout << "tauvec[n_a+i+1] = " << tauvec[n_a+i+1] << std::endl;
      // std::cin.get();
    }
    delete [] apos;
    delete [] aneg;
    delete [] taupos;
    delete [] tauneg;
    return;
  }
}

struct H_LOS_sim {
  int vecsize;
  VecDoub avec, tvec;
  
  H_LOS_sim() {};
  
  void simulate(double lc, double lb, double t0, double eps=1e-8,double astep=0.001,int maxsize=10000) {
    tau_profile(lc, lb, t0, avec, tvec, eps, astep, maxsize);
    vecsize=avec.size();
  }
};


int main() {
  //simulate LOS line profiles for a variety of lc, lb, and t0
  double lc, lb, t0;

  double lcmin=0.5, lcmax=14.5, lcstep=0.5;
  int nlc=ceil((lcmax-lcmin)/lcstep)+1, ilc;

  double lbfracmin=0.05, lbfracmax=1.0, lbfracstep=0.05;
  int nlbfrac=ceil((lbfracmax-lbfracmin)/lbfracstep)+1,ilbfrac;
  
  int nt0=20;
  int it0;
  VecDoub avec(0);
  VecDoub tauvec(0);

  //open file for output
  ofstream outfile;
  outfile.open("H_LOS_prof_coords.dat");
  outfile.width(10);
  outfile << nlc;
  outfile.width(10);
  outfile << nlbfrac;
  outfile.width(10);
  outfile << nt0;
  outfile << endl;
  outfile << endl;

  for (ilc=0;ilc<nlc;ilc++) {
    outfile.width(10);
    outfile << ilc*lcstep+lcmin;
  }
  outfile << endl;

  for (ilbfrac=0;ilbfrac<nlbfrac;ilbfrac++) {
    outfile.width(10);
    outfile << ilbfrac*lbfracstep+lbfracmin;
  }
  outfile << endl;
  
  for (it0=0;it0<nt0;it0++) {
    outfile.width(10);
    outfile << -pi/2+it0*pi/nt0;
  }
  outfile << endl;
    
  for (ilc=0;ilc<nlc;ilc++) {
    lc=ilc*lcstep+lcmin;
    for (ilbfrac=0;ilbfrac<nlbfrac;ilbfrac++) {
      lb=(ilbfrac*lbfracstep+lbfracmin)*lc;
      for (it0=0;it0<nt0;it0++) {
  	t0=-pi/2+it0*pi/nt0;

  	std::cout << "lc = " << lc << ", lb = " << lb << ", t0 = " << t0 << endl;

  	//simulate
  	tau_profile(lc, lb, t0, avec, tauvec);

  	//writeoutput to file
  	outfile.width(10);
  	outfile << lc;
  	outfile.width(10);
  	outfile << lb;
  	outfile.width(10);
  	outfile << t0;
  	outfile.width(10);
  	outfile << avec.size();
  	outfile << endl;
  	for (int itau=0;itau<avec.size();itau++) {
  	  outfile.width(20);
  	  outfile << avec[itau];
  	}
  	outfile << endl;
  	for (int itau=0;itau<avec.size();itau++) {
  	  outfile.width(20);
  	  outfile << tauvec[itau];
  	}
  	outfile << endl;
      }
    }
    //    std::cin.get();
  }
  //close file
  outfile.close();

}



// int main(int argc, char* argv[]) {
//   //get the parameters from the function call:
//   double lc, lb, t0;
//   lc=atof(argv[1]);
//   lb=atof(argv[2]);
//   t0=atof(argv[3]);
//   string outfilename=argv[4];
//   // cout << "lambda_c = " << lc << endl;
//   // cout << "lambda_b = " << lb << endl;
//   // cout << "t0 = " << t0 << endl;

//   VecDoub avec(0);
//   VecDoub tauvec(0);

//   tau_profile(lc,lb,t0,avec,tauvec);

//   ofstream outfile;
//   outfile.open(outfilename.c_str());
//   for (int i=0;i<avec.size();i++) {
//     std::cout << "i = " << i << std::endl;
//     outfile << "i = " << i << std::endl;
//     std::cout << "a = " << avec[i] << std::endl;
//     outfile << "a = " << avec[i] << std::endl;
//     std::cout << "tau = " << tauvec[i] << std::endl;
//     outfile << "tau = " << tauvec[i] << std::endl;
//   }
//   outfile.close();
  
//   return 0;
// }

// int main(int argc, char* argv[]) {
//   //get the parameters from the function call:
//   double lc, lb, v, t, a;
//   lc=atof(argv[1]);
//   lb=atof(argv[2]);
//   v=atof(argv[3]);
//   t=atof(argv[4]);
//   a=atof(argv[5]);

//   cout << "lambda_c = " << lc << endl;
//   cout << "lambda_b = " << lb << endl;
//   cout << "v = " << v << endl;
//   cout << "t = " << t << endl;
//   cout << "a = " << a << endl;

//   VecDoub anglevec(2);
//   double cos_t=cos(t);
//   angles(lc, lb*cos_t, v, t, cos_t, a, anglevec);
//   double losintval=v_integrand(lc, lb*cos_t, v, t, cos_t, a);

//   cout << "angles = [" << anglevec[0] << ", " << anglevec[1] << "]\n";
//   cout << "los_integrand = " << losintval << endl;

//   double totalint = 0.0;
//   double v0=a, vf=10;
//   int nsteps = 100000;
//   double dv=(vf-v0)/(nsteps-1);
//   double val0 = v_integrand(lc,lb*cos_t,v0,t,cos_t,a);
//   double val1;
//   for (int i=1; i<nsteps; i++) {
//     //    std::cout << "v = " << v0+i*dv << std::endl;
//     val1=v_integrand(lc,lb*cos_t,v0+i*dv,t,cos_t,a);
//     totalint+=(val0+val1)/2.*dv;
//     //    std::cout << "int = " << totalint << std::endl;
//     //    if (totalint != totalint) {std::cin.get();}
//     val0=val1;
//   }
//   cout << "total int = " << totalint << endl;

//   cout << "inner int = " << t_integrand(lc,lb,t,a) << endl;

//   // double totalint = 0.0;
//   // double t0=-pi/2, tf=pi/2;
//   // int nsteps = 1000;
//   // double dt=(tf-t0)/(nsteps-1);
//   // double val0 = t_integrand(lc,lb,t0,a);
//   // double val1;
//   // for (int i=1; i<nsteps; i++) {
//   //   val1=t_integrand(lc,lb,t0+i*dt,a);
//   //   totalint+=(val0+val1)/2.*dt;
//   //   val0=val1;
//   // }
//   // cout << "total int = " << totalint << endl;
  
//   tint outerint(lc,lb,a);

//   // Midpnt<tint> q(outerint,-pi/2,pi/2);
//   // double autoint = qromo(q,1e-4);
//   // cout << "auto int = " << autoint << endl;

  
//   double deint=qsimp(outerint,-pi/2,pi/2,1e-2);
//   cout << "DE int = " << deint << endl;
  
//   return 0;
// }


// int main() {
//   cout << "Enter values for the parameters: \n";
//   cout << "lambda_c = ";
//   double lc;
//   cin >> lc;
//   cout << "lambda_b = ";
//   double lb;
//   cin >> lb;
//    cout << "t = ";
//   double t;
//   cin >> t;
//   cout << "a = ";
//   double a;
//   cin >> a;  

//   VecDoub anglevec(2);

//   for (double v=a;v<10;v+=0.5) {
//     cout << "v = " << v << endl; 
//     angles(lc, lb, v, t, a, anglevec);
//     double losintval=v_integrand(lc, lb, v, t, a);
//     cout << "angles = [" << anglevec[0] << ", " << anglevec[1] << "]\n";
//     cout << "los_integrand = " << losintval << endl << endl;
//   }

//   return 0;
// }
