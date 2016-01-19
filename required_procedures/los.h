//los.h -- routines to compute intensity from source functions

#ifndef __LOS_H
#define __LOS_H

#include <iostream> // for file output and dialog
#include <fstream>  // ""
#include <cmath>    // for cos and sin
#include <cstdlib>
#include "physical.h" // physical definitions, nH, nCO2, KK, etc.
#include "quadrature.h" // to compute line integrals
#include "definitions.h" // basic parameter definitions
//#include "losdefs.h" must be added seperately in calling program
#include "nr3.h" // type VecDoub
#include "rad_trans.h" // two-point column density, holstein t and g functions
#include "ludcmp.h" // LU decomposition solver
#include "interpgen.h" // runtime hydrogen interpolation
#include "interp_1d.h" // for linear interpolation basis of Sobj class

class doub3d {

public:
  //class for 3d matrix
  double ***mat;
  int d1, d2, d3;
    
  doub3d() {mat=NULL;d1=0;d2=0;d3=0;}
  doub3d(int d11, int d22, int d33) : d1(d11), d2(d22), d3(d33) {
    mat = new double**[d1];
    for (int i1=0; i1<d1; i1++) {
      mat[i1] = new double*[d2];
      for (int i2=0; i2<d2; i2++)
	mat[i1][i2]=new double[d3];
    }
  }
  doub3d(const doub3d &M) : d1(M.d1), d2(M.d2), d3(M.d3) {
    mat = new double**[d1];
    for (int i1=0; i1<d1; i1++) {
      mat[i1] = new double*[d2];
      for (int i2=0; i2<d2; i2++) {
	mat[i1][i2]=new double[d3];
	for (int i3=0; i3<d3; i3++)
	  mat[i1][i2][i3] = M.mat[i1][i2][i3];
      }
    }
  }

  double** operator[](int i) { return mat[i]; }

  doub3d& operator=(const doub3d &M) {
    d1=M.d1;
    d2=M.d2;
    d3=M.d3;
    mat = new double**[d1];
    for (int i1=0; i1<d1; i1++) {
      mat[i1] = new double*[d2];
      for (int i2=0; i2<d2; i2++) {
	mat[i1][i2]=new double[d3];
	for (int i3=0; i3<d3; i3++)
	  mat[i1][i2][i3] = M.mat[i1][i2][i3];
      }
    }
    return *this;
  }
  
  ~doub3d() {
    for (int i1=0; i1<d1; i1++) {
      for (int i2=0; i2<d2; i2++)
	delete [] mat[i1][i2];
      delete [] mat[i1];
    }
    delete [] mat;
  }

private:
  
};











//define the Sobj class for source functions:
//needs to: load file, return interpolated values
struct Sobj {
// object for read-in and bilinear interpolation on a computed source
// function.

  string fname;
  double nH, T; //exobase H density and temperature for this S
  VecDoub rpts,tpts,ppts;
  int nrpts, ntpts, nppts;
  doub3d Smat;
  Linear_interp rterp, tterp, pterp;//objects for BOUNDARY points
  int ir, it, ip;
  bool init;

  //default constructor
  Sobj() {init=0;}
  
  //constructor to load sobj from file
  Sobj(string fnamee) : fname(fnamee) {
    ifstream Sfile;
    Sfile.open(fname.c_str());//fname must be absolute or local.
    if (Sfile.is_open()) {
      //get nH and T for this observation
      Sfile >> nH;
      // std::cout << "nH = " << nH << std::endl;
      Sfile >> T;
      // std::cout << "T = " << T << std::endl;
      
      //load the points:
      Sfile >> nrpts;
      // std::cout << "nrpts = " << nrpts << std::endl;
      rpts.resize(nrpts);
      // std::cout << " rpts = [ "; 
      for (ir = 0; ir < nrpts; ir++) {
	Sfile >> rpts[ir];
	// std::cout << rpts[ir] << ", ";
      }
      // std::cout << "\b\b ]\n";

      Sfile >> ntpts;
      // std::cout << "ntpts = " << ntpts << std::endl;
      tpts.resize(ntpts);
      // std::cout << " tpts = [ "; 
      for (it = 0; it < ntpts; it++) {
	Sfile >> tpts[it];
	// std::cout << tpts[it] << ", "; // 
      }
      // std::cout << "\b\b ]\n";

      Sfile >> nppts;
      // std::cout << "nppts = " << nppts << std::endl;
      ppts.resize(nppts);
      // std::cout << " ppts = [ "; 
      for (ip = 0; ip < nppts; ip++) {
	Sfile >> ppts[ip];
	// std::cout << ppts[ip] << ", ";
      }
      // std::cout << "\b\b ]\n";

      //pass to the interpolation objects
      rterp = Linear_interp(rpts,rpts);
      tterp = Linear_interp(tpts,tpts);
      pterp = Linear_interp(ppts,ppts);

      //allocate the S array:
      Smat=doub3d(nrpts-1,ntpts-1,nppts-1);
      
      // get the tabulated S values from the specified file
      for (ir = 0; ir < nrpts-1; ir++) {
	for (it = 0; it < ntpts-1; it++) {
	  for (ip = 0; ip < nppts-1; ip++) {
	    Sfile >> Smat[ir][it][ip];
	    // std::cout << "Smat[" << ir << "][" << it << "][" << ip << "] = "
	    // 	    <<  Smat[     ir     ][     it     ][     ip     ] << std::endl;
	    // std::cin.get();
 	  }
	}
      }
      init=1;
    } else {
      std::cout << "Sobj: Cannot find file " << fname << ".\n";
      throw("Bad filename in Sobj")
    }
  }

  //constructor to create Sobj from existing data.
  Sobj(double nHH, double TT,
       Linear_interp rterpp, Linear_interp tterpp, Linear_interp pterpp, 
       VecDoub sol) : nH(nHH), T(TT),
		      nrpts(rterpp.n), ntpts(tterpp.n), nppts(pterpp.n), 
		      rterp(rterpp),   tterp(tterpp),   pterp(pterpp) {
    //allocate the S array
    Smat=doub3d(nrpts-1,ntpts-1,nppts-1);
    //read the passed values into the S array
    for (ir = 0; ir < nrpts-1; ir++) {
      for (it = 0; it < ntpts-1; it++) {
	for (ip = 0; ip < nppts-1; ip++) {
	  int col=ir*(ntpts-1)*(nppts-1)+it*(nppts-1)+ip;
	  Smat[ir][it][ip] = sol[col];
	  // std::cout << "Smat[" << ir << "][" << it << "][" << ip << "] = "
	  // 	    <<  Smat[     ir     ][     it     ][     ip     ] << std::endl;
	  //	std::cin.get();
	}
      }
    }
    init=1;
  }

  double operator()(double rpt, double tpt, double ppt) {
    if (!init)
      throw("Sobj not initialized!");
    //    make sure we're not off-grid
    /* std::cout << "rterp.xx[0] = " << rterp.xx[0] << std::endl; */
    /* std::cout << "rterp.xx[nrpts-1] = " << rterp.xx[nrpts-1] << std::endl; */
    /* std::cout << "tterp.xx[0] = " << tterp.xx[0] << std::endl; */
    /* std::cout << "tterp.xx[ntpts-1] = " << tterp.xx[ntpts-1] << std::endl; */
    /* std::cout << "pterp.xx[0] = " << pterp.xx[0] << std::endl; */
    /* std::cout << "pterp.xx[nppts-1] = " << pterp.xx[nppts-1] << std::endl; */
    bool inrgrid = ((rterp.xx[0] >= rpt && rpt >= rterp.xx[nrpts-1])||
		    (rterp.xx[nrpts-1] >= rpt && rpt >= rterp.xx[0]));
    bool intgrid = ((tterp.xx[0] >= tpt && tpt >= tterp.xx[ntpts-1])||
		    (tterp.xx[ntpts-1] >= tpt && tpt >= tterp.xx[0]));
    bool inpgrid = ((pterp.xx[0] >= ppt && ppt >= pterp.xx[nppts-1])||
		    (pterp.xx[nppts-1] >= ppt && ppt >= pterp.xx[0]));
    if (!intgrid || !inrgrid || !inpgrid) {
      std::cout << "Point ( " << rpt << ", " << tpt << ", " << ppt << ") "
		<< "is not inside the S grid!\n";
      std::cout << "Alt = " << rpt*1e-5 - 3395 
		<< ", SZA = " << tpt*180/pi 
		<< ", phi = " << ppt*180/pi << ".\n";
      throw("bad pt in S")
      return 0.0;
    }

    //find the grid square:
    // std::cout << "rpt = " << rpt << " , tpt = " << tpt << ", ppt = "<< ppt <<".\n";
    ir = rterp.index(rpt);
    it = tterp.index(tpt);
    ip = pterp.index(ppt);
    // std::cout << "ir = " << ir << " , it = " << it << ", ip = "<< ip <<".\n";

    /* We want to interpolate the finite element S grid. The above
       indices give us the lower boundary of the regions, but to
       interpolate we need to know what the next closest point is. */
    double tr, tt, tp, tra, tta, tpa;
    int ira, ita, ipa;
    tr=(rterp.xx[ir]+rterp.xx[ir+1])/2;
    if (rterp.xx[nrpts-1] > rterp.xx[0]) {
      ira = (rpt > tr) ? ir+1 : ir-1;
    } else {
      ira = (rpt > tr) ? ir-1 : ir+1;
    }
    if (ira==-1) ira=0;
    if (ira==nrpts-1) ira=nrpts-2;
    tra=(rterp.xx[ira]+rterp.xx[ira+1])/2;

    // std::cout << "ira = " << ira << std::endl;
    // std::cout << "tr = " << tr << ", tra = " << tra << std::endl;
    
    tt=(tterp.xx[it]+tterp.xx[it+1])/2;
    ita = (tpt > tt) ? it+1 : it-1;
    if (ita==-1) {
      ita=0;
      tta=0;//endpoint is azimuthal average at t=0
    } else if (ita==ntpts-1) {
      ita=ntpts-2;
      tta=pi;//endpoint is azimuthal average at t=pi
    } else {
      tta=(tterp.xx[ita]+tterp.xx[ita+1])/2;
    }

    // std::cout << "ita = " << ita << std::endl;
    // std::cout << "tt = " << tt << ", tta = " << tta << std::endl;
    
    tp=(pterp.xx[ip]+pterp.xx[ip+1])/2;
    ipa = (ppt > tp) ? ip+1 : ip-1;
    if (ipa==-1) {//loop to phi<0
      ipa=nppts-2;
      tpa=(pterp.xx[ipa]+pterp.xx[ipa+1])/2-2*pi;//-2pi for coef
                                                 //calculation (b/c
                                                 //pterp.xx
                                                 //e[0,2pi]).
      // For example in the azimuthally symmetric case with nppts=2,
      // tp=pi, tpa=-pi, so that the phi area contribution is 2pi,
      // covering all azimuths.
    } else if (ipa==nppts-1) {//loop to phi>2pi
      ipa=0;
      tpa=(pterp.xx[ipa]+pterp.xx[ipa+1])/2+2*pi;//+2pi for coef calculation
                                                 //(b/c pterp.xx e[0,2pi])
      // For example in the azimuthally symmetric case with nppts=2,
      // tp=pi, tpa=3pi, so that the phi area contribution is 2pi,
      // covering all azimuths.
    } else {
      tpa=(pterp.xx[ipa]+pterp.xx[ipa+1])/2;
    }

    // std::cout << "ipa = " << ipa << std::endl;
    // std::cout << "tp = " << tp << ", tpa = " << tpa << std::endl;
    
    //first the coefficients (representing the degree with which to
    //weight the adjacent points)
    double ar, at, ap;
    if (ira==nrpts-2||ira==0) {
      ar=0.0;//if we're at the top or bottom of the grid, stay put!
    } else {
      ar=(rpt-tr)/(tra-tr);
    }
    at=(tpt-tt)/(tta-tt);

    ap=(ppt-tp)/(tpa-tp);

    // std::cout << "ar = " << ar << ", at = " << at << ", ap = " << ap << std::endl; 
    
    //now compute the values:
    //S121 = (this, a, this) in (r, t, p)
    double S111, S112, S121, S122, S211, S212, S221, S222;
    //t is tricky because going past the endpoint means we are
    //approaching the pole, where the value should be the azimuthal
    //average of the adjacent elements.
    S111 = Smat[ir ][it ][ip ];
    S112 = Smat[ir ][it ][ipa];
    S211 = Smat[ira][it ][ip ];
    S212 = Smat[ira][it ][ipa];
    if (ita==0||ita==ntpts-2) {
      S121 = S122 = S221 = S222 = 0.0;
      for (int iip=0;iip<nppts-1;iip++) {
	S121 += Smat[ir ][ita][iip];
	S122 += Smat[ir ][ita][iip];
	S221 += Smat[ira][ita][iip];
	S222 += Smat[ira][ita][iip];
      }
      S121 /= (nppts-1);
      S122 /= (nppts-1);
      S221 /= (nppts-1);
      S222 /= (nppts-1);
    } else { 
      S121 = Smat[ir ][ita][ip ];
      S122 = Smat[ir ][ita][ipa];
      S221 = Smat[ira][ita][ip ];
      S222 = Smat[ira][ita][ipa];
    }

    // std::cout << " [ S111 S112 ] = [" << S111 << ", " << S112 << "];\n";
    // std::cout << " [ S121 S122 ] = [" << S121 << ", " << S122 << "];\n";
    // std::cout << std::endl;
    // std::cout << " [ S211 S212 ] = [" << S211 << ", " << S212 << "];\n";
    // std::cout << " [ S221 S222 ] = [" << S221 << ", " << S222 << "];\n";
    // std::cout << std::endl;
        
    //now we can interpolate
    double Sterp;
    Sterp = (1-ar)*(1-at)*(1-ap)*S111
          + (1-ar)*(1-at)*   ap *S112
          + (1-ar)*   at *(1-ap)*S121
          + (1-ar)*   at *   ap *S122
          +    ar *(1-at)*(1-ap)*S211
          +    ar *(1-at)*   ap *S212
          +    ar *   at *(1-ap)*S221
          +    ar *   at *   ap *S222;


    // std::cout << "Sterp = " << Sterp;
    /* std::cin.get(); */

    return Sterp;
  }

};

struct Holinterp
{
  /* Structure to load tabular Holstein functions */
  string Holfilename;
  int npts;
  VecDoub tau_vec, Hol_vec;
  Linear_interp Hol_interp;

Holinterp(string Holfilenamee) 
  : Holfilename(Holfilenamee) {
    std::cout << "Reading interpolated Hol values from " << Holfilename << std::endl;
    //load the data from the file
    ifstream Holfile;
    Holfile.open(Holfilename.c_str());
    //read number of points
    Holfile >> npts;
    //    std::cout << "npts = " << npts << std::endl;
    tau_vec.resize(npts);
    Hol_vec.resize(npts);
    // clear the header
    char dumline[100];
    Holfile.getline(dumline,100);
    Holfile.getline(dumline,100);
    //    std::cout << dumline << std::endl;
    for (int i = 0; i < npts; i++) {
      Holfile >> tau_vec[i];
      //      std::cout << "tau_vec[" << i << "] = " << tau_vec[i] << std::endl;
      Holfile >> Hol_vec[i];
      //      std::cout << "Hol_vec[" << i << "] = " << Hol_vec[i] << std::endl;
      //      std::cin.get();
    }
    Holfile.close();
    Hol_interp = Linear_interp(tau_vec,Hol_vec);
  }

  double interp(const double tau) {
    return Hol_interp.interp(tau);
  }
};


struct LOS_integrator
{
  /* structure to compute intensity along LOS using saved, tabulated
     values for HolT */
  string HolTfname;
  Holinterp HolTinterp;
  
  LOS_integrator(string HolTfnamee=HolTfilename) 
    : HolTfname(HolTfnamee), HolTinterp(HolTfname) {
    std::cout << "Reading interpolated HolT values from " << HolTfname << std::endl;
  }

  double integrate(Sobj &S, atmointerp &thisatmointerp,
		   const VecDoub &pos, const VecDoub &dir,
		   const double &lineintcoef, double &tauHout, double &tauCO2out)
  {
    //get coordinates for the start point
    double x0 = pos[0];
    double y0 = pos[1];
    double z0 = pos[2];
    double r0 = x0*x0 + y0*y0 + z0*z0;
    r0 = sqrt(r0);

    /* std::cout << "r0  = { " << x0 << " , " << y0 << " , " << z0 << " }, \n"; */
    
    //get the unit vector components along the LOS
    double line_x = dir[0];
    double line_y = dir[1];
    double line_z = dir[2];

    /* std::cout << "line = { " << line_x << " , " << line_y << " , " << line_z << "}\n"; */
    /* std::cin.get(); */
	    
    // initialize the variables
    double ds = 0.0;// ds = 0.0; this changes soon
    double s = 0.0;// s = 0 at start of integration
    double rrpt = r0;// first radial point
    double tpt = std::acos(x0/r0);// first theta point
    double ppt = std::atan2(y0, z0);// first phi point
    ppt = ppt < 0 ? ppt+2*pi : ppt;// put phi in the right domain
    double xpt = x0;
    double ypt = y0;
    double zpt = z0;
    double tauH = 0.0;
    double tauCO2 = 0.0;
    double intensity = 0.0;// integral starts out at zero
    int iter = 0;
    /* std::cout << "rrpt = " << rrpt << ", rminatm = " << rminatm << ", rmax = " << rmax << std::endl; */

    //if we're outside the atmosphere to begin, advance to the edge with
    //no penalty:
    while (rrpt > rmax) {
      std::cout << "Moving to start of atmosphere:\n";
      std::cout << "Before moving: rrpt = " << rrpt << " > rmax = " << rmax <<std::endl; 
      s = xpt*line_x + ypt*line_y + zpt*line_z;
      s = -s - sqrt(s*s + rmax*rmax - rrpt*rrpt);// I know it looks negative
      // but actually s=r0.line is negative (basic quad. eq.), because
      // we are outside the atmosphere looking in.
      std::cout << "Moving forward by s = " << s << " cm." << std::endl;
      xpt += s*line_x;
      ypt += s*line_y;
      zpt += s*line_z;
      // compute the new radius and angle
      rrpt = sqrt(xpt*xpt + ypt*ypt + zpt*zpt);
      tpt = std::acos(xpt/rrpt);
      ppt = std::atan2(ypt, zpt);
      ppt = ppt < 0 ? ppt+2*pi : ppt;// put phi in the right domain
      std::cout << "Now rrpt = " << rrpt << ", tpt = " << tpt*180/pi << " degrees.\n";
      std::cin.get();
    }

    // now proceed with integration  
    while (iter == 0 || (rrpt >= rminatm && rrpt <= rmax && iter < maxit)) {
      //while still inside the atmosphere
      iter++;
      // if too many iterations, report error, but don't interrupt
      if (iter == maxit) {
	std::cout << "With r0  = { " << x0 << " , " << y0 << " , " << z0 << " }, \n";
	std::cout << "     dir = { " << line_x << " , " << line_y << " , " << line_z << " }, \n";
	std::cout << "iter > maxit!" << std::endl;
      }

      /* std::cout << "tauH = " << tauH << std::endl; */
      /* std::cout << "HolTinterp(tauH) = " << HolTinterp.interp(tauH) << std::endl; */
      /* std::cout << "tauCO2 = " << tauCO2 << std::endl; */
      /* std::cout << "exp(-tauCO2) = " << exp(-tauCO2) << std::endl; */
      /* std::cout << "S(rrpt="<< rrpt */
      /* 	        <<  ", tpt=" << tpt  */
      /*           <<  ", ppt=" << ppt */
      /* 		<<") = " << S(rrpt,tpt,ppt) << std::endl; */
      /* std::cout << "dI = " << HolTinterp.interp(tauH)*exp(-tauCO2)*S(rrpt,tpt,ppt)*dtau  */
      /* 		<< std::endl; */
      /* std::cout << "dtau = " << dtau << std::endl; */

      //compute new step size
      double ainv;
      ainv = dtau_H(rrpt,tpt,ppt,thisatmointerp); // mean scattering length
      
      ds = dtau/ainv; // new distance computed from dtau and
      // local scattering cross section
      /* std::cout << "ds = " << ds << std::endl; */
      /* std::cout << "dsmin = " << dsmin << std::endl; */
      /* std::cout << "dsmax = " << dsmax << std::endl; */
      ds = ds < dsmin ? dsmin : ds; // keep the distance
      ds = ds > dsmax ? dsmax : ds; // within the bounds
      /* std::cout << "ds = " << ds << std::endl; */
      /* std::cin.get(); */
	    
      //compute the optical depth
      tauH += ainv*ds;
      tauCO2 += dtau_CO2(rrpt,tpt,ppt,thisatmointerp)*ds;

      intensity += HolTinterp.interp(tauH)*exp(-tauCO2)*S(rrpt,tpt,ppt)*ds;

      // move out along the unit vector:	      
      s += ds;
      xpt += line_x*ds;
      ypt += line_y*ds;
      zpt += line_z*ds;

      // compute the new radius and angle
      rrpt = sqrt(xpt*xpt + ypt*ypt + zpt*zpt);
      tpt = std::acos(xpt/rrpt);
      ppt = std::atan2(ypt, zpt);
      ppt = ppt < 0 ? ppt+2*pi : ppt;// put phi in the right domain
    }

    //multiply by factors which are constant outside the integral
    intensity *= lineintcoef;

    //echo the total optical depth to the calling process
    tauHout = tauH;
    tauCO2out = tauCO2;
  
    return intensity;
  }

  //overloaded function in case taus are not desired output
  double integrate(Sobj &S, atmointerp &thisatmointerp,
		   const VecDoub &pos, const VecDoub &dir,
		   const double &g)
  {
    double taudum=0.0;

    return integrate(S, thisatmointerp, pos, dir, g, taudum, taudum);
  }
};


#endif
