//iph_sim.h -- routines for simulating the IPH


#ifndef __IPH_SIM_H
#define __IPH_SIM_H

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include "definitions.h"
#include "interp_1d.h"

using std::log;
using std::exp;
using std::pow;

/* first, a routine to simulate the IPH line shape, using parameters
   from Quemerais' 1999 paper using SWAN/SOHO data. */
void iph_shape(double ra, double dec, VecDoub &iphparms) {
  /* using an input RA and dec, compute the shape and center of the
     IPH line */
  /* double elon0 = pi/180.*252.3; */
  /* double elat0 = pi/180.*8.7; */
  double ra0 = pi/180.*251.98; //origin of the flow in RA, Dec
  double dec0 = pi/180.*-13.697;
  
  ra=pi/180.*ra;
  dec=pi/180.*dec;

  /* the angle between these two directions is given by this formula */
  /* (take the dot product of the two and use a trig formula.) */
  double angle=180./pi*acos(cos(ra - ra0)*cos(dec)*cos(dec0) + sin(dec)*sin(dec0));
  // std::cout << "angle = " << angle << std::endl;

  /* all of the IPH properties can be expressed in terms of this angle
     to a high degree of accuracy. */
  /* for example, the relative velocity, from Quemerais 1999: */
  double iphvelcoords[37]={0., 4.86006, 9.89382, 14.9276, 19.9613, 25.1688, 30.0289, 35.0626, 40.0964, 45.1301, 50.3373, 55.1977, 60.2314, 65.0915, 70.1252, 75.159, 80.0191, 85.0528, 90.4337, 94.9469, 99.9807, 105.014, 110.048, 114.735, 119.942, 130.01, 124.976, 134.87, 139.903, 144.764, 149.971, 155.005, 160.038, 164.899, 170.106, 174.966, 180.};
  double iphvel[37]={-25.3333, -25.1667, -25., -24.6667, -24.1667, -23.3334, -22.4167, -21.25, -20.25, -19.3333, -17.9166, -16.4167, -14.5834, -12.8334, -10.9999, -8.83333, -6.83337, -4.83329, -2.66668, -0.666594, 1.41676, 3.33331, 5.33339, 7.16669, 9.08337, 12.6667, 11.0001, 14.3333, 15.6667, 16.8334, 18., 19.1667, 20.0834, 21.0833, 21.25, 21.4167, 21.5001};
  /* also the temperature, which has fewer reported values. */
  double iphtempcoords[3]={0.,90.,180.};
  double iphtemp[3]={14417,13920,19614};

  /* Now find the interpolated parameters for this ra, dec */
  int ivel=0;
  while (iphvelcoords[ivel]<angle)
    ivel++;
  double thisiphvel= (iphvelcoords[ivel]-angle)*iphvel[ivel-1]
               +(angle-iphvelcoords[ivel-1])*iphvel[ivel];
  thisiphvel/=(iphvelcoords[ivel]-iphvelcoords[ivel-1]);

  int itemp=0;
  while (iphtempcoords[itemp]<angle)
    itemp++;
  double thisiphtemp= (iphtempcoords[itemp]-angle)*iphtemp[itemp-1]
               +(angle-iphtempcoords[itemp-1])*iphtemp[itemp];
  thisiphtemp/=(iphtempcoords[itemp]-iphtempcoords[itemp-1]);

  iphparms.resize(2);
  iphparms[0]=thisiphvel;
  iphparms[1]=thisiphtemp;

  return;
}

struct LOSprofinterp {
  //read simulated exospheric line profiles from file.
  int nlc, nlbfrac, nt0;
  int ilc,ilbfrac,it0,ia;
  string filename;
  VecDoub lc;
  VecDoub lbfrac;
  VecDoub t0;
  double *****profvec;//array of line shapes
  int ***na;//array of line shape sizes
  Linear_interp lc_terp;
  Linear_interp lbfrac_terp;
  Linear_interp t0_terp;
  bool init;

  LOSprofinterp() {init = 0;}
  LOSprofinterp(string filenamee) : filename(filenamee) {

    double dum;
  
    ifstream infile;
    infile.open(filename.c_str());

    //read in the number of gridpoints
    infile >> nlc;
    lc.resize(nlc);
    infile >> nlbfrac;
    lbfrac.resize(nlbfrac);
    infile >> nt0;
    t0.resize(nt0);

    infile.get();
    //now read in the coordinates
    for (ilc=0;ilc<nlc;ilc++) {
      infile >> lc[ilc];
      /* std::cout << "lc["<<ilc<<"] = "<<lc[ilc]<<std::endl; */
      /* std::cin.get(); */
    }
    for (ilbfrac=0;ilbfrac<nlbfrac;ilbfrac++) {
      infile >> lbfrac[ilbfrac];
      /* std::cout << "lbfrac["<<ilbfrac<<"] = "<<lbfrac[ilbfrac]<<std::endl; */
      /* std::cin.get(); */
    }
    for (it0=0;it0<nt0;it0++) {
      infile >> t0[it0];
      /* std::cout << "t0["<<it0<<"] = "<<t0[it0]<<std::endl; */
      /* std::cin.get(); */
    }

    //now we can allocate the pointer arrays for a and for tau
    infile.get();

    profvec=new double****[nlc];
    na=new int**[nlc];
    for (ilc=0;ilc<nlc;ilc++) {
      profvec[ilc]=new double***[nlbfrac];
      na[ilc]=new int*[nlbfrac];
      for (ilbfrac=0;ilbfrac<nlbfrac;ilbfrac++) {
	profvec[ilc][ilbfrac]=new double**[nt0];
	na[ilc][ilbfrac]=new int[nt0]; 
	for (it0=0;it0<nt0;it0++) {
	  /* clear the parameters */
	  infile >> dum;
	  infile >> dum;
	  infile >> dum;
	  infile >> na[ilc][ilbfrac][it0];
	  /* std::cout << "na = " << na[ilc][ilbfrac][it0] << std::endl; */
	  profvec[ilc][ilbfrac][it0] = new double*[na[ilc][ilbfrac][it0]];
	  /* read in the velocities */
	  for (ia = 0; ia < na[ilc][ilbfrac][it0]; ia++) {
	    profvec[ilc][ilbfrac][it0][ia] = new double[2];
	    infile >> profvec[ilc][ilbfrac][it0][ia][0];
	  }
	  /* Now the taus */
	  for (ia = 0; ia < na[ilc][ilbfrac][it0]; ia++) {
	    infile >> profvec[ilc][ilbfrac][it0][ia][1];
	    //profvec[ilc][ilbfrac][it0][ia][1]*=lc[ilc]*lc[ilc]/lbfrac[ilbfrac]*exp(-lc[ilc]);
	    //from an old, obsolete version of the code.
	    /* std::cout << "a = " << profvec[ilc][ilbfrac][it0][ia][0] << std::endl; */
	    /* std::cout << "tau = " << profvec[ilc][ilbfrac][it0][ia][1] << std::endl; */
	    /* std::cin.get(); */
	  }
	}
      }
    }

    infile.close();

    //now create the interpolation objects
    lc_terp = Linear_interp(lc,lc);
    lbfrac_terp = Linear_interp(lbfrac,lbfrac);
    t0_terp = Linear_interp(t0,t0);

    init=1;
  }

  LOSprofinterp(const LOSprofinterp &LOS) : nlc(LOS.nlc), 
					    nlbfrac(LOS.nlbfrac), 
					    nt0(LOS.nt0),
					    ilc(LOS.ilc),
					    ilbfrac(LOS.ilbfrac),
					    it0(LOS.it0),
					    ia(LOS.ia),
					    filename(LOS.filename),
					    lc(LOS.lc),
					    lbfrac(LOS.lbfrac),
					    t0(LOS.t0),
					    lc_terp(LOS.lc_terp),
					    lbfrac_terp(LOS.lbfrac_terp),
					    t0_terp(LOS.t0_terp),
					    init(LOS.init) {
    profvec=new double****[nlc];
    na=new int**[nlc];
    for (ilc=0;ilc<nlc;ilc++) {
      profvec[ilc]=new double***[nlbfrac];
      na[ilc]=new int*[nlbfrac];
      for (ilbfrac=0;ilbfrac<nlbfrac;ilbfrac++) {
	profvec[ilc][ilbfrac]=new double**[nt0];
	na[ilc][ilbfrac]=new int[nt0]; 
	for (it0=0;it0<nt0;it0++) {
	  na[ilc][ilbfrac][it0] = LOS.na[ilc][ilbfrac][it0];
	  /* std::cout << "na = " << na[ilc][ilbfrac][it0] << std::endl; */
	  profvec[ilc][ilbfrac][it0] = new double*[na[ilc][ilbfrac][it0]];
	  /* read in the velocities */
	  for (ia = 0; ia < na[ilc][ilbfrac][it0]; ia++) {
	    profvec[ilc][ilbfrac][it0][ia] = new double[2];
	    profvec[ilc][ilbfrac][it0][ia][0]=LOS.profvec[ilc][ilbfrac][it0][ia][0];
	    profvec[ilc][ilbfrac][it0][ia][1]=LOS.profvec[ilc][ilbfrac][it0][ia][1];
	  }
	}
      }
    }
  }

  LOSprofinterp&operator=(const LOSprofinterp &LOS) {
    nlc=LOS.nlc;
    nlbfrac=LOS.nlbfrac;
    nt0=LOS.nt0;
    ilc=LOS.ilc;
    ilbfrac=LOS.ilbfrac;
    it0=LOS.it0;
    ia=LOS.ia;
    filename=LOS.filename;
    lc=LOS.lc;
    lbfrac=LOS.lbfrac;
    t0=LOS.t0;
    lc_terp=LOS.lc_terp;
    lbfrac_terp=LOS.lbfrac_terp;
    t0_terp=LOS.t0_terp;
    init=LOS.init;
    profvec=new double****[nlc];
    na=new int**[nlc];
    for (ilc=0;ilc<nlc;ilc++) {
      profvec[ilc]=new double***[nlbfrac];
      na[ilc]=new int*[nlbfrac];
      for (ilbfrac=0;ilbfrac<nlbfrac;ilbfrac++) {
	profvec[ilc][ilbfrac]=new double**[nt0];
	na[ilc][ilbfrac]=new int[nt0]; 
	for (it0=0;it0<nt0;it0++) {
	  na[ilc][ilbfrac][it0] = LOS.na[ilc][ilbfrac][it0];
	  /* std::cout << "na = " << na[ilc][ilbfrac][it0] << std::endl; */
	  profvec[ilc][ilbfrac][it0] = new double*[na[ilc][ilbfrac][it0]];
	  /* read in the velocities */
	  for (ia = 0; ia < na[ilc][ilbfrac][it0]; ia++) {
	    profvec[ilc][ilbfrac][it0][ia] = new double[2];
	    profvec[ilc][ilbfrac][it0][ia][0]=LOS.profvec[ilc][ilbfrac][it0][ia][0];
	    profvec[ilc][ilbfrac][it0][ia][1]=LOS.profvec[ilc][ilbfrac][it0][ia][1];
	  }
	}
      }
    }
    return *this;
  }



  
  ~LOSprofinterp() {
    if (init) {
      //      std::cout << "calling destructor" << std::endl;
      for (ilc=0;ilc<nlc;ilc++) {
	for (ilbfrac=0;ilbfrac<nlbfrac;ilbfrac++) {
	  for (it0=0;it0<nt0;it0++) {
	    for (ia = 0; ia < na[ilc][ilbfrac][it0]; ia++) {
	      delete [] profvec[ilc][ilbfrac][it0][ia];
	    }
	    delete [] profvec[ilc][ilbfrac][it0];
	  }
	  delete [] profvec[ilc][ilbfrac];
	  delete [] na[ilc][ilbfrac];
	}
	delete [] profvec[ilc];
	delete [] na[ilc];
      }
      delete [] profvec;
      delete [] na;
    }
  }

  void interp(double tlc, double tlb, double tt0, VecDoub & avec, VecDoub & tauvec) {
    /* returns an interpolated atmospheric absorption vector,
       including all factors except exobase density. Return vectors
       avec and tauvec have variable length depending on when the
       optical depth drops below threshold. */

    if (!init)
      throw("LOS prof interp not initialized!");

    double tlbfrac = tlb/tlc;

    int itlc = lc_terp.index(tlc);
    /* std::cout << "itlc = " << itlc << std::endl; */
    int itlbfrac = lbfrac_terp.index(tlbfrac);
        /* std::cout << "itlbfrac = " << itlbfrac << std::endl; */
    int itt0 = t0_terp.index(tt0);
    /* std::cout << "itt0 = " << itt0 << std::endl; */
    
    double alc = (tlc-lc_terp.xx[itlc])/(lc_terp.xx[itlc+1]-lc_terp.xx[itlc]);
    /* std::cout << "lc["<<itlc<<"] = " << lc_terp.xx[itlc] << std::endl; */
    /* std::cout << "lc["<<itlc+1<<"] = " << lc_terp.xx[itlc+1] << std::endl; */
    /* std::cout << "alc = " << alc << std::endl; */

    double albfrac = (tlbfrac-lbfrac_terp.xx[itlbfrac])/(lbfrac_terp.xx[itlbfrac+1]
							 -lbfrac_terp.xx[itlbfrac]);
    /* std::cout << "lbfrac["<<itlbfrac<<"] = " */
    /* 	      << lbfrac_terp.xx[itlbfrac] << std::endl; */
    /* std::cout << "lbfrac["<<itlbfrac+1<<"] = " */
    /* 	      << lbfrac_terp.xx[itlbfrac+1] << std::endl; */
    /* std::cout << "albfrac = " << albfrac << std::endl; */

    double at0 = (tt0-t0_terp.xx[itt0])/(t0_terp.xx[itt0+1]-t0_terp.xx[itt0]);
    /* std::cout << "t0["<<itt0<<"] = " << t0_terp.xx[itt0] << std::endl; */
    /* std::cout << "t0["<<itt0+1<<"] = " << t0_terp.xx[itt0+1] << std::endl; */
    /* std::cout << "at0 = " << at0 << std::endl; */

    double coef[2][2][2] = {{{ (1-alc)*(1-albfrac)*(1-at0), (1-alc)*(1-albfrac)*at0},
			     { (1-alc)*albfrac*(1-at0), (1-alc)*albfrac*at0}},
			    {{ alc*(1-albfrac)*(1-at0), alc*(1-albfrac)*at0},
			     { alc*albfrac*(1-at0), alc*albfrac*at0}}};

    
    // get the max width of the profiles, and the corresponding a values
    int namax=0;
    for (ilc=itlc;ilc<=itlc+1;ilc++) {
      for (ilbfrac=itlbfrac;ilbfrac<=itlbfrac+1;ilbfrac++) {
	for (it0=itt0;it0<=itt0+1;it0++) {
	  if (na[ilc][ilbfrac][it0] > namax) {
	    namax=na[ilc][ilbfrac][it0];
	    avec.resize(namax);
	    for (ia=0;ia<namax;ia++) {
	      avec[ia]=profvec[ilc][ilbfrac][it0][ia][0];
	    }
	  }
	}
      }
    }

    //now we can interpolate the taus

    //first, set up for output of tau
    tauvec.assign(namax,0.);

    int retmidpos=(namax-1)/2;
    int tmidpos;
    //interpolate
    for (ilc=0;ilc<=1;ilc++) {
      for (ilbfrac=0;ilbfrac<=1;ilbfrac++) {
	for (it0=0;it0<=1;it0++) {
	  tmidpos=na[ilc+itlc][ilbfrac+itlbfrac][it0+itt0];
	  tmidpos=(tmidpos-1)/2;
	  /* std::cout << "tmidpos = " << tmidpos << std::endl; */
	  /* std::cout << "lc = " << lc_terp.xx[ilc+itlc] << std::endl; */
	  /* std::cout << "lbfrac = " << lbfrac_terp.xx[ilbfrac+itlbfrac] << std::endl; */
	  /* std::cout << "t0 = " << t0_terp.xx[it0+itt0] << std::endl; */
	  /* std::cout << "coef = " <<  coef[ilc][ilbfrac][it0] << std::endl << std::endl; */
	  //work from the inside out
	  tauvec[retmidpos]+=(coef[ilc][ilbfrac][it0]*
			      profvec[ilc+itlc][ilbfrac+itlbfrac][it0+itt0][tmidpos][1]);
	  ia = 1;
	  while (ia<=tmidpos) {
	    /* std::cout << "ia = " << ia << std::endl; */
	    tauvec[retmidpos-ia]+=(
				   coef[ilc][ilbfrac][it0]*
				   profvec[ilc+itlc][ilbfrac+itlbfrac][it0+itt0][tmidpos-ia][1]
				   );
	    tauvec[retmidpos+ia]+=(
				   coef[ilc][ilbfrac][it0]*
				   profvec[ilc+itlc][ilbfrac+itlbfrac][it0+itt0][tmidpos+ia][1]
				   );
	    ia++;
	  }
	  //	  std::cin.get();
	}
      }
    }

    return;
  }

  double iph_absorbed(double iphoffset, double iphsigma, double nexo, 
		      double lc,double lb,double t0) 
  {
    /* Given an iph velocity offset and sigma in units of alpha, where
       alpha is the thermal velocity of H atoms at the Mars exobase,
       compute the fraction of the IPH absorbed by the mars
       atmosphere. */
    
    VecDoub avec, tvec;
    interp(lc,lb,t0,avec,tvec);

    if (avec.size() < 2) {//absorption is already below threshold
      return 0.0;
    } else {
      double absfrac=0.0;
      double ipharg;
      /* integrate across the Mars line using the trapezoidal rule.*/
      ipharg=(avec[0]-iphoffset)/iphsigma;
      ipharg*=ipharg;
      for (ia=0;ia<avec.size();ia++) {
	ipharg=(avec[ia]-iphoffset)/iphsigma;
	ipharg*=ipharg;
	absfrac+=exp(-ipharg)*(1.0-exp(-nexo*tvec[ia+1]));
      }
      absfrac*=(avec[1]-avec[0]);
      absfrac/=iphsigma*sqrt(pi);
      
      return absfrac;
    }
  }

  double iph_absorbed_radec(double iphra, double iphdec, double scvel,
			    double nexo, double Texo,
			    double lb, double t0) {
    /* Calculate the extinction of the IPH by the Mars corona, using a
       model to compute the width and doppler shift of the IPH line,
       and the specified spacecraft velocity to obtain the relative
       shape of the IPH line relative to that of the Mars corona. Then
       extinct by the correct amount along the line of sight
       specified by the impace parameter lambda and initial angle. */

    // std::cout << " ra = " << iphra << std::endl;
    // std::cout << " dec = " << iphdec << std::endl;
    
    double lc=G*mMars*mH/(kB*Texo*rexo);
    // std::cout << " lc = " << lc << std::endl;

    
    VecDoub avec, tvec;
    interp(lc,lb,t0,avec,tvec);

    VecDoub iphparms;
    iph_shape(iphra,iphdec,iphparms);

    //now we have the iph velocity offset and temperature width.
    //convert this to alpha space using the temperature at the exobase.
    double alpha=sqrt(2*kB*Texo/mH);//cm/s
    alpha/=1e5;//km/s

    double iphoffset=iphparms[0];
    iphoffset+=scvel;
    iphoffset/=alpha;
    // std::cout << "iphoffset = " << iphoffset << std::endl;
    double iphsigma=iphparms[1];
    iphsigma=sqrt(iphsigma/Texo);
    // std::cout << "iphsigma = " << iphsigma << std::endl;
    
    
    if (avec.size() < 2) {//absorption is already below threshold
      return 0.0;
    } else {
      double absfrac=0.0;
      double ipharg;
      /* integrate across the Mars line.*/
      ipharg=(avec[0]-iphoffset)/iphsigma;
      ipharg*=ipharg;
      absfrac+=0.5*exp(-ipharg)*(1.0-exp(-nexo*tvec[0]));
      for (ia=1;ia<avec.size()-1;ia++) {
	ipharg=(avec[ia]-iphoffset)/iphsigma;
	ipharg*=ipharg;
	absfrac+=exp(-ipharg)*(1.0-exp(-nexo*tvec[ia]));
      }
      ipharg=(avec[avec.size()-1]-iphoffset)/iphsigma;
      ipharg*=ipharg;
      absfrac+=0.5*exp(-ipharg)*(1.0-exp(-nexo*tvec[avec.size()-1]));

      absfrac*=(avec[1]-avec[0]);//assumes uniform grid
      absfrac/=iphsigma*sqrt(pi);//normalize to full width of line to get absorption fraction
      
      return absfrac;
    }
  }

  double iph_absorbed_radec_scvec(double iphra, double iphdec, double scvel,
				  double nexo, double Texo,
				  double* pos, double* dir) {
    /* Calculate the extinction of the IPH by the Mars corona,
       using the cartesian position and look direction of the
       spacecraft. */
    double lb;
    double t0;
    double s = 0;
    double mrh;
    double scr=0.0;
    double tr=0.0, td=0.0;
    for (int j=0;j<3;j++) {
      tr=pos[j];
      td=dir[j];
      s += -tr*td;
      scr+=tr*tr;
      // std::cout << "tr = " << tr << std::endl;
    }
    scr=sqrt(scr);
    mrh=0.0;
    for (int j=0;j<3;j++) {
      tr=pos[j]+s*dir[j];
      mrh+=tr*tr;
    }
    mrh=sqrt(mrh);
    
    lb=G*mMars*mH/(kB*Texo*mrh);
    t0=asin(-s/scr);

    return iph_absorbed_radec(iphra, iphdec, scvel, nexo, Texo, lb, t0);
  }

};

struct IPHsim {
  LOSprofinterp LOS_atm;
  int nobs;
  VecDoub ra, dec, scvel;
  MatDoub pos, dir;
  VecDoub lbtt;//lbtan, Jeans parameter at MRH point
  VecDoub t0;//theta0, initial angle through corona
  bool init;

  //default constructor
  IPHsim() {init=0;}
  
  //constructor
  IPHsim(int nobss, 
	 VecDoub &raa, 
	 VecDoub &decc, 
	 VecDoub &scvell, 
	 MatDoub &poss, 
	 MatDoub &dirr) : nobs(nobss), 
			  ra(raa), 
			  dec(decc), 
			  scvel(scvell), 
			  pos(poss), 
			  dir(dirr)
  {
    std::cout << "Loading IPH simulator coronal line profiles..." << std::endl;
    //load the los data
    LOS_atm=LOSprofinterp(losproffname);
			    
    std::cout << "...line profiles loaded." << std::endl;

    //resize the vectors
    lbtt.resize(nobs);
    t0.resize(nobs);

    //copy the vectors over, getting the MRH radius and LOS initial
    //angle
    double s;
    double mrh;
    for (int iobs = 0; iobs<nobs; iobs++) {
      s=0.0;
      double scr=0.0;
      double tr=0.0, td=0.0;
      for (int j=0;j<3;j++) {
	tr=pos[iobs][j];
	td=dir[iobs][j];
	s += -tr*td;
	scr+=tr*tr;
	// std::cout << "tr = " << tr << std::endl;
      }
      scr=sqrt(scr);
      mrh=0.0;
      for (int j=0;j<3;j++) {
	tr=pos[iobs][j]+s*dir[iobs][j];
	mrh+=tr*tr;
      }
      mrh=sqrt(mrh);

      lbtt[iobs]=G*mMars*mH/(kB*mrh);
      t0[iobs]=asin(-s/scr);

      // std::cout << "For iobs = " << iobs << std::endl;
      // std::cout << "  s = " << s << std::endl;
      // std::cout << "  scr = " << scr << std::endl;
      // std::cout << "  mrh = " << mrh << std::endl;
      // std::cout << "  lbtt[iobs] = " << lbtt[iobs] << std::endl;
      // std::cout << "  t0[iobs] = " << t0[iobs] << std::endl;
      // std::cin.get();

    }
    init = 1;
  }

  double sim(double IPH_brightness, double nexo, double Texo, int iobs) {
    if (!init)
      throw("IPH simulator must be initialized with observation data!");
    // std::cout << "For iobs = " << iobs << std::endl;
    // std::cout << " IPH_brightness = " << IPH_brightness << std::endl;
    // std::cout << " nexo = " << nexo << std::endl;
    // std::cout << " Texo = " << Texo << std::endl;
    // std::cout << " ra = " << ra[iobs] << std::endl;
    // std::cout << " dec = " << dec[iobs] << std::endl;
    // VecDoub iphparms(2);
    // iph_shape(ra[iobs],dec[iobs],iphparms);
    // std::cout << "For this RA, Dec:\n";
    // std::cout << "IPH vel = " << iphparms[0] << std::endl;
    // std::cout << "IPH temp = " << iphparms[1] << std::endl;
    // std::cout << " scvel = " << scvel[iobs] << std::endl;
    // std::cout << " lb = " << lbtt[iobs]/Texo << std::endl;
    // std::cout << " t0 = " << t0[iobs] << std::endl;

    //if we are below the atmo minimum or above the maximum, don't
    //perform the calculation.
    double b = G*mMars*mH/kB/lbtt[iobs]; 
    double absfrac;
    if (b < rminatm) 
      absfrac = 1.0;
    else if (b > rmax)
      absfrac = 0.0;
    else
      absfrac=LOS_atm.iph_absorbed_radec(ra[iobs], 
					 dec[iobs],
					 scvel[iobs],
					 nexo, 
					 Texo,
					 lbtt[iobs]/Texo,
					 t0[iobs]);
    
  // std::cout << " absfrac = " << absfrac << std::endl;
  // std::cin.get();
  

  return IPH_brightness*(1.0-absfrac);
  }
};


  
#endif
