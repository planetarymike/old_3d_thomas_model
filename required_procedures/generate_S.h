//generate_S.h -- code to generate H source functions according to the
//                formalism of Anderson & Hord (1976).

#ifndef __generate_S_h
#define __generate_S_h

#include <ctime> // timing
#include <iostream> // for file output and dialog
#include <fstream>  // ""
#include <cmath>    // for cos and sin
#include <cstdlib>
#include "physical.h" // physical definitions, nH, nCO2, KK, etc.
#include "quadrature.h" // to compute line integrals
#include "definitions.h" // basic parameter definitions
#include "nr3.h" // type VecDoub
#include "rad_trans.h" // two-point column density, holstein t and g functions
#include "los.h" //Holstein lookups, column intergrals
#include "ludcmp.h" // LU decomposition solver
#include "interpgen.h" // runtime hydrogen interpolation

using std::cos;
using std::sin;

string Sfilename(double nH, double T) {
  string fname;
  char cfname[200];
  sprintf(cfname,"S_nH%G_T%G_%dx%dx%d.dat", nH, T, nrpts, nthetapts, nphipts);
  fname=srcfnsloc+cfname;
  return fname;
}

void generate_S(const double nexo, 
		const double Texo, 
		Sobj &Sout,
		atmointerp &atmoout,
		bool printout = 1, 
		bool filesout = 1) {

  std::cout << "#####################################################" << std::endl;
  std::cout << "#####################################################" << std::endl;
  std::cout << "##########KERNEL GENERATION NOW IN PROGRESS##########" << std::endl;
  std::cout << "#####################################################" << std::endl;
  std::cout << "#####################################################" << std::endl;


  //get nexo and Texo from the command line
  std::cout << "nexo = " << nexo << std::endl;
  std::cout << "Texo = " << Texo << std::endl;
  double sH0=sH(Texo);//line center H cross section

  //get the Holstein T function
  //from pre-tabulated file
  Holinterp HolT_lookup(HolTfilename);
  
  //start timing
  clock_t start = clock(); // start whole task timer
    
  //define the nH interpolation function
  atmointerp thisatmointerp(nexo, Texo, nphyspts);

  //define vectors for the quadrature regions and weights
  VecDoub rpts(nrpts), rwts(nrpts);
  VecDoub thetapts(nthetapts), thetawts(nthetapts);
  VecDoub phipts(nphipts), phiwts(nphipts);
  //and the ray directions and weights
  VecDoub raytpts(ntrays), raytwts(ntrays);
  VecDoub rayppts(nprays), raypwts(nprays);
  //now obtain the quadrature points in each dimension:
  getquadpts(rpts,rwts,thetapts,thetawts,phipts,phiwts,printout,thisatmointerp,rmethod);
  getraypts(raytpts,raytwts,rayppts,raypwts,printout,raymethod);

  //  std::cin.get();

  Linear_interp rterp(rpts,rpts);
  Linear_interp tterp(thetapts,thetapts);
  Linear_interp pterp(phipts,phipts);

  // what's the smallest altitude in the grid?
  double rmin = rpts[nrpts-1];
  //rmax is defined in definitions.h as an input parameter
  //rmin is a consequence of rminatm
  
  // variables for matrix manipulation
  // source function is computed for regions, not region boundaries, hence the -1
  const int nrows = (nrpts-1)*(nthetapts-1)*(nphipts-1);
  const int ncols = nrows;
  std::cout << "nrows = " << nrows 
	    << ", ncols = " << ncols << std::endl;

  double** kvals;
  kvals = new double*[nrows];
  for (int i = 0; i < nrows; i++) {
    kvals[i] = new double[ncols];
    //    std::cout << "kvals[" << i << "] = " << kvals[i] << std::endl;
  }
  double yvec[nrows]; // single scattering values

  //create file objects for output
  ofstream solfile;
  string solfname;
  solfname = Sfilename(nexo,Texo);
  solfile.open(solfname.c_str());

  // auxiliary variables used in computing integrals
  double x1, y1, z1, xpt, ypt, zpt; // ephemeral cartesian coordinates
  double rr1, t1, p1; //ephemeral spherical coordinates
  double tray, pray; // working angles for each point
  double line_x,line_y,line_z; //unit vector along the line of interest
  double rrpt, tpt, ppt; // spherical coordinates of the auxiliary grid point
  double s, ds, si, sf; // ephemeral distance and distance across volume
  double tauH, tauCO2; // optical depth
  double domega; // differential solid angle
  double omega; // check parameter the make sure sum(domega) = 4*pi
  double coef; // influence coeffecient computed between atmospheric
	       // point and auxiliary grid points

  // variables for single scattering function
  VecDoub r1(3), r2(3); // vectors for line integration
  double tauHcol, tauCO2col;


  //lots of loop and index variables
  int row,col,irw,itw,ipw,irs,its,ips;//loop variables
  double dr, dt, dp;//size of the current grid box
  int iter; // number of iterations in radial expansion
  int ordx, otdx, opdx, rdx, tdx, pdx; // volume indices
  
  for (row = 0; row < nrows; row++) {
    //initialize the kernel matrix:
    for(col = 0; col < ncols; col++)
      kvals[row][col] = 0.0;
    //    std::cout << "kvals[" << row << "] = " << kvals[row] << std::endl;

    //get the r, theta, phi point corresponding to this row
    irw = row/((nthetapts-1)*(nphipts-1));
    itw = (row-irw*(nthetapts-1)*(nphipts-1))/(nphipts-1);
    ipw = row%(nphipts-1);

    // keep count of the computation
    time_t rawtime;
    struct tm * timeinfo;
    time (&rawtime);
    timeinfo = localtime(&rawtime);
    std::cout << irw << ": "
    	      << itw << ": "
    	      << ipw << " starting at " << asctime(timeinfo);
    //    std::cin.get();

    //get coordinates for the wedge point under consideration
    rr1 = (rpts[irw]+rpts[irw+1])/2;
    t1 = (thetapts[itw]+thetapts[itw+1])/2;
    p1 = (phipts[ipw]+phipts[ipw+1])/2;
    x1 = rr1*sin(t1)*cos(p1);
    y1 = rr1*sin(t1)*sin(p1);
    z1 = rr1*cos(t1);

    // std::cout << "r1 = " << rr1 << std::endl;
    // std::cout << "t1 = " << t1 << std::endl;
    // std::cout << "p1 = " << p1 << std::endl;

    //initialize omega
    //    omega = 0.0;

    //now move out from this point along the sightline grid:
    for (its = 0; its < ntrays; its++) {
	for (ips = 0; ips < nprays; ips++) {
	    //obtain the angles that we're working outward along:
	    tray = raytpts[its];
	    pray = rayppts[ips];

	    // std::cout << "tray = " << tray << std::endl;
	    // std::cout << "pray = " << pray << std::endl;

	    //get the solid angle defined by these points:
	    domega = sin(tray)*raytwts[its]*raypwts[ips]/(4*pi);
	    //	    omega += domega;

	    
	    //get the unit vector components along the angle outward
	    line_x = cos(p1)*cos(pray)*cos(t1)*sin(tray)
	              +cos(p1)*cos(tray)*sin(t1)
	              -sin(tray)*sin(pray)*sin(p1);
	    line_y = cos(tray)*sin(p1)*sin(t1)
	              +sin(tray)*cos(pray)*sin(p1)*cos(t1)
	              +sin(tray)*sin(pray)*cos(p1);
	    line_z = cos(tray)*cos(t1)-cos(pray)*sin(tray)*sin(t1);

	    /* std::cout << "line = { " << line_x */
	    /* 	      << " , " << line_y << " , " */
	    /* 	      << line_z << "}\n"; */
	    // std::cin.get();
	    
	    // initialize the variables
	    s = 0.0;// s = 0 at start of integration
	    rrpt = rr1;// first point is grid point itself
	    tpt = t1;
	    ppt = p1;

	    coef = 0.0;// this changes soon
	    iter = 0;

	    xpt = x1;
	    ypt = y1;
	    zpt = z1;
	    /* std::cout << "rrpt = " << rrpt  */
	    /* 	      << ", rmin = " << rmin  */
	    /* 	      << ", rmax = " << rmax << std::endl; */

	    while ((rrpt >= rmin && rrpt <= rmax && iter < maxit) || iter == 0) {
	      //while still inside the atmosphere
	      iter++;
	      // if too many iterations, report error, but proceed
	      // with calculation
	      if (iter == maxit) {
		std::cout << "On row " << row
			  << " iter > maxit!" << std::endl;
		throw("Max iterations exceeded!")
	      }

	      //get the starting box coordinates
	      ordx=rterp.index(rrpt);
	      rdx=ordx;
	      otdx=tterp.index(tpt);
	      tdx=otdx;
	      opdx=pterp.index(ppt);
	      pdx=opdx;

	      // std::cout << "Ray tray = " << tray << ", pray = " << pray << ":\n"
	      // 		<< "Starting in box {" << ordx << ", " 
	      // 		<< otdx << ", " 
	      // 		<< opdx << "}\n";
	      
	      /* std::cout <<"Here, \n" */
	      /* 	        <<" r=[" << rpts[rdx]  */
	      /* 	        << ", " << rpts[rdx+1]  */
	      /* 	        << "]" << std::endl; */
	      /* std::cout <<" t=[" << thetapts[tdx]  */
	      /* 		<< ", " << thetapts[rdx+1]  */
	      /* 		<< "]" <<std::endl; */
	      /* std::cout <<" p=[" << phipts[pdx]  */
	      /* 		<< ", " << phipts[pdx+1]  */
	      /* 		<< "]" << std::endl; */

	      //get the step size for this box
	      dr=abs(rpts[rdx+1]-rpts[rdx]);
	      dt=rpts[rdx]*abs(thetapts[tdx+1]-thetapts[tdx]);
	      dp=rpts[rdx]*sin((thetapts[tdx+1]+thetapts[tdx])/2.)
		          *abs(phipts[pdx+1]-phipts[pdx]);

	      ds = dsfrac*dr;
	      ds = ds < dsfrac*dt ? ds : dsfrac*dt;
	      ds = ds < dsfrac*dp ? ds : dsfrac*dp;
	      
	      // std::cout << "dr = " << dr << std::endl;
	      // std::cout << "dt = " << dt << std::endl;
	      // std::cout << "dp = " << dp << std::endl;
	      // std::cout << "ds = " << ds << std::endl;

	      si=s;
	      while (rdx==ordx&&tdx==otdx&&pdx==opdx) {
		//while still in the same box
		//move out along the unit vector:	      
		s += ds;
		xpt += line_x*ds;
		ypt += line_y*ds;
		zpt += line_z*ds;
		
		// compute the new radius and angle
		rrpt = sqrt(xpt*xpt + ypt*ypt + zpt*zpt);
		tpt = std::acos(zpt/rrpt);
		ppt = std::atan2(ypt,xpt);
		
		rdx=rterp.index(rrpt);
		tdx=tterp.index(tpt);
		pdx=pterp.index(ppt);
		
		if (rrpt>rmax||rrpt<rmin)  break;
	      }
	      sf=s;

	      // std::cout << "si = " << si << std::endl;
	      // std::cout << "sf = " << sf << std::endl;
	      
	      // std::cout << "Ending in box: " << rdx << ", " 
	      // 		<< tdx << ", " 
	      // 		<< pdx << ".\n";
	      // std::cin.get();
	      
	      
	      //now get the optical depths
	      coef=sH0*thisatmointerp.nH((rpts[ordx]+rpts[ordx+1])/2,
					   (thetapts[otdx]+thetapts[otdx+1])/2,
					   (phipts[opdx]+phipts[opdx+1])/2);
	      // std::cout << "coef (sHtot*nH) = " << coef << std::endl;

	      coef=HolT_lookup.interp(si*coef)-HolT_lookup.interp(sf*coef);
	      // std::cout << "coef (HolTi-HolTf) = " << coef << std::endl;


	      coef*=domega;
	      //need to add CO2 absorption
	      // std::cout << "coef (HolTi-HolTf)*domega = " << coef << std::endl;
	      
	      coef*=exp(-(sf-si)*sCO2*thisatmointerp.nCO2((rpts[ordx]+rpts[ordx+1])/2,
	      					  (thetapts[otdx]+thetapts[otdx+1])/2,
	      					  (phipts[opdx]+phipts[opdx+1])/2));

	      // std::cout << "coef (HolTi-HolTf)*domega*exp(-tCO2) = " << coef << std::endl;
		
	      col=ordx*(nthetapts-1)*(nphipts-1)+otdx*(nphipts-1)+opdx;
	      kvals[row][col]+=coef;
	      // std::cout << "now kvals[" << row << "][" << col <<"] = " 
	      //  		<< kvals[row][col] << std::endl;
	      // std::cin.get();
	    }
	}
    }

    //    std::cout << irw << ": " << itw << ": " << ipw 
    //              << " omega/4*pi = " << omega << std::endl;

    // keep count of the computation
    time (&rawtime);
    timeinfo = localtime(&rawtime);
    std::cout << irw << ": " << itw << ": "<< ipw << " ending at " 
	      << asctime(timeinfo);


    //now compute the single scattering function:
    std::cout << "Getting single scattering at point " << irw << ": "
	      << itw << ": " << ipw << std::endl;

    //see if the point is behind the planet
    if (t1>pi/2&&rr1*sin(t1)<rmin) {
      std::cout << "Point " << irw << ": "
		<< itw << ": " << ipw << " is behind the limb." << std::endl;
      yvec[row] = 0.0;
    } else {
      //compute
      //      std::cout << "r1 = {" << x1 << ", " << y1 << ", " << z1 << "}\n";
      // construct the point at infinity:
      xpt = x1;
      ypt = y1;
      zpt = 2*rmax;
      //      std::cout << "rpt = {" << xpt << ", " << ypt << ", " << zpt << "}\n";

      // store the values in the vector object r2
      rr1 = (rpts[irw]+rpts[irw+1])/2;
      r1[0] = x1; r1[1] = y1; r1[2] = z1;
      r2[0] = xpt; r2[1] = ypt; r2[2] = zpt;

      // get the column densities
      std::cout << "Getting CO2 optical depth...\n";
      dtau_CO2_int dtauCO2(thisatmointerp);
      tauCO2col = qlinetrap_2pt(dtauCO2, r1, r2, 1e-3);
      std::cout << "  ... CO2 optical depth = "<< tauCO2col << ".\n";

      std::cout << "Getting H optical depth...\n";
      dtau_H_int dtauH(thisatmointerp);
      tauHcol = qlinetrap_2pt(dtauH, r1, r2, 1e-3);
      std::cout << "  ... H optical depth = "<< tauHcol << ".\n";

      // put them together into the y-vec
      //std::cout << "irw = " << irw << std::endl;
      yvec[row] = thisatmointerp.nH(rr1)*HolT_lookup.interp(tauHcol)*std::exp(-tauCO2col);
      //    yvec[row] = 0.0;
    }
  }

  // solve using LU decomposition
  std::cout << "Solving for source function by matrix inversion: " << std::endl;
  MatDoub preinv(nrows,ncols);
  for (row = 0; row < nrows; row++)
    for (col = 0; col < ncols/*square*/; col++) 
      preinv[row][col] = row==col ? 1.0-kvals[row][col] : -kvals[row][col];

  VecDoub ynrvec(nrows,yvec);
  LUdcmp mat(preinv);
  VecDoub sol(nrows);
  mat.solve(ynrvec, sol);
  
  //check solution:
  double diff = 0.0;
  double comp, temp;
  for (row = 0; row < nrows; row++) {
    comp = 0.0;
    for(col = 0; col < nrows/*square*/; col++)
      comp += kvals[row][col]*sol[col];
    comp += yvec[row];
    temp = sol[row] - comp;
    diff += temp*temp;
  }
  std::cout << "|S_inv - ( y + K.S_inv )| = diff = " << sqrt(diff) << std::endl;

  //write solution to file:
  std::cout << "Writing solution to file." << std::endl;
  //nH and T
  solfile.width(20);
  solfile << nexo;
  solfile.width(20);
  solfile << Texo << std::endl;
  //rpts, tpts, ppts
  solfile.width(20);
  solfile << nrpts << std::endl;
  for (irw = 0; irw < nrpts; irw++) {
    solfile.width(20);
    solfile << rpts[irw];
  }
  solfile << std::endl;
  solfile.width(20);
  solfile << nthetapts << std::endl;
  for (itw = 0; itw < nthetapts; itw++) {
    solfile.width(20);
    solfile << thetapts[itw];
  }
  solfile << std::endl;
  solfile.width(20);
  solfile << nphipts << std::endl;
  for (ipw = 0; ipw < nphipts; ipw++) {
    solfile.width(20);
    solfile << phipts[ipw];
  }
  solfile << std::endl;
  //finally, the values
  for (row = 0; row < nrows; row++) {
    solfile.width(20);
    solfile << sol[row] << std::endl;
  }
    
  // print time elapsed
  clock_t finish = clock();
  double secs = (finish-start)*1.0/CLOCKS_PER_SEC;
  int hrs = secs/3600;
  int mins = (secs - hrs*3600)/60;
  secs = secs - hrs*3600 - mins*60;
  std::cout << "Elapsed time is " 
	    << hrs << " hours, " 
	    << mins << " minutes, and " 
	    << secs << " seconds.\n";

  //return values to caller
  Sout = Sobj(nexo, Texo, rterp, tterp, pterp, sol);
  atmoout = thisatmointerp;

  //Cleanup:
  solfile.close();

  //free allocated kernel matrix
  for (int i = 0; i < nrows; i++) {
    delete [] kvals[i];
  }
  delete [] kvals;

  return;
}


#endif
