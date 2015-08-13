//generate_S.cpp -- get influence matrix coefficents and compute the source function

// code adapted from: (but modified to match methods of chaufray's Calcul_Matrice.f90)
// opticaldepth_parallel_noinput.cpp -- get point to point optical depths for a grid of points around Mars
// altered to run in parallel using openmp
// altered to remove all inputs from the command line for server run 10 Jul 2011
// altered to run in parallel using openmpi 12 Jul 2011
// altered to sum over phi in the code without wasting memory 13 Jul 2011
// altered to compute the single scattering term



//compile: mpicxx generate_S.cpp `pkg-config --cflags --libs gsl` -I./required_procedures -o generate_S.x

//call: mpirun -n 4 ./generate_S.x 100000 350 //or whatever

#include <mpi.h> // parallelization
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
#include "pipeline_los_defs.h" // definitions for LOS integration
#include "los.h" //Holstein lookups, column intergrals
#include "ludcmp.h" // LU decomposition solver
#include "interpgen.h" // runtime hydrogen interpolation

using std::cos;
using std::sin;


int main(int argc, char* argv[]) {
  int rank,size;

  MPI::Init(argc,argv);
  rank=MPI::COMM_WORLD.Get_rank();
  size=MPI::COMM_WORLD.Get_size();

  if (rank == master) {
    std::cout << "#####################################################" << std::endl;
    std::cout << "#####################################################" << std::endl;
    std::cout << "##########KERNEL GENERATION NOW IN PROGRESS##########" << std::endl;
    std::cout << "#####################################################" << std::endl;
    std::cout << "#####################################################" << std::endl;
  }

  //get nexo and Texo from the command line
  double nexo = atof(argv[1]);
  if (rank == master) std::cout << "nexo = " << nexo << std::endl;
  double Texo = atof(argv[2]);
  if (rank == master) std::cout << "Texo = " << Texo << std::endl;
  double sH0=sH(Texo);//line center H cross section


  
  //start timing
  clock_t start = clock(); // start whole task timer
  int secwait = CLOCKS_PER_SEC;


  //______________________________________________________
  //------------------CONTROL PARAMETERS------------------
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
  //number of interpolation points in the hydrogen atmosphere
  const int nphyspts = 10000;//> 366!! or segfault!

  //smallest resolvable length
  const double dsfrac = 0.05;//fraction of the smallest box dimension

  //max and min path length
  const double dsmin = 0.1e5; // 0.1km
  const double dsmax = 100e5; // 100km
  
  //maximum number of iterations
  const int maxit = 100000;

  //define the number of boundaries in each dimension
  const int nrpts = 15;//includes boundaries at rmin and rmax
  const int nthetapts = 37;//includes boundaries at 0 and 180
  const int nphipts = 2;//includes boundaries at 0 and 360
  
  //number of rays
  const int ntrays = 6;
  const int nprays = 12;

  //use chaufray hydrogen density?
  bool usechauH = FALSE;

  //select the method for distributing the radial points in
  //altitude. Options are: taupts, charpts, logpts, loglinpts.
  string rmethod = "taupts";

  //select the method for distributing the rays in angle.
  //default is gauss-jordan quadrature in theta, uniform in phi
  string raymethod = "gaussian";
  
  //define the HolG interpolation function
  Holinterp HolG_lookup("./tabulated_data/HolGinterp.dat");
  Holinterp HolT_lookup("./tabulated_data/HolTinterp.dat");
  
  
  //some radiative transfer and atmospheric parameters are defined in definitions.h

  //______________________________________________________
  //----------------END CONTROL PARAMETERS----------------
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 

  //wait 1 sec
  clock_t timer = clock(); // start timing
  while (clock() - timer < secwait)
    ;
  MPI_Barrier(MPI_COMM_WORLD);

  //define the nH interpolation function
  static const atmointerp thisatmointerp(nexo, Texo, nphyspts);

  //define vectors for the quadrature regions and weights
  VecDoub rpts(nrpts), rwts(nrpts);
  VecDoub thetapts(nthetapts), thetawts(nthetapts);
  VecDoub phipts(nphipts), phiwts(nphipts);
  VecDoub raytpts(ntrays), raytwts(ntrays);
  VecDoub rayppts(nprays), raypwts(nprays);
  //now obtain the quadrature points in each dimension:
  bool printout = 1;
  //only one process needs to calculate the quadrature points
  //swap variables
  double rptr[nrpts],     rwtptr[nrpts],
         tptr[nthetapts], twtptr[nthetapts], 
         pptr[nphipts],   pwtptr[nphipts],
         rtptr[ntrays],   rtwtptr[ntrays],
         rpptr[nprays],   rpwtptr[nprays];
  if (rank == master) {
    //get the points
    getquadpts(rpts,rwts,thetapts,thetawts,phipts,phiwts,printout,thisatmointerp,rmethod);
    getraypts(raytpts,raytwts,rayppts,raypwts,printout,raymethod);
    //assign them to the swap variables
    for (int n = 0; n < nrpts; n++) { rptr[n] = rpts[n]; rwtptr[n] = rwts[n]; }
    for (int n = 0; n < nthetapts; n++) { tptr[n] = thetapts[n]; twtptr[n] = thetawts[n]; }
    for (int n = 0; n < nphipts; n++) { pptr[n] = phipts[n]; pwtptr[n] = phiwts[n]; }
    for (int n = 0; n < ntrays; n++) { rtptr[n] = raytpts[n]; rtwtptr[n] = raytwts[n]; }
    for (int n = 0; n < nprays; n++) { rpptr[n] = rayppts[n]; rpwtptr[n] = raypwts[n]; }
    //send to the other processes
    for (int proc = 0; proc < size; proc++) {
      if (proc != master) {
	MPI_Send(rptr,nrpts,MPI_DOUBLE,proc,0,MPI_COMM_WORLD);
	MPI_Send(rwtptr,nrpts,MPI_DOUBLE,proc,1,MPI_COMM_WORLD);	
	MPI_Send(tptr,nthetapts,MPI_DOUBLE,proc,2,MPI_COMM_WORLD);
	MPI_Send(twtptr,nthetapts,MPI_DOUBLE,proc,3,MPI_COMM_WORLD);
	MPI_Send(pptr,nphipts,MPI_DOUBLE,proc,4,MPI_COMM_WORLD);
	MPI_Send(pwtptr,nphipts,MPI_DOUBLE,proc,5,MPI_COMM_WORLD);	
      	MPI_Send(rtptr,ntrays,MPI_DOUBLE,proc,4,MPI_COMM_WORLD);
	MPI_Send(rtwtptr,ntrays,MPI_DOUBLE,proc,5,MPI_COMM_WORLD);	
	MPI_Send(rpptr,nprays,MPI_DOUBLE,proc,4,MPI_COMM_WORLD);
	MPI_Send(rpwtptr,nprays,MPI_DOUBLE,proc,5,MPI_COMM_WORLD);	
      }
    }
  }
  else { /*rank != master*/
    MPI_Status status;
    //receive from the swap variables
    MPI_Recv(rptr,nrpts,MPI_DOUBLE,master,0,MPI_COMM_WORLD,&status);
    MPI_Recv(rwtptr,nrpts,MPI_DOUBLE,master,1,MPI_COMM_WORLD,&status);
    MPI_Recv(tptr,nthetapts,MPI_DOUBLE,master,2,MPI_COMM_WORLD,&status);
    MPI_Recv(twtptr,nthetapts,MPI_DOUBLE,master,3,MPI_COMM_WORLD,&status);
    MPI_Recv(pptr,nphipts,MPI_DOUBLE,master,4,MPI_COMM_WORLD,&status);
    MPI_Recv(pwtptr,nphipts,MPI_DOUBLE,master,5,MPI_COMM_WORLD,&status);
    MPI_Recv(rtptr,ntrays,MPI_DOUBLE,master,4,MPI_COMM_WORLD,&status);
    MPI_Recv(rtwtptr,ntrays,MPI_DOUBLE,master,5,MPI_COMM_WORLD,&status);
    MPI_Recv(rpptr,nprays,MPI_DOUBLE,master,4,MPI_COMM_WORLD,&status);
    MPI_Recv(rpwtptr,nprays,MPI_DOUBLE,master,5,MPI_COMM_WORLD,&status);

    // assign to the vector objects
    for (int n = 0; n < nrpts; n++) { rpts[n] = rptr[n]; rwts[n] = rwtptr[n]; }
    for (int n = 0; n < nthetapts; n++) { thetapts[n] = tptr[n]; thetawts[n] = twtptr[n]; }
    for (int n = 0; n < nphipts; n++) { phipts[n] = pptr[n]; phiwts[n] = pwtptr[n]; }
    for (int n = 0; n < ntrays; n++) { raytpts[n] = rtptr[n]; raytwts[n] = rtwtptr[n]; }
    for (int n = 0; n < nprays; n++) { rayppts[n] = rpptr[n]; raypwts[n] = rpwtptr[n]; }
  }

  Linear_interp rterp(rpts,rpts);
  Linear_interp tterp(thetapts,thetapts);
  Linear_interp pterp(phipts,phipts);

  // //check
  // std::cout << "Proc rank = " << rank << std::endl;
  // std::cout << "These are the altitudes in use for this run:" << std::endl;
  // std::cout.width(20);
  // std::cout<<"alt(km)";
  // std::cout.width(20);
  // std::cout<<"nH cm^-3"<<std::endl;
  // for (int i = 0; i < nrpts; i++) {
  //   std::cout.width(20);
  //   std::cout << (rpts[i]-rMars)/1e5;
  //   std::cout.width(20);
  //   std::cout << nHinterp(rpts[i])<<std::endl;
  // }




  // what's the smallest altitude in the grid?
  double rmin = rpts[nrpts-1];
  //rmax is defined in definitions.h
  

  //wait 1 sec
  timer = clock(); // start timing
  while (clock() - timer < secwait)
    ;
  MPI_Barrier(MPI_COMM_WORLD);


  // variables for matrix manipulation
  // now we need the number of regions not region boundaries
  const int nrows = (nrpts-1)*(nthetapts-1)*(nphipts-1);
  const int ncols = nrows;
  if (rank == master)  std::cout << "nrows = " << nrows << ", ncols = " << ncols << std::endl;
  //  double nHvals[nrows][ncols];
  //  double nCO2vals[nrows][ncols];  

  // //Allocate the kernel array:
  // double* kernptr;
  // kernptr = (double *) malloc(nrows * ncols * sizeof(double));
  // if (kernptr == NULL) throw("Failure to allocate room for the array");
  // //Get pointers to the rows:
  // double** kvals;
  // kvals = (double **) malloc(nrows * sizeof(double*));
  // if (kvals == NULL) throw("Failure to allocate room for pointers");
  // //Point the pointers:
  // for (int k = 0; k < nrows; k++)
  //   kvals[k] = kernptr + (k * ncols);

  double** kvals;
  kvals = new double*[nrows];
  for (int i = 0; i < nrows; i++) {
    kvals[i] = new double[ncols];
    //    std::cout << "kvals[" << i << "] = " << kvals[i] << std::endl;
  }
  //  std::cin.get();
  //  double* oldptr = 0;
  double yvec[nrows]; // single scattering values

  // std::cout << "kvals = " << kvals << std::endl;
  // for(int r = 0; r<nrows; r++) {
  //   std::cout << "kvals + " << r*ncols << " = " << kvals + r << std::endl;
  //   std::cout << "kvals[" << r << "] = " << kvals[r] << std::endl;
  //   std::cin.get();
  // }
  

  //create file objects for output
  //ofstream Hfile, CO2file;
  ofstream trcheckfile, kfile, yfile, invsolfile, itersolfile;
  if (rank == master) {
    //    getfname(Hfile,"nH",nrpts,nthetapts,nphipts);
    //    getfname(CO2file,"nCO2",nrpts,nthetapts,nphipts);
    getfname(trcheckfile,"trcheck",nrpts,nthetapts,nphipts);
    getfname(kfile,"kern",nrpts,nthetapts,nphipts);
    getfname(yfile,"y",nrpts,nthetapts,nphipts);
    getfname(invsolfile,"inv_sol",nrpts,nthetapts,nphipts);
    getfname(itersolfile,"iter_sol",nrpts,nthetapts,nphipts);
  }

  //chop up the work so that each processor has its own bit:
  int sec = nrows/size;
  int extra = nrows%size;
  int low = 0, high;
  for (int proc = 0; proc < rank; proc++)
    low += proc < (size-extra) ? sec : sec+1;
  high = low + (rank < (size-extra) ? sec : sec+1);

  // echo the parts being done by each processor to terminal
  std::cout << "Process " << rank << " of " << size-1
	    << ": low = " << low << ", high = " << high << std::endl;


  //wait 1 sec
  timer = clock(); // start timing
  while (clock() - timer < secwait)
    ;
  MPI_Barrier(MPI_COMM_WORLD);

  
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



  
  int row,col,irw,itw,ipw,irs,its,ips;//loop variables
  double dr, dt, dp;//size of the current grid box
  int iter; // number of iterations in radial expansion
  int ordx, otdx, opdx, rdx, tdx, pdx; // volume indices
  
  //delay threads
  timer = clock(); // start timing
  while (clock() - timer < rank*secwait/50)
    ;

  for (row = low; row < high; row++) {
    //initialize the kernel matrix:
    for(col = 0; col < ncols; col++)
      kvals[row][col] = 0.0;
    //    std::cout << "kvals[" << row << "] = " << kvals[row] << std::endl;
    //get the r and theta point corresponding to this row
    irw = row/((nthetapts-1)*(nphipts-1));
    itw = (row-irw*(nthetapts-1)*(nphipts-1))/(nphipts-1);
    ipw = row%(nphipts-1);

    // keep count of the computation
    time_t rawtime;
    struct tm * timeinfo;
    time (&rawtime);
    timeinfo = localtime(&rawtime);
    std::cout << irw << ": " << itw << ": " << ipw << " starting at " << asctime(timeinfo);
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
    omega = 0.0;

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
	    omega += domega;

	    
	    //get the unit vector components along the angle outward
	    //(see rotation unit vectors.nb for derivation):
	    line_x = cos(p1)*cos(pray)*cos(t1)*sin(tray)
	              +cos(p1)*cos(tray)*sin(t1)
	              -sin(tray)*sin(pray)*sin(p1);
	    line_y = cos(tray)*sin(p1)*sin(t1)
	              +sin(tray)*cos(pray)*sin(p1)*cos(t1)
	              +sin(tray)*sin(pray)*cos(p1);
	    line_z = cos(tray)*cos(t1)-cos(pray)*sin(tray)*sin(t1);

	    // std::cout << "line = { " << line_x << " , " << line_y << " , " << line_z << "}\n";
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
	    // std::cout << "rrpt = " << rrpt << ", rmin = " << rmin << ", rmax = " << rmax << std::endl;

	    while ((rrpt >= rmin && rrpt <= rmax && iter < maxit) || iter == 0) {
	      //while still inside the atmosphere
	      iter++;
	      // if too many iterations, report error, but proceed
	      // with calculation
	      if (iter == maxit) {
		std::cout << "On row " << row
			  << " iter > maxit!" << std::endl;
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
	      
	      // std::cout<<"Here, \n"
	      // 	       <<" r=[" << rpts[rdx] << ", " << rpts[rdx+1] << "]" << std::endl;
	      // std::cout<<" t=[" << thetapts[tdx] << ", " << thetapts[rdx+1] << "]" <<std::endl;
	      // std::cout<<" p=[" << phipts[pdx] << ", " << phipts[pdx+1] << "]" << std::endl;

	      //get the step size for this box
	      dr=(rpts[rdx+1]-rpts[rdx]);
	      dt=rpts[rdx]*(thetapts[tdx+1]-thetapts[tdx]);
	      dp=rpts[rdx]*sin((thetapts[tdx+1]+thetapts[tdx])/2.)*(phipts[pdx+1]-phipts[pdx]);

	      ds = dsfrac*abs(dr);
	      ds = ds < dsfrac*abs(dt) ? ds : dsfrac*abs(dt);
	      ds = ds < dsfrac*abs(dp) ? ds : dsfrac*abs(dp);
	      
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
		if (rrpt>rmax||rrpt<rmin)
		  break;
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

    std::cout << irw << ": " << itw << ": " << ipw << " omega/4*pi = " << omega << std::endl;

    // This has been moved to the master process for the purpose of checking transposition
    // //multiply by the factors which are constant across the row
    // //(outside the integral):
    // for (col = 0; col < ncols; col++) {
    //   kvals[row][col] *= dtau_tot(rr1,nHinterp,Texo)/(4*pi);
    // }


    
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
      zpt = rmax;
      //      std::cout << "rpt = {" << xpt << ", " << ypt << ", " << zpt << "}\n";

      // store the values in the vector object r2
      r1[0] = x1; r1[1] = y1; r1[2] = z1;
      r2[0] = xpt; r2[1] = ypt; r2[2] = zpt;

      // get the column densities
      std::cout << "Getting CO2 optical depth...\n";
      dtau_CO2_int dtauCO2(&thisatmointerp);
      tauCO2col = qlinetrap_2pt(dtauCO2, r1, r2, 1e-3);
      std::cout << "  ... CO2 optical depth = "<< tauCO2col << ".\n";

      std::cout << "Getting H optical depth...\n";
      dtau_H_int dtauH(&thisatmointerp);
      tauHcol = qlinetrap_2pt(dtauH, r1, r2, 1e-3);
      std::cout << "  ... H optical depth = "<< tauHcol << ".\n";

      // put them together into the y-vec
      //std::cout << "irw = " << irw << std::endl;
      yvec[row] = HolT_lookup.interp(tauHcol)*std::exp(-tauCO2col);
      //    yvec[row] = 0.0;
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  //wait 1 sec
  timer = clock(); // start timing
  while (clock() - timer < secwait)
    ;

  // send computed values to master process
  if (rank != master) {
    for (row = low; row < high; row++) {
      MPI_Send(kvals[row],ncols,MPI_DOUBLE,master,row,MPI_COMM_WORLD);
      MPI_Send(&yvec[row], 1,MPI_DOUBLE,master,row,MPI_COMM_WORLD);
    }
  }
  else { /*rank == master*/
    MPI_Status status;
    for (int proc = 0; proc < size; proc++) {
      if (proc!=master) {
   	int proclow = 0, prochigh;
   	for (int procd = 0; procd < proc; procd++)
   	  proclow += procd < (size-extra) ? sec : sec+1;
   	prochigh = proclow + (proc < (size-extra) ? sec : sec+1);

	int count = (prochigh-proclow);
   	std::cout << "Receiving " << count << " rows from process " << proc << std::endl;

	//	std::cout << "proclow = " << proclow << ", prochigh = " << prochigh << std::endl;
	for (row = proclow; row < prochigh; row++) {
	  MPI_Recv(kvals[row],ncols,MPI_DOUBLE,proc,row,MPI_COMM_WORLD,&status);
	  MPI_Recv(&yvec[row],1,MPI_DOUBLE,proc,row,MPI_COMM_WORLD,&status);
	}
      }
    }
  }
  
  
  MPI_Barrier(MPI_COMM_WORLD);
  //wait 1 sec
  timer = clock(); // start timing
  while (clock() - timer < secwait)
    ;


  // the solution only needs to be found by the master node:
  if (rank == master) {
  
    // //what we've got right now is a kernel matrix which should equal its own transpose.
    // //let's see if it does:
    // double k1, k2, tdiff = 0.0;
    // for (row = 0; row < nrows; row++) {
    //   for (col = row; col < ncols; col++) {

    // 	k1 = kvals[row][col];
    // 	k2 = kvals[col][row];
    // 	if (k1 != 0.0 || k2 != 0.0) {
    // 	  k1 = k1 > k2 ? 2.0*(k1-k2)/(k1+k2) : 2.0*(k2-k1)/(k1+k2);
    // 	  tdiff += k1*k1;
    // 	}
    //   }
    // }
    // std::cout << "norm difference between kernel and its transpose = " << sqrt(tdiff) << std::endl;

    //write this transpose check matrix to file:

    // std::cout << "Writing trcheck to file." << std::endl;
    // for (row = 0; row < nrows; row++) {
    //   for (col = 0; col < ncols; col++) {
    // 	//	Hfile.width(20);
    // 	//	Hfile << nHvals[row][col];
    // 	//	CO2file.width(20);
    // 	//	CO2file << nCO2vals[row][col];
    // 	trcheckfile.width(20);
    // 	trcheckfile << kvals[row][col];
    //   }
    //   //      Hfile << std::endl;
    //   //      CO2file << std::endl;
    //   trcheckfile << std::endl;
    // }
  

    //Now multiply by the rowwise factors that destroy the transposition symmetry:

    // for (row = 0; row < nrows; row++) {
    //   irw = row/nthetapts;
    //   rr1 = rpts[irw];
    //   //multiply by the factors which are constant across the row
    //   //(outside the integral):
    //   for (col = 0; col < ncols; col++) {
    // 	kvals[row][col] *= dtau_H(rr1,thisatmointerp)/(4*pi);
    //   }
    // }
  

  
    //    write kernel and yvec to file;
    std::cout << "Writing kernel and single scattering vector to file." << std::endl;
    for (row = 0; row < nrows; row++) {
      for (col = 0; col < ncols; col++) {
    	//	Hfile.width(20);
    	//	Hfile << nHvals[row][col];
    	//	CO2file.width(20);
    	//	CO2file << nCO2vals[row][col];
    	kfile.width(20);
    	kfile << kvals[row][col];
      }
      //      Hfile << std::endl;
      //      CO2file << std::endl;
      kfile << std::endl;
      yfile.width(20);
      yfile << yvec[row] << std::endl;
    }
  
    // solve using LU decomposition
    std::cout << "Solving for source function: " << std::endl << "By matrix inversion: " << std::endl;
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

    //compute solution by iterative multiplication
    std::cout << "Solving by iteration:" << std::endl;
    double sol_it[nrows];
    double smax = 0.0, osmax = 0.0;
    //initialize
    for (row = 0; row < nrows; row++) {
      sol_it[row] = yvec[row];
      osmax = sol_it[row] > osmax ? sol_it[row] : osmax;
      //      std::cout << "osmax = " << osmax << std::endl;
    }
      

    // keep multiplying until convergence
    diff = 1.0;
    iter = 1;
    while (diff > 1e-5 && iter < maxit) {
      diff = 0.0;
      iter++;
      smax = 0.0;
      if (iter == maxit) std::cout << "iter > maxit in iterative solution!" << std::endl;
      for (row = 0; row < nrows; row++) {
	comp = 0.0;
	for (col = 0; col < ncols; col++)
	  comp += kvals[row][col]*sol_it[col];
	comp += yvec[row];
	smax = comp > smax ? comp : smax;
	// how much have we moved since last time?
	temp = comp - sol_it[row];
	temp = temp > 0 ? temp : -temp;
	temp = temp/osmax;
	sol_it[row] = comp;
	diff = temp > diff ? temp : diff;
      }
      osmax = smax;
      if (iter%10 == 0) //print diff and iter each N iterations
	std::cout << "iter = " << iter << ", diff = " << diff << ", osmax = " << osmax << std::endl;
    }
    // print iter and diff at end of computation
    std::cout << "iter = " << iter << ", diff = " << diff << ", osmax = " << osmax << std::endl;
    
    //compare two solutions:
    diff = 0.0;
    for (row = 0; row < nrows; row++) {
      temp = sol_it[row]/sol[row]-1.0;
      diff += temp*temp;
    }
    std::cout << "|S_it/S_inv - 1| = " << sqrt(diff) << std::endl;
    
    
    //write solution to file:
    std::cout << "Writing solution to file." << std::endl;
    for (row = 0; row < nrows; row++) {
      invsolfile.width(20);
      invsolfile << sol[row] << std::endl;
      itersolfile.width(20);
      itersolfile << sol_it[row] << std::endl;
    }
  }


  //Cleanup:

  
  //    Hfile.close();
  //    CO2file.close();
  // trcheckfile.close();
  kfile.close();
  yfile.close();
  invsolfile.close();
  itersolfile.close();    

  // //free allocated kernel matrix
  // free (kvals);
  for (int i = 0; i < nrows; i++) {
    //    std::cout << "kvals[" << i << "] = " << kvals[i] << std::endl;
    delete [] kvals[i];
  }
  delete [] kvals;

  // print time elapsed
  if (rank == master) {
    clock_t finish = clock();
    double secs = (finish-start)*1.0/CLOCKS_PER_SEC;
    int hrs = secs/3600;
    int mins = (secs - hrs*3600)/60;
    secs = secs - hrs*3600 - mins*60;
    std::cout << "Elapsed time is " << hrs << " hours, " << mins << " minutes, and " << secs << " seconds.\n";
  }

  MPI::Finalize();

  return 0;
}
