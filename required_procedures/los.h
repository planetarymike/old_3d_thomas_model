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
#include "ray_intersect.h"
#include <Eigen/Dense>

using std::exp;
using std::log;

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

struct atmo_ray;

//structure to hold the atmosphere grid
struct atmo_grid {
  //define vectors for the quadrature regions and weights
  VecDoub rpts, rwts;
  VecDoub thetapts, thetawts;
  VecDoub costhetapts;
  VecDoub phipts, phiwts;

  //and the ray directions and weights
  VecDoub raytpts, raytwts;
  VecDoub rayppts, raypwts;

  // source function is computed for regions, not region boundaries, hence the -1
  int nvoxels;

  atmo_grid(int nrpts, int ntpts, int nppts,
	    int ntrays, int nprays,
	    atmointerp atm,
	    bool printout=true) {
    //set up the vectors
    rpts.resize(nrpts);
    rwts.resize(nrpts);
    thetapts.resize(nthetapts);
    costhetapts.resize(nthetapts);
    thetawts.resize(nthetapts);
    phipts.resize(nphipts);
    phiwts.resize(nphipts);

    raytpts.resize(ntrays);
    raytwts.resize(ntrays);
    rayppts.resize(nprays);
    raypwts.resize(nprays);
    
    nvoxels = (nrpts-1)*(nthetapts-1)*(nphipts-1);

    //now obtain the quadrature points in each dimension:
    getquadpts(rpts,rwts,
	       thetapts,thetawts,
	       phipts,phiwts,
	       printout,
	       atm,
	       rmethod,
	       tmethod);

    for (int i=0;i<nthetapts;i++)
      costhetapts[i]=cos(thetapts[i]);

    getraypts(raytpts,raytwts,
	      rayppts,raypwts,
	      printout,
	      raymethod);
  }

  template <typename C>
  void voxel_function(C &obj, double (C::*function)(const double &,const double &,const double &), double* v) {
    double r,t,p;
    int row, ir, it, ip;
    for (ir=0;ir<nrpts-1;ir++) {
      for (it=0;it<nthetapts-1;it++) {
	for (ip=0;ip<nphipts-1;ip++) {
	  row=ir*(nthetapts-1)*(nphipts-1)+it*(nphipts-1)+ip;
	  get_point_coords(ir, it, ip, r, t, p);
	  
	  v[row]=(obj.*function)(r,t,p);
	}
      }
    }
  }
  
  void get_point_coords(   int ir,    int it,    int ip,
			   double &r, double &t, double &p) {
    r = (rpts[ir]+rpts[ir+1])/2;
    t = (thetapts[it]+thetapts[it+1])/2;
    p = (phipts[ip]+phipts[ip+1])/2;
  }
  
  void get_ray_coords(int itray,   int ipray,
		      double &t,   double &p) {
    t = raytpts[itray];
    p = rayppts[ipray];
  }

  template <typename C>
  void voxel_traverse(atmo_ray &v, C &obj, void (C::*function)(const int /*row*/,
							       const int /*ir*/,
							       const int /*it*/,
							       const int /*ip*/,
							       const double /*pathlength*/,
							       const double /*tau_species_initial*/,
							       double &/*tau_species_final*/,
							       const double /*tau_absorber_initial*/,
							       double &/*tau_absorber_final*/,
							       const double /*domega*/));

};

struct atmo_point {
  //atmosphere point
  double x, y, z;
  double r, t, p;

  int ir, it, ip; //if a grid point, these are populated with the indices of the point
  int row;
  
  atmo_point() {}

  atmo_point(double rr, double tt, double pp) : r(rr), t(tt), p(pp) {
    x = r*sin(t)*cos(p);
    y = r*sin(t)*sin(p);
    z = r*cos(t);
    row=ir=ip=it=-1;
  }

  atmo_point(atmo_grid &grid, int irr, int itt, int ipp) : ir(irr), it(itt), ip(ipp) {
    grid.get_point_coords(ir, it, ip,
			   r,  t,  p);

    x = r*sin(t)*cos(p);
    y = r*sin(t)*sin(p);
    z = r*cos(t);

    row=-1;
  }

  atmo_point(atmo_grid &grid, int roww) : row(roww) {
    //get the r, theta, phi point corresponding to this row
    ir = row/((nthetapts-1)*(nphipts-1));
    it = (row-ir*(nthetapts-1)*(nphipts-1))/(nphipts-1);
    ip = row%(nphipts-1);
    
    grid.get_point_coords(ir, it, ip,
			   r,  t,  p);
    
    x = r*sin(t)*cos(p);
    y = r*sin(t)*sin(p);
    z = r*cos(t);
  }
  
};

struct atmo_ray {
  double line_x, line_y, line_z;//cartesian vector elements
  double tray, pray; // theta and phi from a given point
  double costray, sintray;//needed for sphere intersections later on, might as well get them now
  int itray, ipray; //indices of ray point if on grid
  atmo_point pt;
  double domega;
  
  atmo_ray() {}

  atmo_ray(atmo_point ptt, double trayy, double prayy) : pt(ptt), tray(trayy), pray(prayy) {
    costray=cos(tray);
    sintray=sin(tray);
    
    line_x = cos(pt.p)*cos(pray)*cos(pt.t)*sintray
              +cos(pt.p)*costray*sin(pt.t)
              -sintray*sin(pray)*sin(pt.p);
    line_y = costray*sin(pt.p)*sin(pt.t)
              +sintray*cos(pray)*sin(pt.p)*cos(pt.t)
              +sintray*sin(pray)*cos(pt.p);
    line_z = costray*cos(pt.t)-cos(pray)*sintray*sin(pt.t);
    domega=1.0;
  }

  atmo_ray(atmo_point ptt, double line_xx, double line_yy, double line_zz) :
    pt(ptt), line_x(line_xx), line_y(line_yy), line_z(line_zz) {

    costray=(line_x*pt.x + line_y*pt.y + line_z*pt.z)/pt.r;
    sintray=sqrt(1-costray*costray);

    tray=-1;//we don't actually need to compute these for the ray trace
    pray=-1;
    itray=-1;
    ipray=-1;
    domega=1.0;
  }

  atmo_ray(atmo_point ptt, atmo_grid &grid, int itt, int ipp) : pt(ptt), itray(itt), ipray(ipp) {
    tray = grid.raytpts[itray];
    pray = grid.rayppts[ipray];

    costray=cos(tray);
    sintray=sin(tray);
    
    line_x = cos(pt.p)*cos(pray)*cos(pt.t)*sintray
              +cos(pt.p)*costray*sin(pt.t)
              -sintray*sin(pray)*sin(pt.p);
    line_y = costray*sin(pt.p)*sin(pt.t)
              +sintray*cos(pray)*sin(pt.p)*cos(pt.t)
              +sintray*sin(pray)*cos(pt.p);
    line_z = costray*cos(pt.t)-cos(pray)*sintray*sin(pt.t);

    domega = sintray*grid.raytwts[itray]*grid.raypwts[ipray]/(4*pi);
  }
					       

};

template <typename C>
void atmo_grid::voxel_traverse(atmo_ray &v,
			       C &obj,
			       void (C::*function)(const int /*row*/,
						   const int /*ir*/,
						   const int /*it*/,
						   const int /*ip*/,
						   const double /*pathlength*/,
						   const double /*tau_species_initial*/,
						   double &/*tau_species_final*/,
						   const double /*tau_absorber_initial*/,
						   double &/*tau_absorber_final*/,
						   const double /*domega*/)) {
  
  //define objects for ray-voxel intersection
  vector<boundary_intersection> rsec(2*nrpts);//two intersections are possible per sphere
  vector<boundary_intersection> tsec(2*nthetapts);//two intersections are possible per cone
  vector<boundary_intersection> psec(nphipts-1);//only one intersection is possible with a plane
  //-1 accounts for the fact that the
  //   half-plane at 0 and 2pi are the
  //   same
    
  //get ray-voxel boundary intersections
  //these are arrays of the path length s and the index of the boundary intersected
  sphere_intersections(v.pt.r,
		       v.costray,v.sintray,
		       rpts,rsec);
  cone_intersections(v.pt.x,v.pt.y,v.pt.z,
		     v.line_x,v.line_y,v.line_z,
		     v.pt.r,v.costray,
		     costhetapts,tsec);
  halfplane_intersections(v.pt.x,v.pt.y,
			  v.line_x,v.line_y,
			  phipts,psec);
  //now (r,t,p)sec contains the ray/voxel boundary intersection information
    
    
    
  //these include all intersections, even negative values,
  //so we must advance the table variables to the initial point
  int rbi, tbi, pbi;
  bool notsec, nopsec;
  rbi=0;
  while (rbi<2*nrpts&&rsec[rbi].distance<0)
    rbi++;
  if (rbi==2*nrpts)//we get to 2*nrpts if 2*nrpts-1, the last r pt, has no intersections
    toss("no radial intersections, something is borked.")
	     tbi=0;
  while (tbi<2*nthetapts&&tsec[tbi].distance<0)
    tbi++;
  notsec = tbi==2*nthetapts ? true : false;//no intersections with theta coordinate boundaries
  pbi=0;
  while (pbi<nphipts-1&&psec[pbi].distance<0)
    pbi++;
  nopsec = pbi==nphipts-1 ? true : false;//no intersections with phi coordinate boundaries

  // initialize the variables
  double si, sf;
  si = sf = 0.0;// s = 0 at start of integration
  double tau_species_i, tau_species_f;
  tau_species_i = tau_species_f = 0.0;
  double tau_absorber_i, tau_absorber_f;
  tau_absorber_i = tau_absorber_f = 0.0;

  double rpt, tpt, ppt;
  rpt = v.pt.r;// first point is grid point itself
  tpt = v.pt.t;
  ppt = v.pt.p;


  //step outward along the ray, adding influence in each voxel
  //this is a while loop that moves through the boundary intersections until the ray terminates
  bool inside = true;// by definition we begin inside the atmosphere
    
  int rdx = v.pt.ir;// current volume indicies along the ray
  int tdx = v.pt.it;
  int pdx = v.pt.ip;
  int next_rdx = rdx;// next volume indicies along the ray 
  int next_tdx = tdx; 
  int next_pdx = pdx;
    
  int nextb;//tracker for whether next crossing is r, t, or p

  int iter = 0;//we should limit the number of voxel crossings to at
  //most the number of voxels (a generous upper limit)
  while (inside&&iter<nvoxels) {
    iter++;

    //update path lengths
    si=sf;
    tau_species_i=tau_species_f;
    tau_absorber_i=tau_absorber_f;
    //update the index variables to the next gridpoint
    rdx=next_rdx;tdx=next_tdx;pdx=next_pdx;
    //get the new point coordinates
    get_point_coords(rdx,tdx,pdx,
		     rpt,tpt,ppt);


    //look at the entries in each coordinate table to
    //determine which boundary we cross next
    //nextb = 1,2,3 depending on if next boundary is r, theta, phi
    if (notsec) {
      //no intersections with t
      if (nopsec) {
	//no intersections with t or p, r is always next
	nextb=1;
      } else {
	//no intersections with t, but intersections with p are possible
	if (rsec[rbi] <= psec[pbi]) {
	  //          <= here allows rays to pass through
	  //             the voxel edges/corners without
	  //             breaking the code.
	  nextb=1;
	} else {
	  nextb=3;
	}
      }
    } else {
      if (nopsec) {
	//intersections with t, but no intersections with p
	if (rsec[rbi] <= tsec[tbi]) {
	  nextb=1;
	} else {
	  nextb=2;
	}
      } else {
	//intersections with r, t, and p are possible
	if (rsec[rbi] <= tsec[tbi] &&
	    rsec[rbi] <= psec[pbi]) {
	  //distance to next r boundary is smallest
	  nextb=1;
	} else if (tsec[tbi] <= psec[pbi]) {
	  //distance to next t boundary is smallest
	  nextb=2;
	} else {
	  nextb=3;
	}
      }
    }
	    
    //now we step over the next boundary
    //how we take the step depends on what coordinate boundary we're crossing
    switch (nextb) {
    case 1:
      //crossing r boundary
      sf = rsec[rbi].distance;//distance to next boundary
      next_rdx = rdx < rsec[rbi].boundary ? rdx+1 : rdx-1;
      //rdx matches the smaller boundary index
		
      if (rsec[rbi].boundary==0||rsec[rbi].boundary==nrpts-1)
	inside=false;//whenever we cross rmin or rmax the integration is over
		
      rbi++;//we have crossed this r boundary
      break;
    case 2:
      //crossing t boundary

      //horrible things happen at the singularity along the
      //symmetry axis, so we set up the ray and voxel grid
      //to be offset; that way no intersections occur where
      //theta=0,pi

      sf=tsec[tbi].distance;//distance to next boundary
      next_tdx = tdx < tsec[tbi].boundary ? tdx+1 : tdx-1;
      //tdx matches the smaller boundary index
		
      tbi++;//we have crossed this t boundary
      if (tbi == 2*nthetapts)
	notsec=true;
      break;
    case 3:
      //crossing p boundary 
      sf=psec[pbi].distance;//distance to next boundary
      next_pdx = pdx < psec[pbi].boundary ? pdx+1 : pdx-1;
      //pdx matches the smaller boundary index,
      //unless we cross the boundary at 0, 2pi:
      if (pdx==0 && psec[pbi].boundary==nphipts-2)
	next_pdx=nphipts-2;
      if (pdx==nphipts-2 && psec[pbi].boundary==0)
	next_pdx=0;
		
      pbi++;//we have crossed this p boundary
      if (pbi == nphipts-1)
	nopsec=true;
      break;
    }

    //call the desired function in each voxel
    (obj.*function)(v.pt.row,
		    rdx, tdx, pdx,
		    sf-si,
		    tau_species_i,
		    tau_species_f,//output
		    tau_absorber_i,
		    tau_absorber_f,//output
		    v.domega);
    
  }
  
}

struct Holinterp_abs
{
  /* Structure to load and compute
     integrated Holstein functions 
     with absorption. */
  string Holfilename;
  bool silent;
  int ntaupts, nabspts;
  VecDoub tau_vec, abs_vec;
  MatDoub Hol_vals;
  Linear_interp tau_terp, abs_terp;

  Holinterp_abs(string Holfilenamee, bool silent=FALSE) 
    : Holfilename(Holfilenamee) {
    if (!silent)
      std::cout << "Reading interpolated Hol values from " << Holfilename << std::endl;
    //load the data from the file
    ifstream Holfile;
    Holfile.open(Holfilename.c_str());
    
    char dumline[100];//dummy char string to read headers
    
    //read coordinate vectors
    Holfile.getline(dumline,100);
    Holfile.getline(dumline,100);
    Holfile >> ntaupts;
    tau_vec.resize(ntaupts);
    for (int i = 0; i < ntaupts; i++)
      Holfile >> tau_vec[i];
    tau_terp=Linear_interp(tau_vec,tau_vec);
    

    Holfile.getline(dumline,100);
    Holfile.getline(dumline,100);
    Holfile >> nabspts;
    abs_vec.resize(nabspts);
    for (int i = 0; i < nabspts; i++)
      Holfile >> abs_vec[i];
    abs_terp=Linear_interp(abs_vec,abs_vec);
    
    //now read in the data:
    Hol_vals.resize(ntaupts,nabspts);
    Holfile.getline(dumline,100);
    Holfile.getline(dumline,100);
    for (int itau = 0; itau < ntaupts; itau++) {
      for (int iabs = 0; iabs < nabspts; iabs++) {
	Holfile >> Hol_vals[itau][iabs];
      }
    }
    
    Holfile.close();
  }
  
  double operator()(const double tau, const double abs) {
    double logtau=log(tau);
    double logabs=log(abs);
    
    if (logabs < abs_vec[0])
      logabs = abs_vec[0];
    if (logtau < tau_vec[0]) 
      logtau = tau_vec[0];
    if (logabs > abs_vec[nabspts-1])
      toss("abs > max absorption in Holinterp");
    if (logtau > tau_vec[ntaupts-1])
      toss("tau > max optical depth in Holinterp");
    
    int i, j;
    double yy, t, u;
    //find the grid square:
    i = tau_terp.cor ? tau_terp.hunt(logtau) : tau_terp.locate(logtau);
    j = abs_terp.cor ? abs_terp.hunt(logabs) : abs_terp.locate(logabs);
    
    //interpolate:
    t = (logtau-tau_terp.xx[i])/(tau_terp.xx[i+1]-tau_terp.xx[i]);
    u = (logabs-abs_terp.xx[j])/(abs_terp.xx[j+1]-abs_terp.xx[j]);
    yy = (1.-t)*(1.-u)*Hol_vals[i  ][j  ]
        +    t *(1.-u)*Hol_vals[i+1][j  ]
        +(1.-t)*    u *Hol_vals[i  ][j+1]
        +    t *    u *Hol_vals[i+1][j+1];
    
    return exp(yy);
  }
};


struct Holinterp
{
  /* Structure to load tabular Holstein functions */
  string Holfilename;
  bool silent;
  int npts;
  VecDoub tau_vec, Hol_vec;
  Linear_interp Hol_interp;

  Holinterp(string Holfilenamee, bool silent=FALSE) 
  : Holfilename(Holfilenamee) {
    if (!silent)
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
    double logtau=log(tau);
    if (logtau > tau_vec[npts-1]) {
      std::cout << "tau = " << tau << std::endl;
      std::cout << "logtau = " << logtau << std::endl;
      std::cout << "tau_vec[npts-1] = " << tau_vec[npts-1] << std::endl;
      toss("tau > max optical depth in Holinterp");
    }

    double retval;
    if (logtau<tau_vec[0]) {
      retval=Hol_vec[0];
    } else if (logtau>tau_vec[npts-1]) {
      retval=Hol_vec[npts-1];
    } else {
      retval=Hol_interp.interp(logtau);
    }
    
    return exp(retval);
  }

  double operator()(const double tau) {
    return interp(tau);
  }
};


//an emission structure to hold emission relevant data, compute source
//functions, etc
struct emission {
  atmo_grid grid;
  int nvoxels;
  
  double* species_n;//densities of scatterers and absorbers on the tabulated grid
  double* absorber_n; 
  double* species_sigma;//scatterer cross section on the tabulated grid
  double* absorber_sigma;//absorber cross section on the tabulated grid
  double* dtau_species;
  double* dtau_absorber;
  double* abs; //ratio of dtau_abs to dtau_spe
  double max_tau_species;

  
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXd;
  MatrixXd influence;
  //  MatrixXd kernel;
  VecDoub tau_species;
  VecDoub tau_absorber;
  VectorXd singlescat;
  VectorXd sourcefn;
  VecDoub sourcefn_out;

  //integrated Holstein G and T with absorption
  //Holinterp_abs HolGint;
  //Holinterp_abs HolTint;
  // Holinterp HolG;
  Holinterp HolT;
  Holinterp HolTint;
  
  bool species_init;//are the density and sigma values populated?
  bool influence_init, singlescat_init;//are the RT properties initialized?
  bool solved;//have we solved for the source function?



  
  //constructor
  emission(atmo_grid &gridd) : grid(gridd),
			       // HolG(HolGfilename),
			       HolT(HolTfilename),
			       HolTint(HolTint_noabs_filename)
			       // HolGint(HolGintfilename),
			       // HolTint(HolTintfilename)
  {
    nvoxels = grid.nvoxels;
    
    influence=MatrixXd::Zero(nvoxels,nvoxels);
    //    kernel=MatrixXd::Zero(nvoxels,nvoxels);

    species_n = new double[nvoxels];
    absorber_n = new double[nvoxels]; 
    species_sigma = new double[nvoxels];
    absorber_sigma = new double[nvoxels];
    dtau_species = new double[nvoxels];
    dtau_absorber = new double[nvoxels];
    abs = new double[nvoxels];
    tau_species.assign(nvoxels,0.0);
    tau_absorber.assign(nvoxels,0.0);
    singlescat=VectorXd(nvoxels);
    sourcefn=VectorXd(nvoxels);
    sourcefn_out.resize(nvoxels);
    
    species_init=false;
    influence_init=false;
    singlescat_init=false;
    solved=false;
  }

  ~emission() {
    delete [] species_n;
    delete [] absorber_n; 
    delete [] species_sigma;
    delete [] absorber_sigma;
    delete [] dtau_species;
    delete [] dtau_absorber;
    delete [] abs;
  }





  //function to populate density grid 
  template<typename C>
  void init(C &obj,
	    double(C::*species_density_function)(const double &,const double &,const double &),
	    double(C::*species_sigma_function)(const double &,const double &,const double &),
	    double(C::*absorber_density_function)(const double &,const double &,const double &),
	    double(C::*absorber_sigma_function)(const double &,const double &,const double &))
  {
    
    grid.voxel_function(obj, species_density_function,  species_n);
    grid.voxel_function(obj, species_sigma_function,    species_sigma);
    grid.voxel_function(obj, absorber_density_function, absorber_n);
    grid.voxel_function(obj, absorber_sigma_function,   absorber_sigma);

    for (int row = 0; row < nvoxels; row++) { 
      dtau_species[row] = species_sigma[row]*species_n[row];
      dtau_absorber[row] = absorber_sigma[row]*absorber_n[row];
      abs[row] = dtau_absorber[row]/dtau_species[row];
    }
    
    max_tau_species=0.0;
    
    species_init=true;
  }






  //functions to update influence matrix and single scattering
  void update_influence(const int row,
			const int ir, const int it, const int ip,
			const double pathlength,
			const double tau_species_initial,
			double &tau_species_final,
			const double tau_absorber_initial,
			double &tau_absorber_final,
			const double domega) {
    if (!species_init)
      toss("initialize species densities before computing influence!");
    //update the influence matrix at this point
    int col=ir*(nthetapts-1)*(nphipts-1)+it*(nphipts-1)+ip;

    tau_species_final = tau_species_initial + dtau_species[col]*pathlength;
    tau_absorber_final = tau_absorber_initial + dtau_absorber[col]*pathlength; 
    
    double coef = domega;
    coef *= HolT(tau_species_initial) - HolT(tau_species_final) 
      - abs[col]*( HolTint(tau_species_final) - HolTint(tau_species_initial) );
    /* coef *= (HolGint(tau_species_final,abs)-HolGint(tau_species_initial,abs)); */
    
    influence(row,col) += coef;
    
    double tscomp = tau_species_final*exp(-tau_absorber_final);
    max_tau_species = tscomp > max_tau_species ? tscomp : max_tau_species;
  }
  //after all the updates, the calling routine sets influence_init=true
  











  void update_singlescat(int row,
			 int ir, int it, int ip,
			 double pathlength) {
    if (!species_init)
      toss("initialize species densities before computing single scattering!");
    //update the single scattering vector at this point
    int col=ir*(nthetapts-1)*(nphipts-1)+it*(nphipts-1)+ip;
    tau_species[row]+=dtau_species[col]*pathlength;
    tau_absorber[row]+=dtau_absorber[col]*pathlength;

    double tscomp = tau_species[row]*exp(-tau_absorber[row]);
    max_tau_species = tscomp > max_tau_species ? tscomp : max_tau_species;
  }
  //after all the updates, the calling routine sets singlescat_init=true
  
  //overload so the calling routine can use the same signature for
  //update_singlescat and update_influence
  void update_singlescat(const int row,
			 const int ir, const int it, const int ip,
			 const double pathlength,
			 const double tau_species_initial,
			 double &tau_species_final,
			 const double tau_absorber_initial,
			 double &tau_absorber_final,
			 const double domega) {
    update_singlescat(row,ir,it,ip,pathlength);
  }
  
  
  
  
  
  
  
  //function to solve for the source function
  void solve() {
    if (!species_init || !influence_init || !singlescat_init)
      toss("initialize before solving!");
    //do the matrix stuff to solve the equation
    
    std::cout << "max_tau_species = " << max_tau_species << std::endl;

    for (int row=0;row<nvoxels;row++) {
      //compute single scattering from optical depths
      if (tau_species[row]==-1.0||tau_absorber[row]==-1.0)
	singlescat(row)=0.0; //point is behind the limb
      else 
	singlescat(row)=HolT(tau_species[row])*exp(-tau_absorber[row])/species_sigma[row];
      // singlescat(row)=HolT.interp(tau_species[row])*exp(-tau_absorber[row])/species_sigma[row];
      
      
      for (int col=0;col<nvoxels;col++) {
	//convert influence matrix to kernel matrix
	//kernel(row,col) = 0.5*(influence(row,col)+influence(col,row));
	influence(row,col) = row==col ? 1.0-influence(row,col) : -influence(row,col);
      }
    }
    
    sourcefn=influence.partialPivLu().solve(singlescat);//partialPivLu has multithreading support

    //check solution:
    double diff = (influence*sourcefn - singlescat).norm()/singlescat.norm();
    std::cout << "relative error is = " << diff << std::endl;
    
    //transform variables back to Chaufray formulation of S so that we
    //can use existing code to compute intensities
    for (int row = 0; row < nvoxels; row++) {
      sourcefn(row)*=species_n[row]*species_sigma[row];
      sourcefn_out[row]=sourcefn(row);
    } 
    solved=true;

  }
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
  Sobj(string fnamee, bool silent=FALSE) : fname(fnamee) {
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
      toss("Bad filename in Sobj")
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

  //constructor to create Sobj from existing data.
  Sobj(double nHH, double TT,
       atmo_grid &grid,
       VecDoub sol) : nH(nHH), T(TT) {

    rterp=Linear_interp(grid.rpts,grid.rpts);
    tterp=Linear_interp(grid.thetapts,grid.thetapts);
    pterp=Linear_interp(grid.phipts,grid.phipts);

    nrpts=rterp.n;
    ntpts=tterp.n;
    nppts=pterp.n;
    
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
      toss("Sobj not initialized!");
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
      toss("bad pt in S")
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




struct LOS_integrator
{
  /* structure to compute intensity along LOS using saved, tabulated
     values for HolT */
  string HolTfname;
  Holinterp HolTinterp;
  
  LOS_integrator(string HolTfnamee=HolTfilename, bool silent=FALSE) 
    : HolTfname(HolTfnamee), HolTinterp(HolTfname,silent) {
    if (!silent)
      std::cout << "Reading interpolated Hol values from " << HolTfname << std::endl;
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
    double rmaxpad=0.99*rmax;
    while (rrpt > rmaxpad) {
      std::cout << "Moving to start of atmosphere:\n";
      std::cout << "Before moving: rrpt = " << rrpt << " > rmax = " << rmax <<std::endl; 
      s = xpt*line_x + ypt*line_y + zpt*line_z;
      s = -s - sqrt(s*s + 0.99*rmaxpad*rmaxpad - rrpt*rrpt);// I know it looks negative
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
      /* std::cin.get(); */
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
      if (tauH>taumax)
	toss("maximum optical depth exceeded!\n")

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
    tauCO2out = rrpt<=rminatm ? 1e99 : tauCO2;
  
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
