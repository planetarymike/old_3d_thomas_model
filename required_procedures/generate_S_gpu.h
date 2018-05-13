//generate_S_gpu.h -- code to generate H source functions according to the
//                    formalism of Anderson & Hord (1976).

// Now using efficient ray voxel intersections based on solving
// quadratics instead of line integrals.

//tabulated versions of holstein G and T functions with absorption are
//used to speed execution



#ifndef __generate_S_gpu_h
#define __generate_S_gpu_h

#include <ctime> // timing
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
#include <vector>
#include "ray_intersect.h" //ray-voxel boundary intersections
#include <limits> //sizes of double
typedef std::numeric_limits< double > dbl;
using std::cos;
using std::sin;
using std::vector;

string Sfilename(double nH, double T) {
  string fname;
  char cfname[200];
  sprintf(cfname,"S_gpu_nH%G_T%G_%dx%dx%d.dat", nH, T, nrpts, nthetapts, nphipts);
  fname=srcfnsloc+cfname;
  return fname;
}

string auxfilename(string tag,double nH, double T) {
  string fname;
  char cfname[200];
  sprintf(cfname,"%s_nH%G_T%G_%dx%dx%d.dat", tag.c_str(), nH, T, nrpts, nthetapts, nphipts);
  fname=srcfnsloc+cfname;
  return fname;
}


void generate_S(const double nexo, 
		const double Texo, 
		Sobj &Sout,
		atmointerp &atmoout,
		bool printout = FALSE, 
		bool filesout = FALSE) {

  std::cout << "#####################################################" << std::endl;
  std::cout << "#####################################################" << std::endl;
  std::cout << "##########KERNEL GENERATION NOW IN PROGRESS##########" << std::endl;
  std::cout << "#####################################################" << std::endl;
  std::cout << "#####################################################" << std::endl;


  //get nexo and Texo
  std::cout << "nexo = " << nexo << std::endl;
  std::cout << "Texo = " << Texo << std::endl;

  //start timing
  clock_t start = clock(); // start whole task timer

  //define the physical atmosphere
  atmointerp atm(nexo, Texo, nphyspts,"",printout);
  //physical atm(nexo,Texo);  
  
  //get the quadrature points
  atmo_grid grid(nrpts,nthetapts,nphipts,ntrays,nprays,atm,printout);

  //set up the emission object
  emission lyman_alpha(grid);
  lyman_alpha.init(atm,
  		   &atmointerp::nH,&atmointerp::lya_sH,
  		   &atmointerp::nCO2,&atmointerp::lya_sCO2);
  /* lyman_alpha.init(atm, */
  /* 		   &physical::nH,&physical::lya_sH, */
  /* 		   &physical::nCO2,&physical::lya_sCO2); */

  atmo_point pt;
  atmo_ray ray;

#pragma omp parallel for firstprivate(pt,ray) shared(grid,atm,lyman_alpha,ntrays,nprays,pi,std::cout) default(none)
  for (int row = 0; row < grid.nvoxels; row++) {

    // keep count of the computation
    /* time_t rawtime; */
    /* struct tm * timeinfo; */
    /* time (&rawtime); */
    /* timeinfo = localtime(&rawtime); */
    
    //get the voxel center coordinates
    pt = atmo_point(grid, row);
    
    double omega = 0.0;// check parameter to make sure sum(domega) = 4*pi
    
    //now integrate outward along the ray grid:
    for (int its = 0; its < ntrays; its++) {
      for (int ips = 0; ips < nprays; ips++) {
	ray = atmo_ray(pt, grid, its, ips);
	omega += ray.domega;
	
	grid.voxel_traverse(ray, lyman_alpha, &emission::update_influence);
      }
    }
    
    // std::cout << irw << ": " << itw << ": " << ipw 
    //           << " omega/4*pi = " << omega << std::endl;
    
    
    
    //now compute the single scattering function:
    //see if the point is behind the planet
    if (pt.z<0&&pt.x*pt.x+pt.y*pt.y<rminatm*rminatm) {
      lyman_alpha.tau_species[row]=-1.0;
      lyman_alpha.tau_absorber[row]=-1.0;
    } else {
      //the ray points toward +z, towards the sun
      ray = atmo_ray(pt, 0.0, 0.0, 1.0);
      
      grid.voxel_traverse(ray, lyman_alpha, &emission::update_singlescat);
    }
    
    // keep count of the computation
    /* time (&rawtime); */
    /* timeinfo = localtime(&rawtime); */
    /* std::cout << row << " ending at " */
    /* 	      << asctime(timeinfo); */
  }

  //solve for the source function
  lyman_alpha.influence_init=true;
  lyman_alpha.singlescat_init=true;
  lyman_alpha.solve();

  
  // print time elapsed
  clock_t finish = clock();
  double secs = (finish-start)*1.0/CLOCKS_PER_SEC;
  int hrs = secs/3600;
  int mins = (secs - hrs*3600)/60;
  secs = secs - hrs*3600 - mins*60;
  std::cout << "For nH = " << nexo
	    << " , T = " << Texo << " , "
            << "elapsed time is " 
	    << hrs << " hours, " 
	    << mins << " minutes, and " 
	    << secs << " seconds.\n";

  //return values to caller
  Sout = Sobj(nexo, Texo, grid, lyman_alpha.sourcefn_out);
  atmoout = atm;


  return;
}


#endif
