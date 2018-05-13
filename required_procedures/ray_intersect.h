//ray_intersect.h -- routines to compute ray/conic intersections for radiative transfer

#ifndef __ray_intersect_h
#define __ray_intersect_h

#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>
using std::sort;
using std::sqrt;
using std::copysign;
using std::tan;
using std::cos;
using std::sin;
using std::vector;

//we define a stucture to hold the intersections and integers
struct boundary_intersection {
  int boundary;//index of the boundary we're crossing
  double distance;//distance to the boundary crossing
  bool operator<(const boundary_intersection &rhs) const { return distance < rhs.distance; }
  bool operator<=(const boundary_intersection &rhs) const { return distance < rhs.distance; }
};

//return intersections of ray with sphere boundaries
void sphere_intersections(const double &r1,
			  const double &costray,
			  const double &sintray,
			  VecDoub_I & sphereradii,
			  vector<boundary_intersection> &sphereintersections) {

  //loop through the sphere boundaries and compute the intersections
  double d;
  for (int ib=0;ib<sphereradii.size();ib++) {
    //check if an intersection exists
    d=r1*sintray;
    d=sphereradii[ib]*sphereradii[ib]-d*d;//this is the discriminant
					  //of the quadratic
    if (d>0) {
      //non-tangent intersections exist, we should add them
      sphereintersections[2*ib].boundary=ib;
      sphereintersections[2*ib].distance=-r1*costray-sqrt(d);

      sphereintersections[2*ib+1].boundary=ib;
      sphereintersections[2*ib+1].distance=-r1*costray+sqrt(d);
    } else {
      //no non-tangent intersections exist, let's return a negative
      //distance ignored by the wrapper code
      sphereintersections[2*ib].boundary=ib;
      sphereintersections[2*ib].distance=-1.0;

      sphereintersections[2*ib+1].boundary=ib;
      sphereintersections[2*ib+1].distance=-1.0;
    }
  }

  //sort the intersections by distance and return
  sort(sphereintersections.begin(), sphereintersections.end());
}


//return intersections of ray with cone boundaries
void cone_intersections(const double &p_0,
			const double &p_1,
			const double &p_2,
			const double &l_0,
			const double &l_1,
			const double &l_2,
			const double &p,
			const double &costl,
			VecDoub_I & cosconeangles,//cosines of cone angles
			vector<boundary_intersection> &coneintersections) {
  
  //loop through the cone boundaries and compute the intersections
  double A, B, C;
  double d;
  double costhetab;
  double costhetab2;
  double t;
  for (int ib=0;ib<cosconeangles.size();ib++) {
    costhetab=cosconeangles[ib];
    costhetab2=costhetab*costhetab;
    
    A=l_2*l_2;
    A=A-costhetab2;
    
    B=2.0*(p_2*l_2-costhetab2*p*costl);
    
    C=p_2*p_2-costhetab2*p*p;
    
    //check if an intersection exists
    if (A==0) {
      //the equation is not quadratic (there is only one crossing of this boundary)
      //we should solve the linear system
      if (B==0) {
	//either the ray is embedded in the cone or the math above is incorrect
	if (C==0) {
	  toss("ray embedded in cone error");
	} else {
	  //A==0 and B==0, which means -C/B->infinity, so there are no intersections
	  coneintersections[2*ib].boundary=ib;
	  coneintersections[2*ib].distance=-1.0;
	  coneintersections[2*ib+1].boundary=ib;
	  coneintersections[2*ib+1].distance=-1.0;
	}
      } else {
	coneintersections[2*ib].boundary=ib;
	coneintersections[2*ib].distance=-C/B;
	//only one intersection exists,
	//let's put the other one at a negative distance
	//ignored by the wrapper code
	coneintersections[2*ib+1].boundary=ib;
	coneintersections[2*ib+1].distance=-1.0;
      }
    } else {
      //the equation is quadratic and we need to use the quadratic
      //formula to solve it
      d=B*B-4.0*A*C;//this is the discriminant of the quadratic
      if (d>0) {
	//non-tangent intersections exist, we should add them
	coneintersections[2*ib].boundary=ib;
	t = (-B-sqrt(d))/(2.0*A);
	//check if this intersects the right branch of the cone
	t = (copysign(costhetab,p_2+t*l_2)==costhetab) ? t : -1.0;
	coneintersections[2*ib].distance=t;	  

	coneintersections[2*ib+1].boundary=ib;
	t = (-B+sqrt(d))/(2.0*A);
	//check if this intersects the right branch of the cone
	t = (copysign(costhetab,p_2+t*l_2)==costhetab) ? t : -1.0;
	coneintersections[2*ib+1].distance=t;
      } else {
	//no non-tangent intersections exist, let's return a negative distance
	coneintersections[2*ib].boundary=ib;
	coneintersections[2*ib].distance=-1.0;
	  
	coneintersections[2*ib+1].boundary=ib;
	coneintersections[2*ib+1].distance=-1.0;
      }
    }
  }
  
  //sort the intersections by distance and return
  sort(coneintersections.begin(), coneintersections.end());
}


//return intersections of ray with halfplane boundaries
void halfplane_intersections(const double &p_1,
			     const double &p_2,
			     const double &l_1,
			     const double &l_2,
			     VecDoub_I & halfplaneangles,
			     vector<boundary_intersection> &halfplaneintersections) {

  //loop through the halfplane boundaries and compute the intersections
  //loop through the sphere boundaries and compute the intersections
  double d;
  double tanpb;
  double sinpb;
  double t;

  for (int ib=0;ib<halfplaneangles.size()-1;ib++) {
    tanpb=tan(halfplaneangles[ib]);
    sinpb=sin(halfplaneangles[ib]);
    //check if an intersection exists
    d=l_1-tanpb*l_2;
    if (d!=0) {
      //non-tangent intersections exist, we should add them
      halfplaneintersections[ib].boundary=ib;
      t = (tanpb*p_2-p_1)/d;
      //check if intersection is on this half of the plane
      t = (copysign(sinpb,p_2+t*l_1)==sinpb) ? t : -1.0;
      halfplaneintersections[ib].distance=t;
    } else {
      if (p_1/p_2==tanpb) {//need a quadrant-specific test here...
	//the ray is embedded in the plane, something has gone horribly wrong!
	toss("ray in plane error!");
      } else {
	//the ray is parallel to the plane but outside; no intersections.
	halfplaneintersections[ib].boundary=ib;
	halfplaneintersections[ib].distance=-1.0;
      }
    }
  }

  //sort the intersections by distance and return
  sort(halfplaneintersections.begin(), halfplaneintersections.end());

}








#endif
