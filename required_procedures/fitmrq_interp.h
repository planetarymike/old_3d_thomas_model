//fitmrq_interp.h -- Levenberg-Marquardt chi-squared minimization and fitting routines
// see NR p.803 and surroundings for detailed information about how this works

//modified to take an interpolation object as input

#ifndef __FIT_MRQ_H
#define __FIT_MRQ_H

#include "nr3.h"
#include "gaussj.h"
//#include <iostream>
#include "corona_simulator.h"

struct Fitmrq {
  //this object performs non-linear least-squares fitting using the LM
  //method, and as an added bonus allows for some parameters to be
  //held fixed as user-supplied values. Call the constructor to bind
  //the data vectors and input an initial parameter guess. Then call
  //any combination of hold, free, and fit as often as desired. fit
  //sets the output quantities a (parameter list), covar, alpha, and
  //chisq.
  static const int NDONE=40, ITMAX=1000; //convergence parameters
  int ndat, ma, mfit; //number of data points, parameters
  VecDoub_I &x, &y, &sig;//point coordinates and uncertainities
  double tol;//optional convergence parameter, possibly user-specified
  corona_simulator sim;
  //^^^this function represents the user-supplied model to be fitted
  VecBool ia;//boolean wheter to fit each parameter
  VecDoub a;//model parameters
  MatDoub covar;//model parameter covariance matrix
  MatDoub alpha;//model curvature matrix
  double chisq;//model chi-square value


  //constructor: takes data to be fitted, uncertainties, initial
  //parameter guess, and model function as inputs
  Fitmrq(VecDoub_I &xx, VecDoub_I &yy, VecDoub_I &ssig,
	 VecDoub_I &aa,
	 corona_simulator &simm,
	 const double TOL=1.e-3) : ndat(xx.size()), ma(aa.size()),
				   x(xx), y(yy), sig(ssig), sim(simm),
				   tol(TOL), ia(ma), alpha(ma,ma),
				   a(aa), covar(ma,ma) {
    for (int i = 0; i < ma; i++) ia[i] = true;
    //initalize all parameters as free
  }

  
  void hold(const int i, const double val) { ia[i] = false; a[i] = val; }
  void free(const int i) { ia[i] = true; }
  // call these functions to hold or free parameters as needed, before or after fitting

  void fit() {
    //This routine looks for a chi2 minimum given the data points
    //x[0..ndat-1], y[0..ndat-1], sig[0..ndat-1], and the nonlinear
    //model function funcs that depends on the model parameters
    //a[0..ma-1]. It searches until chi2 is no longer decreasing, then
    //returns the ch2 value, parameters, parameter covariance matrix,
    //and curvature matrix.

    int j, k, l, iter, done=0;
    double alambda=0.001,ochisq;
    VecDoub atry(ma), beta(ma), da(ma);
    mfit=0;
    for (j = 0; j < ma; j++) if (ia[j]) mfit++;
    MatDoub oneda(mfit,1), temp(mfit,mfit);
    mrqcof(a,alpha,beta); //initialize

    
    //    for (int irow = 0; irow < alpha.nrows(); irow++) {
    //      for (int icol = 0; icol < alpha.ncols(); icol++) {
    // 	std::cout << "alpha[" << irow << "][" << icol << "] = ";
    // 	std::cout << alpha[irow][icol] << "\t";
    //      }
    //      std::cout << std::endl;
    //      std::cin.get();
    //    }

    for (j = 0; j < ma; j++) atry[j] = a[j];
    ochisq=chisq;
    for (iter = 0; iter < ITMAX; iter++) {
      if (done == NDONE) alambda = 0.; //if we've found the minimum,
				       //set lambda=0 to return the
				       //full curvature matrix
      for (j = 0; j < mfit; j++) { //construct the L-M descent matrix
	for (k = 0; k < mfit; k++) {
	  covar[j][k] = alpha[j][k];
	}
	covar[j][j] = alpha[j][j]*(1.0+alambda);
	for (k = 0; k < mfit; k++) {
	  temp[j][k] = covar[j][k];
	}
	oneda[j][0] = beta[j];
      }
      //      for (int irow = 0; irow < covar.nrows(); irow++) {
      //       	for (int icol = 0; icol < covar.ncols(); icol++) {
      //       	  std::cout << "covar[" << irow << "][" << icol << "] = ";
      //       	  std::cout << covar[irow][icol] << "\t";
      //       	}
      //       	std::cout << std::endl;
      //      }
      //      for (int irow = 0; irow < temp.nrows(); irow++) {
      //       	for (int icol = 0; icol < temp.ncols(); icol++) {
      //       	  std::cout << "temp[" << irow << "][" << icol << "] = ";
      //       	  std::cout << temp[irow][icol] << "\t";
      //       	}
      //       	std::cout << std::endl;
      //      }

      gaussj(temp,oneda);//replaces temp with its inverse, oneda with
			 //the solution vector
      for (j = 0; j < mfit; j++) {
	for (k = 0; k < mfit; k++) {
	  covar[j][k] = temp[j][k];
	}
	da[j] = oneda[j][0];
      }
      if (done == NDONE) {//done, so insert zeros back in for
			  //non-fitted parameters
	covsrt(covar);
	covsrt(alpha);
	return;
      }
      for (j = 0, l = 0; l < ma; l++) //increment (non-constant) parameters by da
	if (ia[l]) atry[l] = a[l] + da[j++];
      mrqcof(atry, covar, da); //determine new curvature matrix, gradient, and chi2
      if (abs(chisq-ochisq) < MAX(tol, tol*chisq)) done++; //if the
							   //step is
							   //too
							   //small,
							   //increment
							   //the
							   //doomsday
							   //counter
      if (chisq < ochisq) { //if the new chi2 is smaller, accept
	alambda *= 0.1;
	ochisq = chisq;
	for (j = 0; j < mfit; j++) {
	  for (k = 0; k < mfit; k++) alpha[j][k] = covar[j][k];
	  beta[j] = da[j];
	}
	for (l = 0; l < ma; l++) a[l] = atry[l];
      } else { // if new chi2 is larger, increase alambda and try again
	alambda *= 10.0;
	chisq = ochisq;
      }
    }
    throw("Fitmrq:: too many iterations");
  }

  void mrqcof(VecDoub_I &a, MatDoub_O &alpha, VecDoub_O &beta) {
    //Computes the local curvature matrix and gradient at parameters a[0..ma-1]
    int i, j, k, l, m;
    double ymod, wt, sig2i, dy;
    VecDoub dyda(ma);
    for (j = 0; j < mfit; j++) {//initialize alpha and beta
      for (k = 0; k <= j; k++) {
	alpha[j][k] = 0.0; //alpha is symmetric, only initialize half
      }
      beta[j] = 0.0;      
    }
    chisq = 0.0; //initialize chisq
    for (i = 0; i < ndat; i++) { // loop over all data
      //      std::cout << "\n\nIn mrqcof, i = " << i << std::endl; 
      sim(x[i],a,ymod,dyda); //get the local fit value and derivs
      //      std::cout << "\t x[i] = " << x[i] << std::endl; 
      //      std::cout << "\t ymod = " << ymod << std::endl; 
      sig2i = 1.0/(sig[i]*sig[i]); //inverse variance
      dy = y[i]-ymod;//deviation
      chisq += dy*dy*sig2i; // find chisq
      for (j = 0, l = 0; l < ma; l++) {
	if (ia[l]) {
	  wt = dyda[l] * sig2i;
	  //	  std::cout << "dyda[" << l << "] = " << dyda[l] << "\n";
	  for (k = 0, m = 0; m < l+1; m++)
	    if (ia[m]) {
	      //	      std::cout << "\n j = " << j 
	      //	       		<< "\n k = " << k
	      //			<< ", wt*dyda[m] = " << wt*dyda[m] << std::endl;
	      //	      std::cin.get();
	      alpha[j][k++] += wt*dyda[m];
	    }
	       //^^^and the paramater curvature matrix (see NR
	       //p. 800-801, or take deriv of chi2 analytically wrt
	       //params)
	  beta[j++] += dy*wt;//and the gradient
	}
      }
    }
    for (j = 1; j < mfit; j++) //fill in the symmetric side of the curva mat
      for (k = 0; k < j; k++) alpha[k][j] = alpha[j][k];
  }

  void covsrt(MatDoub_IO &covar) {
    //take the covariance matrix and clean up for output to the user,
    //inserting zeroes where the user has specifed that a parameter
    //take on a fixed value
    int i, j, k;
    for (i = mfit; i < ma; i++) 
      for (j = 0; j < i+1; j++) covar[i][j]=covar[j][i]=0.0;
    k = mfit - 1;
    for (j = ma - 1; j >= 0; j--) {
      if (ia[j]) {
	for (i = 0; i < ma; i++) SWAP(covar[i][k],covar[i][j]);
	for (i = 0; i < ma; i++) SWAP(covar[k][i],covar[j][i]);
	k--;
      }
    }
  }

};















#endif	     
