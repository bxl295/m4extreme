// MLS.cpp: implementation of the MLS class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
////////////////////////////////////////////////////////////////////////////////

#include <cassert>
#include "./MLS.h"

#include <cstdio>
#include <ctime>
#include <sys/time.h>

namespace Element
{
  namespace Interpolation
  {
    namespace MLS
    {

      // Class Shape<0>

      Shape<0>::Shape(const Base *B_) : B(B_) {}
			
      Shape<0>::~Shape() {}

      map<Set::Manifold::Point *, double> 
      Shape<0>::operator () (const Set::VectorSpace::Vector &y) const 
      {
	struct timeval start_reset, stop_reset;
	gettimeofday(&start_reset, NULL);
	int i,j,l;
 
	// dimension
	unsigned int dim = B->GetDim();
	
	
	
	/*	std::cout << "mPt = " << std::endl;
		for (i=0;i<dim;i++){
		std::cout << y[i] << " ";
		}
		std::cout << std::endl << std::endl;
	*/
 

 
	// get nodes
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> & x = *(B->GetNodes());
 
	// get polynomial order
	int m = B->GetOrder();
 
	// get radius
	double beta = B->GetRadius();
 
 	// 	dilatation parameter 
	double rho = sqrt(1/beta);
	
 	// calculation of the number of term in the polynomial base
	int k;
	if (dim == 1) k = m+1;
	if (dim == 2) k = (m+1)*(m+2)/2;
	if (dim == 3) k = (m+1)*(m+2)*(m+3)/6;


	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator px;
 
	// number of node
	int nNode = x.size();
 
 
 
 	/////////////////////////////////////////////////////////
	// definition of polynomial basis 
	
	// Polynomial basis for the mPt
	double Px[k];
 
	// 1D
	if (dim == 1) {
	  Px[0] = 1;
	  for (i=1;i<k;i++) {
	    Px[i] = Px[i-1]*y[0];
	  }
	}
 
	// 2D
	if (dim == 2) {
	  double tabx[m+1];
	  double taby[m+1];
		
	  // define 1st line and first column
	  tabx[0] = 1;
	  taby[0] = 1;
	  for (i=1;i<m+1;i++) {
	    tabx[i] = tabx[i-1]*y[0];
	    taby[i] = taby[i-1]*y[1];
	  }
		
	  /*		for (i=1;i<m+1;i++) {
			for (j=1;j<m+1;j++)
			tab[j][i] = tab[j][0]*tab[0][i];
			}
		
			int inc = 0;
			for (i=0;i<m+1;i++) {
			for (j=0;j<i+1;j++) {
			Px[inc] = tab[i-j][j];
			inc++;
			}
			}*/
	  int inc = 0;
	  for (i=0;i<m+1;i++) {
	    for (j=0;j<i+1;j++) {
	      Px[inc] = tabx[i-j]*taby[j];
	      inc++;
	    }
	  }
		
	}	
	
	// 3D
	if (dim == 3) {
	  if (m==1){
	    Px[0] = 1;
	    Px[1] = y[0];
	    Px[2] = y[1];
	    Px[3] = y[2];
	  }
	  if (m == 2) {
	    Px[0] = 1;
	    Px[1] = y[0];
	    Px[2] = y[1];
	    Px[3] = y[2];
	    Px[4] = y[0]*y[0];
	    Px[5] = y[0]*y[1];
	    Px[6] = y[1]*y[1];
	    Px[7] = y[1]*y[2];
	    Px[8] = y[2]*y[2];
	    Px[9] = y[0]*y[2];
	  }
	}
	
	// Polynomial basis for the nodes
	//Set::VectorSpace::Hom P(nNode,k);
	double* P = new double[nNode*k];
	
	l=0;
	for (px=x.begin(); px!=x.end(); px++) {
	  // 1D
	  if (dim == 1) {
			
	    P[l*k+0] = 1;
	    for (i=1;i<k;i++) {
	      P[l*k+i] = P[l*k+(i-1)]*px->second[0];
	    }
	  }
	
	  // 2D
	  if (dim == 2) {
			
	    double tab[m+1][m+1];
	    // define 1st line and first column
	    tab[0][0] = 1.;
	    for (i=1;i<m+1;i++) {
	      tab[0][i] = tab[0][i-1]*px->second[0];
	      tab[i][0] = tab[i-1][0]*px->second[1];
	    }
			
	    if (m>1){
	      for (i=1;i<m+1;i++) {
		for (j=1;j<m+1;j++) {
		  tab[i][j] = tab[0][j]*tab[i][0];
		}
	      }
	    }
			
	    int inc = 0;
	    for (i=0;i<m+1;i++) {
	      for (j=0;j<i+1;j++) {
		P[l*k+inc] = tab[j][i-j];
		inc++;
	      }
	    }
	  }
	  if (dim == 3) {
			
	    if (m==1){
	      P[l*k] = 1;
	      P[l*k+1] = px->second[0];
	      P[l*k+2] = px->second[1];
	      P[l*k+3] = px->second[2];
	    }
	    if (m == 2) {
	      P[l*k] = 1;
	      P[l*k+1] = px->second[0];
	      P[l*k+2] = px->second[1];
	      P[l*k+3] = px->second[2];
	      P[l*k+4] = px->second[0]*px->second[0];
	      P[l*k+5] = px->second[0]*px->second[1];
	      P[l*k+6] = px->second[1]*px->second[1];
	      P[l*k+7] = px->second[1]*px->second[2];
	      P[l*k+8] = px->second[2]*px->second[2];
	      P[l*k+9] = px->second[0]*px->second[2];
	    }
	  }
	  l++;
	}

	//////////////////////////////////////////////////
	// Definition of the weight W
 
	double W[nNode];
	
	
	l=0;
	for (px=x.begin(); px!=x.end(); px++) {
	  // distance between mPt and nodes
	  double dx=0;
	  for (i=0;i<dim;i++){
	    dx+= (y[i] - px->second[i])*(y[i] - px->second[i]);
	  }
		
	  double q = beta*dx;
	  W[l]=exp(-q);
		
	  /*		double q = sqrt(dx)/rho;
		
			if (q > 1.){
			W[l]=0.;
			}	
			else {
			W[l]= 1-6*q*q+8*q*q*q-3*q*q*q*q;
			}
	  */
	  l++;
	}
	 
	/////////////////////////////////////////////////////////
	// Definition of B
	
	//Set::VectorSpace::Hom B(k,nNode);	
	double* B = new double[k*nNode];
	
	for (i=0;i<k;i++) {
	  for (j=0;j<nNode;j++) {
	    B[i*nNode+j] = P[j*k+i]*W[j];
	  }
	}
 
 	/////////////////////////////////////////////////////////
	// Definition of A
	
	Set::VectorSpace::Hom A(k,k);	
 
	// Initialization of A
 	for (i=0;i<k;i++) {
	  for (j=0;j<k;j++) {
	    A[j][i] = 0;
	  }
	}
	
	for (i=0;i<k;i++) {
	  for (j=0;j<k;j++) {
	    //A[j][i] = 0;
	    for (l=0;l<nNode;l++)
	      A[j][i] += W[l]*P[l*k+i]*P[l*k+j];
	  }
	}	
	delete[] P;
 	/////////////////////////////////////////////////////////
	// Definition of AI
	
	Set::VectorSpace::Hom AI(k,k);
		
	AI = Inverse(A);
	/////////////////////////////////////////////////////////
	// Definition of C= AI*B
	Set::VectorSpace::Vector D(nNode);
 
	for (j=0;j<nNode;j++) {
	  D[j]=0;
	  for (i=0;i<k;i++) {
	    for (l=0;l<k;l++) {
	      D[j] += Px[i]*AI[l][i]*B[l*nNode+j];
	    }
	  }
	}

	delete[] B;
	/////////////////////////////////////////////////////////
	// Build of N
	
	map<Set::Manifold::Point *, double> N;
	
	// loop on nodes
	j=0;
	for (px=x.begin(); px!=x.end(); px++)
	  {
	    N.insert(make_pair(px->first,D[j]));
	    j++;
	  }
	return N;

      }

      /////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////
      // 								Class Shape<1>
      /////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////


      Shape<1>::Shape(const Base *B0) : B(B0) {}
	
      Shape<1>::~Shape() {}


      //////////////////////////////////////////////////////////////////////////////
      //		INSERT HERE THE COMPUTATION OF THE DIFFUSE DERIVATIVE OF THE MLS SHAPE FUNCTIONS

      map<Set::Manifold::Point *, Set::VectorSpace::Vector> 
      Shape<1>::operator () (const Set::VectorSpace::Vector &y) const
      {

	int i,j,l;
 
	// dimension
	unsigned int dim = B->GetDim();
	
	// get nodes
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> & x = *(B->GetNodes());
 
	// get polynomial order
	int m = B->GetOrder();
 
	// get radius
	double beta = B->GetRadius();
	
 	// 	dilatation parameter 
	double rho = sqrt(1/beta);
 
 	// calculation of the number of term in the polynomial base
	int k;
	if (dim == 1) k = m+1;
	if (dim == 2) k = (m+1)*(m+2)/2;
	if (dim == 3) k = (m+1)*(m+2)*(m+3)/6;


	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator px;
 
	// number of node
	int nNode = x.size();

	
	
	
	/////////////////////////////////////////////////////////
	// definition of polynomial basis 
	
	// Polynomial basis for the mPt
	double Px[k];
 
	// 1D
	if (dim == 1) {
	  Px[0] = 1;
	  for (i=1;i<k;i++) {
	    Px[i] = Px[i-1]*y[0];
	  }
	}
 
	// 2D
	if (dim == 2) {
	  double tabx[m+1];
	  double taby[m+1];
		
	  // define 1st line and first column
	  tabx[0] = 1;
	  taby[0] = 1;
	  for (i=1;i<m+1;i++) {
	    tabx[i] = tabx[i-1]*y[0];
	    taby[i] = taby[i-1]*y[1];
	  }
		
	  /*		for (i=1;i<m+1;i++) {
			for (j=1;j<m+1;j++)
			tab[j][i] = tab[j][0]*tab[0][i];
			}
		
			int inc = 0;
			for (i=0;i<m+1;i++) {
			for (j=0;j<i+1;j++) {
			Px[inc] = tab[i-j][j];
			inc++;
			}
			}*/
	  int inc = 0;
	  for (i=0;i<m+1;i++) {
	    for (j=0;j<i+1;j++) {
	      Px[inc] = tabx[i-j]*taby[j];
	      inc++;
	    }
	  }
		
	}	
	
	// 3D
	if (dim == 3) {
	  if (m==1){
	    Px[0] = 1;
	    Px[1] = y[0];
	    Px[2] = y[1];
	    Px[3] = y[2];
	  }
	  if (m == 2) {
	    Px[0] = 1;
	    Px[1] = y[0];
	    Px[2] = y[1];
	    Px[3] = y[2];
	    Px[4] = y[0]*y[0];
	    Px[5] = y[0]*y[1];
	    Px[6] = y[1]*y[1];
	    Px[7] = y[1]*y[2];
	    Px[8] = y[2]*y[2];
	    Px[9] = y[0]*y[2];
	  }
	}
	
	// Polynomial basis for the mPt

	double DPx[k], DPy[k], DPz[k];

	// 1D
	if (dim == 1) {
	  if (m==1){
	    DPx[0]=0;
	    DPx[1]=1;
	  }
	  if (m>=2){
	    DPx[0]=0;
	    DPx[1]=1;
	    //DPx[3]=3*y[0]*y[0];
	    for (i=2;i<m+1;i++){
	      DPx[i] = y[0]*DPx[i-1]*i/(i-1);
	    }
			
	  }	
	}
	
	// 2D
	if (dim == 2) {
	  // define polynomial table
	  if (m==1){
	    DPx[0]=0;
	    DPx[1]=1;
	    DPx[2]=0;
	    DPy[0]=0;
	    DPy[1]=0;
	    DPy[2]=1;
	  }
	  if (m>=2){
	    double tabx[m+1][m+1];
	    double taby[m+1][m+1];
	    // define 1st line and first column
	    tabx[0][0] = 0;
	    tabx[0][1] = 1;
	    tabx[1][0] = y[1];
	    taby[0][0] = 0;
	    taby[0][1] = y[0];
	    taby[1][0] = 1;
			
	    for (i=2;i<m+1;i++) {
	      tabx[0][i] = y[0]*tabx[0][i-1]*i/(i-1);
	      tabx[i][0] = y[1]*tabx[i-1][0];
	      taby[0][i] = y[0]*taby[0][i-1];
	      taby[i][0] = y[1]*taby[i-1][0]*i/(i-1);
	    }
	    for (i=1;i<m+1;i++) {
	      for (j=1;j<m+1;j++){
		tabx[i][j] = tabx[0][j]*tabx[i][0];
		taby[i][j] = taby[0][j]*taby[i][0];
	      }
	    }
	    int inc = 0;
	    for (i=0;i<m+1;i++) {
	      for (j=0;j<i+1;j++) {
		if (i==j){
		  DPx[inc] = 0;
		}
		else {
		  DPx[inc] = tabx[j][i-j];
		}
		if (j==0){
		  DPy[inc] = 0;
		}
		else {
		  DPy[inc] = taby[j][i-j];
		}
		inc++;
	      }
	    }
	  }		
	}

	// 3D
	if (dim == 3) {
	  // define polynomial table
	  if (m==1){
	    DPx[0]=0;
	    DPx[1]=1;
	    DPx[2]=0;
	    DPx[3]=0;
	    DPy[0]=0;
	    DPy[1]=0;
	    DPy[2]=1;
	    DPy[3]=0;
	    DPz[0]=0;
	    DPz[1]=0;
	    DPz[2]=0;
	    DPz[3]=1;			
	  }
	  if (m==2){
	    DPx[0] = 0;
	    DPx[1] = 1;
	    DPx[2] = 0;
	    DPx[3] = 0;
	    DPx[4] = 2*y[0];
	    DPx[5] = y[1];
	    DPx[6] = 0;
	    DPx[7] = 0;
	    DPx[8] = 0;
	    DPx[9] = y[2];
			
	    DPy[0] = 0;
	    DPy[1] = 0;
	    DPy[2] = 1;
	    DPy[3] = 0;
	    DPy[4] = 0;
	    DPy[5] = y[0];
	    DPy[6] = 2*y[1];
	    DPy[7] = y[2];
	    DPy[8] = 0;
	    DPy[9] = 0;
			
	    DPz[0] = 0;
	    DPz[1] = 0;
	    DPz[2] = 0;
	    DPz[3] = 1;
	    DPz[4] = 0;
	    DPz[5] = 0;
	    DPz[6] = 0;
	    DPz[7] = y[1];
	    DPz[8] = 2*y[2];
	    DPz[9] = y[0];			
	  }
	}
			
	// Polynomial basis for the nodes
	double* P = new double[nNode*k];
	l=0;
	for (px=x.begin(); px!=x.end(); px++) {
	  // 1D
	  if (dim == 1) {
			
	    P[l*k+0] = 1;
	    for (i=1;i<k;i++) {
	      P[l*k+i] = P[l*k+(i-1)]*px->second[0];
	    }
	  }
	  // 2D
	  if (dim == 2) {
			
	    double tab[m+1][m+1];
	    // define 1st line and first column
	    tab[0][0] = 1.;
	    for (i=1;i<m+1;i++) {
	      tab[0][i] = tab[0][i-1]*px->second[0];
	      tab[i][0] = tab[i-1][0]*px->second[1];
	    }
			
	    if (m>1){
	      for (i=1;i<m+1;i++) {
		for (j=1;j<m+1;j++) {
		  tab[i][j] = tab[0][j]*tab[i][0];
		}
	      }
	    }
			
	    int inc = 0;
	    for (i=0;i<m+1;i++) {
	      for (j=0;j<i+1;j++) {
		P[l*k+inc] = tab[j][i-j];
		inc++;
	      }
	    }
	  }
	  if (dim == 3) {
			
	    if (m==1){
	      P[l*k] = 1;
	      P[l*k+1] = px->second[0];
	      P[l*k+2] = px->second[1];
	      P[l*k+3] = px->second[2];
	    }
	    if (m == 2) {
	      P[l*k] = 1;
	      P[l*k+1] = px->second[0];
	      P[l*k+2] = px->second[1];
	      P[l*k+3] = px->second[2];
	      P[l*k+4] = px->second[0]*px->second[0];
	      P[l*k+5] = px->second[0]*px->second[1];
	      P[l*k+6] = px->second[1]*px->second[1];
	      P[l*k+7] = px->second[1]*px->second[2];
	      P[l*k+8] = px->second[2]*px->second[2];
	      P[l*k+9] = px->second[0]*px->second[2];
	    }
	  }
	  l++;
	}
	//////////////////////////////////////////////////
	// Definition of the weight W
 
	double W[nNode], Wx[nNode], Wy[nNode], Wz[nNode];
	
	
	l=0;
	for (px=x.begin(); px!=x.end(); px++) {
	  // distance between mPt and nodes
	  double dx=0;
	  for (i=0;i<dim;i++){
	    dx+= (y[i] - px->second[i])*(y[i] - px->second[i]);
	  }
		
	  double q = beta*dx;
	  W[l]=exp(-q);
	  Wx[l]=-2*beta*(y[0]-px->second[0])*W[l];
	  if (dim >= 2) {Wy[l]=-2*beta*(y[1]-px->second[1])*W[l];}
	  if (dim >= 3) {Wz[l]=-2*beta*(y[2]-px->second[2])*W[l];}
		
		
	  l++;
	}

	/////////////////////////////////////////////////////////
	// Definition of B, Bx, By
	
 	/////////////////////////////////////////////////////////
	// Definition of A, Ax, Ay
	
	Set::VectorSpace::Hom A(k,k);
	double Ax[k][k], Ay[k][k], Az[k][k];
 
	// Initialization of A
 	for (i=0;i<k;i++) {
	  for (j=0;j<k;j++) {
	    A[j][i] = 0;
	    Ax[i][j] = 0;
	    if (dim >= 2) {Ay[i][j] = 0;}
	    if (dim >= 3) {Az[i][j] = 0;}
	  }
	}

	for (i=0;i<k;i++) {
	  for (j=0;j<k;j++) {
	    for (l=0;l<nNode;l++){
	      A[j][i] += W[l]*P[l*k+i]*P[l*k+j];
	      Ax[i][j] += Wx[l]*P[l*k+i]*P[l*k+j];
	      if (dim >= 2) {Ay[i][j] += Wy[l]*P[l*k+i]*P[l*k+j];}
	      if (dim >= 3) {Az[i][j] += Wz[l]*P[l*k+i]*P[l*k+j];}
	    }
	  }
	}	
	
	delete[] P;
 	/////////////////////////////////////////////////////////
	// Definition of AI
	
	Set::VectorSpace::Hom AI(k,k);

	/*for (i=0;i<k;i++) {
	  for (j=0;j<k;j++){
	  std::cout << A[j][i] << " ";
	  }
	  std::cout << "; " << std::endl;
	  }
	*/
	
	AI = Inverse(A);
	//////////////////////////////////////////////////////////////////////
	//		Computation of the derivatives
	

	double alpha[k];
	
	for (i=0;i<k;i++) {
	  alpha[i]=0;
	  for (j=0;j<k;j++){
	    alpha[i]+=AI[j][i]*Px[j];
	  }
	} 

	double alphaX[k],alphaY[k],alphaZ[k];
	double Axalpha[k],Ayalpha[k],Azalpha[k];

	for (i=0;i<k;i++) {
	  Axalpha[i]=0;Ayalpha[i]=0;Azalpha[i]=0;
	  for (j=0;j<k;j++){
	    if (dim >= 1) Axalpha[i] += Ax[i][j]*alpha[j];
	    if (dim >= 2) Ayalpha[i] += Ay[i][j]*alpha[j];
	    if (dim >= 3) Azalpha[i] += Az[i][j]*alpha[j];
	  }
	}

	for (i=0;i<k;i++) {
	  alphaX[i]=0;alphaY[i]=0;alphaZ[i]=0;
	  for (j=0;j<k;j++){
	    if (dim >= 1) alphaX[i] += AI[j][i]*(DPx[j]-Axalpha[j]);
	    if (dim >= 2) alphaY[i] += AI[j][i]*(DPy[j]-Ayalpha[j]);
	    if (dim >= 3) alphaZ[i] += AI[j][i]*(DPz[j]-Azalpha[j]);
	  }
	}	
	
	
	double* D = new double[dim*nNode];
	double DP[k];
	
	l=0;
	for (px=x.begin(); px!=x.end(); px++){
	
	  //computation of DP
	  DP[0]=1;
	  if (dim == 1){
	    DP[1]=px->second[0];
	    if (m>=2) DP[2]=(px->second[0])*(px->second[0]);
	    if (m>=3) DP[3]=(px->second[0])*(px->second[0])*(px->second[0]);		
	  }
	  if (dim == 2){
	    DP[1]=px->second[0];
	    DP[2]=px->second[1];
	    if (m>=2) {
	      DP[3]=(px->second[0])*(px->second[0]);
	      DP[4]=(px->second[0])*(px->second[1]);
	      DP[5]=(px->second[1])*(px->second[1]);
	    }
	    if (m>=3) {
	      DP[6]=(px->second[0])*(px->second[0])*(px->second[0]);
	      DP[7]=(px->second[0])*(px->second[0])*(px->second[1]);
	      DP[8]=(px->second[0])*(px->second[1])*(px->second[1]);
	      DP[9]=(px->second[1])*(px->second[1])*(px->second[1]);
	    }	
	  }
	  if (dim == 3){
	    DP[1]=px->second[0];
	    DP[2]=px->second[1];
	    DP[3]=px->second[2];
	    if (m>=2) {
	      DP[4]=(px->second[0])*(px->second[0]);
	      DP[5]=(px->second[0])*(px->second[1]);
	      DP[6]=(px->second[1])*(px->second[1]);
	      DP[7]=(px->second[1])*(px->second[2]);
	      DP[8]=(px->second[2])*(px->second[2]);
	      DP[9]=(px->second[0])*(px->second[2]);
	    }
	  }
		
	  double alphaDP=0,alphaxDP=0,alphayDP=0,alphazDP=0;
	  for (i=0;i<k;i++) {
	    alphaDP += alpha[i]*DP[i];
	    alphaxDP += alphaX[i]*DP[i];
	    if (dim>=2) alphayDP += alphaY[i]*DP[i];
	    if (dim>=3) alphazDP += alphaZ[i]*DP[i];
	  }
		
	  D[l] = alphaDP*Wx[l] + alphaxDP*W[l];
	  if (dim>=2) D[nNode+l] = alphaDP*Wy[l] + alphayDP*W[l];
	  if (dim>=3) D[2*nNode+l] = alphaDP*Wz[l] + alphazDP*W[l];
				
	  l++;
	}

	/////////////////////////////////////////////////////////
	// Build of DN
	// loop on nodes
	map<Set::Manifold::Point *, Set::VectorSpace::Vector> DN;
	Set::VectorSpace::Vector Dv(dim);
	j=0;
	for (px=x.begin(); px!=x.end(); px++)
	  {
	    for (i=0;i<dim;i++){
	      Dv[i]=D[i*nNode+j];
	    }
	    DN.insert(make_pair(px->first,Dv));
	    j++;
	  }
	delete[] D;
	return DN;

      }		
      ////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////
      //							 Class Jet<0>
      ////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////

      Jet<0>::Jet(const Element::Interpolation::MLS::Base * B0) 
	: B(B0){}
			
      Jet<0>::~Jet() {}
	
      void Jet < 0 > ::operator () (const domain_type & y, shape_type & N, shape_derivative_type & DN) const {
  
	if ( !N.empty() ) N.clear();
	if ( !DN.empty() ) DN.clear();
	
	int i,j,l;
	int ii,ij,ik,il,im;
 
	// dimension
	unsigned int dim = B->GetDim();
	
	// get nodes
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> & x = *(B->GetNodes());
 
	// get polynomial order
	int m = B->GetOrder();
 
	// get radius
	double beta = B->GetRadius();
	
 	// 	dilatation parameter 
	double rho = sqrt(1/beta);
 
 	// calculation of the number of term in the polynomial base
	int k = B->GetNumberOfTerms();

	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator px;
 
	// number of node
	int nNode = x.size();

	/////////////////////////////////////////////////////////
	// definition of polynomial basis 
	
	// Polynomial basis for the mPt
	double Px[k];
 
	// 1D
	if (dim == 1) {
	  Px[0] = 1;
	  for (i=1;i<k;i++) {
	    Px[i] = Px[i-1]*y[0];
	  }
	}
 
	// 2D
	if (dim == 2) {
	  double tabx[m+1];
	  double taby[m+1];
		
	  // define 1st line and first column
	  tabx[0] = 1;
	  taby[0] = 1;
	  for (i=1;i<m+1;i++) {
	    tabx[i] = tabx[i-1]*y[0];
	    taby[i] = taby[i-1]*y[1];
	  }
		
	  /*		for (i=1;i<m+1;i++) {
			for (j=1;j<m+1;j++)
			tab[j][i] = tab[j][0]*tab[0][i];
			}
		
			int inc = 0;
			for (i=0;i<m+1;i++) {
			for (j=0;j<i+1;j++) {
			Px[inc] = tab[i-j][j];
			inc++;
			}
			}*/
	  int inc = 0;
	  for (i=0;i<m+1;i++) {
	    for (j=0;j<i+1;j++) {
	      Px[inc] = tabx[i-j]*taby[j];
	      inc++;
	    }
	  }
		
	}	
	
	// 3D
	if (dim == 3) {
	  if (m==1){
	    Px[0] = 1;
	    Px[1] = y[0];
	    Px[2] = y[1];
	    Px[3] = y[2];
	  }
	  if (m == 2) {
	    Px[0] = 1;
	    Px[1] = y[0];
	    Px[2] = y[1];
	    Px[3] = y[2];
	    Px[4] = y[0]*y[0];
	    Px[5] = y[0]*y[1];
	    Px[6] = y[1]*y[1];
	    Px[7] = y[1]*y[2];
	    Px[8] = y[2]*y[2];
	    Px[9] = y[0]*y[2];
	  }
	}
	
	// Polynomial basis for the mPt

	double DPx[k], DPy[k], DPz[k];

	// 1D
	if (dim == 1) {
	  if (m==1){
	    DPx[0]=0;
	    DPx[1]=1;
	  }
	  if (m>=2){
	    DPx[0]=0;
	    DPx[1]=1;
	    //DPx[3]=3*y[0]*y[0];
	    for (i=2;i<m+1;i++){
	      DPx[i] = y[0]*DPx[i-1]*i/(i-1);
	    }
			
	  }	
	}
	
	// 2D
	if (dim == 2) {
	  // define polynomial table
	  if (m==1){
	    DPx[0]=0;
	    DPx[1]=1;
	    DPx[2]=0;
	    DPy[0]=0;
	    DPy[1]=0;
	    DPy[2]=1;
	  }
	  if (m>=2){
	    double tabx[m+1][m+1];
	    double taby[m+1][m+1];
	    // define 1st line and first column
	    tabx[0][0] = 0;
	    tabx[0][1] = 1;
	    tabx[1][0] = y[1];
	    taby[0][0] = 0;
	    taby[0][1] = y[0];
	    taby[1][0] = 1;
			
	    for (i=2;i<m+1;i++) {
	      tabx[0][i] = y[0]*tabx[0][i-1]*i/(i-1);
	      tabx[i][0] = y[1]*tabx[i-1][0];
	      taby[0][i] = y[0]*taby[0][i-1];
	      taby[i][0] = y[1]*taby[i-1][0]*i/(i-1);
	    }
	    for (i=1;i<m+1;i++) {
	      for (j=1;j<m+1;j++){
		tabx[i][j] = tabx[0][j]*tabx[i][0];
		taby[i][j] = taby[0][j]*taby[i][0];
	      }
	    }
	    int inc = 0;
	    for (i=0;i<m+1;i++) {
	      for (j=0;j<i+1;j++) {
		if (i==j){
		  DPx[inc] = 0;
		}
		else {
		  DPx[inc] = tabx[j][i-j];
		}
		if (j==0){
		  DPy[inc] = 0;
		}
		else {
		  DPy[inc] = taby[j][i-j];
		}
		inc++;
	      }
	    }
	  }		
	}

	// 3D
	if (dim == 3) {
	  // define polynomial table
	  if (m==1){
	    DPx[0]=0;
	    DPx[1]=1;
	    DPx[2]=0;
	    DPx[3]=0;
	    DPy[0]=0;
	    DPy[1]=0;
	    DPy[2]=1;
	    DPy[3]=0;
	    DPz[0]=0;
	    DPz[1]=0;
	    DPz[2]=0;
	    DPz[3]=1;			
	  }
	  if (m==2){
	    DPx[0] = 0;
	    DPx[1] = 1;
	    DPx[2] = 0;
	    DPx[3] = 0;
	    DPx[4] = 2*y[0];
	    DPx[5] = y[1];
	    DPx[6] = 0;
	    DPx[7] = 0;
	    DPx[8] = 0;
	    DPx[9] = y[2];
			
	    DPy[0] = 0;
	    DPy[1] = 0;
	    DPy[2] = 1;
	    DPy[3] = 0;
	    DPy[4] = 0;
	    DPy[5] = y[0];
	    DPy[6] = 2*y[1];
	    DPy[7] = y[2];
	    DPy[8] = 0;
	    DPy[9] = 0;
			
	    DPz[0] = 0;
	    DPz[1] = 0;
	    DPz[2] = 0;
	    DPz[3] = 1;
	    DPz[4] = 0;
	    DPz[5] = 0;
	    DPz[6] = 0;
	    DPz[7] = y[1];
	    DPz[8] = 2*y[2];
	    DPz[9] = y[0];			
	  }
	}
			
	// Polynomial basis for the nodes
	double* P = new double[nNode*k];
	l=0;
	for (px=x.begin(); px!=x.end(); px++) {
	  // 1D
	  if (dim == 1) {
			
	    P[l*k+0] = 1;
	    for (i=1;i<k;i++) {
	      P[l*k+i] = P[l*k+(i-1)]*px->second[0];
	    }
	  }
	  // 2D
	  if (dim == 2) {
			
	    double tab[m+1][m+1];
	    // define 1st line and first column
	    tab[0][0] = 1.;
	    for (i=1;i<m+1;i++) {
	      tab[0][i] = tab[0][i-1]*px->second[0];
	      tab[i][0] = tab[i-1][0]*px->second[1];
	    }
			
	    if (m>1){
	      for (i=1;i<m+1;i++) {
		for (j=1;j<m+1;j++) {
		  tab[i][j] = tab[0][j]*tab[i][0];
		}
	      }
	    }
			
	    int inc = 0;
	    for (i=0;i<m+1;i++) {
	      for (j=0;j<i+1;j++) {
		P[l*k+inc] = tab[j][i-j];
		inc++;
	      }
	    }
	  }
	  if (dim == 3) {
			
	    if (m==1){
	      P[l*k] = 1;
	      P[l*k+1] = px->second[0];
	      P[l*k+2] = px->second[1];
	      P[l*k+3] = px->second[2];
	    }
	    if (m == 2) {
	      P[l*k] = 1;
	      P[l*k+1] = px->second[0];
	      P[l*k+2] = px->second[1];
	      P[l*k+3] = px->second[2];
	      P[l*k+4] = px->second[0]*px->second[0];
	      P[l*k+5] = px->second[0]*px->second[1];
	      P[l*k+6] = px->second[1]*px->second[1];
	      P[l*k+7] = px->second[1]*px->second[2];
	      P[l*k+8] = px->second[2]*px->second[2];
	      P[l*k+9] = px->second[0]*px->second[2];
	    }
	  }
	  l++;
	}
	//////////////////////////////////////////////////
	// Definition of the weight W
 
	double W[nNode], Wx[nNode], Wy[nNode], Wz[nNode];
	
	
	l=0;
	for (px=x.begin(); px!=x.end(); px++) {
	  // distance between mPt and nodes
	  double dx=0;
	  for (i=0;i<dim;i++){
	    dx+= (y[i] - px->second[i])*(y[i] - px->second[i]);
	  }
		
	  double q = beta*dx;
	  W[l]=exp(-q);
	  Wx[l]=-2*beta*(y[0]-px->second[0])*W[l];
	  if (dim >= 2) {Wy[l]=-2*beta*(y[1]-px->second[1])*W[l];}
	  if (dim >= 3) {Wz[l]=-2*beta*(y[2]-px->second[2])*W[l];}
		
		
	  /*double q = sqrt(dx)/rho;
		
	    if (q > 1.){
	    W[l]=0.;
	    Wx[l]=0.;
	    Wy[l]=0.;
	    Wz[l]=0.;
	    }	
	    else {
	    W[l]= 1-6*q*q+8*q*q*q-3*q*q*q*q;
	    Wx[l]= (-12+24*q-12*q*q)*(y[0]-px->second[0])/(rho*rho);
	    if (dim >= 2) {Wy[l]= (-12+24*q-12*q*q)*(y[1]-px->second[1])/(rho*rho);}
	    if (dim >= 3) {Wz[l]= (-12+24*q-12*q*q)*(y[2]-px->second[2])/(rho*rho);}
	    }
	  */
	  l++;
	}

	/////////////////////////////////////////////////////////
	// Definition of B, Bx, By
	
	
	
	double* B0 = new double[k*nNode];	
	for (i=0;i<k;i++) {
	  for (j=0;j<nNode;j++) {
	    B0[i*nNode+j] = P[j*k+i]*W[j];
	  }
	}
	
 	/////////////////////////////////////////////////////////
	// Definition of A, Ax, Ay
	
	Set::VectorSpace::Hom A(k,k);
	double Ax[k][k], Ay[k][k], Az[k][k];
 
	// Initialization of A
 	for (i=0;i<k;i++) {
	  for (j=0;j<k;j++) {
	    A[j][i] = 0;
	    Ax[i][j] = 0;
	    if (dim >= 2) {Ay[i][j] = 0;}
	    if (dim >= 3) {Az[i][j] = 0;}
	  }
	}

	for (i=0;i<k;i++) {
	  for (j=0;j<k;j++) {
	    for (l=0;l<nNode;l++){
	      A[j][i] += W[l]*P[l*k+i]*P[l*k+j];
	      Ax[i][j] += Wx[l]*P[l*k+i]*P[l*k+j];
	      if (dim >= 2) {Ay[i][j] += Wy[l]*P[l*k+i]*P[l*k+j];}
	      if (dim >= 3) {Az[i][j] += Wz[l]*P[l*k+i]*P[l*k+j];}
	    }
	  }
	}	
	
	delete[] P;
 	/////////////////////////////////////////////////////////
	// Definition of AI
	
	Set::VectorSpace::Hom AI(k,k);

	/*for (i=0;i<k;i++) {
	  for (j=0;j<k;j++){
	  std::cout << A[j][i] << " ";
	  }
	  std::cout << "; " << std::endl;
	  }
	*/
	
	AI = Inverse(A);       
	
	//////////////////////////////////////////////////////////////////////
	//		Computation of the derivatives
	
	/////////////////////////////////////////////////////////
	// 1st term
	/////////////////////////////////////////////////////////

	double alpha[k];
	
	for (i=0;i<k;i++) {
	  alpha[i]=0;
	  for (j=0;j<k;j++){
	    alpha[i]+=AI[j][i]*Px[j];
	  }
	} 

	double alphaX[k],alphaY[k],alphaZ[k];
	double Axalpha[k],Ayalpha[k],Azalpha[k];

	for (i=0;i<k;i++) {
	  Axalpha[i]=0;Ayalpha[i]=0;Azalpha[i]=0;
	  for (j=0;j<k;j++){
	    if (dim >= 1) Axalpha[i] += Ax[i][j]*alpha[j];
	    if (dim >= 2) Ayalpha[i] += Ay[i][j]*alpha[j];
	    if (dim >= 3) Azalpha[i] += Az[i][j]*alpha[j];
	  }
	}

	for (i=0;i<k;i++) {
	  alphaX[i]=0;alphaY[i]=0;alphaZ[i]=0;
	  for (j=0;j<k;j++){
	    if (dim >= 1) alphaX[i] += AI[j][i]*(DPx[j]-Axalpha[j]);
	    if (dim >= 2) alphaY[i] += AI[j][i]*(DPy[j]-Ayalpha[j]);
	    if (dim >= 3) alphaZ[i] += AI[j][i]*(DPz[j]-Azalpha[j]);
	  }
	}	
	
	
	double* D = new double[dim*nNode];
	double DP[k];
	
	l=0;
	for (px=x.begin(); px!=x.end(); px++){
	
	  //computation of DP
	  DP[0]=1;
	  if (dim == 1){
	    DP[1]=px->second[0];
	    if (m>=2) DP[2]=(px->second[0])*(px->second[0]);
	    if (m>=3) DP[3]=(px->second[0])*(px->second[0])*(px->second[0]);		
	  }
	  if (dim == 2){
	    DP[1]=px->second[0];
	    DP[2]=px->second[1];
	    if (m>=2) {
	      DP[3]=(px->second[0])*(px->second[0]);
	      DP[4]=(px->second[0])*(px->second[1]);
	      DP[5]=(px->second[1])*(px->second[1]);
	    }
	    if (m>=3) {
	      DP[6]=(px->second[0])*(px->second[0])*(px->second[0]);
	      DP[7]=(px->second[0])*(px->second[0])*(px->second[1]);
	      DP[8]=(px->second[0])*(px->second[1])*(px->second[1]);
	      DP[9]=(px->second[1])*(px->second[1])*(px->second[1]);
	    }	
	  }
	  if (dim == 3){
	    DP[1]=px->second[0];
	    DP[2]=px->second[1];
	    DP[3]=px->second[2];
	    if (m>=2) {
	      DP[4]=(px->second[0])*(px->second[0]);
	      DP[5]=(px->second[0])*(px->second[1]);
	      DP[6]=(px->second[1])*(px->second[1]);
	      DP[7]=(px->second[1])*(px->second[2]);
	      DP[8]=(px->second[2])*(px->second[2]);
	      DP[9]=(px->second[0])*(px->second[2]);
	    }
	  }
		
	  double alphaDP=0,alphaxDP=0,alphayDP=0,alphazDP=0;
	  for (i=0;i<k;i++) {
	    alphaDP += alpha[i]*DP[i];
	    alphaxDP += alphaX[i]*DP[i];
	    if (dim>=2) alphayDP += alphaY[i]*DP[i];
	    if (dim>=3) alphazDP += alphaZ[i]*DP[i];
	  }
		
	  D[l] = alphaDP*Wx[l] + alphaxDP*W[l];
	  if (dim>=2) D[nNode+l] = alphaDP*Wy[l] + alphayDP*W[l];
	  if (dim>=3) D[2*nNode+l] = alphaDP*Wz[l] + alphazDP*W[l];
				
	  l++;
	}
	
	/////////////////////////////////////////////////////////
	// Build of DN
	// loop on nodes
	Set::VectorSpace::Vector Dv(dim);
	j=0;
	for (px=x.begin(); px!=x.end(); px++)
	  {
	    for (i=0;i<dim;i++){
	      Dv[i]=D[i*nNode+j];
	    }
	    DN.insert(make_pair(px->first,Dv));
	    j++;
	  }
	delete[] D;
		
	

	// Definition of C= AI*B 
	double* C1 = new double[k*nNode];
	//double C1[k][nNode];

	for (i=0;i<k;i++) {
	  for (j=0;j<nNode;j++) {
	    C1[i*nNode+j]=0;
	    for (l=0;l<k;l++) {
	      C1[i*nNode+j] += AI[l][i]*B0[l*nNode+j];
	    }
	  }
	}

	delete[] B0;
	
	double D0[nNode];
 
	for (i=0;i<nNode;i++){
	  D0[i]=0;
	}
	for (i=0;i<nNode;i++){
	  for (j=0;j<k;j++){
	    D0[i]+=Px[j]*C1[j*nNode+i];
	  }
	}
	
	delete[] C1;
	/////////////////////////////////////////////////////////
	// Build of N
	
	
	// loop on nodes
	j=0;
	for (px=x.begin(); px!=x.end(); px++)
	  {
	    N.insert(make_pair(px->first,D0[j]));
	    j++;
	  }	


	
	//Element::Interpolation::MLS::Shape<0> S(B);
	//N = S(y);
		
	//Element::Interpolation::MLS::Shape<1> DS(B);
	//DN = DS(y);

	return;
      }

      void Jet < 0 > ::operator () (const domain_type & y, shape_type & N, shape_derivative_type & DN,
				    map<Set::Manifold::Point *, Set::VectorSpace::Vector> & AIPW,
				    map<Set::Manifold::Point *, Set::VectorSpace::Hom> & AIxPW,
				    map<Set::Manifold::Point *, Set::VectorSpace::Hom> & AIPWx) {
  
	if ( !N.empty() ) N.clear();
	if ( !DN.empty() ) DN.clear();
	
	AIPW.clear();
	AIxPW.clear();
	AIPWx.clear();


	int i,j,l;
	int ii,ij,ik,il,im;
 
	// dimension
	unsigned int dim = B->GetDim();
	
	// get nodes
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> & x = *(B->GetNodes());
 
	// get polynomial order
	int m = B->GetOrder();
 
	// get radius
	double beta = B->GetRadius();
	
	// 	dilatation parameter 
	double rho = sqrt(1/beta);
 
	// calculation of the number of term in the polynomial base
	int k = B->GetNumberOfTerms();

	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator px;
 
	// number of node
	int nNode = x.size();

	/////////////////////////////////////////////////////////
	// definition of polynomial basis 
	
	// Polynomial basis for the mPt
	double Px[k];
 
	// 1D
	if (dim == 1) {
	  Px[0] = 1;
	  for (i=1;i<k;i++) {
	    Px[i] = Px[i-1]*y[0];
	  }
	}
 
	// 2D
	if (dim == 2) {
	  double tabx[m+1];
	  double taby[m+1];
		
	  // define 1st line and first column
	  tabx[0] = 1;
	  taby[0] = 1;
	  for (i=1;i<m+1;i++) {
	    tabx[i] = tabx[i-1]*y[0];
	    taby[i] = taby[i-1]*y[1];
	  }
		
	  /*		for (i=1;i<m+1;i++) {
			for (j=1;j<m+1;j++)
			tab[j][i] = tab[j][0]*tab[0][i];
			}
		
			int inc = 0;
			for (i=0;i<m+1;i++) {
			for (j=0;j<i+1;j++) {
			Px[inc] = tab[i-j][j];
			inc++;
			}
			}*/
	  int inc = 0;
	  for (i=0;i<m+1;i++) {
	    for (j=0;j<i+1;j++) {
	      Px[inc] = tabx[i-j]*taby[j];
	      inc++;
	    }
	  }
		
	}	
	
	// 3D
	if (dim == 3) {
	  if (m==1){
	    Px[0] = 1;
	    Px[1] = y[0];
	    Px[2] = y[1];
	    Px[3] = y[2];
	  }
	  if (m == 2) {
	    Px[0] = 1;
	    Px[1] = y[0];
	    Px[2] = y[1];
	    Px[3] = y[2];
	    Px[4] = y[0]*y[0];
	    Px[5] = y[0]*y[1];
	    Px[6] = y[1]*y[1];
	    Px[7] = y[1]*y[2];
	    Px[8] = y[2]*y[2];
	    Px[9] = y[0]*y[2];
	  }
	}
	
	// Polynomial basis for the mPt

	double DPx[k], DPy[k], DPz[k];

	// 1D
	if (dim == 1) {
	  if (m==1){
	    DPx[0]=0;
	    DPx[1]=1;
	  }
	  if (m>=2){
	    DPx[0]=0;
	    DPx[1]=1;
	    //DPx[3]=3*y[0]*y[0];
	    for (i=2;i<m+1;i++){
	      DPx[i] = y[0]*DPx[i-1]*i/(i-1);
	    }
			
	  }	
	}
	
	// 2D
	if (dim == 2) {
	  // define polynomial table
	  if (m==1){
	    DPx[0]=0;
	    DPx[1]=1;
	    DPx[2]=0;
	    DPy[0]=0;
	    DPy[1]=0;
	    DPy[2]=1;
	  }
	  if (m>=2){
	    double tabx[m+1][m+1];
	    double taby[m+1][m+1];
	    // define 1st line and first column
	    tabx[0][0] = 0;
	    tabx[0][1] = 1;
	    tabx[1][0] = y[1];
	    taby[0][0] = 0;
	    taby[0][1] = y[0];
	    taby[1][0] = 1;
			
	    for (i=2;i<m+1;i++) {
	      tabx[0][i] = y[0]*tabx[0][i-1]*i/(i-1);
	      tabx[i][0] = y[1]*tabx[i-1][0];
	      taby[0][i] = y[0]*taby[0][i-1];
	      taby[i][0] = y[1]*taby[i-1][0]*i/(i-1);
	    }
	    for (i=1;i<m+1;i++) {
	      for (j=1;j<m+1;j++){
		tabx[i][j] = tabx[0][j]*tabx[i][0];
		taby[i][j] = taby[0][j]*taby[i][0];
	      }
	    }
	    int inc = 0;
	    for (i=0;i<m+1;i++) {
	      for (j=0;j<i+1;j++) {
		if (i==j){
		  DPx[inc] = 0;
		}
		else {
		  DPx[inc] = tabx[j][i-j];
		}
		if (j==0){
		  DPy[inc] = 0;
		}
		else {
		  DPy[inc] = taby[j][i-j];
		}
		inc++;
	      }
	    }
	  }		
	}

	// 3D
	if (dim == 3) {
	  // define polynomial table
	  if (m==1){
	    DPx[0]=0;
	    DPx[1]=1;
	    DPx[2]=0;
	    DPx[3]=0;
	    DPy[0]=0;
	    DPy[1]=0;
	    DPy[2]=1;
	    DPy[3]=0;
	    DPz[0]=0;
	    DPz[1]=0;
	    DPz[2]=0;
	    DPz[3]=1;			
	  }
	  if (m==2){
	    DPx[0] = 0;
	    DPx[1] = 1;
	    DPx[2] = 0;
	    DPx[3] = 0;
	    DPx[4] = 2*y[0];
	    DPx[5] = y[1];
	    DPx[6] = 0;
	    DPx[7] = 0;
	    DPx[8] = 0;
	    DPx[9] = y[2];
			
	    DPy[0] = 0;
	    DPy[1] = 0;
	    DPy[2] = 1;
	    DPy[3] = 0;
	    DPy[4] = 0;
	    DPy[5] = y[0];
	    DPy[6] = 2*y[1];
	    DPy[7] = y[2];
	    DPy[8] = 0;
	    DPy[9] = 0;
			
	    DPz[0] = 0;
	    DPz[1] = 0;
	    DPz[2] = 0;
	    DPz[3] = 1;
	    DPz[4] = 0;
	    DPz[5] = 0;
	    DPz[6] = 0;
	    DPz[7] = y[1];
	    DPz[8] = 2*y[2];
	    DPz[9] = y[0];			
	  }
	}
			
	// Polynomial basis for the nodes
	double* P = new double[nNode*k];
	l=0;
	for (px=x.begin(); px!=x.end(); px++) {
	  // 1D
	  if (dim == 1) {
			
	    P[l*k+0] = 1;
	    for (i=1;i<k;i++) {
	      P[l*k+i] = P[l*k+(i-1)]*px->second[0];
	    }
	  }
	  // 2D
	  if (dim == 2) {
			
	    double tab[m+1][m+1];
	    // define 1st line and first column
	    tab[0][0] = 1.;
	    for (i=1;i<m+1;i++) {
	      tab[0][i] = tab[0][i-1]*px->second[0];
	      tab[i][0] = tab[i-1][0]*px->second[1];
	    }
			
	    if (m>1){
	      for (i=1;i<m+1;i++) {
		for (j=1;j<m+1;j++) {
		  tab[i][j] = tab[0][j]*tab[i][0];
		}
	      }
	    }
			
	    int inc = 0;
	    for (i=0;i<m+1;i++) {
	      for (j=0;j<i+1;j++) {
		P[l*k+inc] = tab[j][i-j];
		inc++;
	      }
	    }
	  }
	  if (dim == 3) {
			
	    if (m==1){
	      P[l*k] = 1;
	      P[l*k+1] = px->second[0];
	      P[l*k+2] = px->second[1];
	      P[l*k+3] = px->second[2];
	    }
	    if (m == 2) {
	      P[l*k] = 1;
	      P[l*k+1] = px->second[0];
	      P[l*k+2] = px->second[1];
	      P[l*k+3] = px->second[2];
	      P[l*k+4] = px->second[0]*px->second[0];
	      P[l*k+5] = px->second[0]*px->second[1];
	      P[l*k+6] = px->second[1]*px->second[1];
	      P[l*k+7] = px->second[1]*px->second[2];
	      P[l*k+8] = px->second[2]*px->second[2];
	      P[l*k+9] = px->second[0]*px->second[2];
	    }
	  }
	  l++;
	}
	//////////////////////////////////////////////////
	// Definition of the weight W
 
	double W[nNode], Wx[nNode], Wy[nNode], Wz[nNode];
	
	
	l=0;
	for (px=x.begin(); px!=x.end(); px++) {
	  // distance between mPt and nodes
	  double dx=0;
	  for (i=0;i<dim;i++){
	    dx+= (y[i] - px->second[i])*(y[i] - px->second[i]);
	  }
		
	  double q = beta*dx;
	  W[l]=exp(-q);
	  Wx[l]=-2*beta*(y[0]-px->second[0])*W[l];
	  if (dim >= 2) {Wy[l]=-2*beta*(y[1]-px->second[1])*W[l];}
	  if (dim >= 3) {Wz[l]=-2*beta*(y[2]-px->second[2])*W[l];}
		
		
	  /*double q = sqrt(dx)/rho;
		
	    if (q > 1.){
	    W[l]=0.;
	    Wx[l]=0.;
	    Wy[l]=0.;
	    Wz[l]=0.;
	    }	
	    else {
	    W[l]= 1-6*q*q+8*q*q*q-3*q*q*q*q;
	    Wx[l]= (-12+24*q-12*q*q)*(y[0]-px->second[0])/(rho*rho);
	    if (dim >= 2) {Wy[l]= (-12+24*q-12*q*q)*(y[1]-px->second[1])/(rho*rho);}
	    if (dim >= 3) {Wz[l]= (-12+24*q-12*q*q)*(y[2]-px->second[2])/(rho*rho);}
	    }
	  */
	  l++;
	}

	/////////////////////////////////////////////////////////
	// Definition of B, Bx, By
	
	
	
	double* B0 = new double[k*nNode];	
	for (i=0;i<k;i++) {
	  for (j=0;j<nNode;j++) {
	    B0[i*nNode+j] = P[j*k+i]*W[j];
	  }
	}
	
	/////////////////////////////////////////////////////////
	// Definition of A, Ax, Ay
	
	Set::VectorSpace::Hom A(k,k);
	double Ax[k][k], Ay[k][k], Az[k][k];
 
	// Initialization of A
	for (i=0;i<k;i++) {
	  for (j=0;j<k;j++) {
	    A[j][i] = 0;
	    Ax[i][j] = 0;
	    if (dim >= 2) {Ay[i][j] = 0;}
	    if (dim >= 3) {Az[i][j] = 0;}
	  }
	}

	for (i=0;i<k;i++) {
	  for (j=0;j<k;j++) {
	    for (l=0;l<nNode;l++){
	      A[j][i] += W[l]*P[l*k+i]*P[l*k+j];
	      Ax[i][j] += Wx[l]*P[l*k+i]*P[l*k+j];
	      if (dim >= 2) {Ay[i][j] += Wy[l]*P[l*k+i]*P[l*k+j];}
	      if (dim >= 3) {Az[i][j] += Wz[l]*P[l*k+i]*P[l*k+j];}
	    }
	  }
	}	
	
	
	/////////////////////////////////////////////////////////
	// Definition of AI
	
	Set::VectorSpace::Hom AI(k,k);

	/*for (i=0;i<k;i++) {
	  for (j=0;j<k;j++){
	  std::cout << A[j][i] << " ";
	  }
	  std::cout << "; " << std::endl;
	  }
	*/
	
	AI = Inverse(A);
	
	
	j=0;
	for (px=x.begin(); px!=x.end(); px++)
	  {
	    Set::VectorSpace::Hom AIxPWloc(dim, k);
	    //AIxPWloc=-AI.Ax.AI.B
	    for (ii=0;ii<k;ii++){
	      for (ij=0;ij<k;ij++){
		for (ik=0;ik<k;ik++){
		  for (il=0;il<k;il++){
		      AIxPWloc[ii][0] -= AI[ij][ii]*Ax[ij][ik]*AI[il][ik]*B0[il*nNode+j];
		      if (dim > 1) AIxPWloc[ii][1] -= AI[ij][ii]*Ay[ij][ik]*AI[il][ik]*B0[il*nNode+j];
		      if (dim > 2) AIxPWloc[ii][2] -= AI[ij][ii]*Az[ij][ik]*AI[il][ik]*B0[il*nNode+j];
		  }
		}
	      }	
	    }
	    // map<Set::Manifold::Point *, Set::VectorSpace::Hom (k,dim) > AIxPW
	    AIxPW.insert(make_pair(px->first,AIxPWloc));
	    j++;
	  }			
		
	j=0;
	for (px=x.begin(); px!=x.end(); px++)
	  {
	    Set::VectorSpace::Hom AIPWxloc(dim, k);
	    for (ii=0;ii<k;ii++){
	      for (ij=0;ij<k;ij++){
		AIPWxloc[ii][0] += AI[ij][ii]*P[j*k+ij]*Wx[j];
		if (dim > 1) AIPWxloc[ii][1] += AI[ij][ii]*P[j*k+ij]*Wy[j];
		if (dim > 2) AIPWxloc[ii][2] += AI[ij][ii]*P[j*k+ij]*Wz[j];
	      }
	    }
	    // map<Set::Manifold::Point *, Set::VectorSpace::Hom (k,dim) > AIPWx
	    AIPWx.insert(make_pair(px->first,AIPWxloc));
	    j++;
	  }
	delete[] P;
	//////////////////////////////////////////////////////////////////////
	//		Computation of the derivatives
	
	/////////////////////////////////////////////////////////
	// 1st term
	/////////////////////////////////////////////////////////

	double alpha[k];
	
	for (i=0;i<k;i++) {
	  alpha[i]=0;
	  for (j=0;j<k;j++){
	    alpha[i]+=AI[j][i]*Px[j];
	  }
	} 

	double alphaX[k],alphaY[k],alphaZ[k];
	double Axalpha[k],Ayalpha[k],Azalpha[k];

	for (i=0;i<k;i++) {
	  Axalpha[i]=0;Ayalpha[i]=0;Azalpha[i]=0;
	  for (j=0;j<k;j++){
	    if (dim >= 1) Axalpha[i] += Ax[i][j]*alpha[j];
	    if (dim >= 2) Ayalpha[i] += Ay[i][j]*alpha[j];
	    if (dim >= 3) Azalpha[i] += Az[i][j]*alpha[j];
	  }
	}

	for (i=0;i<k;i++) {
	  alphaX[i]=0;alphaY[i]=0;alphaZ[i]=0;
	  for (j=0;j<k;j++){
	    if (dim >= 1) alphaX[i] += AI[j][i]*(DPx[j]-Axalpha[j]);
	    if (dim >= 2) alphaY[i] += AI[j][i]*(DPy[j]-Ayalpha[j]);
	    if (dim >= 3) alphaZ[i] += AI[j][i]*(DPz[j]-Azalpha[j]);
	  }
	}	
	
	
	double* D = new double[dim*nNode];
	double DP[k];
	
	l=0;
	for (px=x.begin(); px!=x.end(); px++){
	
	  //computation of DP
	  DP[0]=1;
	  if (dim == 1){
	    DP[1]=px->second[0];
	    if (m>=2) DP[2]=(px->second[0])*(px->second[0]);
	    if (m>=3) DP[3]=(px->second[0])*(px->second[0])*(px->second[0]);		
	  }
	  if (dim == 2){
	    DP[1]=px->second[0];
	    DP[2]=px->second[1];
	    if (m>=2) {
	      DP[3]=(px->second[0])*(px->second[0]);
	      DP[4]=(px->second[0])*(px->second[1]);
	      DP[5]=(px->second[1])*(px->second[1]);
	    }
	    if (m>=3) {
	      DP[6]=(px->second[0])*(px->second[0])*(px->second[0]);
	      DP[7]=(px->second[0])*(px->second[0])*(px->second[1]);
	      DP[8]=(px->second[0])*(px->second[1])*(px->second[1]);
	      DP[9]=(px->second[1])*(px->second[1])*(px->second[1]);
	    }	
	  }
	  if (dim == 3){
	    DP[1]=px->second[0];
	    DP[2]=px->second[1];
	    DP[3]=px->second[2];
	    if (m>=2) {
	      DP[4]=(px->second[0])*(px->second[0]);
	      DP[5]=(px->second[0])*(px->second[1]);
	      DP[6]=(px->second[1])*(px->second[1]);
	      DP[7]=(px->second[1])*(px->second[2]);
	      DP[8]=(px->second[2])*(px->second[2]);
	      DP[9]=(px->second[0])*(px->second[2]);
	    }
	  }
		
	  double alphaDP=0,alphaxDP=0,alphayDP=0,alphazDP=0;
	  for (i=0;i<k;i++) {
	    alphaDP += alpha[i]*DP[i];
	    alphaxDP += alphaX[i]*DP[i];
	    if (dim>=2) alphayDP += alphaY[i]*DP[i];
	    if (dim>=3) alphazDP += alphaZ[i]*DP[i];
	  }
		
	  D[l] = alphaDP*Wx[l] + alphaxDP*W[l];
	  if (dim>=2) D[nNode+l] = alphaDP*Wy[l] + alphayDP*W[l];
	  if (dim>=3) D[2*nNode+l] = alphaDP*Wz[l] + alphazDP*W[l];
				
	  l++;
	}
	
	/////////////////////////////////////////////////////////
	// Build of DN
	// loop on nodes
	Set::VectorSpace::Vector Dv(dim);
	j=0;
	for (px=x.begin(); px!=x.end(); px++)
	  {
	    for (i=0;i<dim;i++){
	      Dv[i]=D[i*nNode+j];
	    }
	    DN.insert(make_pair(px->first,Dv));
	    j++;
	  }
	delete[] D;
		
	

	// Definition of C= AI*B 
	double* C1 = new double[k*nNode];
	//double C1[k][nNode];

	for (i=0;i<k;i++) {
	  for (j=0;j<nNode;j++) {
	    C1[i*nNode+j]=0;
	    for (l=0;l<k;l++) {
	      C1[i*nNode+j] += AI[l][i]*B0[l*nNode+j];
	    }
	  }
	}
	
	j=0;
	for (px=x.begin(); px!=x.end(); px++)
	  {
	    Set::VectorSpace::Vector AIPWloc(k);
	    for (i=0;i<k;i++){
	      AIPWloc[i]=C1[i*nNode+j];
	    }
	    // map<Set::Manifold::Point *, Set::VectorSpace::Vector (k) > AIPW
	    AIPW.insert(make_pair(px->first,AIPWloc));
	    j++;
	  }	

	delete[] B0;
	
	double D0[nNode];
 
	for (i=0;i<nNode;i++){
	  D0[i]=0;
	}
	for (i=0;i<nNode;i++){
	  for (j=0;j<k;j++){
	    D0[i]+=Px[j]*C1[j*nNode+i];
	  }
	}
	
	delete[] C1;
	/////////////////////////////////////////////////////////
	// Build of N
	
	
	// loop on nodes
	j=0;
	for (px=x.begin(); px!=x.end(); px++)
	  {
	    N.insert(make_pair(px->first,D0[j]));
	    j++;
	  }	


	
	//Element::Interpolation::MLS::Shape<0> S(B);
	//N = S(y);
		
	//Element::Interpolation::MLS::Shape<1> DS(B);
	//DN = DS(y);

	return;
      }

    }

  }

}
