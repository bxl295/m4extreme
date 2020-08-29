// EigenSym.cpp: Implementation of the EigenSym class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <cassert>
#include "EigenSym.h"

namespace LinearAlgebra  
{
EigenSym::EigenSym() {}

EigenSym::~EigenSym() {}

void 
EigenSym::operator () (
		const Set::VectorSpace::Sym &E,
		Set::VectorSpace::Diagonal &Lambda,
		Set::VectorSpace::Hom &V)
{
	unsigned int i, j, n = E.size1();
	assert(n == Lambda.size());
	assert(n == V.size1());
	Set::VectorSpace::Hom A = E.Embed();
	TNT::Array1D<double> TNT_real(n);
	TNT::Array2D<double> TNT_A(n,n);
	TNT::Array2D<double> TNT_V(n,n);
	for (i=0; i<n; i++)
		for (j=0; j<n; j++)
			TNT_A[i][j] = A[j][i];
	JAMA::Eigenvalue<double> JAMA_Eig(TNT_A);
	JAMA_Eig.getRealEigenvalues(TNT_real);
	JAMA_Eig.getV(TNT_V);
	for (i=0; i<n; i++) 
		Lambda[i] = TNT_real[i];
	for (i=0; i<n; i++)
		for (j=0; j<n; j++)
			V[j][i] = TNT_V[i][j];
}

}