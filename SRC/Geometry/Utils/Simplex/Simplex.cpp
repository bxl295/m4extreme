// Simplex.cpp: Implementation for the Simplex class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include "./Simplex.h"

namespace Geometry
{
namespace Utils
{
Simplex::Simplex(){}

Simplex::Simplex(
	const vector<Set::VectorSpace::Vector> &rhs_x) 
	: n((unsigned int)rhs_x.size()-1), npp((unsigned int)rhs_x.size()), x(rhs_x) 
{
	assert(n > 0);
}

Simplex::~Simplex(){}

double 
Simplex::CircumRadius()
{
	Set::VectorSpace::Vector xc = CircumCenter();
	Set::VectorSpace::Vector u = x[0] - xc;
	return sqrt(u(u));
}

Set::VectorSpace::Vector 
Simplex::CircumCenter()
{
	unsigned int i, a;
	Set::VectorSpace::Hom A(npp);
	for (a=0; a<npp; a++)
	{
		for (i=0; i<n; i++) A[a][i] = x[a][i]; 
		A[a][n] = 1.0;
	}
	Set::VectorSpace::Vector b(npp);
	for (a=0; a<npp; a++)
	{
		Set::VectorSpace::Vector u(x[a]);
		b[a] = 0.5*u(u);
	}
	Set::VectorSpace::Vector c = Adjoint(Inverse(A))*b;
	Set::VectorSpace::Vector xc(n);
	for (i=0; i<n; i++) xc[i] = c[i];
	return xc;
}

void 
Simplex::CircumSphere(
	Set::VectorSpace::Vector &xc, 
	double &r)
{
	xc = CircumCenter();
	Set::VectorSpace::Vector u = x[0] - xc;
	r = sqrt(u(u));
}

double 
Simplex::Volume()
{
	unsigned int i, a;
	Set::VectorSpace::Hom A(npp);
	for (a=0; a<npp; a++)
	{
		for (i=0; i<n; i++) A[a][i] = x[a][i]; 
		A[a][n] = 1.0;
	}
	return fabs(Jacobian(A))/(double)Factorial(n);
}

Set::VectorSpace::Vector 
Simplex::BaryCenter()
{
	unsigned int a;
	Set::VectorSpace::Vector O(n);
	Set::VectorSpace::Vector A(n);
	for (a=0; a<npp; a++)
	{
		Set::VectorSpace::Vector B(x[a]);
		A += B;
	}
	return O + A/(double)npp;
}

vector<double>
Simplex::Lambda(const point_type &P)
{
	unsigned int i, a;
	vector<double> N(npp,0.0);
	Set::VectorSpace::Hom A(npp);
	for (a=0; a<npp; a++)
	{
		for (i=0; i<n; i++) A[a][i] = x[a][i]; 
		A[a][n] = 1.0;
	}
	double J = Jacobian(A);
	for (a=0; a<npp; a++)
	{
		Set::VectorSpace::Hom B=A;
		for (i=0; i<n; i++) B[a][i] = P[i]; 
		N[a] = Jacobian(B)/J;
	}

	return N;
}

vector<Set::VectorSpace::Vector>
Simplex::DLambda()
{
	unsigned int i, a;
	Set::VectorSpace::Vector O(n);
	vector<Set::VectorSpace::Vector> DN(npp,O);
	Set::VectorSpace::Vector P(n), Q(n);
	vector<double> NP = Lambda(P);
	for (i=0; i<n; i++)
	{
		Q = P; Q[i]=1.0; 
		vector<double> NQ = Lambda(Q);
		for (a=0; a<npp; a++) DN[a][i] = NQ[a] - NP[a];
	}
	return DN;
}

unsigned int 
Simplex::Factorial(const unsigned int &k)
{
	unsigned int i, fact=1;
	if (k <= 1) return fact;
	for (i=2; i<=k; i++) fact *= i;
	return fact;
}

}

}

