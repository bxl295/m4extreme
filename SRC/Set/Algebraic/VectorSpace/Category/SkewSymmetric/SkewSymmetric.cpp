// SkewSkwmetric.cpp: Implementation of the Skw class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <cassert>
#include "./SkewSymmetric.h"
#include "../../../../../Utils/Indexing/Indexing.h"
#include "tnt_array1d.h"
#include "tnt_array2d.h"
#include "jama_eig.h"

namespace Set
{
namespace VectorSpace
{
//////////////////////////////////////////////////////////////////////
// Class Skw
//////////////////////////////////////////////////////////////////////

Skw::Skw() : Vector(), n1(0){}

Skw::Skw(const unsigned int &m1) : Vector(m1*(m1-1)/2), n1(m1){}

Skw::Skw(const unsigned int &m1, double * const &u) 
: Vector(m1*(m1-1)/2,u), n1(m1){}

Skw::~Skw(){}

Skw::Skw(const Vector &A) : Vector(A)
{
	unsigned int n = A.size();
	float x = (float)(1 + 8*n); 
	float y = sqrt(x);
	n1 = (1 + (unsigned int)y)/2;
	assert (n == n1*(n1-1)/2);
}

Skw::Skw(const Skw &A) : Vector(A), n1(A.n1){}

Skw::Skw(const Hom &A) 
: Vector(A.size1()*(A.size1()-1)/2), n1(A.size1())
{
	assert (A.size1() == A.size2());
	unsigned int i, j; 
	double *p, *q;
	Hom B=A;
	for (j=0, p=this->head; j<n1-1; j++)
		for (i=j+1, q=B.begin()+j*n1+i; i<n1; i++, p++, q++) *p = -(*q);
	Hom C=Adjoint(A);
	for (j=0, p=this->head; j<n1-1; j++)
		for (i=j+1, q=C.begin()+j*n1+i; i<n1; i++, p++, q++) *p += *q;
	for (p=this->head; p<this->tail; p++) *p *= 0.5;
}

Skw & Skw::operator = (const Vector &A)
{
	assert(Vector::n == A.size()); 
	if (this == &A) return *this;
	Vector::operator = (A); 
	return *this;
}

Skw & Skw::operator = (const Skw &A)
{
	assert(n1 == A.n1); 
	if (this == &A) return *this;
	Vector::operator = (A);
	return *this;
}

Skw & Skw::operator = (const Hom &A)
{
	assert (n1 = A.size1());
	assert (n1 = A.size2());
	unsigned int i, j; 
	double *p, *q; Hom B=A;
	for (j=0, p=this->head; j<n1-1; j++)
		for (i=j+1, q=B.begin()+j*n1+i; i<n1; i++, p++, q++) *p = -(*q);
	Hom C=Adjoint(A);
	for (j=0, p=this->head; j<n1-1; j++)
		for (i=j+1, q=C.begin()+j*n1+i; i<n1; i++, p++, q++) *p += *q;
	for (p=this->head; p<this->tail; p++) *p *= 0.5;
	return *this;
}

double Skw::operator () (const Vector &A) const
{
	assert (Vector::n == A.size());
	Hom B = Embed();
	Skw C(A);
	return B(C.Embed());
}

double Skw::operator () (const Skw &A) const
{
	assert (n1 == A.n1);
	Hom B = Embed();
	Hom C = A.Embed();
	return B(C);
}

double Skw::operator () (const Hom &A) const
{
	assert (n1 == A.size1());
	assert (n1 == A.size2());
	Hom B = Embed();
	return B(A);
}

Hom Skw::Embed() const
{
	unsigned int i, j; 
	double *p, *q;
	Hom A(n1);
	for (j=0, p=this->head; j<n1-1; j++)
		for (i=j+1, q=A.begin()+j*n1+i; i<n1; i++, p++, q++) *q = *p;
	Hom B=Adjoint(A);
	for (j=0, p=this->head; j<n1-1; j++)
		for (i=j+1, q=B.begin()+j*n1+i; i<n1; i++, p++, q++) *q = -(*p);
	return B;
}

unsigned int Skw::size1() const{return n1;}

unsigned int Skw::size2() const{return n1;}	

//////////////////////////////////////////////////////////////////////
// Class SkwZero
//////////////////////////////////////////////////////////////////////

SkwZero::SkwZero() : VectorZero(), n1(0){}

SkwZero::SkwZero(const unsigned int &m1) 
: VectorZero(m1*(m1-1)/2), n1(m1){}

SkwZero::SkwZero(const unsigned int &m1, double * const &u) 
: VectorZero(m1*(m1-1)/2,u), n1(m1){}

SkwZero::~SkwZero(){}

SkwZero::SkwZero(const VectorZero &A) : VectorZero(A)
{
	unsigned int n = A.size();
	float x = (float)(1 + 8*n); 
	float y = sqrt(x);
	n1 = (1 + (unsigned int)y)/2;
	assert (n == n1*(n1-1)/2);
}

SkwZero::SkwZero(const SkwZero &A) : VectorZero(A), n1(A.n1){}

SkwZero::SkwZero(const HomZero &A) 
: VectorZero(A.size1()*(A.size1()-1)/2), n1(A.size1()){}

SkwZero & SkwZero::operator = (const VectorZero &A)
{
	assert(Vector::n == A.size()); 
	if (this == &A) return *this;
	return *this;
}

SkwZero & SkwZero::operator = (const SkwZero &A)
{
	assert(n1 == A.n1); 
	return *this;
}

SkwZero & SkwZero::operator = (const HomZero &A)
{
	assert (n1 = A.size1());
	assert (n1 = A.size2());
	return *this;
}

double SkwZero::operator () (const Vector &A) const
{
	assert (Vector::n == A.size());
	return 0.0;
}

double SkwZero::operator () (const Skw &A) const
{
	assert (n1 == A.size1());
	return 0.0;
}

double SkwZero::operator () (const Hom &A) const
{
	assert (n1 == A.size1());
	assert (n1 == A.size2());
	return 0.0;
}

HomZero SkwZero::Embed() const
{
	HomZero A(n1);
	return A;
}

unsigned int SkwZero::size1() const{return n1;}

unsigned int SkwZero::size2() const{return n1;}	

namespace SkewSymmetric
{
//////////////////////////////////////////////////////////////////////
// Class Exp and Log classes
//////////////////////////////////////////////////////////////////////

void Copy(
	const Set::VectorSpace::Hom &A,
	TNT::Array2D<double> &B)
{
	unsigned int i1, n1 = A.size1(); int j1;
	unsigned int i2, n2 = A.size2(); int j2;
	for (i2=0, j2=0; i2<n2; i2++, j2++)
		for (i1=0, j1=0; i1<n1; i1++, j1++)
			B[j1][j2] = A[i2][i1];
}

void Copy(
	const TNT::Array2D<double> &B, 
	Set::VectorSpace::Hom &A)
{
	unsigned int i1, n1 = A.size1(); int j1;
	unsigned int i2, n2 = A.size2(); int j2;
	for (i2=0, j2=0; i2<n2; i2++, j2++)
		for (i1=0, j1=0; i1<n1; i1++, j1++)
			A[i2][i1] = B[j1][j2];
}

Set::VectorSpace::Hom
DiagExpEven(const Set::VectorSpace::Hom &A)
{
	unsigned int n = A.size1(); 
	Set::VectorSpace::Hom B(n);
	unsigned int i, ipp;
	for (i=0; i<n; i+=2)
	{
		ipp = i+1;
		double u = A[i  ][i  ];
		double v = A[ipp][i  ];
		double r=exp(u);
		double a = r*cos(v);
		double b = r*sin(v);
		B[i  ][i  ] =  a;
		B[ipp][i  ] =  b;
		B[i  ][ipp] = -b;
		B[ipp][ipp] =  a;
	}
	return B;
}

Set::VectorSpace::Hom
DiagExp(const Set::VectorSpace::Hom &A)
{
	unsigned int n = A.size1(); 
	unsigned int nh = n/2;
	unsigned int m = 2*nh;
	if (m == n) return DiagExpEven(A);
	unsigned int nmm = n-1;
	unsigned int i, j, jpp, k, kpp, l;
	double Tol = 1.0e-12;
	l = nmm;
	for (i=0; i<nmm; i+=2)
		if (fabs(A[i][i+1]) < Tol) {l=i; break;}
	Set::VectorSpace::Hom B(nmm); 
	if (l > 0) 
		for (j=0; j<l; j+=2)
		{
			jpp = j+1;
			B[j  ][j  ] = A[j  ][j  ];
			B[jpp][j  ] = A[jpp][j  ];
			B[j  ][jpp] = A[j  ][jpp];
			B[jpp][jpp] = A[jpp][jpp];
		}
	if (l < nmm) 
		for (j=l; j<nmm; j+=2)
		{
			k = j+1;
			kpp = k+1;
			jpp = j+1;
			B[j  ][j  ] = A[k  ][k  ];
			B[jpp][j  ] = A[kpp][k  ];
			B[j  ][jpp] = A[k  ][kpp];
			B[jpp][jpp] = A[kpp][kpp];
		}
	Set::VectorSpace::Hom C = DiagExpEven(B);
	Set::VectorSpace::Hom D(n); 
	D[l][l] = 1.0;
	if (l > 0) 
		for (j=0; j<l; j+=2)
		{
			jpp = j+1;
			D[j  ][j  ] = C[j  ][j  ];
			D[jpp][j  ] = C[jpp][j  ];
			D[j  ][jpp] = C[j  ][jpp];
			D[jpp][jpp] = C[jpp][jpp];
		}
	if (l < nmm) 
		for (j=l; j<nmm; j+=2)
		{
			k = j+1;
			kpp = k+1;
			jpp = j+1;
			D[k  ][k  ] = C[j  ][j  ];
			D[kpp][k  ] = C[jpp][j  ];
			D[k  ][kpp] = C[j  ][jpp];
			D[kpp][kpp] = C[jpp][jpp];
		}
	return D;
}

Set::VectorSpace::Hom
DiagLogEven(const Set::VectorSpace::Hom &A)
{
	unsigned int n = A.size1(); 
	Set::VectorSpace::Hom B(n);
	unsigned int i, ipp;
	for (i=0; i<n; i+=2)
	{
		ipp = i+1;
		double a, b, x, y, t, r;
		x = A[i  ][i  ];
		y = A[ipp][i  ];
		r = sqrt(x*x + y*y);
		if (y >= 0.0) {t = acos(x/r);}
		else {t = - acos(x/r);}
		a = log(r); b = t;
		B[i  ][i  ] =  a;
		B[ipp][i  ] =  b;
		B[i  ][ipp] = -b;
		B[ipp][ipp] =  a;
	}
	return B;
}

Set::VectorSpace::Hom
DiagLog(const Set::VectorSpace::Hom &A)
{
	unsigned int n = A.size1(); 
	unsigned int nh = n/2;
	unsigned int m = 2*nh;
	if (m == n) return DiagLogEven(A);
	unsigned int nmm = n-1;
	unsigned int i, j, jpp, k, kpp, l;
	double Tol = 1.0e-12;
	l = nmm;
	for (i=0; i<nmm; i+=2)
		if (fabs(A[i][i+1]) < Tol) {l=i; break;}
	Set::VectorSpace::Hom B(nmm); 
	if (l > 0) 
		for (j=0; j<l; j+=2)
		{
			jpp = j+1;
			B[j  ][j  ] = A[j  ][j  ];
			B[jpp][j  ] = A[jpp][j  ];
			B[j  ][jpp] = A[j  ][jpp];
			B[jpp][jpp] = A[jpp][jpp];
		}
	if (l < nmm) 
		for (j=l; j<nmm; j+=2)
		{
			k = j+1;
			kpp = k+1;
			jpp = j+1;
			B[j  ][j  ] = A[k  ][k  ];
			B[jpp][j  ] = A[kpp][k  ];
			B[j  ][jpp] = A[k  ][kpp];
			B[jpp][jpp] = A[kpp][kpp];
		}
	Set::VectorSpace::Hom C = DiagLogEven(B);
	Set::VectorSpace::Hom D(n); 
	D[l][l] = 1.0;
	if (l > 0) 
		for (j=0; j<l; j+=2)
		{
			jpp = j+1;
			D[j  ][j  ] = C[j  ][j  ];
			D[jpp][j  ] = C[jpp][j  ];
			D[j  ][jpp] = C[j  ][jpp];
			D[jpp][jpp] = C[jpp][jpp];
		}
	if (l < nmm) 
		for (j=l; j<nmm; j+=2)
		{
			k = j+1;
			kpp = k+1;
			jpp = j+1;
			D[k  ][k  ] = C[j  ][j  ];
			D[kpp][k  ] = C[jpp][j  ];
			D[k  ][kpp] = C[j  ][jpp];
			D[kpp][kpp] = C[jpp][jpp];
		}
	return D;
}

//////////////////////////////////////////////////////////////////////
// Class Exp
//////////////////////////////////////////////////////////////////////

Exp::Exp() {}

Exp::~Exp() {}

Exp::Exp(const Exp &rhs)  {}

Exp & 
Exp::operator = (const Exp &rhs ) 
{
	return *this;
}

Set::VectorSpace::Hom 
Exp::operator () (const Set::VectorSpace::Skw &W)
{
	unsigned int n = W.size1();
	TNT::Array2D<double> A(n,n); Copy(W.Embed(),A);
	JAMA::Eigenvalue<double> Eigen(A);
	TNT::Array2D<double> D(n,n); Eigen.getD(D); 
	TNT::Array2D<double> V(n,n); Eigen.getV(V);
	Set::VectorSpace::Hom C(n); Copy(D,C);
	Set::VectorSpace::Hom U(n); Copy(V,U);
	return U*DiagExp(C)*Inverse(U);
}

//////////////////////////////////////////////////////////////////////
// Class Log
//////////////////////////////////////////////////////////////////////

Log::Log() {}

Log::~Log() {}

Log::Log(const Log &rhs)  {}

Log & 
Log::operator = (const Log &rhs ) 
{
	return *this;
}

Set::VectorSpace::Skw 
Log::operator () (const Set::VectorSpace::Hom &R)
{
	unsigned int n = R.size1();
	TNT::Array2D<double> A(n,n); Copy(R,A);
	JAMA::Eigenvalue<double> Eigen(A);
	TNT::Array2D<double> D(n,n); Eigen.getD(D); 
	TNT::Array2D<double> V(n,n); Eigen.getV(V);
	Set::VectorSpace::Hom C(n); Copy(D,C);
	Set::VectorSpace::Hom U(n); Copy(V,U);
	return U*DiagLog(C)*Inverse(U);
}

}

}

}
