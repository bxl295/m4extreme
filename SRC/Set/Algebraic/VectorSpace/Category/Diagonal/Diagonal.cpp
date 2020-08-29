// Diagonal.cpp: implementation for the Diagonal class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <cassert>

#include "./Diagonal.h"

namespace Set
{
namespace VectorSpace
{
//////////////////////////////////////////////////////////////////////
// Class Diagonal
//////////////////////////////////////////////////////////////////////

Diagonal::Diagonal() {}
	
Diagonal::Diagonal(const unsigned int &n) :
	Set::VectorSpace::Vector(n) {}
	
Diagonal::Diagonal(const unsigned int &n, double * const &u) :	
	Set::VectorSpace::Vector(n,u) {}

Diagonal::~Diagonal() {}
	
Diagonal::Diagonal(const Diagonal &A) :
	Set::VectorSpace::Vector(A) {}

Diagonal::Diagonal(const Set::VectorSpace::Vector &A) :
	Set::VectorSpace::Vector(A) {}

Diagonal::Diagonal(const Set::VectorSpace::Hom &A)
{
	unsigned int i; double *p, *q;
	for (i=0, p=this->head, q=A.begin(); i<n; p++, q+=n+1, i++) *p = *q;
}

Diagonal & 
Diagonal::operator = (const Diagonal &A)
{
	if (this == &A) return *this;
	Set::VectorSpace::Vector::operator = (A);
	return *this;
}

Diagonal & 
Diagonal::operator = (const Set::VectorSpace::Vector &A)
{
	Set::VectorSpace::Vector::operator = (A);
	return *this;
}

Set::VectorSpace::Vector 
Diagonal::operator () (const Set::VectorSpace::Vector &A) const
{
	assert (n == A.size());
	double *p, *q, *r; Set::VectorSpace::Vector B(n);
	for (p=this->head, q=A.begin(), r=B.begin(); p<this->tail; p++, q++, r++) 
		*r = (*p)*(*q); return B;
}

double 
Diagonal::operator () (
	const Set::VectorSpace::Vector &A, 
	const Set::VectorSpace::Vector &B) const
{
	assert (n == A.size()); assert (n == B.size());
	double *p, *q, *r; Set::VectorSpace::Vector C(n); 
	for (p=this->head, q=A.begin(), r=C.begin(); p<this->tail; p++, q++, r++) 
		*r = (*p)*(*q); return C(B);
}

double 
Diagonal::operator () (const Diagonal &A) const
{
	double a = 0.0; double *p, *q; 
	for (p=this->head, q=A.begin(); p<this->tail; p++, q++) a += (*q)*(*p);
	return a;
}
double 
Diagonal::operator () (const Set::VectorSpace::Hom &A) const
{
	assert (n == A.size1()); assert (n == A.size2());
	unsigned int n1 = n+1; double a = 0.0; double *p, *q;
	for (p=this->head, q=A.begin(); p<this->tail; p++, q+=n1) 
		a += (*q)*(*p); return a;
}

Set::VectorSpace::Hom 
Diagonal::Embed() const
{
	unsigned int n1 = n+1; 
	double *p, *q; Set::VectorSpace::Hom B(n); 
	for (p=this->begin(), q=B.begin(); p<this->end(); p++, q+=n1) *q = *p;
	return B;
}

}

}

//////////////////////////////////////////////////////////////////////
// Operators
//////////////////////////////////////////////////////////////////////
/*
bool operator != (const Set::VectorSpace::Diagonal &A, 
				  const Set::VectorSpace::Hom &B)
{
	return (A.Embed() != B);
}

bool operator != (const Set::VectorSpace::Hom &B, 
				  const Set::VectorSpace::Diagonal &A)
{
	return (A.Embed() != B);
}

bool operator == (const Set::VectorSpace::Diagonal &A, 
				  const Set::VectorSpace::Hom &B)
{
	return (A.Embed() == B);
}

bool operator == (const Set::VectorSpace::Hom &B,
				  const Set::VectorSpace::Diagonal &A)
{
	return (A.Embed() == B);
}
*/
void Random(Set::VectorSpace::Diagonal &A)
{
	double *p;
	for (p=A.begin(); p<A.end(); p++) 
		*p = 2.0*((double)rand()/(double)RAND_MAX) - 1.0;
}

Set::VectorSpace::Diagonal
operator + (const Set::VectorSpace::Diagonal &A, 
			const Set::VectorSpace::Diagonal &B)
{
	Set::VectorSpace::Diagonal C=A; C+=B; return C;
}

Set::VectorSpace::Hom
operator + (const Set::VectorSpace::Diagonal &A, 
			const Set::VectorSpace::Hom &B)
{
	unsigned int n = A.size(), n1 = n+1;
	assert (n == B.size1()); assert (n == B.size2()); 
	Set::VectorSpace::Hom C=B; double *p, *q;
	for (p=A.begin(), q=C.begin(); p<A.end(); p++, q+=n1) 
		*q += *p; return C;
}

Set::VectorSpace::Hom
operator + (const Set::VectorSpace::Hom &B,
			const Set::VectorSpace::Diagonal &A)
{
	unsigned int n = A.size(), n1 = n+1;
	assert (n == B.size1()); assert (n == B.size2()); 
	Set::VectorSpace::Hom C=B; double *p, *q;
	for (p=A.begin(), q=C.begin(); p<A.end(); p++, q+=n1) 
		*q += *p; return C;
}

Set::VectorSpace::Diagonal
operator - (const Set::VectorSpace::Diagonal &A, 
			const Set::VectorSpace::Diagonal &B)
{
	Set::VectorSpace::Diagonal C=A; C-=B; return C;
}

Set::VectorSpace::Hom
operator - (const Set::VectorSpace::Diagonal &A, 
			const Set::VectorSpace::Hom &B)
{
	unsigned int n = A.size(), n1 = n+1;
	assert (n == B.size1()); assert (n == B.size2()); 
	Set::VectorSpace::Hom C=B; double *p, *q;
	for (p=A.begin(), q=C.begin(); p<A.end(); p++, q+=n1) 
		*q -= *p; return C;
}

Set::VectorSpace::Hom
operator - (const Set::VectorSpace::Hom &B,
			const Set::VectorSpace::Diagonal &A)
{
	unsigned int n = A.size(), n1 = n+1;
	assert (n == B.size1()); assert (n == B.size2()); 
	Set::VectorSpace::Hom C=B; double *p, *q;
	for (p=A.begin(), q=C.begin(); p<A.end(); p++, q+=n1) 
		*q -= *p; return C;
}
 
Set::VectorSpace::Diagonal 
operator - (const Set::VectorSpace::Diagonal &A)
{
	Set::VectorSpace::Diagonal B=A; double *p;
	for (p=B.begin(); p<B.end(); p++) *p = -(*p);
	return B;
}

void Null(Set::VectorSpace::Diagonal &A)
{
	double *p;
	for (p=A.begin(); p<A.end(); p++) *p = 0.0;
}
 
Set::VectorSpace::Diagonal 
operator * (const Set::VectorSpace::Diagonal &A, 
			const double &b)
{
	Set::VectorSpace::Diagonal B=A; double *p;
	for (p=B.begin(); p<B.end(); p++) *p *= b;
	return B;
}

Set::VectorSpace::Diagonal 
operator * (const double &b, 
			const Set::VectorSpace::Diagonal &A)
{
	Set::VectorSpace::Diagonal B=A; double *p;
	for (p=B.begin(); p<B.end(); p++) *p *= b;
	return B;
}

Set::VectorSpace::Diagonal 
operator / (const Set::VectorSpace::Diagonal &A, 
			const double &b)
{
	Set::VectorSpace::Diagonal B=A; double *p;
	for (p=B.begin(); p<B.end(); p++) *p /= b;
	return B;
}

Set::VectorSpace::Vector
operator * (const Set::VectorSpace::Diagonal &A, 
			const Set::VectorSpace::Vector &B)
{
	return A(B);
}

 
Set::VectorSpace::Diagonal 
operator * (const Set::VectorSpace::Diagonal &A, 
			const Set::VectorSpace::Diagonal &B)
{
	unsigned int n = A.size(); assert(n == B.size());
	Set::VectorSpace::Diagonal C(n); double *p, *q, *r;
	for (p=A.begin(), q=B.begin(), r=C.begin(); p<A.end(); p++, q++, r++) 
		*r = (*p)*(*q); return C;
}

 
Set::VectorSpace::Hom
operator * (const Set::VectorSpace::Diagonal &A, 
			const Set::VectorSpace::Hom &B)
{
	unsigned int n = A.size(); 
	assert(n == B.size1());
	assert(n == B.size2());
	unsigned int i, j; 
	double *p, *q, *r;
	Set::VectorSpace::Hom C(n);
	for (j=0, q=B.begin(), r=C.begin(); j<n; j++)
		for (i=0, p=A.begin(); i<n; i++, p++, q++, r++)
			*r = (*p)*(*q); return C;
}

Set::VectorSpace::Hom
operator * (const Set::VectorSpace::Hom &B, 
			const Set::VectorSpace::Diagonal &A)
{
	unsigned int n = A.size(); 
	assert(n == B.size1());
	assert(n == B.size2());
	unsigned int i, j; 
	double *p, *q, *r;
	Set::VectorSpace::Hom C(n);
	for (j=0, p=A.begin(), q=B.begin(), r=C.begin(); j<n; j++, p++)
		for (i=0; i<n; i++, q++, r++)
			*r = (*p)*(*q); return C;
}

Set::VectorSpace::Diagonal
Adjoint(const Set::VectorSpace::Diagonal &A)
{
	return A;
}

 
double Norm(const Set::VectorSpace::Diagonal &A)
{
	double a=0.0; double *p;
	for (p=A.begin(); p<A.end(); p++) 
		if (fabs(*p) > a) a = fabs(*p);
	return a;
}
 
Set::VectorSpace::Diagonal
Exp(const Set::VectorSpace::Diagonal &A)
{
	double *p, *q; Set::VectorSpace::Diagonal B(A.size());
	for (p=A.begin(), q=B.begin(); p<A.end(); p++, q++)
		*q = exp(*p); return B;
}

 
Set::VectorSpace::Diagonal
Log(const Set::VectorSpace::Diagonal &A)
{
	double *p, *q; Set::VectorSpace::Diagonal B(A.size());
	for (p=A.begin(), q=B.begin(); p<A.end(); p++, q++)
		*q = log(*p); return B;
}

Set::VectorSpace::Diagonal
Inverse(const Set::VectorSpace::Diagonal &A)
{
	double *p, *q; Set::VectorSpace::Diagonal B(A.size());
	for (p=A.begin(), q=B.begin(); p<A.end(); p++, q++)
		*q = 1.0/(*p); return B;
}

void Identity(Set::VectorSpace::Diagonal &A)
{
	double *p;
	for (p=A.begin(); p<A.end(); p++) *p = 1.0;
}

double 
Jacobian(const Set::VectorSpace::Diagonal &A)
{
	double a=1.0; double *p;
	for (p=A.begin(); p<A.end(); p++) a *= *p;
	return a;
}

 
double 
Trace(const Set::VectorSpace::Diagonal &A)
{
	double a=0.0; double *p;
	for (p=A.begin(); p<A.end(); p++) a += *p;
	return a;
}
