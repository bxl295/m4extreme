// Orthonormal.cpp: Implementation of the Orthonormal class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <cassert>
#include "./Orthonormal.h"

namespace Set
{
namespace Euclidean
{	
namespace Orthonormal
{
//////////////////////////////////////////////////////////////////////
// Class Point
//////////////////////////////////////////////////////////////////////

Point::Point() {}

Point::Point(const unsigned int &n) : Set::Array(n) {}

Point::Point(const unsigned int &n, double * const &u) : Set::Array(n,u) {}

Set::Manifold::Point *Point::Clone() const
{
	return new Point(*this);
}

Point::~Point(){}

Point::Point(const Point &A) : Set::Array(A) {}

Point & 
Point::operator = (const Point &A)
{
	if (this == &A) return *this; 
	Set::Array::operator = (A);
	return *this;
}

void
Point::Randomize()
{
	Set::Array::Randomize();
}

unsigned int
Point::size() const
{
	return Set::Array::size();
}

void 
Point::print(ostream *os)
{
	Set::Array::print(os);
}

void
Point::operator += (const Vector &A)
{
	assert(n == A.size()); double *p, *q; 
	for (p=head, q=A.begin(); p<tail; p++, q++) *p += *q;
}

void 
Point::operator -= (const Vector &A)
{
	assert(n == A.size()); double *p, *q; 
	for (p=head, q=A.begin(); p<tail; p++, q++) *p -= *q;
}

void 
Point::operator *= (const double &a)
{
	double *p; 
	for (p=this->head; p<this->tail; p++) *p *= a;
}

void 
Point::operator /= (const double &a)
{
	double *p; 
	for (p=this->head; p<this->tail; p++) *p /= a;
}

Set::Manifold::Point & 
Point::operator = (const Set::Manifold::Point &A)
{
	if (this == &A) return *this; 
	operator = ((Point &)A);
	return *this;
}

bool 
Point::operator != (const Set::Manifold::Point &P) const
{
	Point &Q = (Point &)P;
	return ((Set::Array &)(*this) != (Set::Array &)Q); 
}

bool 
Point::operator == (const Set::Manifold::Point &P) const
{
	Point &Q = (Point &)P;
	return ((Set::Array &)(*this) == (Set::Array &)Q); 
}

void 
Point::operator += (const Set::VectorSpace::Vector &A)
{
	operator += (Vector(A));
}

void 
Point::operator -= (const Set::VectorSpace::Vector &A)
{
	operator -= (Vector(A));
}

//////////////////////////////////////////////////////////////////////
// Class Vector
//////////////////////////////////////////////////////////////////////

Vector::Vector() {}

Vector::Vector(const unsigned int &n) 
: Set::VectorSpace::Vector(n) {}

Vector::Vector(const unsigned int &n, double * const &u)
: Set::VectorSpace::Vector(n,u) {}

Set::Manifold::Point *Vector::Clone() const
{
	return new Vector(*this);
}

Vector::~Vector() {}

Vector::Vector(const Vector &A) 
: Set::VectorSpace::Vector(A) {}

Vector::Vector(const Set::VectorSpace::Vector &A) 
: Set::VectorSpace::Vector(A) {}

Vector::Vector(const Covector &A) 
: Set::VectorSpace::Vector(A) {}

Vector::Vector(const Set::Euclidean::Orthonormal::Point &A) 
: Set::VectorSpace::Vector(A) {}

Vector & 
Vector::operator = (const Vector &A)
{
	if (this == &A) return *this;
	Set::VectorSpace::Vector::operator = (A);
	return *this;
}

Vector & 
Vector::operator = (const Set::VectorSpace::Vector &A)
{
	if (this == &A) return *this;
	Set::VectorSpace::Vector::operator = (A);
	return *this;
}

Vector & 
Vector::operator = (const Set::Euclidean::Orthonormal::Point &A)
{
	Set::VectorSpace::Vector::operator = (A);
	return *this;
}

void
Vector::Randomize()
{
	Set::VectorSpace::Vector::Randomize();
}

unsigned int
Vector::size() const
{
	return n;
}

void
Vector::print(ostream *os)
{
	Set::VectorSpace::Vector::print(os);
}

Set::Manifold::Point & 
Vector::operator = (const Set::Manifold::Point &A)
{
	if (this == &A) return *this; 
	operator = ((Vector &)A);
	return *this;
}

bool 
Vector::operator != (const Set::Manifold::Point &P) const
{
	return (*this != (Vector &)P); 
}

bool 
Vector::operator == (const Set::Manifold::Point &P) const
{
	return (*this == (Vector &)P); 
}

void 
Vector::operator += (const Set::VectorSpace::Vector &A)
{
	Set::VectorSpace::Vector::operator += (A);
}

void 
Vector::operator -= (const Set::VectorSpace::Vector &A)
{
	Set::VectorSpace::Vector::operator -= (A);
}

//////////////////////////////////////////////////////////////////////
// Class VectorZero
//////////////////////////////////////////////////////////////////////

VectorZero::VectorZero() 
: Set::VectorSpace::VectorZero() {}

VectorZero::VectorZero(const unsigned int &n) 
: Set::VectorSpace::VectorZero(n) {}

VectorZero::VectorZero(const unsigned int &n, double * const &u) 
: Set::VectorSpace::VectorZero(n) {}

VectorZero::~VectorZero(){}

VectorZero::VectorZero(const VectorZero &A) 
: Set::VectorSpace::VectorZero(A) {}

VectorZero::VectorZero(const Set::VectorSpace::VectorZero &A) 
: Set::VectorSpace::VectorZero(A) {}

VectorZero::VectorZero(const CovectorZero &A) 
: Set::VectorSpace::VectorZero(A) {}

VectorZero & 
VectorZero::operator = (const VectorZero &A)
{
	assert(n == A.n); 
	if (this == &A) return *this;
	return *this;
}

VectorZero & 
VectorZero::operator = (const Set::VectorSpace::VectorZero &A)
{
	assert(n == A.size()); 
	return *this;
}

//////////////////////////////////////////////////////////////////////
// Class Covector
//////////////////////////////////////////////////////////////////////

Covector::Covector() {}

Covector::~Covector() {}

Covector::Covector(const unsigned int &n) 
: Set::VectorSpace::Vector(n) {}

Covector::Covector(const unsigned int &n, double * const &u)
: Set::VectorSpace::Vector(n,u) {}

Covector::Covector(const Covector &A) 
: Set::VectorSpace::Vector(A) {}

Covector::Covector(const Set::VectorSpace::Vector &A) 
: Set::VectorSpace::Vector(A) {}

Covector::Covector(const Set::Euclidean::Orthonormal::Vector &A) 
: Set::VectorSpace::Vector(A) {}

Covector & Covector::operator = (const Set::VectorSpace::Vector &A)
{
	Set::VectorSpace::Vector::operator = (A);
	return *this;
}

Covector & Covector::operator = (const Covector &A)
{
	if (this == &A) return *this;
	Set::VectorSpace::Vector::operator = (A);
	return *this;
}

//////////////////////////////////////////////////////////////////////
// Class CovectorZero
//////////////////////////////////////////////////////////////////////

CovectorZero::CovectorZero() 
: Set::VectorSpace::VectorZero() {}

CovectorZero::CovectorZero(const unsigned int &n) 
: Set::VectorSpace::VectorZero(n) {}

CovectorZero::CovectorZero(const unsigned int &n, double * const &u) 
: Set::VectorSpace::VectorZero(n) {}

CovectorZero::~CovectorZero(){}

CovectorZero::CovectorZero(const Set::VectorSpace::VectorZero &A) 
: Set::VectorSpace::VectorZero(A) {}

CovectorZero::CovectorZero(const CovectorZero &A) 
: Set::VectorSpace::VectorZero(A) {}

CovectorZero::CovectorZero(const Set::Euclidean::Orthonormal::VectorZero &A) 
: Set::VectorSpace::VectorZero(A) {}

CovectorZero & 
CovectorZero::operator = (const CovectorZero &A)
{
	assert(n == A.n); 
	if (this == &A) return *this;
	return *this;
}

CovectorZero & 
CovectorZero::operator = (const Set::VectorSpace::VectorZero &A)
{
	assert(n == A.size()); 
	return *this;
}

}

}

}

Set::Euclidean::Orthonormal::Point
operator + (const Set::Euclidean::Orthonormal::Point &P,
			const Set::Euclidean::Orthonormal::Vector &A)
{
	Set::Euclidean::Orthonormal::Point Q = P; 
	Q += A; return Q;
}

Set::Euclidean::Orthonormal::Point
operator - (const Set::Euclidean::Orthonormal::Point &P,
			const Set::Euclidean::Orthonormal::Vector &A)
{
	Set::Euclidean::Orthonormal::Point Q = P; 
	Q -= A; return Q;
}

Set::Euclidean::Orthonormal::Vector
operator - (const Set::Euclidean::Orthonormal::Point &Q,
			const Set::Euclidean::Orthonormal::Point &P)
{
	Set::Euclidean::Orthonormal::Vector B(Q);
	Set::Euclidean::Orthonormal::Vector A(P);
	return B - A;
}

Set::Euclidean::Orthonormal::Point
operator + (const Set::Euclidean::Orthonormal::Point &P,
			const Set::VectorSpace::Vector &A)
{
	Set::Euclidean::Orthonormal::Point Q = P; 
	Q += A; return Q;
}

Set::Euclidean::Orthonormal::Point
operator - (const Set::Euclidean::Orthonormal::Point &P,
			const Set::VectorSpace::Vector &A)
{
	Set::Euclidean::Orthonormal::Point Q = P; 
	Q -= A; return Q;
}

void Random(Set::Euclidean::Orthonormal::Point &P){P.Randomize();}
