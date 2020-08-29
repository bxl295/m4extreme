// LinearMapping.cpp: Implementation of the LinearMapping class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <cassert>
#include "LinearMapping.h"

namespace Set
{
namespace LinearMapping
{
//////////////////////////////////////////////////////////////////////
// Class Point
//////////////////////////////////////////////////////////////////////

Point::Point() {}

Point::Point(
	const unsigned int &n1, 
	const unsigned int &n2) : 
	Set::VectorSpace::Hom(n1,n2) {}

Point::Point(
	const unsigned int &n1, 
	const unsigned int &n2, 
	double * const &u) : 
	Set::VectorSpace::Hom(n1,n2,u) {}

Point::Point(const unsigned int &n) : 
	Set::VectorSpace::Hom(n) {}

Point::Point(const unsigned int &n, double * const &u) : 
	Set::VectorSpace::Hom(n,u) {}

Set::Manifold::Point *Point::Clone() const
{
	return new Point(*this);
}

Point::~Point(){}

Point::Point(const Point &A) : 
	Set::VectorSpace::Hom(A) {}

Point & 
Point::operator = (const Point &A)
{
	if (this == &A) return *this; 
	Set::VectorSpace::Hom::operator = (A);
	return *this;
}

Point & 
Point::operator = (const Set::VectorSpace::Hom &A)
{
	if (this == &A) return *this; 
	Set::VectorSpace::Hom::operator = (A);
	return *this;
}

void 
Point::Randomize()
{
	Set::VectorSpace::Hom::Randomize();
}

unsigned int 
Point::size() const
{
	return n;
}

void 
Point::print(ostream *os)
{
	Set::VectorSpace::Hom::print(os);
}

void 
Point::operator += (const Set::LinearMapping::Vector &A)
{
	assert(n == A.size()); double *p, *q; 
	for (p=head, q=A.begin(); p<tail; p++, q++) *p += *q;
}

void 
Point::operator -= (const Set::LinearMapping::Vector &A)
{
	assert(n == A.size()); double *p, *q; 
	for (p=head, q=A.begin(); p<tail; p++, q++) *p -= *q;
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
	return ((Set::VectorSpace::Hom &)(*this) != (Set::VectorSpace::Hom &)Q); 
}

bool 
Point::operator == (const Set::Manifold::Point &P) const
{
	Point &Q = (Point &)P;
	return ((Set::VectorSpace::Hom &)(*this) == (Set::VectorSpace::Hom &)Q); 
}

void 
Point::operator += (const Set::VectorSpace::Vector &A)
{
	Set::VectorSpace::Hom::operator += (A);
}

void 
Point::operator -= (const Set::VectorSpace::Vector &A)
{
	Set::VectorSpace::Hom::operator -= (A);
}

void 
Point::operator += (const Set::VectorSpace::Hom &A)
{
//	operator += (Vector(A));
	Set::VectorSpace::Hom::operator += (A);
}

void 
Point::operator -= (const Set::VectorSpace::Hom &A)
{
//	operator -= (Vector(A));
	Set::VectorSpace::Hom::operator -= (A);
}

//////////////////////////////////////////////////////////////////////
// Class Vector
//////////////////////////////////////////////////////////////////////

Vector::Vector() {}

Vector::Vector(
	const unsigned int &n1, 
	const unsigned int &n2) : 
	Set::VectorSpace::Hom(n1,n2) {}

Vector::Vector(
	const unsigned int &n1, 
	const unsigned int &n2, 
	double * const &u) : 
	Set::VectorSpace::Hom(n1,n2,u) {}

Vector::Vector(const unsigned int &n) 
: Set::VectorSpace::Hom(n) {}

Vector::Vector(const unsigned int &n, double * const &u)
: Set::VectorSpace::Hom(n,u) {}

Set::Manifold::Point *Vector::Clone() const
{
	return new Vector(*this);
}

Vector::~Vector() {}

Vector::Vector(const Vector &A) 
: Set::VectorSpace::Hom(A) {}

Vector::Vector(const Set::VectorSpace::Hom &A) 
: Set::VectorSpace::Hom(A) {}

Vector::Vector(const Covector &A) 
: Set::VectorSpace::Hom(A) {}

Vector::Vector(const Set::LinearMapping::Point &A) 
: Set::VectorSpace::Hom(A) {}

Vector & 
Vector::operator = (const Vector &A)
{
	if (this == &A) return *this;
	Set::VectorSpace::Hom::operator = (A);
	return *this;
}

Vector & 
Vector::operator = (const Set::VectorSpace::Vector &A)
{
	if (this == &A) return *this;
	Set::VectorSpace::Hom::operator = (A);
	return *this;
}

Vector & 
Vector::operator = (const Set::VectorSpace::Hom &A)
{
	if (this == &A) return *this;
	Set::VectorSpace::Hom::operator = (A);
	return *this;
}

Vector & 
Vector::operator = (const Set::LinearMapping::Point &A)
{
	Set::VectorSpace::Hom::operator = (A);
	return *this;
}

void 
Vector::Randomize()
{
	Set::VectorSpace::Hom::Randomize();
}

unsigned int 
Vector::size() const
{
	return n;
}

void 
Vector::print(ostream *os)
{
	Set::VectorSpace::Hom::print(os);
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
	Set::VectorSpace::Hom::operator += (A);
}

void 
Vector::operator -= (const Set::VectorSpace::Vector &A)
{
	Set::VectorSpace::Hom::operator -= (A);
}

void 
Vector::operator += (const Set::VectorSpace::Hom &A)
{
	Set::VectorSpace::Hom::operator += (A);
}

void 
Vector::operator -= (const Set::VectorSpace::Hom &A)
{
	Set::VectorSpace::Hom::operator -= (A);
}

//////////////////////////////////////////////////////////////////////
// Class VectorZero
//////////////////////////////////////////////////////////////////////

VectorZero::VectorZero() 
: Set::VectorSpace::HomZero() {}

VectorZero::VectorZero(
	const unsigned int &n1, 
	const unsigned int &n2) : 
	Set::VectorSpace::HomZero(n1,n2) {}

VectorZero::VectorZero(
	const unsigned int &n1, 
	const unsigned int &n2, 
	double * const &u) : 
	Set::VectorSpace::HomZero(n1,n2,u) {}

VectorZero::VectorZero(const unsigned int &n) 
: Set::VectorSpace::HomZero(n) {}

VectorZero::VectorZero(const unsigned int &n, double * const &u) 
: Set::VectorSpace::HomZero(n) {}

VectorZero::~VectorZero(){}

VectorZero::VectorZero(const VectorZero &A) 
: Set::VectorSpace::HomZero(A) {}

VectorZero::VectorZero(const Set::VectorSpace::HomZero &A) 
: Set::VectorSpace::HomZero(A) {}

VectorZero::VectorZero(const CovectorZero &A) 
: Set::VectorSpace::HomZero(A) {}

VectorZero & 
VectorZero::operator = (const VectorZero &A)
{
	assert(n == A.n); 
	if (this == &A) return *this;
	return *this;
}

VectorZero & 
VectorZero::operator = (const Set::VectorSpace::HomZero &A)
{
	assert(n == A.size()); 
	return *this;
}

//////////////////////////////////////////////////////////////////////
// Class Covector
//////////////////////////////////////////////////////////////////////

Covector::Covector() {}

Covector::Covector(
	const unsigned int &n1, 
	const unsigned int &n2) : 
	Set::VectorSpace::Hom(n1,n2) {}

Covector::Covector(
	const unsigned int &n1, 
	const unsigned int &n2, 
	double * const &u) : 
	Set::VectorSpace::Hom(n1,n2,u) {}

Covector::Covector(const unsigned int &n) 
: Set::VectorSpace::Hom(n) {}

Covector::Covector(const unsigned int &n, double * const &u)
: Set::VectorSpace::Hom(n,u) {}

Covector::~Covector() {}

Covector::Covector(const Covector &A) 
: Set::VectorSpace::Hom(A) {}

Covector::Covector(const Set::VectorSpace::Hom &A) 
: Set::VectorSpace::Hom(A) {}

Covector::Covector(const Set::LinearMapping::Vector &A) 
: Set::VectorSpace::Hom(A) {}

Covector & Covector::operator = (const Set::VectorSpace::Vector &A)
{
	Set::VectorSpace::Hom::operator = (A);
	return *this;
}

Covector & Covector::operator = (const Set::VectorSpace::Hom &A)
{
	Set::VectorSpace::Hom::operator = (A);
	return *this;
}

Covector & Covector::operator = (const Covector &A)
{
	if (this == &A) return *this;
	Set::VectorSpace::Hom::operator = (A);
	return *this;
}

//////////////////////////////////////////////////////////////////////
// Class CovectorZero
//////////////////////////////////////////////////////////////////////

CovectorZero::CovectorZero() 
: Set::VectorSpace::HomZero() {}

CovectorZero::CovectorZero(
	const unsigned int &n1, 
	const unsigned int &n2) : 
	Set::VectorSpace::HomZero(n1,n2) {}

CovectorZero::CovectorZero(
	const unsigned int &n1, 
	const unsigned int &n2, 
	double * const &u) : 
	Set::VectorSpace::HomZero(n1,n2,u) {}

CovectorZero::CovectorZero(const unsigned int &n) 
: Set::VectorSpace::HomZero(n) {}

CovectorZero::CovectorZero(const unsigned int &n, double * const &u) 
: Set::VectorSpace::HomZero(n) {}

CovectorZero::~CovectorZero(){}

CovectorZero::CovectorZero(const Set::VectorSpace::HomZero &A) 
: Set::VectorSpace::HomZero(A) {}

CovectorZero::CovectorZero(const CovectorZero &A) 
: Set::VectorSpace::HomZero(A) {}

CovectorZero::CovectorZero(const Set::LinearMapping::VectorZero &A) 
: Set::VectorSpace::HomZero(A) {}

CovectorZero & 
CovectorZero::operator = (const CovectorZero &A)
{
	assert(n == A.n); 
	if (this == &A) return *this;
	return *this;
}

CovectorZero & 
CovectorZero::operator = (const Set::VectorSpace::HomZero &A)
{
	assert(n == A.size()); 
	return *this;
}

}

}

Set::LinearMapping::Point
operator + (const Set::LinearMapping::Point &P,
			const Set::LinearMapping::Vector &A)
{
	Set::LinearMapping::Point Q = P; 
	Q += A; return Q;
}

Set::LinearMapping::Point
operator - (const Set::LinearMapping::Point &P,
			const Set::LinearMapping::Vector &A)
{
	Set::LinearMapping::Point Q = P; 
	Q -= A; return Q;
}

Set::LinearMapping::Vector
operator - (const Set::LinearMapping::Point &Q,
			const Set::LinearMapping::Point &P)
{
	Set::LinearMapping::Vector B(Q);
	Set::LinearMapping::Vector A(P);
	return B - A;
}

Set::LinearMapping::Point
operator + (const Set::LinearMapping::Point &P,
			const Set::VectorSpace::Vector &A)
{
	Set::LinearMapping::Point Q = P; 
	Q += A; return Q;
}

Set::LinearMapping::Point
operator - (const Set::LinearMapping::Point &P,
			const Set::VectorSpace::Vector &A)
{
	Set::LinearMapping::Point Q = P; 
	Q -= A; return Q;
}

Set::LinearMapping::Point
operator + (const Set::LinearMapping::Point &P,
			const Set::VectorSpace::Hom &A)
{
	Set::LinearMapping::Point Q = P; 
	Q += A; return Q;
}

Set::LinearMapping::Point
operator - (const Set::LinearMapping::Point &P,
			const Set::VectorSpace::Hom &A)
{
	Set::LinearMapping::Point Q = P; 
	Q -= A; return Q;
}

void Random(Set::LinearMapping::Point &P){P.Randomize();}
