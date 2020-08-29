// Cartesian.cpp: Implementation of the Cartesian class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <cassert>
#include "Cartesian.h"

namespace Set
{
namespace Euclidean
{	
namespace Cartesian
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
	return n;
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

Vector::Vector(const Set::Euclidean::Cartesian::Point &A) 
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
Vector::operator = (const Set::Euclidean::Cartesian::Point &A)
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

Covector::Covector(const Set::Euclidean::Cartesian::Vector &A) 
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

CovectorZero::CovectorZero(const Set::Euclidean::Cartesian::VectorZero &A) 
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

//////////////////////////////////////////////////////////////////////
// Class Embedding<0>
//////////////////////////////////////////////////////////////////////

Embedding<0>::Embedding() {}

Embedding<0>::Embedding(
	const unsigned int &m1, 
	const unsigned int &m2)
: n1(m1), n2(m2), x0(m1), y0(m2), y(m2), A(m2,m1) {}

Embedding<0>::Embedding(
	const Set::Euclidean::Cartesian::Point &a0,
	const Set::Euclidean::Orthonormal::Point &b0,
	const Set::VectorSpace::Hom &B) : 
	x0(a0), y0(b0), y(B.size1()), A(B)
{
	n1 = x0.size(); n2 = y0.size();
	assert (A.size1() == n2);
	assert (A.size2() == n1);
}

Embedding<0>::Embedding(const Embedding<0> &f) 
: n1(f.n1), n2(f.n2), x0(f.x0), y0(f.y0), y(f.y), A(f.A){}

Embedding<0>::~Embedding(){}

Embedding<0> & 
Embedding<0>::operator = (const Embedding<0> &f)
{
	assert(f.n1 == n1); assert(f.n2 == n2);
	if (this == &f) return *this;
	x0 = f.x0; y0 = f.y0; A = f.A;
	return *this;
}

Set::Euclidean::Cartesian::Point 
Embedding<0>::DomainOrigin() const
{
	return x0;
}

Set::Euclidean::Cartesian::Point & 
Embedding<0>::DomainOrigin()
{
	return x0;
}

Set::Euclidean::Orthonormal::Point 
Embedding<0>::RangeOrigin() const
{
	return y0;
}

Set::Euclidean::Orthonormal::Point & 
Embedding<0>::RangeOrigin()
{
	return y0;
}

Set::VectorSpace::Hom 
Embedding<0>::LinearMapping() const
{
	return A;
}

Set::VectorSpace::Hom & 
Embedding<0>::LinearMapping()
{
	return A;
}

Set::Manifold::Map * 
Embedding<0>::Clone()
{
	return new Embedding<0>(*this);
}
	
Set::Manifold::TMap *
Embedding<0>::Diff()
{
	return new Embedding<1>(*this);
}

void 
Embedding<0>::Randomize()
{
	Random(x0); Random(y0); Random(A); 
}

const Set::Euclidean::Orthonormal::Point &
Embedding<0>::operator () (
	const Set::Euclidean::Cartesian::Point &x)
{
	y = y0 + A*(x-x0); return y;
}

unsigned int 
Embedding<0>::size1() const
{
	return n1;

}
unsigned int 
Embedding<0>::size2() const
{
	return n2;
}

const Set::Manifold::Point & 
Embedding<0>::operator () (const Set::Manifold::Point &P)
{
	return operator () ((Point &)P);
}

//////////////////////////////////////////////////////////////////////
// Class Embedding<1>
//////////////////////////////////////////////////////////////////////

Embedding<1>::Embedding() {}

Embedding<1>::Embedding(const Set::VectorSpace::Hom &B) : A(B) {}

Embedding<1>::Embedding(const Embedding<0> &f) 
: n1(f.size1()), A(f.LinearMapping())
{
	n2 = A.size();
}

Embedding<1>::Embedding(const Embedding<1> &f)
: n1(f.n1), n2(f.n2), A(f.A) {}

Embedding<1>::~Embedding(){}

Embedding<1> & 
Embedding<1>::operator = (const Embedding<1> &f)
{
	assert(f.n1 == n1); assert(f.n2 == n2);
	if (this == &f) return *this; 
	A = f.A; return *this;
}

Set::VectorSpace::Hom 
Embedding<1>::LinearMapping() const
{
	return A;
}

Set::VectorSpace::Hom & 
Embedding<1>::LinearMapping()
{
	return A;
}

Set::Manifold::TMap * 
Embedding<1>::Clone()
{
	return new Embedding<1>(*this);
}
	
Set::Manifold::TMap *
Embedding<1>::Diff()
{
	return new Embedding<2>(*this);
}

void 
Embedding<1>::Randomize()
{
	Random(A);
}

const Set::VectorSpace::Hom &
Embedding<1>::operator () (
	const Set::Euclidean::Cartesian::Point &)
{
	return A;
}

unsigned int 
Embedding<1>::size1() const
{
	return n1;

}
unsigned int 
Embedding<1>::size2() const
{
	return n2;
}

const Set::VectorSpace::Hom & 
Embedding<1>::operator () (const Set::Manifold::Point &P)
{
	return operator () ((Point &)P);
}

//////////////////////////////////////////////////////////////////////
// Class Embedding<2>
//////////////////////////////////////////////////////////////////////

Embedding<2>::Embedding() {}

Embedding<2>::Embedding(const Embedding<1> &f) 
: n1(f.size1()), A(f.size2(),f.size1())
{
	n2 = A.size();
}

Embedding<2>::Embedding(const Embedding<2> &f)
: n1(f.n1), n2(f.n2), A(f.A) {}

Embedding<2>::~Embedding(){}

Embedding<2> & 
Embedding<2>::operator = (const Embedding<2> &f)
{
	assert(f.n1 == n1); assert(f.n2 == n2);
	if (this == &f) return *this; 
	A = f.A; return *this;
}

Set::VectorSpace::Hom 
Embedding<2>::LinearMapping() const
{
	return A;
}

Set::VectorSpace::Hom & 
Embedding<2>::LinearMapping()
{
	return A;
}

Set::Manifold::TMap * 
Embedding<2>::Clone()
{
	return new Embedding<2>(*this);
}

void 
Embedding<2>::Randomize() {}

const Set::VectorSpace::HomZero &
Embedding<2>::operator () (
	const Set::Euclidean::Cartesian::Point &)
{
	return A;
}

unsigned int 
Embedding<2>::size1() const
{
	return n1;

}
unsigned int 
Embedding<2>::size2() const
{
	return n2;
}

const Set::VectorSpace::Hom & 
Embedding<2>::operator () (const Set::Manifold::Point &P)
{
	return operator () ((Point &)P);
}

//////////////////////////////////////////////////////////////////////
// Class Submersion<0>
//////////////////////////////////////////////////////////////////////

Submersion<0>::Submersion() {}

Submersion<0>::Submersion(
	const unsigned int &m1, 
	const unsigned int &m2)
: n1(m1), n2(m2), x0(m1), y0(m2), y(m2), A(m2,m1) {}

Submersion<0>::Submersion(
	const Set::Euclidean::Orthonormal::Point &a0,
	const Set::Euclidean::Cartesian::Point &b0,
	const Set::VectorSpace::Hom &B) : 
	x0(a0), y0(b0), y(B.size1()), A(B)
{
	n1 = x0.size(); n2 = y0.size();
	assert (A.size1() == n2);
	assert (A.size2() == n1);
}

Submersion<0>::Submersion(const Submersion<0> &f) 
: n1(f.n1), n2(f.n2), x0(f.x0), y0(f.y0), y(f.y), A(f.A){}

Submersion<0>::~Submersion(){}

Submersion<0> & 
Submersion<0>::operator = (const Submersion<0> &f)
{
	assert(f.n1 == n1); assert(f.n2 == n2);
	if (this == &f) return *this;
	x0 = f.x0; y0 = f.y0; A = f.A;
	return *this;
}

Set::Euclidean::Orthonormal::Point 
Submersion<0>::DomainOrigin() const
{
	return x0;
}

Set::Euclidean::Orthonormal::Point & 
Submersion<0>::DomainOrigin()
{
	return x0;
}

Set::Euclidean::Cartesian::Point 
Submersion<0>::RangeOrigin() const
{
	return y0;
}

Set::Euclidean::Cartesian::Point & 
Submersion<0>::RangeOrigin()
{
	return y0;
}

Set::VectorSpace::Hom 
Submersion<0>::LinearMapping() const
{
	return A;
}

Set::VectorSpace::Hom & 
Submersion<0>::LinearMapping()
{
	return A;
}

Set::Manifold::Map * 
Submersion<0>::Clone()
{
	return new Submersion<0>(*this);
}
	
Set::Manifold::TMap *
Submersion<0>::Diff()
{
	return new Submersion<1>(*this);
}

void 
Submersion<0>::Randomize()
{
	Random(x0); Random(y0); Random(A); 
}

const Set::Euclidean::Cartesian::Point &
Submersion<0>::operator () (
	const Set::Euclidean::Orthonormal::Point &x)
{
	y = y0 + A*(x-x0); return y;
}

unsigned int 
Submersion<0>::size1() const
{
	return n1;

}
unsigned int 
Submersion<0>::size2() const
{
	return n2;
}

const Set::Manifold::Point & 
Submersion<0>::operator () (const Set::Manifold::Point &P)
{
	return operator () ((Set::Euclidean::Orthonormal::Point &)P);
}

//////////////////////////////////////////////////////////////////////
// Class Submersion<1>
//////////////////////////////////////////////////////////////////////

Submersion<1>::Submersion() {}

Submersion<1>::Submersion(const Set::VectorSpace::Hom &B) : A(B) {}

Submersion<1>::Submersion(const Submersion<0> &f) 
: n1(f.size1()), A(f.LinearMapping())
{
	n2 = A.size();
}

Submersion<1>::Submersion(const Submersion<1> &f)
: n1(f.n1), n2(f.n2), A(f.A) {}

Submersion<1>::~Submersion(){}

Submersion<1> & 
Submersion<1>::operator = (const Submersion<1> &f)
{
	assert(f.n1 == n1); assert(f.n2 == n2);
	if (this == &f) return *this; 
	A = f.A; return *this;
}

Set::VectorSpace::Hom 
Submersion<1>::LinearMapping() const
{
	return A;
}

Set::VectorSpace::Hom & 
Submersion<1>::LinearMapping()
{
	return A;
}

Set::Manifold::TMap * 
Submersion<1>::Clone()
{
	return new Submersion<1>(*this);
}
	
Set::Manifold::TMap *
Submersion<1>::Diff()
{
	return new Submersion<2>(*this);
}

void 
Submersion<1>::Randomize()
{
	Random(A);
}

const Set::VectorSpace::Hom &
Submersion<1>::operator () (
	const Set::Euclidean::Orthonormal::Point &)
{
	return A;
}

unsigned int 
Submersion<1>::size1() const
{
	return n1;

}
unsigned int 
Submersion<1>::size2() const
{
	return n2;
}

const Set::VectorSpace::Hom & 
Submersion<1>::operator () (const Set::Manifold::Point &P)
{
	return operator () ((Set::Euclidean::Orthonormal::Point &)P);
}

//////////////////////////////////////////////////////////////////////
// Class Submersion<2>
//////////////////////////////////////////////////////////////////////

Submersion<2>::Submersion() {}

Submersion<2>::Submersion(const Submersion<1> &f) 
: n1(f.size1()), A(f.size2(),f.size1())
{
	n2 = A.size();
}

Submersion<2>::Submersion(const Submersion<2> &f)
: n1(f.n1), n2(f.n2), A(f.A) {}

Submersion<2>::~Submersion(){}

Submersion<2> & 
Submersion<2>::operator = (const Submersion<2> &f)
{
	assert(f.n1 == n1); assert(f.n2 == n2);
	if (this == &f) return *this; 
	A = f.A; return *this;
}

Set::VectorSpace::Hom 
Submersion<2>::LinearMapping() const
{
	return A;
}

Set::VectorSpace::Hom & 
Submersion<2>::LinearMapping()
{
	return A;
}

Set::Manifold::TMap * 
Submersion<2>::Clone()
{
	return new Submersion<2>(*this);
}

void 
Submersion<2>::Randomize() {}

const Set::VectorSpace::HomZero &
Submersion<2>::operator () (
	const Set::Euclidean::Orthonormal::Point &)
{
	return A;
}

unsigned int 
Submersion<2>::size1() const
{
	return n1;

}
unsigned int 
Submersion<2>::size2() const
{
	return n2;
}

const Set::VectorSpace::Hom & 
Submersion<2>::operator () (const Set::Manifold::Point &P)
{
	return operator () ((Set::Euclidean::Orthonormal::Point &)P);
}

}

}

}

Set::Euclidean::Cartesian::Point
operator + (const Set::Euclidean::Cartesian::Point &P,
			const Set::Euclidean::Cartesian::Vector &A)
{
	Set::Euclidean::Cartesian::Point Q = P; 
	Q += A; return Q;
}

Set::Euclidean::Cartesian::Point
operator - (const Set::Euclidean::Cartesian::Point &P,
			const Set::Euclidean::Cartesian::Vector &A)
{
	Set::Euclidean::Cartesian::Point Q = P; 
	Q -= A; return Q;
}

Set::Euclidean::Cartesian::Vector
operator - (const Set::Euclidean::Cartesian::Point &Q,
			const Set::Euclidean::Cartesian::Point &P)
{
	Set::Euclidean::Cartesian::Vector B(Q);
	Set::Euclidean::Cartesian::Vector A(P);
	return B - A;
}

Set::Euclidean::Cartesian::Point
operator + (const Set::Euclidean::Cartesian::Point &P,
			const Set::VectorSpace::Vector &A)
{
	Set::Euclidean::Cartesian::Point Q = P; 
	Q += A; return Q;
}

Set::Euclidean::Cartesian::Point
operator - (const Set::Euclidean::Cartesian::Point &P,
			const Set::VectorSpace::Vector &A)
{
	Set::Euclidean::Cartesian::Point Q = P; 
	Q -= A; return Q;
}

void Random(Set::Euclidean::Cartesian::Point &P){P.Randomize();}
