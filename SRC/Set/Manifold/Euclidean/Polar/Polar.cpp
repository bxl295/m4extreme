// Polar.cpp: Implementation of the Polar class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <cassert>
#include "./Polar.h"

namespace Set
{
namespace Euclidean
{	
namespace Polar
{
//////////////////////////////////////////////////////////////////////
// Class Point
//////////////////////////////////////////////////////////////////////

Point::Point()  : Set::Array(2) {}

Point::Point(double * const &u) : Set::Array(2,u) {}

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
	double &r = *head; double &t = *(head+1);
	r = (double)rand()/(double)RAND_MAX;
	t = PI2*(double)rand()/(double)RAND_MAX;
}

unsigned int 
Point::size() const
{
	return 2;
}

void 
Point::print(ostream *os)
{
	Set::Array::print(os);
}

void 
Point::operator += (const Vector &A)
{
	double &r = *head; double &t = *(head+1);
	double &dr = *A.begin(); double &dt = *(A.begin()+1);
	double x, y, dx, dy;
	double c=cos(t), s=sin(t);
	dx = dr*c - r*s*dt;
	dy = dr*s + r*c*dt;
	x = r*c + dx;
	y = r*s + dy;
	r = sqrt(x*x + y*y);
	if (y >= 0.0){t = acos(x/r);}
	else   {t = PI2 - acos(x/r);}
}

void 
Point::operator -= (const Vector &A)
{
	double &r = *head; double &t = *(head+1);
	double &dr = *A.begin(); double &dt = *(A.begin()+1);
	double x, y, dx, dy;
	double c=cos(t), s=sin(t);
	dx = dr*c - r*s*dt;
	dy = dr*s + r*c*dt;
	x = r*c - dx;
	y = r*s - dy;
	r = sqrt(x*x + y*y);
	if (y >= 0.0){t = acos(x/r);}
	else   {t = PI2 - acos(x/r);}
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

Vector::Vector() : Set::VectorSpace::Vector(2) {}

Vector::Vector(double * const &u) : Set::VectorSpace::Vector(2,u) {}

Set::Manifold::Point *Vector::Clone() const
{
	return new Vector(*this);
}

Vector::~Vector() {}

Vector::Vector(const Vector &A) : 
	Set::VectorSpace::Vector(A) {}

Vector::Vector(const Set::VectorSpace::Vector &A) : 
	Set::VectorSpace::Vector(A) {}

Vector::Vector(const Covector &A) : 
	Set::VectorSpace::Vector(A) {}

Vector::Vector(const Set::Euclidean::Polar::Point &A) : 
	Set::VectorSpace::Vector(A) {}

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
Vector::operator = (const Set::Euclidean::Polar::Point &A)
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
	return 2;
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

VectorZero::VectorZero() : 
	Set::VectorSpace::VectorZero(2) {}

VectorZero::VectorZero(double * const &u) : 
	Set::VectorSpace::VectorZero(2) {}

VectorZero::~VectorZero(){}

VectorZero::VectorZero(const VectorZero &A) : 
	Set::VectorSpace::VectorZero(A) {}

VectorZero::VectorZero(const Set::VectorSpace::VectorZero &A) : 
	Set::VectorSpace::VectorZero(A) {}

VectorZero::VectorZero(const CovectorZero &A) : 
	Set::VectorSpace::VectorZero(A) {}

VectorZero & 
VectorZero::operator = (const VectorZero &A)
{
	assert(2 == A.n); 
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

Covector::~Covector() {}

Covector::Covector() : 
	Set::VectorSpace::Vector(2) {}

Covector::Covector(double * const &u) : 
	Set::VectorSpace::Vector(2,u) {}

Covector::Covector(const Covector &A) : 
	Set::VectorSpace::Vector(A) {}

Covector::Covector(const Set::VectorSpace::Vector &A) : 
	Set::VectorSpace::Vector(A) {}

Covector::Covector(const Set::Euclidean::Polar::Vector &A) : 
	Set::VectorSpace::Vector(A) {}

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

CovectorZero::CovectorZero() : 
	Set::VectorSpace::VectorZero(2) {}

CovectorZero::CovectorZero(double * const &u) : 
	Set::VectorSpace::VectorZero(2) {}

CovectorZero::~CovectorZero(){}

CovectorZero::CovectorZero(const Set::VectorSpace::VectorZero &A) : 
	Set::VectorSpace::VectorZero(A) {}

CovectorZero::CovectorZero(const CovectorZero &A) : 
	Set::VectorSpace::VectorZero(A) {}

CovectorZero::CovectorZero(const Set::Euclidean::Polar::VectorZero &A) : 
	Set::VectorSpace::VectorZero(A) {}

CovectorZero & 
CovectorZero::operator = (const CovectorZero &A)
{
	assert(2 == A.n); 
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

Embedding<0>::Embedding() : y(2) {}

Embedding<0>::Embedding(const Embedding<0> &f) : y(f.y) {}

Embedding<0>::~Embedding() {}

Embedding<0> & 
Embedding<0>::operator = (const Embedding<0> &f)
{
	if (this == &f) return *this; 
	y = f.y; return *this;
}

Set::Manifold::Map * 
Embedding<0>::Clone()
{
	return new Embedding<0>;
}
	
Set::Manifold::TMap *
Embedding<0>::Diff()
{
	return new Embedding<1>;
}

void 
Embedding<0>::Randomize() {}

const Set::Euclidean::Orthonormal::Point &
Embedding<0>::operator () (const Set::Euclidean::Polar::Point &x)
{
	y[0] = x[0]*cos(x[1]); 
	y[1] = x[0]*sin(x[1]);
	return y;
}

unsigned int 
Embedding<0>::size1() const
{
	return 2;

}
unsigned int 
Embedding<0>::size2() const
{
	return 2;
}

const Set::Manifold::Point & 
Embedding<0>::operator () (const Set::Manifold::Point &P)
{
	return operator () ((Point &)P);
}

//////////////////////////////////////////////////////////////////////
// Class Embedding<1>
//////////////////////////////////////////////////////////////////////

Embedding<1>::Embedding() : A(2) {}

Embedding<1>::Embedding(const Embedding<0> &f) : A(2) {}

Embedding<1>::Embedding(const Embedding<1> &f) : A(f.A) {}

Embedding<1>::~Embedding(){}

Embedding<1> & 
Embedding<1>::operator = (const Embedding<1> &f)
{
	if (this == &f) return *this; 
	A = f.A; return *this;
}

Set::Manifold::TMap * 
Embedding<1>::Clone()
{
	return new Embedding<1>;
}
	
Set::Manifold::TMap *
Embedding<1>::Diff()
{
	return new Embedding<2>;
}

void 
Embedding<1>::Randomize() {}

const Set::VectorSpace::Hom &
Embedding<1>::operator () (const Set::Euclidean::Polar::Point &x)
{
	double c = cos(x[1]); 
	double s = sin(x[1]); 
	A[0][0] =   c; 
	A[0][1] =   s;
	A[1][0] = - x[0]*s; 
	A[1][1] =   x[0]*c;
	return A;
}

unsigned int 
Embedding<1>::size1() const
{
	return 2;

}
unsigned int 
Embedding<1>::size2() const
{
	return 4;
}

const Set::VectorSpace::Hom & 
Embedding<1>::operator () (const Set::Manifold::Point &P)
{
	return operator () ((Point &)P);
}

//////////////////////////////////////////////////////////////////////
// Class Embedding<2>
//////////////////////////////////////////////////////////////////////

Embedding<2>::Embedding() : A(4,2) {}

Embedding<2>::Embedding(const Embedding<1> &f) : A(4,2) {}

Embedding<2>::Embedding(const Embedding<2> &f) : A(f.A) {}

Embedding<2>::~Embedding(){}

Embedding<2> & 
Embedding<2>::operator = (const Embedding<2> &f)
{
	if (this == &f) return *this; 
	A = f.A; return *this;
}

Set::Manifold::TMap * 
Embedding<2>::Clone()
{
	return new Embedding<2>;
}

void 
Embedding<2>::Randomize() {}

const Set::VectorSpace::Hom &
Embedding<2>::operator () (const Set::Euclidean::Polar::Point &x)
{
	double c = cos(x[1]); 
	double s = sin(x[1]); 
	Set::VectorSpace::Hom A0(2,A[0].begin());
	Set::VectorSpace::Hom A1(2,A[1].begin());
	A0[0][0] =   0.0; 
	A0[0][1] =   0.0;
	A0[1][0] = - s; 
	A0[1][1] =   c;
	A1[0][0] = - s; 
	A1[0][1] =   c;
	A1[1][0] = - x[0]*c; 
	A1[1][1] = - x[0]*s;
	return A;
}

unsigned int 
Embedding<2>::size1() const
{
	return 2;

}
unsigned int 
Embedding<2>::size2() const
{
	return 8;
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

Submersion<0>::Submersion(const Submersion<0> &f) : x(f.x) {}

Submersion<0>::~Submersion(){}

Submersion<0> & 
Submersion<0>::operator = (const Submersion<0> &f)
{
	if (this == &f) return *this;
	x = f.x; return *this;
}

Set::Manifold::Map * 
Submersion<0>::Clone()
{
	return new Submersion<0>(*this);
}
	
Set::Manifold::TMap *
Submersion<0>::Diff()
{
	return new Submersion<1>;
}

void 
Submersion<0>::Randomize() {}

const Set::Euclidean::Polar::Point &
Submersion<0>::operator () (const Set::Euclidean::Orthonormal::Point &y)
{
	assert(y.size() == 2); 
	double r = sqrt(y[0]*y[0] + y[1]*y[1]); x[0] = r;
	if (y[1] >= 0.0){x[1] =       acos(y[0]/r);}
	else            {x[1] = PI2 - acos(y[0]/r);}
	return x;
}

unsigned int 
Submersion<0>::size1() const
{
	return 2;

}
unsigned int 
Submersion<0>::size2() const
{
	return 2;
}

const Set::Manifold::Point & 
Submersion<0>::operator () (const Set::Manifold::Point &P)
{
	return operator () ((Set::Euclidean::Orthonormal::Point &)P);
}

//////////////////////////////////////////////////////////////////////
// Class Submersion<1>
//////////////////////////////////////////////////////////////////////

Submersion<1>::Submersion() : A(2) {}

Submersion<1>::Submersion(const Submersion<0> &f) : A(2) {}

Submersion<1>::Submersion(const Submersion<1> &f) : A(f.A) {}

Submersion<1>::~Submersion(){}

Submersion<1> & 
Submersion<1>::operator = (const Submersion<1> &f)
{
	if (this == &f) return *this; 
	A = f.A; return *this;
}

Set::Manifold::TMap * 
Submersion<1>::Clone()
{
	return new Submersion<1>;
}
	
Set::Manifold::TMap *
Submersion<1>::Diff()
{
	return new Submersion<2>;
}

void 
Submersion<1>::Randomize() {}

const Set::VectorSpace::Hom &
Submersion<1>::operator () (const Set::Euclidean::Orthonormal::Point &y)
{
	assert(y.size() == 2); 
	double r2 = y[0]*y[0] + y[1]*y[1];
	double r = sqrt(r2); 
	A[0][0] =   y[0]/r;
	A[0][1] = - y[1]/r2;
	A[1][0] =   y[1]/r;
	A[1][1] =   y[0]/r2;
	return A;
}

unsigned int 
Submersion<1>::size1() const
{
	return 2;

}
unsigned int 
Submersion<1>::size2() const
{
	return 4;
}

const Set::VectorSpace::Hom & 
Submersion<1>::operator () (const Set::Manifold::Point &P)
{
	return operator () ((Set::Euclidean::Orthonormal::Point &)P);
}

//////////////////////////////////////////////////////////////////////
// Class Submersion<2>
//////////////////////////////////////////////////////////////////////

Submersion<2>::Submersion() : A(4,2) {}

Submersion<2>::Submersion(const Submersion<1> &f) : A(4,2) {}

Submersion<2>::Submersion(const Submersion<2> &f) : A(f.A) {}

Submersion<2>::~Submersion(){}

Submersion<2> & 
Submersion<2>::operator = (const Submersion<2> &f)
{
	if (this == &f) return *this; 
	A = f.A; return *this;
}

Set::Manifold::TMap * 
Submersion<2>::Clone()
{
	return new Submersion<2>(*this);
}

void 
Submersion<2>::Randomize() {}

const Set::VectorSpace::Hom &
Submersion<2>::operator () (const Set::Euclidean::Orthonormal::Point &y)
{
	assert(y.size() == 2); 
	double r2 = y[0]*y[0] + y[1]*y[1];
	double r = sqrt(r2); 
	Set::VectorSpace::Hom A0(2,A[0].begin());
	Set::VectorSpace::Hom A1(2,A[1].begin());
	double c = y[0]/r;
	double s = y[1]/r;
	double cc = c*c;
	double cs = c*s;
	double ss = s*s;
	A0[0][0] =   ss/r;
	A0[0][1] =   2.0*cs/r2;
	A0[1][0] = - cs/r;
	A0[1][1] =   (ss - cc)/r2;
	A1[0][0] = - cs/r;
	A1[0][1] =   (ss - cc)/r2;
	A1[1][0] =   cc/r;
	A1[1][1] = - 2.0*cs/r2;
	return A;
}

unsigned int 
Submersion<2>::size1() const
{
	return 2;

}
unsigned int 
Submersion<2>::size2() const
{
	return 8;
}

const Set::VectorSpace::Hom & 
Submersion<2>::operator () (const Set::Manifold::Point &P)
{
	return operator () ((Set::Euclidean::Orthonormal::Point &)P);
}

}

}

}

Set::Euclidean::Polar::Point
operator + (const Set::Euclidean::Polar::Point &P,
			const Set::Euclidean::Polar::Vector &A)
{
	Set::Euclidean::Polar::Point Q = P; 
	Q += A; return Q;
}

Set::Euclidean::Polar::Point
operator - (const Set::Euclidean::Polar::Point &P,
			const Set::Euclidean::Polar::Vector &A)
{
	Set::Euclidean::Polar::Point Q = P; 
	Q -= A; return Q;
}

Set::Euclidean::Polar::Vector
operator - (const Set::Euclidean::Polar::Point &P1,
			const Set::Euclidean::Polar::Point &P2)
{
	Set::Euclidean::Polar::Vector A;
	double &dr = *A.begin();  double &dt = *(A.begin()+1);
	double &r1 = *P1.begin(); double &t1 = *(P1.begin()+1);
	double &r2 = *P2.begin(); double &t2 = *(P2.begin()+1);
	double c1 = cos(t1), s1 = sin(t1);
	double c2 = cos(t2), s2 = sin(t2);
	double x1 = r1*c1, y1 = r1*s1;
	double x2 = r2*c2, y2 = r2*s2;
	double dx = x1 - x2, dy = y1 - y2;
	dr =    dx*c2 + s2*dy;
	dt = (- dx*s2 + c2*dy)/r2;
	return A;
}

Set::Euclidean::Polar::Point
operator + (const Set::Euclidean::Polar::Point &P,
			const Set::VectorSpace::Vector &A)
{
	Set::Euclidean::Polar::Point Q = P; 
	Q += A; return Q;
}

Set::Euclidean::Polar::Point
operator - (const Set::Euclidean::Polar::Point &P,
			const Set::VectorSpace::Vector &A)
{
	Set::Euclidean::Polar::Point Q = P; 
	Q -= A; return Q;
}

void Random(Set::Euclidean::Polar::Point &P){P.Randomize();}
