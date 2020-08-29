// Vector.h: implementation of the Vector class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <cassert>
#include "Vector.h"

using namespace std;

namespace Set
{
namespace VectorSpace
{
//////////////////////////////////////////////////////////////////////
// Vector member methods
//////////////////////////////////////////////////////////////////////
Vector::Vector(const VectorZero & A):Array(A) {}
void Vector::operator += (const VectorZero &O){assert(n == O.n); }
void Vector::operator -= (const VectorZero &O){assert(n == O.n); }

double Vector::operator () (const VectorZero &O) const
{
	assert(n == O.n); return 0.0;
}

//////////////////////////////////////////////////////////////////////
// VectorZero member methods
//////////////////////////////////////////////////////////////////////

VectorZero::VectorZero() : Vector(){}

VectorZero::VectorZero(const unsigned int &n) : Vector(n){}

VectorZero::VectorZero(const unsigned int &n, double * const &u) : Vector(n,u){}

VectorZero::~VectorZero(){}

VectorZero::VectorZero(const Vector &A) : Vector(A){}

VectorZero::VectorZero(const VectorZero &A) : Vector(A){}

VectorZero & VectorZero::operator = (const VectorZero &){return *this;}

void VectorZero::operator += (const VectorZero &O){assert(n == O.n);}

void VectorZero::operator -= (const VectorZero &O){assert(n == O.n);}

void VectorZero::operator *= (const double &){}

void VectorZero::operator /= (const double &){}

double VectorZero::operator () (const Vector &A) const
{
	assert(n == A.size()); return 0.0;
}

}

}

//////////////////////////////////////////////////////////////////////
// Additive group operators
//////////////////////////////////////////////////////////////////////
Set::VectorSpace::VectorZero 
operator + (const Set::VectorSpace::VectorZero &A, 
			const Set::VectorSpace::VectorZero &B)
{
	return A;
}

Set::VectorSpace::Vector 
operator + (const Set::VectorSpace::Vector &A, 
			const Set::VectorSpace::VectorZero &B)
{
	return A;
}

Set::VectorSpace::Vector 
operator + (const Set::VectorSpace::VectorZero &A, 
			const Set::VectorSpace::Vector &B)
{
	return B;
}

Set::VectorSpace::VectorZero 
operator - (const Set::VectorSpace::VectorZero &A, 
			const Set::VectorSpace::VectorZero &B)
{
	return A;
}

Set::VectorSpace::Vector 
operator - (const Set::VectorSpace::Vector &A, 
			const Set::VectorSpace::VectorZero &B)
{
	return A;
}

Set::VectorSpace::Vector 
operator - (const Set::VectorSpace::VectorZero &A, 
			const Set::VectorSpace::Vector &B)
{
	return -B;
}

Set::VectorSpace::VectorZero 
operator - (const Set::VectorSpace::VectorZero &A)
{
	return A;
}

//////////////////////////////////////////////////////////////////////
// Vector space operators
//////////////////////////////////////////////////////////////////////
Set::VectorSpace::VectorZero 
operator * (const Set::VectorSpace::VectorZero &A, 
			const double &p)
{
	return A;
}

Set::VectorSpace::VectorZero 
operator * (const double &p, 
			const Set::VectorSpace::VectorZero &A)
{
	return A;
}

Set::VectorSpace::VectorZero 
operator / (const Set::VectorSpace::VectorZero &A, 
			const double &p)
{
	return A;
}

//////////////////////////////////////////////////////////////////////
// Normed space operators
//////////////////////////////////////////////////////////////////////
double Norm(const Set::VectorSpace::VectorZero &A)
{
	return 0.0;
}
