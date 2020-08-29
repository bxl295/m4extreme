// SparseArray.cpp: implementation the SparseArray class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include "./SparseArray.h"

namespace Geometry
{
namespace Algebraic
{
//////////////////////////////////////////////////////////////////////
// Member methods
//////////////////////////////////////////////////////////////////////

SparseArray::SparseArray() {}

SparseArray::~SparseArray() {}

SparseArray::SparseArray(const SparseArray &A) : 
	map<unsigned int, int>(A) {}

SparseArray & SparseArray::operator = (const SparseArray &A)
{
	if (this == &A) return *this; this->empty();
	map<unsigned int, int>::operator = (A);
	return *this;
}

void 
SparseArray::operator += (const SparseArray &A)
{
	map<unsigned int, int>::const_iterator pA;
	for (pA = A.begin(); pA!=A.end(); pA++)
		(*this)[pA->first] += pA->second;
}

void 
SparseArray::operator -= (const SparseArray &A)
{
	map<unsigned int, int>::const_iterator pA;
	for (pA = A.begin(); pA!=A.end(); pA++)
		(*this)[pA->first] -= pA->second;
}

void 
SparseArray::operator *= (const int &q)
{
	map<unsigned int, int>::iterator pA;
	for (pA = this->begin(); pA!=this->end(); pA++)
		pA->second *= q;
}

void SparseArray::Randomize()
{
	map<unsigned int, int>::iterator pA;
	for (pA = this->begin(); pA!=this->end(); pA++)
		pA->second *= (int)rand();
}

void SparseArray::print(ostream *os)
{
	map<unsigned int, int>::iterator pMap;
	for (pMap=this->begin(); pMap!=this->end(); pMap++)
		*os << pMap->second << ", ";
}

}

}

//////////////////////////////////////////////////////////////////////
// Set operators
//////////////////////////////////////////////////////////////////////

void Random(Geometry::Algebraic::SparseArray &A){A.Randomize();}

//////////////////////////////////////////////////////////////////////
// Algebraic operators
//////////////////////////////////////////////////////////////////////

Geometry::Algebraic::SparseArray 
operator - (const Geometry::Algebraic::SparseArray &A)
{
	Geometry::Algebraic::SparseArray B=A;
	map<unsigned int, int>::iterator pB;
	for (pB = B.begin(); pB!=B.end(); pB++)
		pB->second = -(pB->second);
	return B;
}

Geometry::Algebraic::SparseArray 
operator * (const int &q, const Geometry::Algebraic::SparseArray &A)
{
	Geometry::Algebraic::SparseArray B=A;
	map<unsigned int, int>::iterator pB;
	for (pB = B.begin(); pB!=B.end(); pB++)
		pB->second = q*pB->second;
	return B;
}

Geometry::Algebraic::SparseArray 
operator * (const Geometry::Algebraic::SparseArray &A, const int &q)
{
	Geometry::Algebraic::SparseArray B=A;
	map<unsigned int, int>::iterator pB;
	for (pB = B.begin(); pB!=B.end(); pB++)
		pB->second = q*pB->second;
	return B;
}

//////////////////////////////////////////////////////////////////////
// Printing
//////////////////////////////////////////////////////////////////////

ostream & 
operator<<(ostream &os, Geometry::Algebraic::SparseArray &A)
{
	A.print(&os); return os;
}
