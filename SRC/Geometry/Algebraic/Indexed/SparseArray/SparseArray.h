// SparseArray.h: interface for the SparseArray class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(GEOMETRY_ALGEBRAIC_SPARSEARRAY_H__INCLUDED_)
#define GEOMETRY_ALGEBRAIC_SPARSEARRAY_H__INCLUDED_

#pragma once
#include <stdlib.h>
#include <map>
#include <iostream>

using namespace std;

namespace Geometry
{
namespace Algebraic
{
//////////////////////////////////////////////////////////////////////
// Class SparseArray
//////////////////////////////////////////////////////////////////////

class SparseArray : public map<unsigned int, int>
{
public:

	SparseArray();
	virtual ~SparseArray();
	SparseArray(const SparseArray &);
	SparseArray & operator = (const SparseArray &);
	void operator += (const SparseArray &);
	void operator -= (const SparseArray &);
	void operator *= (const int &);
	void Randomize();
	void print(ostream *);
};

}

}

Geometry::Algebraic::SparseArray operator - (const Geometry::Algebraic::SparseArray &);
Geometry::Algebraic::SparseArray operator * (const int &, const Geometry::Algebraic::SparseArray &);
Geometry::Algebraic::SparseArray operator * (const Geometry::Algebraic::SparseArray &, const int &);
void Random(Geometry::Algebraic::SparseArray &);
ostream & operator<<(ostream &, Geometry::Algebraic::SparseArray &);

#endif // !defined(GEOMETRY_ALGEBRAIC_SPARSEARRAY_H__INCLUDED_)
