// Simplex.h: interface for the Simplex class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(GEOMETRY_UTILS_SIMPLEX_H__INCLUDED_)
#define GEOMETRY_UTILS_SIMPLEX_H__INCLUDED_

#pragma once

#include <cassert>
#include <vector>
#include "../../../Set/Algebraic/AlgLib.h"

namespace Geometry
{
namespace Utils
{
class Simplex
{
public:

	typedef Set::VectorSpace::Vector point_type;

	Simplex();
	Simplex(const vector<point_type> &);
	virtual ~Simplex();
	double CircumRadius();
	point_type CircumCenter();
	void CircumSphere(point_type &, double &);
	double Volume();
	point_type BaryCenter();
	vector<double> Lambda(const point_type &);
	vector<Set::VectorSpace::Vector> DLambda();

private:

	unsigned int n;
	unsigned int npp;
	vector<point_type> x;

	unsigned int Factorial(const unsigned int &);

private:

	Simplex(const Simplex &);
	Simplex & operator = (const Simplex &);	
};

}

}

#endif // !defined(GEOMETRY_UTILS_SIMPLEX_H__INCLUDED_)
