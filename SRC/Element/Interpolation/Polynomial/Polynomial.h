// Polynomial.h: interface for the Polynomial class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(ELEMENT_INTERPOLATION_POLYNOMIAL_H__INCLUDED_)
#define ELEMENT_INTERPOLATION_POLYNOMIAL_H__INCLUDED_

#pragma once

#include <map>
#include <vector>
#include "../Interpolation.h"
#include "../Utils/Polynomial/Polynomial.h"
#include "../Utils/MultiIndex/MultiIndex.h"

using namespace std;

namespace Element
{
namespace Interpolation
{
namespace Polynomial
{
typedef Element::Interpolation::Utils::MultiIndex multiindex_type;
typedef Element::Interpolation::Utils::Polynomial polynomial_type;
typedef Set::Manifold::Point point_type;
typedef Set::VectorSpace::Vector vector_type;

//////////////////////////////////////////////////////////////////////
// Class Shape<p>
//////////////////////////////////////////////////////////////////////

template<unsigned int p> class Shape;

//////////////////////////////////////////////////////////////////////
// Class Shape<0>
//////////////////////////////////////////////////////////////////////

template<> 
class Shape<0> : public Element::Interpolation::Shape<0>
{
friend class Shape<1>;

public: 

	Shape();
	Shape(const set<point_type *> &);
	Shape(const map<point_type *, multiindex_type> &);
	virtual ~Shape();
	map<point_type *, double> 
		operator () (const vector_type &) const;

private:

	unsigned int n;
	map<point_type *, polynomial_type> w;

private:

	Shape(const Shape<0> &);
	Shape<0> & operator = (const Shape<0> &);
};

//////////////////////////////////////////////////////////////////////
// Class Shape<1>
//////////////////////////////////////////////////////////////////////

template<> 
class Shape<1> : public Element::Interpolation::Shape<1>
{
public: 

	Shape();
	Shape(const Shape<0> &);
	virtual ~Shape();
	map<point_type *, vector_type> 
		operator () (const vector_type &) const;

private:

	unsigned int n;
	map<point_type *, vector<polynomial_type> > Dw;

private:

	Shape(const Shape<1> &);
	Shape<1> & operator = (const Shape<1> &);
};

//////////////////////////////////////////////////////////////////////
// Utilities
//////////////////////////////////////////////////////////////////////

vector_type RandomPoint(const unsigned int &);

map<point_type *, vector_type>
NodalCoordinates(const set<point_type *> &);

map<point_type *, vector_type>
NodalCoordinates(const map<point_type *, multiindex_type> &);

set<multiindex_type>
IndexSetReduced(const unsigned int &, const unsigned int &);

set<multiindex_type>
IndexSetInteriorReduced(const unsigned int &, const unsigned int &);

}

}

}

#endif // !defined(ELEMENT_INTERPOLATION_POLYNOMIAL_H__INCLUDED_

