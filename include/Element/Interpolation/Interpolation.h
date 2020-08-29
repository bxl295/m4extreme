// Interpolation.h: interface for the Interpolation class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(ELEMENT_INTERPOLATION_H__INCLUDED_)
#define ELEMENT_INTERPOLATION_H__INCLUDED_

#pragma once

#include <map>
#include "../../Set/Manifold/Manifold.h"
#include "../../Set/Algebraic/VectorSpace/VectorSpace.h"
#include "../../Set/Algebraic/VectorSpace/Category/Category.h"

using namespace std;

namespace Element
{
namespace Interpolation
{
//////////////////////////////////////////////////////////////////////
// Class Shape<p>
//////////////////////////////////////////////////////////////////////

template<unsigned int p> class Shape
{
public: 

	Shape() {}
	virtual ~Shape() {}
	virtual map<Set::Manifold::Point *, Set::VectorSpace::Hom>
		operator () (const Set::VectorSpace::Vector &) const=0;
};

//////////////////////////////////////////////////////////////////////
// Class Shape<0>
//////////////////////////////////////////////////////////////////////

template<> class Shape<0>
{
public: 

	Shape();
	virtual ~Shape();
	virtual map<Set::Manifold::Point *, double> 
		operator () (const Set::VectorSpace::Vector &) const=0;
};

//////////////////////////////////////////////////////////////////////
// Class Shape<1>
//////////////////////////////////////////////////////////////////////

template<> class Shape<1>
{
public: 

	Shape();
	virtual ~Shape();
	virtual map<Set::Manifold::Point *, Set::VectorSpace::Vector>
		operator () (const Set::VectorSpace::Vector &) const=0;
};

}

}

#endif // !defined(ELEMENT_INTERPOLATION_H__INCLUDED_)
