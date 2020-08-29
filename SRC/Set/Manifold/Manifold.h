// Manifold.h: interface for the Manifold class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(SET_MANIFOLD_H__INCLUDED_)
#define SET_MANIFOLD_H__INCLUDED_

#pragma once

#include <iostream>
#include "../Algebraic/VectorSpace/Vector/Vector.h"

using namespace std;

namespace Set
{
namespace Manifold
{
//////////////////////////////////////////////////////////////////////
// Point class
//////////////////////////////////////////////////////////////////////

class Point
{
public: 

	Point() {}
	virtual Point *Clone() const=0;
	virtual ~Point() {}
	virtual Point & operator = (const Point &)=0;
	virtual bool operator != (const Point &) const=0;
	virtual bool operator == (const Point &) const=0;
	virtual void Randomize()=0;
	virtual void operator += (const Set::VectorSpace::Vector &)=0;
	virtual void operator -= (const Set::VectorSpace::Vector &)=0;
	virtual unsigned int size() const=0;
	virtual void print(ostream *)=0;
};

}

}

#endif // !defined(SET_MANIFOLD_H__INCLUDED_)
