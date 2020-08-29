// Search.h: interface for the Search class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(GEOMETRY_SEARCH_H__INCLUDED_)
#define GEOMETRY_SEARCH_H__INCLUDED_

#pragma once

#include <set>
#include "../Algebraic/AlgLib.h"
#include "../../Set/Manifold/Manifold.h"
#include "../../Set/Manifold/Euclidean/Orthonormal/Orthonormal.h"

using namespace std;

namespace Geometry
{
//////////////////////////////////////////////////////////////////////
// Class Search
//////////////////////////////////////////////////////////////////////

template<typename T>
class Search
{
public: 
	typedef Set::Euclidean::Orthonormal::Point  point_type;
        typedef T Record;

	Search(){}
	virtual ~Search(){}

	// return a vector of records which are in the search range of an arbitrary location
	virtual void operator () (const point_type &, double, vector<Record> &) = 0;

	// return a set of records which are in the search range of an input record
        virtual const set<Record> * operator () (Record) const = 0;
};

}

#endif // !defined(GEOMETRY_SEARCH_H__INCLUDED_
