// LumpedMass.h: interface for the LumpedMass class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(ELEMENT_LUMPEDMASS_H__INCLUDED_)
#define ELEMENT_LUMPEDMASS_H__INCLUDED_

#pragma once

#include <set>
#include <map>
#include <cassert>
#include "../Set/Manifold/Manifold.h"

using namespace std;

namespace Element
{
class LumpedMass
{
public: 

	LumpedMass(){}
	virtual ~LumpedMass(){}
	virtual set<Set::Manifold::Point *> GetNodes() const=0;
	virtual map<Set::Manifold::Point *, double> GetMass() const=0;

        virtual void write(ostream & os) const {}
	virtual void read(istream & is) {}

        virtual void GetMass(map<Set::Manifold::Point *, double> &) const { assert(false); }
};

}

#endif // !defined(ELEMENT_LUMPEDMASS_H__INCLUDED_
