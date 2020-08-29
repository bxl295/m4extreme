// LumpedMass.h: interface for the LumpedMass class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(ELEMENT_CONFORMING_LUMPEDMASS_H__INCLUDED_)
#define ELEMENT_CONFORMING_LUMPEDMASS_H__INCLUDED_

#pragma once

#include <set>
#include <map>
#include <vector>
#include "../LumpedMass.h"

using namespace std;

namespace Element
{
namespace Conforming
{
class LumpedMass : public Element::LumpedMass
{
public: 

	LumpedMass();
	LumpedMass(const double &);
	LumpedMass(
		const double &,
		const vector<map<Set::Manifold::Point *, double> > &,
		const vector<double> &);
	virtual ~LumpedMass();
	LumpedMass(const LumpedMass &);
	virtual LumpedMass & operator = (const LumpedMass &);
	virtual set<Set::Manifold::Point *> GetNodes() const;
	virtual map<Set::Manifold::Point *, double> GetMass() const;

        virtual void GetMass(map<Set::Manifold::Point *, double> &) const;

protected:

	double rho;
	vector<map<Set::Manifold::Point *, double> > N;
	vector<double> QW;
};

}

}

#endif // !defined(ELEMENT_CONFORMING_LUMPEDMASS_H__INCLUDED_
