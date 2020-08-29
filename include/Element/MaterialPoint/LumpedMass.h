// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
////////////////////////////////////////////////////////////////////////////

#if !defined(ELEMENT_MATERIALPOINT_LUMPEDMASS_H__INCLUDED_)
#define ELEMENT_MATERIALPOINT_LUMPEDMASS_H__INCLUDED_

#pragma once

#include <set>
#include <map>
#include "./MaterialPoint.h"
#include "Element/LumpedMass.h"

using namespace std;

namespace Element
{
namespace MaterialPoint
{
class LumpedMass : public Element::LumpedMass
{
public: 
    
	LumpedMass();
	LumpedMass(Data *, const double &);
        LumpedMass(LocalState * pLS, const double &);
	virtual ~LumpedMass();
	LumpedMass(const LumpedMass &);
	virtual LumpedMass & operator = (const LumpedMass &);
	virtual set<Set::Manifold::Point *> GetNodes() const;
	virtual map<Set::Manifold::Point *, double> GetMass() const;

        virtual void GetMass(map<Set::Manifold::Point *, double> &) const;
        
private:

	Data *D;
        vector<double> dm_0;
        
        LocalState * _pLS;
};

}

}

#endif // !defined(ELEMENT_MATERIALPOINT_LUMPEDMASS_H__INCLUDED_
