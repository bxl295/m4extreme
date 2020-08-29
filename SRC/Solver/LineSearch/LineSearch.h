// LineSearch.h: interface for the LineSearch class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(SOLVER_LINESEARCH_INCLUDED_)
#define SOLVER_LINESEARCH_INCLUDED_

#pragma once

#include "../../Model/Model.h"
#include "../../Set/Manifold/Manifold.h"

namespace Solver
{
namespace LineSearch
{
class Propagator 
{
public:

	Propagator(){}
	virtual ~Propagator(){}
	virtual void operator () (
		const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &)=0;

private:

	Propagator(Propagator &);
	void operator=(Propagator &);
};

}

}

#endif // !defined(SOLVER_LINESEARCH_INCLUDED_)
