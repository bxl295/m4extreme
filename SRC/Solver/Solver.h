// Solver.h: interface for the Solver class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(SOLVER__INCLUDED_)
#define SOLVER__INCLUDED_

#pragma once

namespace Solver
{
class Propagator 
{
public:

	Propagator(){}
	virtual ~Propagator(){}
	virtual void operator ++ ()=0;

private:

	Propagator(Propagator &);
	void operator=(Propagator &);
};

}

#endif // !defined(SOLVER__INCLUDED_)
