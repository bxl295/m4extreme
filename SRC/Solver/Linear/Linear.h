// Linear.h: interface for the Linear class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(SOLVER_LINEAR__INCLUDED_)
#define SOLVER_LINEAR__INCLUDED_

#pragma once

#include <map>
#include <set>

using namespace std;

namespace Solver
{
namespace Linear
{
template <class key> class System
{
public:

	System(){}
	virtual ~System(){}

	virtual double Get(key)=0;
	virtual double Get(key, key)=0;
	virtual void Set(key, double)=0;
	virtual void Set(key, key, double)=0;

	virtual void SetToZero1()=0;
	virtual void SetToZero2()=0;
	virtual void Add(key, double)=0;
	virtual void Add(key, key, double)=0;
	virtual void Solve()=0;
	virtual double Norm()=0;

private:

	System(System &);
	void operator=(System &);
};

}

}

#endif // !defined(SOLVER_LINEAR__INCLUDED_)
