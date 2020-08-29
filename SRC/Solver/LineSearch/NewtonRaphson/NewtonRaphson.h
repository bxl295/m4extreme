// NewtonRaphson.h: interface for the NewtonRaphson class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(SOLVER_LINESEARCH_NEWTONRAPHSON__INCLUDED_)
#define SOLVER_LINESEARCH_NEWTONRAPHSON__INCLUDED_

#pragma once

#include <set>
#include "../LineSearch.h"
#include "../../../Model/Model.h"
#include "../../../Model/Utils/Operators/Operators.h"
#include "../../../Set/Manifold/Manifold.h"

using namespace std;

namespace Solver
{
namespace LineSearch
{
class NewtonRaphson : public Propagator  
{
public:

	NewtonRaphson();
	virtual ~NewtonRaphson();
	NewtonRaphson(
		Model::LocalState *, 
		Model::Jet<1> *, 
		set<Set::Manifold::Point *> *);
	NewtonRaphson(
		Model::LocalState *, 
		Model::Energy<1> *, 
		Model::Energy<2> *, 
		set<Set::Manifold::Point *> *);

	bool & SetPrint();
	unsigned int & SetNitMax();
	double & SetTol();
	void operator () (const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &);

private:

	bool Print;
	unsigned int NitMax;
	double Tol;
	Model::LocalState *LS;
	Model::Jet<1> *DJ;
	Model::Energy<1> *DE;
	Model::Energy<2> *DDE;
	set<Set::Manifold::Point *> *x;

private:

	NewtonRaphson(NewtonRaphson &);
	void operator=(NewtonRaphson &);
};

}

}

#endif // !defined(SOLVER_LINESEARCH_NEWTONRAPHSON__INCLUDED_)
