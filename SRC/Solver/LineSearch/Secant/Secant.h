// Secant.h: interface for the Secant class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(SOLVER_LINESEARCH_SECANT__INCLUDED_)
#define SOLVER_LINESEARCH_SECANT__INCLUDED_

#pragma once

#include <map>
#include "../LineSearch.h"
#include "../../../Model/Model.h"
#include "../../../Model/Utils/Operators/Operators.h"
#include "../../../Set/Manifold/Manifold.h"

using namespace std;

namespace Solver
{
namespace LineSearch
{
class Secant : public Propagator  
{
public:

	Secant();
	virtual ~Secant();
	Secant(
		Model::LocalState *, 
		Model::Energy<1> *, 
		set<Set::Manifold::Point *> *);

	bool & SetPrint();
	unsigned int & SetNitMax();
	double & SetTol();
	double & SetStep();
	void operator () (
		const map<Set::Manifold::Point *, 
		Set::VectorSpace::Vector> &);

private:

	bool Print;
	unsigned int NitMax;
	double Tol; double Step;
	Model::LocalState *LS;
	Model::Energy<1> *DE;
	set<Set::Manifold::Point *> *x;

private:

	Secant(Secant &);
	void operator=(Secant &);
};

}

}

#endif // !defined(SOLVER_LINESEARCH_SECANT__INCLUDED_)
