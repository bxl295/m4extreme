// NewtonRaphson.h: interface for the NewtonRaphson class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(SOLVER_NEWTONRAPHSON__INCLUDED_)
#define SOLVER_NEWTONRAPHSON__INCLUDED_

#pragma once

#include <set>
#include <map>
#include "Solver/Solver.h"
#include "Solver/Linear/Linear.h"
#include "Solver/LineSearch/LineSearch.h"
#include "Model/Model.h"
#include "Model/Static/Static.h"
#include "Set/Manifold/Manifold.h"
#include "Solver/LineSearch/Secant/Secant.h"

using namespace std;

enum UPDATE_TYPE {DONOT_UPDATE = 0, DO_UPDATE = 1};

namespace Solver
{
class NewtonRaphson : public Propagator  
{
public:

	NewtonRaphson();
	virtual ~NewtonRaphson();
	NewtonRaphson(
		Model::LocalState *, 
		Model::Jet<1> *, 
		Solver::Linear::System<DOF> *, 
		set<Set::Manifold::Point *> *, 
		UPDATE_TYPE = DO_UPDATE);
        NewtonRaphson(
		Model::LocalState *, 
                Model::Energy<1> *, 
		Model::Jet<1> *, 
		Solver::Linear::System<DOF> *, 
		set<Set::Manifold::Point *> *, 
		UPDATE_TYPE = DO_UPDATE);
	NewtonRaphson(
		Model::LocalState *, 
		Model::Energy<1> *, 
		Model::Energy<2> *, 
		Solver::Linear::System<DOF> *, 
		set<Set::Manifold::Point *> *, 
		UPDATE_TYPE = DO_UPDATE);

	bool & SetPrint();
	unsigned int & SetNitMax();
	double & SetTol();
	void operator ++ ();

private:

	bool Print;
	UPDATE_TYPE type;
	unsigned int NitMax;
	double Tol;
	Model::LocalState *LS;
	Model::Jet<1> *J;
	Model::Energy<1> *DE;
	Model::Energy<2> *DDE;
	Solver::Linear::System<DOF> *S;
	set<Set::Manifold::Point *> *x;
        Solver::LineSearch::Secant *LSS;

private:

	NewtonRaphson(NewtonRaphson &);
	void operator=(NewtonRaphson &);
};

}

#endif // !defined(SOLVER_NEWTONRAPHSON__INCLUDED_)
