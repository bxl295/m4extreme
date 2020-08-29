// ExplicitDynamics.h: interface for the ExplicitDynamics class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(SOLVER_EXPLICITDYNAMICS__INCLUDED_)
#define SOLVER_EXPLICITDYNAMICS__INCLUDED_

#pragma once

#include <map>
#include "../Solver.h"
#include "../../Clock/Clock.h"
#include "../../Model/Model.h"
#include "../../Set/Manifold/Manifold.h"
#include "../../Model/LumpedMass/LumpedMass.h"

using namespace std;

namespace Solver
{
class ExplicitDynamics : public Propagator  
{
public:

	ExplicitDynamics();
	virtual ~ExplicitDynamics();

	ExplicitDynamics(
		Clock *,
		Model::LocalState *, 
		Model::Energy<1> *, 
		map<Set::Manifold::Point *, double> *, 
		set<Set::Manifold::Point *> *,
		map<Set::Manifold::Point *, Set::VectorSpace::Vector> *,
		const double=0.5);

	ExplicitDynamics(
		Clock *,
		Model::LocalState *, 
		Model::Energy<1> *, 
		map<Set::Manifold::Point *, double> *, 
		set<Set::Manifold::Point *> *,
		map<Set::Manifold::Point *, Set::VectorSpace::Vector> *,
		const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &,
		const double=0.5);

	bool & SetPrint();
        void UpdateMass(map<Set::Manifold::Point*, double> *);
	map<Set::Manifold::Point *, Set::VectorSpace::Vector> & getA() { return a;}
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> & getA() const { return a;}
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> & getF() const { return f;}

	void operator ++ ();
        void BallisticUpdate();
	void Predictor();
	void Predictor(const set<Set::Manifold::Point *> &);
	void Corrector();
	void PullBack();

private:

	double GamOld;
	double GamNew;
	bool Print; Clock *T;
	Model::LocalState *LS;
	Model::Energy<1> *DE;
	map<Set::Manifold::Point *, double> *m;
	set<Set::Manifold::Point *> *x;
	map<Set::Manifold::Point *, Set::VectorSpace::Vector> *v;
	map<Set::Manifold::Point *, Set::VectorSpace::Vector> a;

	map<Set::Manifold::Point *, Set::VectorSpace::Vector> f;

private:

	ExplicitDynamics(ExplicitDynamics &);
	void operator=(ExplicitDynamics &);       
};

}

#endif // !defined(SOLVER_EXPLICITDYNAMICS__INCLUDED_)
