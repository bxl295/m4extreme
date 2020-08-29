// NewtonRaphson.cpp: implementation of the NewtonRaphson class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include "./NewtonRaphson.h"

namespace Solver
{
namespace LineSearch
{
NewtonRaphson::NewtonRaphson() : 
	LS(0), DJ(0), DE(0), DDE(0), x(0), 
	Print(false), Tol(0.0), NitMax(0) {}

NewtonRaphson::~NewtonRaphson() {}

NewtonRaphson::NewtonRaphson(
	Model::LocalState *LS_, 
	Model::Jet<1> *DJ_, 
	set<Set::Manifold::Point *> *x_) : 
	LS(LS_), DJ(DJ_), DE(0), DDE(0), x(x_), 
	Print(false), Tol(0.0), NitMax(0) {}

NewtonRaphson::NewtonRaphson(
	Model::LocalState *LS_, 
	Model::Energy<1> *rhs_DE, 
	Model::Energy<2> *rhs_DDE, 
	set<Set::Manifold::Point *> *x_) : 
	LS(LS_), DJ(0), DE(rhs_DE), DDE(rhs_DDE), x(x_), 
	Print(false), Tol(0.0), NitMax(0) {}

bool & 
NewtonRaphson::SetPrint()
{
	return Print;
}

unsigned int & 
NewtonRaphson::SetNitMax()
{
	return NitMax;
}

double & 
NewtonRaphson::SetTol()
{
	return Tol;
}

void 
NewtonRaphson::operator () (
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &u)
{
	unsigned int k; pair<double,double> dj;
	double Error, df, ddf, DStep, Step=0.0;
	if (Print) cout << "   Newton-Raphson line-search iteration:" << endl;
	for (k=0; k<NitMax; k++)
	{
		if (DJ != 0)
		{
			dj = (*DJ)(*x,u);
			df = dj.first; 
			ddf = dj.second;
		}
		else
		{
			df = (*DE)(*x,u); 
			ddf = (*DDE)(*x,u); 
		}
		DStep = - df/ddf; 
		*x += DStep*u;
		Error = fabs(DStep);
		Step += DStep;
		if (Print) cout 
			<< "   Iteration = " << k << ", " 	
			<< "Error = " << Error << endl;
		if (Error < Tol) return;
	}
	throw(0);
}

}

}
