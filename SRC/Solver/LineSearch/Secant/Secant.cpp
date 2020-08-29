// Secant.cpp: implementation of the Secant class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include "./Secant.h"

namespace Solver
{
namespace LineSearch
{
Secant::Secant() : 
	LS(0), DE(0), x(0), 
	Print(false), Tol(0.0), Step(0.0), NitMax(0) {}

Secant::~Secant() {}

Secant::Secant(
	Model::LocalState *LS_, 
	Model::Energy<1> *DE_, 
	set<Set::Manifold::Point *> *x_) : 
	LS(LS_), DE(DE_), x(x_), 
	Print(false), Tol(0.0), Step(0.0), NitMax(0) {}

bool & 
Secant::SetPrint()
{
	return Print;
}

unsigned int & 
Secant::SetNitMax()
{
	return NitMax;
}

double & 
Secant::SetTol()
{
	return Tol;
}

double & 
Secant::SetStep()
{
	return Step;
}

void 
Secant::operator () (
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &u)
{
	unsigned int k; 
	double ds = Step, dfold, dfnew, ddf, Error;
	if (Print) cout << "   Secant line-search iteration:" << endl;
	dfold = (*DE)(*x,u); 
	*x += ds*u;
	dfnew = (*DE)(*x,u); 
	for (k=0; k<NitMax; k++)
	{
		ddf = (dfnew - dfold)/ds; //cout << "ddf = " << ddf << endl;
                if (ddf > Tol)
                    ds = - dfnew/ddf; //cout << "ds = "<< ds << endl;
                else
                    ds += Step;
		*x += ds*u;
		dfold = dfnew; 
		dfnew = (*DE)(*x,u); 
		Error = fabs(dfnew);
		if (Print) cout 
			<< "   Iteration = " << k << ", " 	
			<< "Error = " << Error << endl;
		if (Error < Tol) return;
	}
	//throw(0);
}

}

}
