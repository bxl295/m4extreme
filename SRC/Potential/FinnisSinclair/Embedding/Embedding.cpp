// Embedding.cpp: implementation for the Embedding class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <cassert>
#include "./Embedding.h"

namespace Potential
{
namespace FinnisSinclair
{
namespace Embedding
{
//////////////////////////////////////////////////////////////////////
// Class Data
//////////////////////////////////////////////////////////////////////


Data::Data(){}

Data::~Data(){}

Data::Data(const double * const rhs) : A(rhs[0]) {}

Data::Data(const double &rhs_A) : A(rhs_A) {}
	
Data::Data(const Data &rhs) : A(rhs.A) {}

Data & 
Data::operator = (const Data &rhs)
{
	if (this == &rhs) return *this;
	A = rhs.A; return *this;
}

void 
Data::Randomize()
{
	A = (double)rand()/(double)RAND_MAX;
}

//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

LocalState::LocalState() {Prop = 0;}

LocalState::~LocalState() {}

Potential::Radial::LocalState *
LocalState::Clone() const 
{
	return new LocalState(*this);
}

LocalState::LocalState(const Data *rhs_Prop)
	: Prop(rhs_Prop) {}

LocalState::LocalState(const LocalState &) {}

LocalState & 
LocalState::operator = (const LocalState &rhs)
{
	if (this == &rhs) return *this; 
	Prop = rhs.Prop; 
	return *this;
}

void 
LocalState::operator ++ () {}

//////////////////////////////////////////////////////////////////////
// Class Energy<0>
//////////////////////////////////////////////////////////////////////

Energy<0>::Energy() {LS=0;}

Energy<0>::~Energy(){}

Potential::Radial::Energy<0> *
Energy<0>::Clone() const
{
	return new Energy<0>(*this);
}

Energy<0>::Energy(LocalState *rhs_LS) : LS(rhs_LS){}

Energy<0>::Energy(const Energy<0> &rhs) : LS(rhs.LS){}

const Energy<0> & 
Energy<0>::operator = (const Energy<0> &rhs)
{
	if (this == &rhs) return *this; 
	LS = rhs.LS; return *this;
}

double 
Energy<0>::operator () (const double &n) const
{
	assert (n >= 0.0);
	return -LS->Prop->A*sqrt(n);
}

//////////////////////////////////////////////////////////////////////
// Class Energy<1>
//////////////////////////////////////////////////////////////////////

Energy<1>::Energy() {LS=0;}

Energy<1>::~Energy(){}

Potential::Radial::Energy<1> *
Energy<1>::Clone() const
{
	return new Energy<1>(*this);
}

Energy<1>::Energy(LocalState *rhs_LS) : LS(rhs_LS){}

Energy<1>::Energy(const Energy<1> &rhs) : LS(rhs.LS){}

const Energy<1> & 
Energy<1>::operator = (const Energy<1> &rhs)
{
	if (this == &rhs) return *this; 
	LS = rhs.LS; return *this;
}

double 
Energy<1>::operator () (const double &n) const
{
	assert (n >= 0.0);
	if (n == 0) return 0.0;
	return -0.5*LS->Prop->A/sqrt(n);
}

//////////////////////////////////////////////////////////////////////
// Class Energy<2>
//////////////////////////////////////////////////////////////////////

Energy<2>::Energy() {LS=0;}

Energy<2>::~Energy(){}

Potential::Radial::Energy<2> *
Energy<2>::Clone() const
{
	return new Energy<2>(*this);
}

Energy<2>::Energy(LocalState *rhs_LS) : LS(rhs_LS){}

Energy<2>::Energy(const Energy<2> &rhs) : LS(rhs.LS){}

const Energy<2> & 
Energy<2>::operator = (const Energy<2> &rhs)
{
	if (this == &rhs) return *this; 
	LS = rhs.LS; return *this;
}

double 
Energy<2>::operator () (const double &n) const
{
	assert (n >= 0.0);
	if (n == 0) return 0.0;
	return 0.25*LS->Prop->A/pow(n,1.5);
}

//////////////////////////////////////////////////////////////////////
// Class Jet<0>
//////////////////////////////////////////////////////////////////////

Jet<0>::Jet() {LS=0;}

Jet<0>::~Jet() {}

Potential::Radial::Jet<0> *
Jet<0>::Clone() const
{
	return new Jet<0>(*this);
}

Jet<0>::Jet(LocalState *rhs_LS) : LS(rhs_LS) {}

Jet<0>::Jet(const Jet<0> &rhs) : LS(rhs.LS) {}

Jet<0> & 
Jet<0>::operator = (const Jet<0> &rhs)
{
	if (&rhs == this) return *this;
	LS = rhs.LS; return *this;
}

pair<double,double>
Jet<0>::operator () (const double &n) const
{
	assert (n >= 0.0);
	if (n == 0) return make_pair(0.0,0.0);
	return make_pair(-LS->Prop->A*sqrt(n),-0.5*LS->Prop->A/sqrt(n));
}

//////////////////////////////////////////////////////////////////////
// Class Jet<1>
//////////////////////////////////////////////////////////////////////

Jet<1>::Jet() {LS=0;}

Jet<1>::~Jet() {}

Potential::Radial::Jet<1> *
Jet<1>::Clone() const
{
	return new Jet<1>(*this);
}

Jet<1>::Jet(LocalState *rhs_LS) : LS(rhs_LS) {}

Jet<1>::Jet(const Jet<1> &rhs) : LS(rhs.LS) {}

Jet<1> & 
Jet<1>::operator = (const Jet<1> &rhs)
{
	if (&rhs == this) return *this;
	LS = rhs.LS; return *this;
}

pair<double,double> 
Jet<1>::operator () (const double &n) const
{
	assert (n >= 0.0);
	if (n == 0) return make_pair(0.0,0.0);
	return make_pair(-0.5*LS->Prop->A/sqrt(n),0.25*LS->Prop->A/pow(n,1.5));
}

}

}

}
