// PairPotential.cpp: implementation for the PairPotential class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include "./PairPotential.h"

namespace Potential
{
namespace Johnson
{
namespace PairPotential
{
//////////////////////////////////////////////////////////////////////
// Class Data
//////////////////////////////////////////////////////////////////////

Data::Data() {}

Data::~Data() {}

Data::Data(const double * const rhs) 
 : phie(rhs[0]), gamma(rhs[1]), ao(rhs[2]) {}

Data::Data(const double &rhs_phie, 
           const double &rhs_gamma,
           const double &rhs_ao) 
 : phie(rhs_phie), gamma(rhs_gamma), ao(rhs_ao) {}

Data::Data(const Data &rhs) 
 : phie(rhs.phie), gamma(rhs.gamma), ao(rhs.ao) {}

Data & Data::operator = (const Data &rhs)
{
	if (this == &rhs) return *this; 
	phie = rhs.phie; gamma = rhs.gamma; ao = rhs.ao; 
	return *this;
}

void Data::Randomize()
{
	phie  = (double)rand()/(double)RAND_MAX;
	gamma = (double)rand()/(double)RAND_MAX;
	ao    = (double)rand()/(double)RAND_MAX;
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
Energy<0>::operator () (const double &r) const
{
	double V = 0.0;
	double re = LS->Prop->ao/sqrt(2.0);
	V = LS->Prop->phie*exp( -LS->Prop->gamma*(r/re-1.0) ); 
	return 0.5*V;
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
Energy<1>::operator () (const double &r) const
{
	double dV = 0.0;
	double re = LS->Prop->ao/sqrt(2.0);
	dV = -LS->Prop->phie*LS->Prop->gamma/re*
             exp( -LS->Prop->gamma*(r/re-1.0) ); 
	return 0.5*dV;
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
Energy<2>::operator () (const double &r) const
{
	double ddV = 0.0;
        double re = LS->Prop->ao/sqrt(2.0);
	ddV = pow(LS->Prop->gamma,2)*LS->Prop->phie/pow(re,2)*
                exp( -LS->Prop->gamma*(r/re-1.0) ); 
	return 0.5*ddV;
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
Jet<0>::operator () (const double &r) const
{
	double V = 0.0, dV = 0.0;
	double re = LS->Prop->ao/sqrt(2.0);
	V = LS->Prop->phie*exp( -LS->Prop->gamma*(r/re-1.0) ); 
	dV = -LS->Prop->phie*LS->Prop->gamma/re*
		exp( -LS->Prop->gamma*(r/re-1.0) ); 
        return make_pair(0.5*V,0.5*dV);
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
Jet<1>::operator () (const double &r) const
{
	double dV = 0.0, ddV=0.0;
	double re = LS->Prop->ao/sqrt(2.0);
	dV = -LS->Prop->phie*LS->Prop->gamma/re*
		exp( -LS->Prop->gamma*(r/re-1.0) ); 
        ddV = pow(LS->Prop->gamma,2)*LS->Prop->phie/pow(re,2)*
		exp( -LS->Prop->gamma*(r/re-1.0) ); 
	return make_pair(0.5*dV,0.5*ddV);
}

}

}

}
