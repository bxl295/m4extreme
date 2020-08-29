// Ideal.cpp: implementation for the Ideal class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include "./Ideal.h"

namespace Material
{
namespace Gas
{
namespace EoS
{
namespace Ideal
{
//////////////////////////////////////////////////////////////////////
// Class Data
//////////////////////////////////////////////////////////////////////

Data::Data(): J0(1.0), W0(1.0) {}
	
Data::Data(const double &rhs_J0, const double &rhs_W0) 
: J0(rhs_J0), W0(rhs_W0){}
	
Data::Data(const Data &rhs) : J0(rhs.J0), W0(rhs.W0){}

Data::~Data(){}

Data & 
Data::operator = (const Data &rhs)
{
	if (this == &rhs) return *this;
	J0 = rhs.J0; W0 = rhs.W0;
	return *this;
}

const double &
Data::GetJ0() const
{
	return J0;
}

double & 
Data::GetJ0()
{
	return J0;
}

const double &
Data::GetW0() const
{
	return W0;
}

double & 
Data::GetW0()
{
	return W0;
}

void 
Data::Randomize()
{
	J0 = (double)rand()/(double)RAND_MAX;
	W0 = (double)rand()/(double)RAND_MAX;
}

//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

LocalState::~LocalState() {}

Material::Gas::EoS::LocalState *
LocalState::Clone() const 
{
	return new LocalState(*this);
}

LocalState::LocalState(Data *Prop_) : Prop(Prop_) {}

LocalState::LocalState(const LocalState &rhs) : Prop(rhs.Prop) {}

void 
LocalState::operator ++ () {}

//////////////////////////////////////////////////////////////////////
// Class Energy<0>
//////////////////////////////////////////////////////////////////////

Energy<0>::~Energy(){}

Material::Gas::EoS::Energy<0> *
Energy<0>::Clone() const
{
	return new Energy<0>(*this);
}

Energy<0>::Energy(LocalState *LS_) : LS(LS_){}

Energy<0>::Energy(const Energy<0> &rhs) : LS(rhs.LS){}

double 
Energy<0>::operator () (const double &J) const
{
	const double &J0 = LS->Prop->J0;
	const double &W0 = LS->Prop->W0;
	return -W0*log(J/J0);
}

//////////////////////////////////////////////////////////////////////
// Class Energy<1>
//////////////////////////////////////////////////////////////////////

Energy<1>::~Energy(){}

Material::Gas::EoS::Energy<1> *
Energy<1>::Clone() const
{
	return new Energy<1>(*this);
}

Energy<1>::Energy(LocalState *LS_) : LS(LS_){}

Energy<1>::Energy(const Energy<1> &rhs) : LS(rhs.LS){}

double
Energy<1>::operator () (const double &J) const
{
	const double &W0 = LS->Prop->W0;
	double p = -W0/J;
	LS->_pressure = p;
	return p;
}

//////////////////////////////////////////////////////////////////////
// Class Energy<2>
//////////////////////////////////////////////////////////////////////

Energy<2>::~Energy(){}

Material::Gas::EoS::Energy<2> *
Energy<2>::Clone() const
{
	return new Energy<2>(*this);
}

Energy<2>::Energy(LocalState *LS_) : LS(LS_){}

Energy<2>::Energy(const Energy<2> &rhs) : LS(rhs.LS){}

double 
Energy<2>::operator () (const double &J) const
{
	const double &W0 = LS->Prop->W0;
	return W0/(J*J);
}

//////////////////////////////////////////////////////////////////////
// Class Jet<0>
//////////////////////////////////////////////////////////////////////

Jet<0>::~Jet() {}

Material::Gas::EoS::Jet<0> *
Jet<0>::Clone() const
{
	return new Jet<0>(*this);
}

Jet<0>::Jet(LocalState *LS_) : LS(LS_){}

Jet<0>::Jet(const Jet<0> &rhs) : LS(rhs.LS){}

pair<double,double>
Jet<0>::operator () (const double &J) const
{
	const double &J0 = LS->Prop->J0;
	const double &W0 = LS->Prop->W0;
	double W = -W0*log(J/J0);
	double DW = -W0/J;
	return make_pair(W,DW);

}

//////////////////////////////////////////////////////////////////////
// Class Jet<1>
//////////////////////////////////////////////////////////////////////

Jet<1>::~Jet() {}

Material::Gas::EoS::Jet<1> *
Jet<1>::Clone() const
{
	return new Jet<1>(*this);
}

Jet<1>::Jet(LocalState *LS_) : LS(LS_){}

Jet<1>::Jet(const Jet<1> &rhs) : LS(rhs.LS){}

pair<double,double>
Jet<1>::operator () (const double &J) const
{
	const double &W0 = LS->Prop->W0;
	double DW = -W0/J;
	double DDW = W0/(J*J);
	return make_pair(DW,DDW);
}

}

}

}

}
