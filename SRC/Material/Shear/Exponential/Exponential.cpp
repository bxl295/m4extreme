// Exponential.cpp: implementation for the Exponential class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <cassert>
#include "./Exponential.h"

namespace Material
{
namespace Shear
{
namespace Exponential
{
//////////////////////////////////////////////////////////////////////
// Class Data
//////////////////////////////////////////////////////////////////////

Data::Data(): J0(1.0), Mu0(0.0) {}
	
Data::Data(const double &rhs_J0, const double &rhs_Mu0) : 
	J0(rhs_J0), Mu0(rhs_Mu0){}
	
Data::Data(const Data &rhs) : J0(rhs.J0), Mu0(rhs.Mu0){}

Data::~Data(){}

Data & 
Data::operator = (const Data &rhs)
{
	if (this == &rhs) return *this;
	J0 = rhs.J0; Mu0 = rhs.Mu0;
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
Data::GetMu0() const
{
	return Mu0;
}

double & 
Data::GetMu0()
{
	return Mu0;
}

void 
Data::Randomize()
{
	J0 = -(double)rand()/(double)RAND_MAX;
	Mu0 = (double)rand()/(double)RAND_MAX;
}

//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

LocalState::~LocalState() {}

Material::Shear::LocalState *
LocalState::Clone() const 
{
	return new LocalState(*this);
}

LocalState::LocalState(Data *Prop_) : Prop(Prop_) {}

LocalState::LocalState(const LocalState &rhs) : Prop(rhs.Prop) {}

void 
LocalState::operator ++ () {}

//////////////////////////////////////////////////////////////////////
// Class Modulus<0>
//////////////////////////////////////////////////////////////////////

Modulus<0>::~Modulus(){}

Material::Shear::Modulus<0> *
Modulus<0>::Clone() const
{
	return new Modulus<0>(*this);
}

Modulus<0>::Modulus(LocalState *LS_) : LS(LS_){}

Modulus<0>::Modulus(const Modulus<0> &rhs) : LS(rhs.LS){}

double 
Modulus<0>::operator () (const double &J) const
{
	const double &J0 = LS->Prop->J0;
	const double &Mu0 = LS->Prop->Mu0;
	return Mu0*exp(J/J0);
}

//////////////////////////////////////////////////////////////////////
// Class Modulus<1>
//////////////////////////////////////////////////////////////////////

Modulus<1>::~Modulus(){}

Material::Shear::Modulus<1> *
Modulus<1>::Clone() const
{
	return new Modulus<1>(*this);
}

Modulus<1>::Modulus(LocalState *LS_) : LS(LS_){}

Modulus<1>::Modulus(const Modulus<1> &rhs) : LS(rhs.LS){}

double
Modulus<1>::operator () (const double &J) const
{
	const double &J0 = LS->Prop->J0;
	const double &Mu0 = LS->Prop->Mu0;
	return (Mu0/J0)*exp(J/J0);
}

//////////////////////////////////////////////////////////////////////
// Class Modulus<2>
//////////////////////////////////////////////////////////////////////

Modulus<2>::~Modulus(){}

Material::Shear::Modulus<2> *
Modulus<2>::Clone() const
{
	return new Modulus<2>(*this);
}

Modulus<2>::Modulus(LocalState *LS_) : LS(LS_){}

Modulus<2>::Modulus(const Modulus<2> &rhs) : LS(rhs.LS){}

double 
Modulus<2>::operator () (const double &J) const
{
	const double &J0 = LS->Prop->J0;
	const double &Mu0 = LS->Prop->Mu0;
	return (Mu0/(J0*J0))*exp(J/J0);
}

//////////////////////////////////////////////////////////////////////
// Class Jet<0>
//////////////////////////////////////////////////////////////////////

Jet<0>::~Jet() {}

Material::Shear::Jet<0> *
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
	const double &Mu0 = LS->Prop->Mu0;
	double Mu = Mu0*exp(J/J0);
	double DMu = Mu/J0;
	return make_pair(Mu,DMu);
}

//////////////////////////////////////////////////////////////////////
// Class Jet<1>
//////////////////////////////////////////////////////////////////////

Jet<1>::~Jet() {}

Material::Shear::Jet<1> *
Jet<1>::Clone() const
{
	return new Jet<1>(*this);
}

Jet<1>::Jet(LocalState *LS_) : LS(LS_){}

Jet<1>::Jet(const Jet<1> &rhs) : LS(rhs.LS){}

pair<double,double>
Jet<1>::operator () (const double &J) const
{
	const double &J0 = LS->Prop->J0;
	const double &Mu0 = LS->Prop->Mu0;
	double DMu = (Mu0/J0)*exp(J/J0);
	double DDMu = DMu/J0;
	return make_pair(DMu,DDMu);
}

//////////////////////////////////////////////////////////////////////
// Class JetJet<0>
//////////////////////////////////////////////////////////////////////

JetJet<0>::~JetJet() {}

Material::Shear::JetJet<0> *
JetJet<0>::Clone() const
{
	return new JetJet<0>(*this);
}

JetJet<0>::JetJet(LocalState *LS_) : LS(LS_){}

JetJet<0>::JetJet(const JetJet<0> &rhs) : LS(rhs.LS){}

triplet<double,double,double>
JetJet<0>::operator () (const double &J) const 
{
	const double &J0 = LS->Prop->J0;
	const double &Mu0 = LS->Prop->Mu0;
	double Mu = Mu0*exp(J/J0);
	double DMu = (Mu0/J0)*exp(J/J0);
	double DDMu = DMu/J0;
	return make_triplet(Mu,DMu,DDMu);
}

}

}

}
