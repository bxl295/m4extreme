// Constant.cpp: implementation for the Constant class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <cassert>
#include "./Constant.h"

namespace Material
{
namespace Shear
{
namespace Constant
{
//////////////////////////////////////////////////////////////////////
// Class Data
//////////////////////////////////////////////////////////////////////

Data::Data(): Mu(0.0) {}
	
Data::Data(const double &Mu_) : Mu(Mu_)
{
	assert(Mu > 0.0);
}
	
Data::Data(const Data &rhs) : Mu(rhs.Mu) {}

Data::~Data(){}

Data & 
Data::operator = (const Data &rhs)
{
	if (this == &rhs) return *this;
	Mu = rhs.Mu; return *this;
}

const double &
Data::GetMu() const
{
	return Mu;
}

double & 
Data::GetMu()
{
	return Mu;
}

void 
Data::Randomize()
{
	Mu = (double)rand()/(double)RAND_MAX;
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
Modulus<0>::operator () (const double &Theta) const
{
	return LS->Prop->Mu; 
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
Modulus<1>::operator () (const double &Theta) const
{
	return 0.0; 
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
Modulus<2>::operator () (const double &Theta) const
{
	return 0.0; 
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
Jet<0>::operator () (const double &Theta) const 
{
	return make_pair(LS->Prop->Mu,0.0);
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
Jet<1>::operator () (const double &Theta) const
{
	return make_pair(0.0,0.0);
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
JetJet<0>::operator () (const double &Theta) const 
{
	return make_triplet(LS->Prop->Mu,0.0,0.0);
}

}

}

}
