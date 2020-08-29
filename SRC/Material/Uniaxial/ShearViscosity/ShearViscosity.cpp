// ShearViscosity.cpp: implementation for the ShearViscosity class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <cassert>
#include "./ShearViscosity.h"

namespace Material
{
namespace Uniaxial
{
namespace ShearViscosity
{
//////////////////////////////////////////////////////////////////////
// Class Data
//////////////////////////////////////////////////////////////////////

Data::Data() {}
	
Data::Data(const double &eta0_,
	   const double &a1_,
	   const double &a2_)
  : eta0(eta0_), a1(a1_), a2(a2_) {}
	
Data::Data(const Data &rhs) 
  : eta0(rhs.eta0), a1(rhs.a1), a2(rhs.a2) {}

Data::~Data(){}

Data & 
Data::operator = (const Data &rhs)
{
	if (this == &rhs) return *this;
	eta0 = rhs.eta0; a1 = rhs.a1; a2 = rhs.a2;
	return *this;
}

void 
Data::Randomize()
{
        eta0 = (double)rand() / (double)RAND_MAX;
	a1   = (double)rand() / (double)RAND_MAX;
	a2   = (double)rand() / (double)RAND_MAX;
}

//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

LocalState::~LocalState() {}

Material::Uniaxial::LocalState *
LocalState::Clone() const 
{
	return new LocalState(*this);
}

LocalState::LocalState(Data *Prop_) : 
	Prop(Prop_) {}

LocalState::LocalState(const LocalState &rhs) :
	Prop(rhs.Prop) {}

void 
LocalState::operator ++ () {}

//////////////////////////////////////////////////////////////////////
// Class Energy<0>
//////////////////////////////////////////////////////////////////////

Energy<0>::~Energy(){}

Material::Uniaxial::Energy<0> *
Energy<0>::Clone() const
{
	return new Energy<0>(*this);
}

Energy<0>::Energy(LocalState *LS_) : LS(LS_){}

Energy<0>::Energy(const Energy<0> &rhs) : LS(rhs.LS){}

double 
Energy<0>::operator () (const double &T) const
{
    assert( T > 0.0 );
    const Data * Prop = LS->Prop;
    return Prop->eta0 * exp(Prop->a1 + Prop->a2 / T);
}

//////////////////////////////////////////////////////////////////////
// Class Energy<1>
//////////////////////////////////////////////////////////////////////

Energy<1>::~Energy(){}

Material::Uniaxial::Energy<1> *
Energy<1>::Clone() const
{
	return new Energy<1>(*this);
}

Energy<1>::Energy(LocalState *LS_) : LS(LS_){}

Energy<1>::Energy(const Energy<1> &rhs) : LS(rhs.LS){}

double 
Energy<1>::operator () (const double &T) const
{
    assert( T > 0.0 );
    const Data * Prop = LS->Prop;
    return -Prop->a2 * Prop->eta0 * exp(Prop->a1 + Prop->a2 / T) / T / T;
}

//////////////////////////////////////////////////////////////////////
// Class Energy<2>
//////////////////////////////////////////////////////////////////////

Energy<2>::~Energy(){}

Material::Uniaxial::Energy<2> *
Energy<2>::Clone() const
{
	return new Energy<2>(*this);
}

Energy<2>::Energy(LocalState *LS_) : LS(LS_){}

Energy<2>::Energy(const Energy<2> &rhs) : LS(rhs.LS){}

double 
Energy<2>::operator () (const double &T) const
{
    assert( T > 0.0 );
    const Data * Prop = LS->Prop;
    const double & a2 = Prop->a2;
    return a2 * Prop->eta0 * exp(Prop->a1 + a2 / T) * (2.0 + a2 / T) / T / T / T;
}

//////////////////////////////////////////////////////////////////////
// Class Jet<0>
//////////////////////////////////////////////////////////////////////

Jet<0>::~Jet() {}

Material::Uniaxial::Jet<0> *
Jet<0>::Clone() const
{
	return new Jet<0>(*this);
}

Jet<0>::Jet(LocalState *LS_) : LS(LS_) {}

Jet<0>::Jet(const Jet<0> &rhs) : LS(rhs.LS) {}

pair<double,double>
Jet<0>::operator () (const double &T) const
{
    assert( T > 0.0 );

    const double & eta0 = LS->Prop->eta0;    
    const double & a2 = LS->Prop->a2;
    double eta = eta0 * exp(LS->Prop->a1 + a2 / T);

    return make_pair( eta, -a2 * eta / T / T );
}

//////////////////////////////////////////////////////////////////////
// Class Jet<1>
//////////////////////////////////////////////////////////////////////

Jet<1>::~Jet() {}

Material::Uniaxial::Jet<1> *
Jet<1>::Clone() const
{
	return new Jet<1>(*this);
}

Jet<1>::Jet(LocalState *LS_) : LS(LS_) {}

Jet<1>::Jet(const Jet<1> &rhs) : LS(rhs.LS) {}

pair<double,double> 
Jet<1>::operator () (const double &T) const
{
    assert( T > 0.0 );

    const double & eta0 = LS->Prop->eta0;    
    const double & a2 = LS->Prop->a2;
    double eta = eta0 * exp(LS->Prop->a1 + a2 / T);

    return make_pair( -a2 * eta / T / T, a2 * eta * (2.0 + a2/T) / T / T / T );    
}

}

}

}
