// LinearThermalSoftening.cpp: implementation for the LinearThermalSoftening class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <cassert>
#include "./LinearThermalSoftening.h"

namespace Material
{
namespace Uniaxial
{
namespace LinearThermalSoftening
{
//////////////////////////////////////////////////////////////////////
// Class Data
//////////////////////////////////////////////////////////////////////

Data::Data() {}
	
Data::Data(
	const double &T0_,
	const double &Tm_)
	: T0(T0_), Tm_T0(Tm_ - T0_){}
	
Data::Data(const Data &rhs) 
	: T0(rhs.T0), Tm_T0(rhs.Tm_T0){}

Data::~Data(){}

Data & 
Data::operator = (const Data &rhs)
{
	if (this == &rhs) return *this;
	T0 = rhs.T0; Tm_T0 = rhs.Tm_T0;
	return *this;
}

void 
Data::Randomize()
{
	T0 = 278.0 + 20.0 * (double)rand()/(double)RAND_MAX;
	Tm_T0 = 1000.0 * (double)rand()/(double)RAND_MAX;
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
    double T_T0  = T - LS->Prop->T0;
    const double & Tm_T0 = LS->Prop->Tm_T0;
    
    if ( T_T0 >= Tm_T0 ) {
        return 0.0;
    }
    else if ( T_T0 < 0.0 ) {
        return 1.0;
    }
    else {
        return 1.0 - T_T0/Tm_T0;
    }
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
    double T_T0  = T - LS->Prop->T0;
    const double & Tm_T0 = LS->Prop->Tm_T0;

    if ( T_T0 >= Tm_T0 || T_T0 < 0.0 ) {
        return 0.0;
    }
    else {
        return -1.0 / Tm_T0;
    }
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
  return 0.0;
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
    double T_T0  = T - LS->Prop->T0;
    const double & Tm_T0 = LS->Prop->Tm_T0;
    
    if ( T_T0 >= Tm_T0 ) {
        return  make_pair(0.0, 0.0);
    }
    else if ( T_T0 < 0.0 ) {
        return make_pair(1.0, 0.0);
    }
    else {
        return make_pair( 1.0 - T_T0/Tm_T0, -1.0 / Tm_T0 );
    }
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
    double T_T0  = T - LS->Prop->T0;
    const double & Tm_T0 = LS->Prop->Tm_T0;

    if ( T_T0 >= Tm_T0 || T_T0 < 0.0 ) {
        return  make_pair(0.0, 0.0);
    }
    else {
      return make_pair( -1.0 / Tm_T0, 0.0);
    }
}

}

}

}
