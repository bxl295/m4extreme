// ThermalSoftening.cpp: implementation for the ThermalSoftening class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <cassert>
#include "./ThermalSoftening.h"

namespace Material
{
namespace Uniaxial
{
namespace ThermalSoftening
{
//////////////////////////////////////////////////////////////////////
// Class Data
//////////////////////////////////////////////////////////////////////

Data::Data() {}
	
Data::Data(
	const double &T0_,
	const double &Tm_,
	const double &a_)
	: T0(T0_), Tm_T0(Tm_ - T0_), a(a_) {}
	
Data::Data(const Data &rhs) 
	: T0(rhs.T0), Tm_T0(rhs.Tm_T0), a(rhs.a) {}

Data::~Data(){}

Data & 
Data::operator = (const Data &rhs)
{
	if (this == &rhs) return *this;
	T0 = rhs.T0; Tm_T0 = rhs.Tm_T0; a = rhs.a;
	return *this;
}

void 
Data::Randomize()
{
	T0 = 278.0 + 20.0 * (double)rand()/(double)RAND_MAX;
	Tm_T0 = 1000.0 * (double)rand()/(double)RAND_MAX;
	a = (double)rand()/(double)RAND_MAX;
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
    else {
        return pow(1.0 - T_T0/Tm_T0, LS->Prop->a);
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
    const double & a = LS->Prop->a;

    if ( T_T0 >= Tm_T0 ) {
        return 0.0;
    }
    else {
        return -a * pow(1.0 - T_T0/Tm_T0, a-1.0) / Tm_T0;
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
    double T_T0  = T - LS->Prop->T0;
    const double & Tm_T0 = LS->Prop->Tm_T0;
    const double & a = LS->Prop->a;

    if ( T_T0 >= Tm_T0 ) {
        return 0.0;
    }
    else {
        return a * (a-1.0) * pow(1.0 - T_T0/Tm_T0, a-2.0) / Tm_T0 / Tm_T0;
    }
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
    const double & a = LS->Prop->a;
   
    if ( T_T0 >= Tm_T0 ) {
        return  make_pair(0.0, 0.0);
    }
    else {
        return make_pair( pow(1.0 - T_T0/Tm_T0, a),
                -a * pow(1.0 - T_T0/Tm_T0, a-1.0) / Tm_T0);
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
    const double & a = LS->Prop->a;


    if ( T_T0 >= Tm_T0 ) {
        return  make_pair(0.0, 0.0);
    }
    else {
        return make_pair( -a * pow(1.0 - T_T0/Tm_T0, a-1.0) / Tm_T0 ,
                a * (a-1.0) * pow(1.0 - T_T0/Tm_T0, a-2.0) / Tm_T0 / Tm_T0);
    }
}

}

}

}
