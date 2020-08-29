// LogLaw.cpp: implementation for the LogLaw class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <cassert>
#include "./LogLaw.h"

namespace Material
{
namespace Uniaxial
{
namespace LogLaw
{
//////////////////////////////////////////////////////////////////////
// Class Data
//////////////////////////////////////////////////////////////////////

Data::Data() {}
	
Data::Data(const double &Cv_)
  : Cv(Cv_) {}
	
Data::Data(const Data &rhs) 
  : Cv(rhs.Cv) {}

Data::~Data(){}

Data & 
Data::operator = (const Data &rhs)
{
	if (this == &rhs) return *this;
	Cv = rhs.Cv;
	return *this;
}

void 
Data::Randomize()
{
	Cv = (double)rand()/(double)RAND_MAX;
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

LocalState::LocalState(double T0_, Data *Prop_) : 
  Prop(Prop_), TOld(T0_), T(T0_) {}

LocalState::LocalState(const LocalState &rhs) :
  Prop(rhs.Prop), TOld(rhs.TOld), T(rhs.T) {}

void 
LocalState::operator ++ () {
  TOld = T;
}

void 
LocalState::Reset(double T_) {
  T = T_;
}

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
    const double &TOld = LS->TOld;
    if ( fabs(T- TOld)<1.0e-8 ) return 0.0;
    else return LS->Prop->Cv * ( T - TOld - T*log(T/TOld));
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
    const double &TOld = LS->TOld;
    if ( fabs(T- TOld)<1.0e-8 ) return 0.0;
    else return -LS->Prop->Cv * log(T/TOld);
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
    return -LS->Prop->Cv / T;
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
    const double & TOld = LS->TOld;
    const double & Cv = LS->Prop->Cv;

    assert( T > 0.0 );
    if ( fabs(T- TOld)<1.0e-8 ) return make_pair(0.0, 0.0);
    else return make_pair( Cv * T * ( 1.0 - log(T/TOld)),
			   -Cv * log(T/TOld));
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
    const double & TOld = LS->TOld;
    const double & Cv = LS->Prop->Cv;

    assert( T > 0.0 );
    if ( fabs(T- TOld)<1.0e-8 ) return make_pair(0.0, -Cv/T);
    else return make_pair( -Cv * log(T/TOld), -Cv / T );
}

}

}

}
