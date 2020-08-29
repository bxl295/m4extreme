// PowerLaw.cpp: implementation for the PowerLaw class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <cassert>
#include "./PowerLaw.h"

namespace Material
{
namespace Uniaxial
{
namespace PowerLaw
{
//////////////////////////////////////////////////////////////////////
// Class Data
//////////////////////////////////////////////////////////////////////

Data::Data() {}
	
Data::Data(
	const double &f0_, 
	const double &C1_,
	const double &r0_, 
	const double &x_,
	double C3_) 
	: f0(f0_), C1(C1_), r0(r0_), x(x_), C3(C3_)
{
	xpp = x + 1.0;
	xmm = x - 1.0;
	C0 = C1*r0/xpp;
	C2 = C1*x/r0;
	C4 = C0*pow(C3,xpp);
}
	
Data::Data(const Data &rhs) 
	: f0(rhs.f0), r0(rhs.r0), 
	x(rhs.x), xpp(rhs.xpp), xmm(rhs.xmm),
	C0(rhs.C0), C1(rhs.C1), C2(rhs.C2), 
	C3(rhs.C3), C4(rhs.C4) {}

Data::~Data(){}

Data & 
Data::operator = (const Data &rhs)
{
	if (this == &rhs) return *this;
	f0 = rhs.f0; r0 = rhs.r0;
	x = rhs.x; xpp = rhs.xpp; xmm = rhs.xmm;
	C0 = rhs.C0; C1 = rhs.C1; C2 = rhs.C2;
	C3 = rhs.C3; C4 = rhs.C4;
	return *this;
}

void 
Data::Randomize()
{
	f0 = 0.05 + (double)rand()/(double)RAND_MAX;
	C1 = 0.05 + (double)rand()/(double)RAND_MAX;
	r0 = 0.05 + (double)rand()/(double)RAND_MAX;
	x  = 0.05 + (double)rand()/(double)RAND_MAX;
	xpp = x + 1.0;
	xmm = x - 1.0;
	C0 = C1*r0/xpp;
	C2 = C1*x/r0;
	C3 = (double)rand()/(double)RAND_MAX;
	C4 = C0*pow(C3,xpp);
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
Energy<0>::operator () (const double &r) const
{
	assert(r >= 0.0); 
	double s = LS->Prop->C3+r/LS->Prop->r0;
        if ( s != 0.0 ) {
            return LS->Prop->f0*r - LS->Prop->C4 +
		LS->Prop->C0*pow(s,LS->Prop->xpp);
        }
        else {
            return LS->Prop->f0*r - LS->Prop->C4;
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
Energy<1>::operator () (const double &r) const
{
	assert(r >= 0.0); 
	double s = LS->Prop->C3+r/LS->Prop->r0;
        if ( s != 0.0 ) {
            return LS->Prop->f0 + LS->Prop->C1*pow(s,LS->Prop->x);
        }
        else {
            return LS->Prop->f0;
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
Energy<2>::operator () (const double &r) const
{
	assert(r >= 0.0); 
	double s = LS->Prop->C3+r/LS->Prop->r0;
        if ( s != 0.0 ) {
            return LS->Prop->C2*pow(s,LS->Prop->xmm);
        }
        else {
            return 0.0;
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
Jet<0>::operator () (const double &r) const
{
	assert(r >= 0.0); 
	double s = LS->Prop->C3+r/LS->Prop->r0;
        if ( s != 0.0 ) {
            return make_pair(LS->Prop->f0*r - LS->Prop->C4 +
                	LS->Prop->C0*pow(s,LS->Prop->xpp),
                        LS->Prop->f0 + LS->Prop->C1*pow(s,LS->Prop->x));
        }
        else {
            return make_pair(LS->Prop->f0*r - LS->Prop->C4,
                        LS->Prop->f0);
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
Jet<1>::operator () (const double &r) const
{
	assert(r >= 0.0);  
	double s = LS->Prop->C3+r/LS->Prop->r0;
        if ( s != 0.0 ) {
            return make_pair(
                    LS->Prop->f0 + LS->Prop->C1*pow(s,LS->Prop->x),
                    LS->Prop->C2*pow(s,LS->Prop->xmm));
        }
        else {
            return make_pair(LS->Prop->f0, 0.0);
        }
}

}

}

}
