// BlazKo.cpp: implementation for the BlazKo class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <cassert>
#include "./BlatzKo.h"

namespace Material
{
namespace Gas
{
namespace EoS
{
namespace BlatzKo
{
//////////////////////////////////////////////////////////////////////
// Class Data
//////////////////////////////////////////////////////////////////////

Data::Data(): Mu(1.0), Nu(0.25) {}
	
Data::Data(const double &Mu_, const double &Nu_) : Mu(Mu_), Nu(Nu_) 
{
	assert(Mu > 0.0);
	assert((Nu > 0.0) && (Nu < 0.5));

	Alpha = Nu/(1-2*Nu);
	C1 = Mu/(2.0*Alpha); 
	x1 = -2.0*Alpha;
	C2 = 1.5*Mu;
	x2 = 2.0/3.0;
	DC1 = C1*x1;
	dx1 = x1 - 1.0;
	DC2 = C2*x2;
	dx2 = x2 - 1.0;
	DDC1 = DC1*dx1;
	ddx1 = dx1 - 1.0;
	DDC2 = DC2*dx2;
	ddx2 = dx2 - 1.0;
}
	
Data::Data(const Data &rhs) : 
	Mu(rhs.Mu), Nu(rhs.Nu), Alpha(rhs.Alpha),
	C1(rhs.C1), x1(rhs.x1), C2(rhs.C2), x2(rhs.x2),
	DC1(rhs.DC1), dx1(rhs.dx1), DC2(rhs.DC2), dx2(rhs.dx2),
	DDC1(rhs.DDC1), ddx1(rhs.ddx1), DDC2(rhs.DDC2), ddx2(rhs.ddx2) {}

Data::~Data(){}

Data & 
Data::operator = (const Data &rhs)
{
	if (this == &rhs) return *this;
	Mu = rhs.Mu; 
	Nu = rhs.Nu;
	Alpha = rhs.Alpha;
	C1 = rhs.C1; 
	x1 = rhs.x1;
	C2 = rhs.C2;
	x2 = rhs.x2;
	DC1 = rhs.DC1; 
	dx1 = rhs.dx1;
	DC2 = rhs.DC2;
	dx2 = rhs.dx2;
	DDC1 = rhs.DDC1; 
	ddx1 = rhs.ddx1;
	DDC2 = rhs.DDC2;
	ddx2 = rhs.ddx2;
	return *this;
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

const double &
Data::GetNu() const
{
	return Nu;
}

double & 
Data::GetNu()
{
	return Nu;
}

void 
Data::Randomize()
{
	Mu = (double)rand()/(double)RAND_MAX;
	Nu = 0.5*(double)rand()/(double)RAND_MAX;
	Alpha = Nu/(1-2*Nu);
	C1 = Mu/(2.0*Alpha); 
	x1 = - 2.0*Alpha;
	C2 = 1.5*Mu;
	x2 = 2.0/3.0;
	DC1 = C1*x1;
	dx1 = x1 - 1.0;
	DC2 = C2*x2;
	dx2 = x2 - 1.0;
	DDC1 = DC1*dx1;
	ddx1 = dx1 - 1.0;
	DDC2 = DC2*dx2;
	ddx2 = dx2 - 1.0;
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
	const double &C1 = LS->Prop->C1;
	const double &x1 = LS->Prop->x1;
	const double &C2 = LS->Prop->C2;
	const double &x2 = LS->Prop->x2;
	double W = C1*(pow(J,x1) - 1.0) + C2*(pow(J,x2) - 1.0);
	if (W < 0.0) W = 0.0; return W; 
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
	const double &DC1 = LS->Prop->DC1;
	const double &dx1 = LS->Prop->dx1;
	const double &DC2 = LS->Prop->DC2;
	const double &dx2 = LS->Prop->dx2;
	return DC1*pow(J,dx1) + DC2*pow(J,dx2); 
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
	const double &DDC1 = LS->Prop->DDC1;
	const double &ddx1 = LS->Prop->ddx1;
	const double &DDC2 = LS->Prop->DDC2;
	const double &ddx2 = LS->Prop->ddx2;
	return DDC1*pow(J,ddx1) + DDC2*pow(J,ddx2); 
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
	const double &C1 = LS->Prop->C1;
	const double &x1 = LS->Prop->x1;
	const double &C2 = LS->Prop->C2;
	const double &x2 = LS->Prop->x2;
	const double &DC1 = LS->Prop->DC1;
	const double &dx1 = LS->Prop->dx1;
	const double &DC2 = LS->Prop->DC2;
	const double &dx2 = LS->Prop->dx2;
	double W = C1 * ( pow(J,x1) - 1.0 ) + C2 * ( pow(J,x2) - 1.0 ); 
	if (W < 0.0) W = 0.0; double DW = DC1*pow(J,dx1) + DC2*pow(J,dx2);
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
	const double &DC1 = LS->Prop->DC1;
	const double &dx1 = LS->Prop->dx1;
	const double &DC2 = LS->Prop->DC2;
	const double &dx2 = LS->Prop->dx2;
	const double &DDC1 = LS->Prop->DDC1;
	const double &ddx1 = LS->Prop->ddx1;
	const double &DDC2 = LS->Prop->DDC2;
	const double &ddx2 = LS->Prop->ddx2;
	double DW = DC1*pow(J,dx1) + DC2*pow(J,dx2);
	double DDW = DDC1*pow(J,ddx1) + DDC2*pow(J,ddx2); 
	return make_pair(DW,DDW);
}

}

}

}

}
