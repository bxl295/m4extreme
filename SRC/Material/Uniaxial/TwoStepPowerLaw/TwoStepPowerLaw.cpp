// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <cassert>
#include "./TwoStepPowerLaw.h"

namespace Material
{
namespace Uniaxial
{
namespace TwoStepPowerLaw
{
//////////////////////////////////////////////////////////////////////
// Class Data
//////////////////////////////////////////////////////////////////////

Data::Data() {}
	
Data::Data(
	const double &f0_, 
	const double &C1_,
	const double &r0_, 
	const double &rc_,
	const double &x_,
	const double &y_,
	double C3_) 
  : f0(f0_), C1(C1_), r0(r0_), rc(rc_), x(x_), y(y_), C3(C3_)
{}
	
Data::Data(const Data &rhs) 
	: f0(rhs.f0), r0(rhs.r0),rc(rhs.rc), 
	x(rhs.x), y(rhs.y), 
	C1(rhs.C1), C3(rhs.C3) {}

Data::~Data(){}

Data & 
Data::operator = (const Data &rhs)
{
	if (this == &rhs) return *this;
	f0 = rhs.f0; r0 = rhs.r0; rc = rhs.rc;
	x = rhs.x; y = rhs.y; 
	C1 = rhs.C1; C3 = rhs.C3;
	return *this;
}

void 
Data::Randomize()
{
	f0 = 0.05 + (double)rand()/(double)RAND_MAX;
	C1 = 0.05 + (double)rand()/(double)RAND_MAX;
	r0 = 0.05 + (double)rand()/(double)RAND_MAX;
	rc = 0.05 + (double)rand()/(double)RAND_MAX;
	x  = 0.05 + (double)rand()/(double)RAND_MAX;
	y  = 0.05 + (double)rand()/(double)RAND_MAX;
	C3 = (double)rand()/(double)RAND_MAX;	
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

	const double & x  =  LS->Prop->x;
	const double & y  =  LS->Prop->y;
	const double & f0 =  LS->Prop->f0;
	const double & r0 =  LS->Prop->r0;
	const double & rc =  LS->Prop->rc;
	const double & C1 =  LS->Prop->C1;
	const double & C3 =  LS->Prop->C3;

	if ( r > rc ) {
	  double ypp = 1.0 + y;
	  double C2 = C1 * r0 / ypp;
	  return f0 * r + C2 * pow(C3 + rc/r0, x - y ) * (pow(C3 + r/r0, ypp) - pow(C3, ypp));
	}
	else {
	  double xpp = 1.0 + x;
	  double C2 = C1 * r0 / xpp;
	  return f0 * r + C2 * (pow(C3 + r/r0, xpp) - pow(C3, xpp));
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
	const double & x  =  LS->Prop->x;
	const double & y  =  LS->Prop->y;
	const double & f0 =  LS->Prop->f0;
	const double & r0 =  LS->Prop->r0;
	const double & rc =  LS->Prop->rc;
	const double & C1 =  LS->Prop->C1;
	const double & C3 =  LS->Prop->C3;

	if ( r > rc ) {
	  return f0 + C1 * pow(C3 + r/r0, y) * pow(C3 + rc/r0, x - y );
	}
	else {
	  return f0 + C1 * pow(C3 + r/r0, x);
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

	const double & x  =  LS->Prop->x;
	const double & y  =  LS->Prop->y;
	const double & r0 =  LS->Prop->r0;
	const double & rc =  LS->Prop->rc;
	const double & C1 =  LS->Prop->C1;
	const double & C3 =  LS->Prop->C3;

	if ( r > rc ) {
	  return C1 * y * pow(C3 + r/r0, y - 1.0) * pow(C3 + rc/r0, x - y ) / r0;
	}
	else {
	  return C1 * x * pow(C3 + r/r0, x - 1.0) / r0;
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
	
	const double & x  =  LS->Prop->x;
	const double & y  =  LS->Prop->y;
	const double & f0 =  LS->Prop->f0;
	const double & r0 =  LS->Prop->r0;
	const double & rc =  LS->Prop->rc;
	const double & C1 =  LS->Prop->C1;
	const double & C3 =  LS->Prop->C3;

	if ( r > rc ) {	  
	  double ypp = 1.0 + y;
	  double C2 = C1 * r0 / ypp;
	  return make_pair(f0 * r + C2 * pow(C3 + rc/r0, x - y ) * (pow(C3 + r/r0, ypp) - pow(C3, ypp)),
			   f0 + C1 * pow(C3 + r/r0, y) * pow(C3 + rc/r0, x - y ));
	}
	else {
	  double xpp = 1.0 + x;
	  double C2 = C1 * r0 / xpp;
	  return make_pair(f0 * r + C2 * (pow(C3 + r/r0, xpp) - pow(C3, xpp)),
			   f0 + C1 * pow(C3 + r/r0, x) );
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

	const double & x  =  LS->Prop->x;
	const double & y  =  LS->Prop->y;
	const double & f0 =  LS->Prop->f0;
	const double & r0 =  LS->Prop->r0;
	const double & rc =  LS->Prop->rc;
	const double & C1 =  LS->Prop->C1;
	const double & C3 =  LS->Prop->C3;

	double s = LS->Prop->C3+r/LS->Prop->r0;

	if ( r > LS->Prop->rc ) {
	  return make_pair(f0 + C1 * pow(C3 + r/r0, y) * pow(C3 + rc/r0, x - y ),
			   C1 * y * pow(C3 + r/r0, y - 1.0) * pow(C3 + rc/r0, x - y ) / r0);
	}
	else {
	  return make_pair(f0 + C1 * pow(C3 + r/r0, x),
			   C1 * x * pow(C3 + r/r0, x - 1.0) / r0);
	}
}

}

}

}
