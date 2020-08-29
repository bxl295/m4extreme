// LennardJones.cpp: implementation for the LennardJones class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <cassert>
#include <iostream>
#include "./LennardJones.h"

namespace Potential
{
namespace Radial
{
namespace LennardJones
{
//////////////////////////////////////////////////////////////////////
// Class Data
//////////////////////////////////////////////////////////////////////

Data::Data(){}

Data::Data(const double &rhs_EMin, const double &rhs_RMin,
	   const double &rhs_m, const double &rhs_n,
	   bool rhs_IsCompressiveOnly) :
  EMin(rhs_EMin), RMin(rhs_RMin), m(rhs_m), n(rhs_n),
  IsCompressiveOnly(rhs_IsCompressiveOnly) {
  assert(n>0.0 && m>n);
}
	
Data::Data(const Data &rhs) :
  EMin(rhs.EMin), RMin(rhs.RMin), m(rhs.m), n(rhs.n),
  IsCompressiveOnly(rhs.IsCompressiveOnly){}

Data::~Data(){}

Data & 
Data::operator = (const Data &rhs)
{
	if (this == &rhs) return *this;
	EMin = rhs.EMin; RMin = rhs.RMin;
	m = rhs.m; n = rhs.n;
	IsCompressiveOnly = rhs.IsCompressiveOnly;
	return *this;
}

void 
Data::Randomize()
{
        IsCompressiveOnly = false;
        EMin = (double)rand()/(double)RAND_MAX;
	RMin = (double)rand()/(double)RAND_MAX;
	m = 10.0 * (double)rand()/(double)RAND_MAX;
	n = 10.0 * (double)rand()/(double)RAND_MAX;
}

//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

LocalState::LocalState() {Prop = 0;}

LocalState::~LocalState() {}

Potential::Radial::LocalState *
LocalState::Clone() const 
{
	return new LocalState(*this);
}

LocalState::LocalState(const Data *rhs_Prop)
	: Prop(rhs_Prop) {}

LocalState::LocalState(const LocalState &) {}

LocalState & 
LocalState::operator = (const LocalState &rhs)
{
	if (this == &rhs) return *this; 
	Prop = rhs.Prop; 
	return *this;
}

void 
LocalState::operator ++ () {}

//////////////////////////////////////////////////////////////////////
// Class Energy<0>
//////////////////////////////////////////////////////////////////////

Energy<0>::Energy() {LS=0;}

Energy<0>::~Energy(){}

Potential::Radial::Energy<0> *
Energy<0>::Clone() const
{
	return new Energy<0>(*this);
}

Energy<0>::Energy(LocalState *rhs_LS) : LS(rhs_LS){}

Energy<0>::Energy(const Energy<0> &rhs) : LS(rhs.LS){}

const Energy<0> & 
Energy<0>::operator = (const Energy<0> &rhs)
{
	if (this == &rhs) return *this; 
	LS = rhs.LS; return *this;
}

double 
Energy<0>::operator () (const double &r) const
{
        const Data *Prop = LS->Prop;
	double x = Prop->RMin/r;
	bool isValid = true;
	if ( Prop->IsCompressiveOnly && x < 1.0 ) {
	  isValid = false;
	}

	if ( isValid ) {
	  const double & m = Prop->m;
	  const double & n = Prop->n;
	  return Prop->EMin*(pow(x, m) - m*pow(x, n)/n);
	}
	else {
	  return 0.0;
	}
}

//////////////////////////////////////////////////////////////////////
// Class Energy<1>
//////////////////////////////////////////////////////////////////////

Energy<1>::Energy() {LS=0;}

Energy<1>::~Energy(){}

Potential::Radial::Energy<1> *
Energy<1>::Clone() const
{
	return new Energy<1>(*this);
}

Energy<1>::Energy(LocalState *rhs_LS) : LS(rhs_LS){}

Energy<1>::Energy(const Energy<1> &rhs) : LS(rhs.LS){}

const Energy<1> & 
Energy<1>::operator = (const Energy<1> &rhs)
{
	if (this == &rhs) return *this; 
	LS = rhs.LS; return *this;
}

double 
Energy<1>::operator () (const double &r) const
{
        const Data *Prop = LS->Prop;
	double x = Prop->RMin/r;	
	bool isValid = true;
	if ( Prop->IsCompressiveOnly && x < 1.0 ) {
	  isValid = false;
	}

	if ( isValid ) {

	  // std::cout << Prop->m << "\t" << Prop->n << "\t"
	  // 	    << Prop->EMin << "\t" << Prop->RMin << "\t"
	  // 	    << r << "\t" << x << std::endl;
	  
	  const double & m = Prop->m;
	  const double & n = Prop->n;
	  return (m*Prop->EMin/r)*(pow(x,n) - pow(x,m));
	}
	else {
	  return 0.0;
	}	      
}

//////////////////////////////////////////////////////////////////////
// Class Energy<2>
//////////////////////////////////////////////////////////////////////

Energy<2>::Energy() {LS=0;}

Energy<2>::~Energy(){}

Potential::Radial::Energy<2> *
Energy<2>::Clone() const
{
	return new Energy<2>(*this);
}

Energy<2>::Energy(LocalState *rhs_LS) : LS(rhs_LS){}

Energy<2>::Energy(const Energy<2> &rhs) : LS(rhs.LS){}

const Energy<2> & 
Energy<2>::operator = (const Energy<2> &rhs)
{
	if (this == &rhs) return *this; 
	LS = rhs.LS; return *this;
}

double 
Energy<2>::operator () (const double &r) const
{
        const Data *Prop = LS->Prop;
	double x = Prop->RMin/r;
	bool isValid = true;
	if ( Prop->IsCompressiveOnly && x < 1.0 ) {
	  isValid = false;
	}

	if ( isValid ) {
	  const double & m = Prop->m;
	  const double & n = Prop->n;
	  return (m*LS->Prop->EMin/(r*r))*( (m+1)*pow(x,m) - (n+1)*pow(x,n));
	}
	else {
	  return 0.0;
	}
}

//////////////////////////////////////////////////////////////////////
// Class Jet<0>
//////////////////////////////////////////////////////////////////////

Jet<0>::Jet() {LS=0;}

Jet<0>::~Jet() {}

Potential::Radial::Jet<0> *
Jet<0>::Clone() const
{
	return new Jet<0>(*this);
}

Jet<0>::Jet(LocalState *rhs_LS) : LS(rhs_LS) {}

Jet<0>::Jet(const Jet<0> &rhs) : LS(rhs.LS) {}

Jet<0> & 
Jet<0>::operator = (const Jet<0> &rhs)
{
	if (&rhs == this) return *this;
	LS = rhs.LS; return *this;
}

pair<double,double>
Jet<0>::operator () (const double &r) const
{
        const Data *Prop = LS->Prop;
	const double & m = Prop->m;
	const double & n = Prop->n;
	double x = Prop->RMin/r;
	bool isValid = true;
	if ( Prop->IsCompressiveOnly && x < 1.0 ) {
	  isValid = false;
	}

	if ( isValid ) {
	  double xn = pow(x,n);
	  double xm = pow(x,m);
	  double V = Prop->EMin*(xm - m*xn/n);
	  double DV = (m*Prop->EMin/r)*(xn - xm);
	  return make_pair(V,DV);
	}
	else {
	  return make_pair(0.0, 0.0);
	}
}

//////////////////////////////////////////////////////////////////////
// Class Jet<1>
//////////////////////////////////////////////////////////////////////

Jet<1>::Jet() {LS=0;}

Jet<1>::~Jet() {}

Potential::Radial::Jet<1> *
Jet<1>::Clone() const
{
	return new Jet<1>(*this);
}

Jet<1>::Jet(LocalState *rhs_LS) : LS(rhs_LS) {}

Jet<1>::Jet(const Jet<1> &rhs) : LS(rhs.LS) {}

Jet<1> & 
Jet<1>::operator = (const Jet<1> &rhs)
{
	if (&rhs == this) return *this;
	LS = rhs.LS; return *this;
}

pair<double,double> 
Jet<1>::operator () (const double &r) const
{
        const Data *Prop = LS->Prop;
	const double & m = Prop->m;
	const double & n = Prop->n;
	double x = Prop->RMin/r;
	bool isValid = true;
	if ( Prop->IsCompressiveOnly && x < 1.0 ) {
	  isValid = false;
	}

	if ( isValid ) {
	  double xn = pow(x,n);
	  double xm = pow(x,m);
	  double A = m*Prop->EMin/r;
	  double DV = A*(xn - xm);
	  double DDV = (A/r)*((m+1) * xm - (n+1)*xn);
	  return make_pair(DV, DDV);
	}
	else {
	  return make_pair(0.0, 0.0);
	}
}

}

}

}
