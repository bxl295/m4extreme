// LK2FK.cpp: implementation for the LK2FK class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <cassert>
#include "./LK2FK.h"

namespace Material
{
namespace Gas
{
namespace EoS
{
namespace LK2FK
{
//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

LocalState::~LocalState() {}

Material::Gas::EoS::LocalState *
LocalState::Clone() const 
{
	return new LocalState(*this);
}

LocalState::LocalState(Material::Gas::EoS::LocalState *LS_) : LS(LS_) {}

LocalState::LocalState(const LocalState &rhs) : LS(rhs.LS) {}

void 
LocalState::operator ++ () 
{
	++(*LS);
}

//////////////////////////////////////////////////////////////////////
// Class Energy<0>
//////////////////////////////////////////////////////////////////////

Energy<0>::~Energy(){}

Material::Gas::EoS::Energy<0> *
Energy<0>::Clone() const
{
	return new Energy<0>(*this);
}

Energy<0>::Energy(Material::Gas::EoS::Energy<0> *f_) : f(f_){}

Energy<0>::Energy(const Energy<0> &rhs) : f(rhs.f){}

double 
Energy<0>::operator () (const double & J) const
{
  return (*f)(log(J));
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

Energy<1>::Energy(Material::Gas::EoS::Energy<1> *Df_) : Df(Df_){}

Energy<1>::Energy(const Energy<1> &rhs) : Df(rhs.Df){}

double
Energy<1>::operator () (const double & J) const
{
	double theta = log(J);
	return (*Df)(theta) / J;
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

Energy<2>::Energy(Material::Gas::EoS::Jet<1> *Dg_) : Dg(Dg_){}

Energy<2>::Energy(const Energy<2> &rhs) : Dg(rhs.Dg){}

double 
Energy<2>::operator () (const double & J) const
{
  double theta = log(J);
  pair<double,double> DK = (*Dg)(theta);
  return (DK.second - DK.first) / J / J;
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

Jet<0>::Jet(Material::Gas::EoS::Jet<0> *g_) : g(g_){}

Jet<0>::Jet(const Jet<0> &rhs) : g(rhs.g){}

pair<double,double>
Jet<0>::operator () (const double & J) const 
{
	double theta = log(J);
	pair<double,double> K = (*g)(theta);
	return make_pair(K.first, K.second/J);
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

Jet<1>::Jet(Material::Gas::EoS::Jet<1> *Dg_) : Dg(Dg_){}

Jet<1>::Jet(const Jet<1> &rhs) : Dg(rhs.Dg){}

pair<double,double>
Jet<1>::operator () (const double & J) const
{
        double theta = log(J);
	pair<double,double> DK = (*Dg)(theta);
	return make_pair(DK.first/J,(DK.second - DK.first)/J/J);
}

}

}

}

}
