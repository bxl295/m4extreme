// FK2LK.cpp: implementation for the FK2LK class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <cassert>
#include "./FK2LK.h"

namespace Material
{
namespace Shear
{
namespace FK2LK
{
//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

LocalState::~LocalState() {}

Material::Shear::LocalState *
LocalState::Clone() const 
{
	return new LocalState(*this);
}

LocalState::LocalState(Material::Shear::LocalState *LS_) : LS(LS_) {}

LocalState::LocalState(const LocalState &rhs) : LS(rhs.LS) {}

void 
LocalState::operator ++ () 
{
	++(*LS);
}

//////////////////////////////////////////////////////////////////////
// Class Modulus<0>
//////////////////////////////////////////////////////////////////////

Modulus<0>::~Modulus(){}

Material::Shear::Modulus<0> *
Modulus<0>::Clone() const
{
	return new Modulus<0>(*this);
}

Modulus<0>::Modulus(Material::Shear::Modulus<0> *f_) : f(f_){}

Modulus<0>::Modulus(const Modulus<0> &rhs) : f(rhs.f){}

double 
Modulus<0>::operator () (const double &Theta) const
{
	return (*f)(exp(Theta));
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

Modulus<1>::Modulus(Material::Shear::Modulus<1> *Df_) : Df(Df_){}

Modulus<1>::Modulus(const Modulus<1> &rhs) : Df(rhs.Df){}

double
Modulus<1>::operator () (const double &Theta) const
{
	double J = exp(Theta);
	return (*Df)(J)*J;
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

Modulus<2>::Modulus(Material::Shear::Jet<1> *Dg_) : Dg(Dg_){}

Modulus<2>::Modulus(const Modulus<2> &rhs) : Dg(rhs.Dg){}

double 
Modulus<2>::operator () (const double &Theta) const
{
	double J = exp(Theta);
	pair<double,double> DK = (*Dg)(J);
	return (DK.first + DK.second*J)*J;
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

Jet<0>::Jet(Material::Shear::Jet<0> *g_) : g(g_){}

Jet<0>::Jet(const Jet<0> &rhs) : g(rhs.g){}

pair<double,double>
Jet<0>::operator () (const double &Theta) const 
{
	double J = exp(Theta);
	pair<double,double> K = (*g)(J);
	return make_pair(K.first,K.second*J);
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

Jet<1>::Jet(Material::Shear::Jet<1> *Dg_) : Dg(Dg_){}

Jet<1>::Jet(const Jet<1> &rhs) : Dg(rhs.Dg){}

pair<double,double>
Jet<1>::operator () (const double &Theta) const
{
	double J = exp(Theta);
	pair<double,double> DK = (*Dg)(J);
	return make_pair(DK.first*J,(DK.first + DK.second*J)*J);
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

JetJet<0>::JetJet(Material::Shear::JetJet<0> *h_) : h(h_){}

JetJet<0>::JetJet(const JetJet<0> &rhs) : h(rhs.h){}

triplet<double,double,double>
JetJet<0>::operator () (const double &Theta) const 
{
	double J = exp(Theta);
	triplet<double,double,double> KK = (*h)(J);
	return make_triplet(KK.first,KK.second*J,(KK.second + KK.third*J)*J);
}

}

}

}
