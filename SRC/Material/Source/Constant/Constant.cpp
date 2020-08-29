// Constant.cpp: implementation of the Constant class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <cassert>
#include "./Constant.h"

namespace Material
{
namespace Source
{
namespace Constant
{
//////////////////////////////////////////////////////////////////////
// Class Data
//////////////////////////////////////////////////////////////////////

Data::Data()  {}

Data::Data(const unsigned int &n) : B(n) {}

Data::Data(const Set::VectorSpace::Vector &rhs_B) : B(rhs_B) {}

Data::~Data(){}

Data::Data(const Data &rhs) : B(rhs.B) {}

Data &
Data::operator = (const Data &rhs)
{
	if (this == &rhs) return *this;
	B = rhs.B; return *this;
}

const Set::VectorSpace::Vector &
Data::GetDensity() const
{
	return B;
}

Set::VectorSpace::Vector & 
Data::GetDensity()
{
	return B;
}

void 
Data::Randomize()
{
	B.Randomize();
}

//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

LocalState::~LocalState() {}

Material::LocalState *
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

Energy<0>::~Energy() {}

Material::Energy<0> *
Energy<0>::Clone() const
{
	return new Energy<0>(*this);
}

Energy<0>::Energy(LocalState *LS_) : LS(LS_) {}

Energy<0>::Energy(const Energy<0> &rhs) : LS(rhs.LS) {}

double 
Energy<0>::operator () (const Set::VectorSpace::Vector &A) const
{
	assert (A.size() == LS->Prop->B.size());
	return -LS->Prop->B(A);
}

//////////////////////////////////////////////////////////////////////
// Class Energy<1>
//////////////////////////////////////////////////////////////////////

Energy<1>::~Energy() {}

Material::Energy<1> *
Energy<1>::Clone() const
{
	return new Energy<1>(*this);
}

Energy<1>::Energy(LocalState *LS_) : LS(LS_) {}

Energy<1>::Energy(const Energy<1> &rhs) : LS(rhs.LS) {}

Set::VectorSpace::Vector 
Energy<1>::operator () (const Set::VectorSpace::Vector &A) const
{
	assert (A.size() == LS->Prop->B.size());
	return -LS->Prop->B;
}

//////////////////////////////////////////////////////////////////////
// Class Energy<2>
//////////////////////////////////////////////////////////////////////

Energy<2>::~Energy() {}

Material::Energy<2> *
Energy<2>::Clone() const
{
	return new Energy<2>(*this);
}

Energy<2>::Energy(LocalState *LS_) : LS(LS_) {}

Energy<2>::Energy(const Energy<2> &rhs) : LS(rhs.LS) {}

Set::VectorSpace::Hom 
Energy<2>::operator () (const Set::VectorSpace::Vector &A) const
{
	assert (A.size() == LS->Prop->B.size());
	Set::VectorSpace::Hom DDW(LS->Prop->B.size());
	return DDW;
}

//////////////////////////////////////////////////////////////////////
// Class Jet<0>
//////////////////////////////////////////////////////////////////////

Jet<0>::~Jet() {}

Material::Jet<0> *
Jet<0>::Clone() const
{
	return new Jet<0>(*this);
}

Jet<0>::Jet(LocalState *LS_) : LS(LS_) {}

Jet<0>::Jet(const Jet<0> &rhs) : LS(rhs.LS) {}

pair<double,Set::VectorSpace::Vector>
Jet<0>::operator () (const Set::VectorSpace::Vector &A) const
{
	assert (A.size() == LS->Prop->B.size());
	return make_pair(-LS->Prop->B(A),-LS->Prop->B);
}

//////////////////////////////////////////////////////////////////////
// Class Jet<1>
//////////////////////////////////////////////////////////////////////

Jet<1>::~Jet() {}

Material::Jet<1> *
Jet<1>::Clone() const
{
	return new Jet<1>(*this);
}

Jet<1>::Jet(LocalState *LS_) : LS(LS_) {}

Jet<1>::Jet(const Jet<1> &rhs) : LS(rhs.LS) {}

pair<Set::VectorSpace::Vector,Set::VectorSpace::Hom> 
Jet<1>::operator () (const Set::VectorSpace::Vector &A) const
{
	assert (A.size() == LS->Prop->B.size());
	Set::VectorSpace::Hom DDW(LS->Prop->B.size());
	return make_pair(-LS->Prop->B,DDW);
}

}

}

}
