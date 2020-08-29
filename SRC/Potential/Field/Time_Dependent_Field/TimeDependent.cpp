// TimeDependent.cpp: implementation of the TimeDependent class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <cassert>
#include "./TimeDependent.h"

namespace Potential
{
namespace Field
{
namespace TimeDependent
{
//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

LocalState::~LocalState() {}

Potential::Field::LocalState *
LocalState::Clone() const 
{
	return new LocalState(*this);
}

LocalState::LocalState(double QW_, TimeFunction *Prop_) : Prop(Prop_), QW(QW_) {}

LocalState::LocalState(const LocalState &rhs) : Prop(rhs.Prop), QW(rhs.QW) {}

void 
LocalState::operator ++ () {}

//////////////////////////////////////////////////////////////////////
// Class Energy<0>
//////////////////////////////////////////////////////////////////////

Energy<0>::~Energy() {}

Potential::Field::Energy<0> *
Energy<0>::Clone() const
{
	return new Energy<0>(*this);
}

Energy<0>::Energy(LocalState *LS_) : LS(LS_) {}

Energy<0>::Energy(const Energy<0> &rhs) : LS(rhs.LS) {}

double 
Energy<0>::operator () (const Set::VectorSpace::Vector &A) const
{
  const Set::VectorSpace::Vector & B = LS->Prop->GetDensity();
  assert (A.size() == B.size());
  return -B(A) * LS->GetQW();
}

//////////////////////////////////////////////////////////////////////
// Class Energy<1>
//////////////////////////////////////////////////////////////////////

Energy<1>::~Energy() {}

Potential::Field::Energy<1> *
Energy<1>::Clone() const
{
	return new Energy<1>(*this);
}

Energy<1>::Energy(LocalState *LS_) : LS(LS_) {}

Energy<1>::Energy(const Energy<1> &rhs) : LS(rhs.LS) {}

Set::VectorSpace::Vector 
Energy<1>::operator () (const Set::VectorSpace::Vector &A) const
{
  const Set::VectorSpace::Vector & B = LS->Prop->GetDensity();
  assert (A.size() == B.size());
  return -B * LS->GetQW();
}

//////////////////////////////////////////////////////////////////////
// Class Energy<2>
//////////////////////////////////////////////////////////////////////

Energy<2>::~Energy() {}

Potential::Field::Energy<2> *
Energy<2>::Clone() const
{
	return new Energy<2>(*this);
}

Energy<2>::Energy(LocalState *LS_) : LS(LS_) {}

Energy<2>::Energy(const Energy<2> &rhs) : LS(rhs.LS) {}

Set::VectorSpace::Hom 
Energy<2>::operator () (const Set::VectorSpace::Vector &A) const
{
  const Set::VectorSpace::Vector & B = LS->Prop->GetDensity();
  assert (A.size() == B.size());
  Set::VectorSpace::Hom DDW(B.size());
  return DDW;
}

//////////////////////////////////////////////////////////////////////
// Class Jet<0>
//////////////////////////////////////////////////////////////////////

Jet<0>::~Jet() {}

Potential::Field::Jet<0> *
Jet<0>::Clone() const
{
	return new Jet<0>(*this);
}

Jet<0>::Jet(LocalState *LS_) : LS(LS_) {}

Jet<0>::Jet(const Jet<0> &rhs) : LS(rhs.LS) {}

pair<double,Set::VectorSpace::Vector>
Jet<0>::operator () (const Set::VectorSpace::Vector &A) const
{
  const Set::VectorSpace::Vector & B = LS->Prop->GetDensity();
  double qw = LS->GetQW();
  assert (A.size() == B.size());
  return make_pair(-B(A) * qw, -B * qw);
}

//////////////////////////////////////////////////////////////////////
// Class Jet<1>
//////////////////////////////////////////////////////////////////////

Jet<1>::~Jet() {}

Potential::Field::Jet<1> *
Jet<1>::Clone() const
{
	return new Jet<1>(*this);
}

Jet<1>::Jet(LocalState *LS_) : LS(LS_) {}

Jet<1>::Jet(const Jet<1> &rhs) : LS(rhs.LS) {}

pair<Set::VectorSpace::Vector,Set::VectorSpace::Hom> 
Jet<1>::operator () (const Set::VectorSpace::Vector &A) const
{
  const Set::VectorSpace::Vector & B = LS->Prop->GetDensity();
  assert (A.size() == B.size());
  Set::VectorSpace::Hom DDW(B.size());
  return make_pair(-B * LS->GetQW(), DDW);
}

}

}

}
