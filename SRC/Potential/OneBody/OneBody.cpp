// OneBody.cpp: implementation for the OneBody class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include "./OneBody.h"

namespace Potential
{
namespace OneBody
{
//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

LocalState::LocalState() {}

LocalState::~LocalState() {}

Element::LocalState *
LocalState::Clone() const
{
	return new LocalState(*this);
}

LocalState::LocalState(Set::Manifold::Point *e0_) : e0(e0_) {}

LocalState::LocalState(const LocalState &rhs) : e0(rhs.e0) {}

LocalState & 
LocalState::operator = (const LocalState &rhs)
{
	if (this == &rhs) return *this; 
	e0 = rhs.e0; return *this;
}

void 
LocalState::operator ++ () {}

set<Set::Manifold::Point *> 
LocalState::GetNodes() const
{
	set<Set::Manifold::Point *> Nodes;
	Nodes.insert(e0); return Nodes;
}

void
LocalState::GetNodes(set<Set::Manifold::Point *>  & Nodes) const
{
  if ( !Nodes.empty()) Nodes.clear();
  Nodes.insert(e0); 
  return;
}

set<pair<Set::Manifold::Point *, Set::Manifold::Point *> > 
LocalState::GetNodePairs() const
{
	set<pair<Set::Manifold::Point *, Set::Manifold::Point *> > NodePairs;
	NodePairs.insert(make_pair(e0,e0)); return NodePairs;
}

//////////////////////////////////////////////////////////////////////
// Class Energy<0>
//////////////////////////////////////////////////////////////////////

Energy<0>::Energy() {}

Energy<0>::~Energy() {}

Element::Energy<0> *
Energy<0>::Clone() const
{
	return new Energy<0>(*this);
}

Element::LocalState *
Energy<0>::GetLocalState() const
{
	return LS;
}

Energy<0>::Energy(LocalState *LS_, 
	Potential::Field::Energy<0> *V_) : LS(LS_), V(V_) {}

Energy<0>::Energy(const Energy<0> &rhs) : V(rhs.V){}

const Energy<0> & 
Energy<0>::operator = (const Energy<0> &rhs)
{
	if (this == &rhs) return *this; 
	LS = rhs.LS; V = rhs.V; 
	return *this;
}

double 
Energy<0>::operator () (
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x) const 
{
	return (*V)(x.begin()->second);
}

//////////////////////////////////////////////////////////////////////
// Class Energy<1>
//////////////////////////////////////////////////////////////////////

Energy<1>::Energy() {}

Energy<1>::~Energy() {}

Element::Energy<1> *
Energy<1>::Clone() const
{
	return new Energy<1>(*this);
}

Element::LocalState *
Energy<1>::GetLocalState() const
{
	return LS;
}

Energy<1>::Energy(LocalState *LS_,
	Potential::Field::Energy<1> *DV_) : LS(LS_), DV(DV_) {
    Element::Energy<1>::_ELS = LS;
}

Energy<1>::Energy(const Energy<1> &rhs) : DV(rhs.DV){}

const Energy<1> & 
Energy<1>::operator = (const Energy<1> &rhs)
{
	if (this == &rhs) return *this; 
	LS = rhs.LS; DV = rhs.DV; 
	return *this;
}

map<Set::Manifold::Point *, Set::VectorSpace::Vector>
Energy<1>::operator () (
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x) const
{
	map<Set::Manifold::Point *, Set::VectorSpace::Vector> DE;
	DE.insert(make_pair(x.begin()->first,(*DV)(x.begin()->second)));
	return DE;
}

void
Energy<1>::operator() (const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x,
		     map<Set::Manifold::Point *, Set::VectorSpace::Vector> & DE) const
{
  if ( !DE.empty() ) DE.clear();

  DE.insert(make_pair(x.begin()->first,(*DV)(x.begin()->second)));
  return;
}

//////////////////////////////////////////////////////////////////////
// Class Energy<2>
//////////////////////////////////////////////////////////////////////

Energy<2>::Energy() {}

Energy<2>::~Energy() {}

Element::Energy<2> *
Energy<2>::Clone() const
{
	return new Energy<2>(*this);
}

Element::LocalState *
Energy<2>::GetLocalState() const
{
	return LS;
}

Energy<2>::Energy(LocalState *LS_,
	Potential::Field::Energy<2> *DDV_) : LS(LS_), DDV(DDV_){}

Energy<2>::Energy(const Energy<2> &rhs) : DDV(rhs.DDV){}

const Energy<2> & 
Energy<2>::operator = (const Energy<2> &rhs)
{
	if (this == &rhs) return *this; 
	LS = rhs.LS; DDV = rhs.DDV; 
	return *this;
}

map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, Set::VectorSpace::Hom>
Energy<2>::operator () (const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x) const
{
	Set::Manifold::Point *e0 = x.begin()->first;
	map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, Set::VectorSpace::Hom> DDE;
	DDE.insert(make_pair(make_pair(e0,e0),(*DDV)(x.begin()->second)));
	return DDE;
}

//////////////////////////////////////////////////////////////////////
// Class Jet<0>
//////////////////////////////////////////////////////////////////////

Jet<0>::Jet() {}

Jet<0>::~Jet() {} 

Element::Jet<0> *
Jet<0>::Clone() const
{
	return new Jet<0>(*this);
}

Element::LocalState *
Jet<0>::GetLocalState() const
{
	return LS;
}

Jet<0>::Jet(LocalState *LS_,
	Potential::Field::Jet<0> *J_) : LS(LS_), J(J_) {}

Jet<0>::Jet(const Jet<0> &rhs)
	: LS(rhs.LS), J(rhs.J) {}

Jet<0> & 
Jet<0>::operator = (const Jet<0> &rhs)
{
	if (this == &rhs) return *this;
	LS = rhs.LS; J = rhs.J;
	return *this;
}

pair<double, map<Set::Manifold::Point *, Set::VectorSpace::Vector> > 
Jet<0>::operator () (
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x) const
{
	pair<double, Set::VectorSpace::Vector> K = (*J)(x.begin()->second);
	map<Set::Manifold::Point *, Set::VectorSpace::Vector> DE;
	DE.insert(make_pair(x.begin()->first,K.second));
	return make_pair(K.first,DE);
}

//////////////////////////////////////////////////////////////////////
// Class Jet<1>
//////////////////////////////////////////////////////////////////////

Jet<1>::Jet() {}

Jet<1>::~Jet() {}

Element::Jet<1> *
Jet<1>::Clone() const
{
	return new Jet<1>(*this);
}

Element::LocalState *
Jet<1>::GetLocalState() const
{
	return LS;
}

Jet<1>::Jet(LocalState *LS_,
	Potential::Field::Jet<1> *DJ_) : LS(LS_), DJ(DJ_) {}

Jet<1>::Jet(const Jet<1> &rhs)
	: LS(rhs.LS), DJ(rhs.DJ) {}

Jet<1> & 
Jet<1>::operator = (const Jet<1> &rhs)
{
	if (this == &rhs) return *this;
	LS = rhs.LS; DJ = rhs.DJ;
	return *this;
}

pair<map<Set::Manifold::Point *, Set::VectorSpace::Vector>, 
	map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, Set::VectorSpace::Hom> >
Jet<1>::operator () (const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x) const
{
	Set::Manifold::Point *e0 = x.begin()->first;
	pair<Set::VectorSpace::Vector, Set::VectorSpace::Hom> DK = (*DJ)(x.begin()->second);
	map<Set::Manifold::Point *, Set::VectorSpace::Vector> DE;
	DE.insert(make_pair(e0,DK.first));
	map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, Set::VectorSpace::Hom> DDE; 
	DDE.insert(make_pair(make_pair(e0,e0),DK.second));
	return make_pair(DE,DDE);
}

}

}
