// ThreeBody.cpp: implementation for the ThreeBody class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include "./ThreeBody.h"

namespace Potential
{
namespace ThreeBody
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

LocalState::LocalState(Set::Manifold::Point *rhs_e0, 
	Set::Manifold::Point *rhs_e1, Set::Manifold::Point *rhs_e2) : 
	e0(rhs_e0), e1(rhs_e1), e2(rhs_e2) {}

LocalState::LocalState(const set<Set::Manifold::Point *> &N) 
{
	set<Set::Manifold::Point *>::const_iterator pN = N.begin();
	e0 = *pN; pN++; e1 = *pN;  pN++; e2 = *pN; 
}

LocalState::LocalState(const LocalState &rhs)
	: e0(rhs.e0), e1(rhs.e1), e2(rhs.e2) {}

LocalState & 
LocalState::operator = (const LocalState &rhs)
{
	if (this == &rhs) return *this; 
	e0 = rhs.e0; e1 = rhs.e1; e2 = rhs.e2;
	return *this;
}

void 
LocalState::operator ++ () {}

set<Set::Manifold::Point *> 
LocalState::GetNodes() const
{
	set<Set::Manifold::Point *> Nodes;
	Nodes.insert(e0); 
	Nodes.insert(e1); 
	Nodes.insert(e2);
	return Nodes;
}

set<pair<Set::Manifold::Point *, Set::Manifold::Point *> > 
LocalState::GetNodePairs() const
{
	set<pair<Set::Manifold::Point *, Set::Manifold::Point *> > NodePairs;
	NodePairs.insert(make_pair(e0,e0));
	NodePairs.insert(make_pair(e0,e1));
	NodePairs.insert(make_pair(e0,e2));
	NodePairs.insert(make_pair(e1,e1));
	NodePairs.insert(make_pair(e1,e0));
	NodePairs.insert(make_pair(e1,e2));
	NodePairs.insert(make_pair(e2,e2));
	NodePairs.insert(make_pair(e2,e0));
	NodePairs.insert(make_pair(e2,e1));
	return NodePairs;
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

Energy<0>::Energy(
	LocalState *rhs_LS,
	Potential::Angular::Energy<0> *rhs_V) 
	: LS(rhs_LS), V(rhs_V) {}

Energy<0>::Energy(const Energy<0> &rhs) : V(rhs.V){}

const Energy<0> & 
Energy<0>::operator = (const Energy<0> &rhs)
{
	if (this == &rhs) return *this; 
	LS = rhs.LS; V = rhs.V; 
	return *this;
}

double 
Energy<0>::operator () (const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x0) const
{
	typedef Set::VectorSpace::Vector vector_type;
	typedef map<Set::Manifold::Point *, vector_type> map_type;
	map_type x = x0;
	Set::Manifold::Point *e0 = LS->e0; 
	Set::Manifold::Point *e1 = LS->e1; 
	Set::Manifold::Point *e2 = LS->e2;
	vector_type dx1 = x[e1] - x[e0];
	vector_type dx2 = x[e2] - x[e0];
	double r1 = Norm(dx1);
	double r2 = Norm(dx2);
	vector_type n1 = dx1/r1;
	vector_type n2 = dx2/r2;
	double c = n1(n2);
	return (*V)(c);
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

Energy<1>::Energy(
	LocalState *rhs_LS,
	Potential::Angular::Energy<1> *rhs_DV) 
	: LS(rhs_LS), DV(rhs_DV) {}

Energy<1>::Energy(const Energy<1> &rhs) : DV(rhs.DV){}

const Energy<1> & 
Energy<1>::operator = (const Energy<1> &rhs)
{
	if (this == &rhs) return *this; 
	LS = rhs.LS; DV = rhs.DV; 
	return *this;
}

map<Set::Manifold::Point *, Set::VectorSpace::Vector>
Energy<1>::operator () (const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x0) const
{
	typedef Set::VectorSpace::Vector vector_type;
	typedef map<Set::Manifold::Point *, vector_type> map_type;
	map_type x = x0;
	Set::Manifold::Point *e0 = LS->e0; 
	Set::Manifold::Point *e1 = LS->e1; 
	Set::Manifold::Point *e2 = LS->e2;
	vector_type dx1 = x[e1] - x[e0];
	vector_type dx2 = x[e2] - x[e0];
	double r1 = Norm(dx1);
	double r2 = Norm(dx2);
	vector_type n1 = dx1/r1;
	vector_type n2 = dx2/r2;
	double c = n1(n2);
	vector_type dc1 = (n2 - c*n1)/r1;
	vector_type dc2 = (n1 - c*n2)/r2;
	vector_type dc0 = -(dc1 + dc2);
	double dv = (*DV)(c); map_type DE;
	DE.insert(map_type::value_type(e0,dv*dc0));
	DE.insert(map_type::value_type(e1,dv*dc1));
	DE.insert(map_type::value_type(e2,dv*dc2));
	return DE;
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

Energy<2>::Energy(
	LocalState *rhs_LS,
	Potential::Angular::Energy<1> *rhs_DV, 
	Potential::Angular::Energy<2> *rhs_DDV)
	: LS(rhs_LS), DV(rhs_DV), DDV(rhs_DDV){}

Energy<2>::Energy(const Energy<2> &rhs) 
	: DV(rhs.DV), DDV(rhs.DDV){}

const Energy<2> & 
Energy<2>::operator = (const Energy<2> &rhs)
{
	if (this == &rhs) return *this; 
	LS = rhs.LS; DV = rhs.DV; DDV = rhs.DDV; 
	return *this;
}

map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, Set::VectorSpace::Hom>
Energy<2>::operator () (const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x0) const
{
	typedef Set::VectorSpace::Vector vector_type;
	typedef Set::VectorSpace::Hom hom_type;
	typedef map<Set::Manifold::Point *, vector_type> map_type;
	typedef map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, hom_type> dmap_type;
	map_type x = x0;
	Set::Manifold::Point *e0 = LS->e0; 
	Set::Manifold::Point *e1 = LS->e1; 
	Set::Manifold::Point *e2 = LS->e2;
	vector_type dx1 = x[e1] - x[e0];
	vector_type dx2 = x[e2] - x[e0];
	double r1 = Norm(dx1);
	double r2 = Norm(dx2);
	vector_type n1 = dx1/r1;
	vector_type n2 = dx2/r2;
	double c = n1(n2), c3 = 3.0*c;
	vector_type dc1 = (n2 - c*n1)/r1;
	vector_type dc2 = (n1 - c*n2)/r2;
	vector_type dc0 = -(dc1 + dc2);
	Set::VectorSpace::Hom n11 = Dyadic(n1,n1);
	Set::VectorSpace::Hom n12 = Dyadic(n1,n2);
	Set::VectorSpace::Hom n21 = Dyadic(n2,n1);
	Set::VectorSpace::Hom n22 = Dyadic(n2,n2);
	double r11 = r1*r1, r12 = r1*r2, r21 = r12, r22 = r2*r2; 
	Set::VectorSpace::Hom ddc11 = c3*n11 - n12 - n21;
	ddc11[0][0] -= c; ddc11[1][1] -= c; ddc11[2][2] -= c; ddc11 /= r11;
	Set::VectorSpace::Hom ddc22 = c3*n22 - n21 - n12;
	ddc22[0][0] -= c; ddc22[1][1] -= c; ddc22[2][2] -= c; ddc22 /= r22;
	Set::VectorSpace::Hom ddc12 = c*n12 - n11 - n22;
	ddc12[0][0] += 1.0; ddc12[1][1] += 1.0; ddc12[2][2] += 1.0; ddc12 /= r12;
	Set::VectorSpace::Hom ddc21 = c*n21 - n22 - n11;
	ddc21[0][0] += 1.0; ddc21[1][1] += 1.0; ddc21[2][2] += 1.0; ddc21 /= r21;
	Set::VectorSpace::Hom ddc01 = -(ddc11 + ddc21);
	Set::VectorSpace::Hom ddc10 = -(ddc11 + ddc12);
	Set::VectorSpace::Hom ddc02 = -(ddc12 + ddc22);
	Set::VectorSpace::Hom ddc20 = -(ddc21 + ddc22);
	Set::VectorSpace::Hom ddc00 =   ddc11 + ddc21 + ddc12 + ddc22;
	double dv = (*DV)(c), ddv = (*DDV)(c); dmap_type DDE;
	DDE.insert(dmap_type::value_type(make_pair(e0,e0),ddv*Dyadic(dc0,dc0)+dv*ddc00));
	DDE.insert(dmap_type::value_type(make_pair(e1,e0),ddv*Dyadic(dc0,dc1)+dv*ddc01));
	DDE.insert(dmap_type::value_type(make_pair(e2,e0),ddv*Dyadic(dc0,dc2)+dv*ddc02));
	DDE.insert(dmap_type::value_type(make_pair(e1,e1),ddv*Dyadic(dc1,dc1)+dv*ddc11));
	DDE.insert(dmap_type::value_type(make_pair(e0,e1),ddv*Dyadic(dc1,dc0)+dv*ddc10));
	DDE.insert(dmap_type::value_type(make_pair(e2,e1),ddv*Dyadic(dc1,dc2)+dv*ddc12));
	DDE.insert(dmap_type::value_type(make_pair(e2,e2),ddv*Dyadic(dc2,dc2)+dv*ddc22));
	DDE.insert(dmap_type::value_type(make_pair(e0,e2),ddv*Dyadic(dc2,dc0)+dv*ddc20));
	DDE.insert(dmap_type::value_type(make_pair(e1,e2),ddv*Dyadic(dc2,dc1)+dv*ddc21));
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

Jet<0>::Jet(
	LocalState *rhs_LS,
	Potential::Angular::Jet<0> *rhs_J)
	: LS(rhs_LS), J(rhs_J) {}

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
Jet<0>::operator () (const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x0) const
{
	typedef Set::VectorSpace::Vector vector_type;
	typedef map<Set::Manifold::Point *, vector_type> map_type;
	map_type x = x0;
	Set::Manifold::Point *e0 = LS->e0; 
	Set::Manifold::Point *e1 = LS->e1; 
	Set::Manifold::Point *e2 = LS->e2;
	vector_type dx1 = x[e1] - x[e0];
	vector_type dx2 = x[e2] - x[e0];
	double r1 = Norm(dx1);
	double r2 = Norm(dx2);
	vector_type n1 = dx1/r1;
	vector_type n2 = dx2/r2;
	double c = n1(n2);
	vector_type dc1 = (n2 - c*n1)/r1;
	vector_type dc2 = (n1 - c*n2)/r2;
	vector_type dc0 = -(dc1 + dc2);
	pair<double,double> j = (*J)(c);
	double dv = j.second; map_type DE;
	DE.insert(map_type::value_type(e0,dv*dc0));
	DE.insert(map_type::value_type(e1,dv*dc1));
	DE.insert(map_type::value_type(e2,dv*dc2));
	return make_pair(j.first,DE);
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

Jet<1>::Jet(
	LocalState *rhs_LS,
	Potential::Angular::Jet<1> *rhs_DJ)
	: LS(rhs_LS), DJ(rhs_DJ) {}

Jet<1>::Jet(const Jet<1> &rhs)
	: LS(rhs.LS), DJ(rhs.DJ) {}

Jet<1> & 
Jet<1>::operator = (const Jet<1> &rhs)
{
	if (this == &rhs) return *this;
	LS = rhs.LS; DJ = rhs.DJ;
	return *this;
}

pair<
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>, 
	map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, Set::VectorSpace::Hom> >
Jet<1>::operator () (const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x0) const
{
	typedef Set::VectorSpace::Vector vector_type;
	typedef Set::VectorSpace::Hom hom_type;
	typedef map<Set::Manifold::Point *, vector_type> map_type;
	typedef map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, hom_type> dmap_type;
	map_type x = x0;
	Set::Manifold::Point *e0 = LS->e0; 
	Set::Manifold::Point *e1 = LS->e1; 
	Set::Manifold::Point *e2 = LS->e2;
	vector_type dx1 = x[e1] - x[e0];
	vector_type dx2 = x[e2] - x[e0];
	double r1 = Norm(dx1);
	double r2 = Norm(dx2);
	vector_type n1 = dx1/r1;
	vector_type n2 = dx2/r2;
	double c = n1(n2), c3 = 3.0*c;
	vector_type dc1 = (n2 - c*n1)/r1;
	vector_type dc2 = (n1 - c*n2)/r2;
	vector_type dc0 = -(dc1 + dc2);
	Set::VectorSpace::Hom n11 = Dyadic(n1,n1);
	Set::VectorSpace::Hom n12 = Dyadic(n1,n2);
	Set::VectorSpace::Hom n21 = Dyadic(n2,n1);
	Set::VectorSpace::Hom n22 = Dyadic(n2,n2);
	double r11 = r1*r1, r12 = r1*r2, r21 = r12, r22 = r2*r2; 
	Set::VectorSpace::Hom ddc11 = c3*n11 - n12 - n21;
	ddc11[0][0] -= c; ddc11[1][1] -= c; ddc11[2][2] -= c; ddc11 /= r11;
	Set::VectorSpace::Hom ddc22 = c3*n22 - n21 - n12;
	ddc22[0][0] -= c; ddc22[1][1] -= c; ddc22[2][2] -= c; ddc22 /= r22;
	Set::VectorSpace::Hom ddc12 = c*n12 - n11 - n22;
	ddc12[0][0] += 1.0; ddc12[1][1] += 1.0; ddc12[2][2] += 1.0; ddc12 /= r12;
	Set::VectorSpace::Hom ddc21 = c*n21 - n22 - n11;
	ddc21[0][0] += 1.0; ddc21[1][1] += 1.0; ddc21[2][2] += 1.0; ddc21 /= r21;
	Set::VectorSpace::Hom ddc01 = -(ddc11 + ddc21);
	Set::VectorSpace::Hom ddc10 = -(ddc11 + ddc12);
	Set::VectorSpace::Hom ddc02 = -(ddc12 + ddc22);
	Set::VectorSpace::Hom ddc20 = -(ddc21 + ddc22);
	Set::VectorSpace::Hom ddc00 = -(ddc11 + ddc21 + ddc12 + ddc22);
	pair<double,double> dj = (*DJ)(c);
	double dv = dj.first; map_type DE;
	double ddv = dj.second; dmap_type DDE;
	DE.insert(map_type::value_type(e0,dv*dc0));
	DE.insert(map_type::value_type(e1,dv*dc1));
	DE.insert(map_type::value_type(e2,dv*dc2));
	DDE.insert(dmap_type::value_type(make_pair(e0,e0),ddv*Dyadic(dc0,dc0)+dv*ddc00));
	DDE.insert(dmap_type::value_type(make_pair(e1,e0),ddv*Dyadic(dc0,dc1)+dv*ddc01));
	DDE.insert(dmap_type::value_type(make_pair(e2,e0),ddv*Dyadic(dc0,dc2)+dv*ddc02));
	DDE.insert(dmap_type::value_type(make_pair(e1,e1),ddv*Dyadic(dc1,dc1)+dv*ddc11));
	DDE.insert(dmap_type::value_type(make_pair(e0,e1),ddv*Dyadic(dc1,dc0)+dv*ddc10));
	DDE.insert(dmap_type::value_type(make_pair(e2,e1),ddv*Dyadic(dc1,dc2)+dv*ddc12));
	DDE.insert(dmap_type::value_type(make_pair(e2,e2),ddv*Dyadic(dc2,dc2)+dv*ddc22));
	DDE.insert(dmap_type::value_type(make_pair(e0,e2),ddv*Dyadic(dc2,dc0)+dv*ddc20));
	DDE.insert(dmap_type::value_type(make_pair(e1,e2),ddv*Dyadic(dc2,dc1)+dv*ddc21));
	return make_pair(DE,DDE);
}

}

}
