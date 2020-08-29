// TwoBody.cpp: implementation for the TwoBody class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include "./TwoBody.h"

namespace Potential
{
namespace TwoBody
{
//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

LocalState::LocalState():LS(NULL) {}

LocalState::~LocalState() {}

Element::LocalState *
LocalState::Clone() const
{
	return new LocalState(*this);
}

LocalState::LocalState(Set::Manifold::Point *rhs_e1, Set::Manifold::Point *rhs_e2,
		       Potential::Radial::LocalState *rhs_LS) 
  : e1(rhs_e1), e2(rhs_e2), LS(rhs_LS) {}

LocalState::LocalState(const set<Set::Manifold::Point *> &N,
		       Potential::Radial::LocalState *rhs_LS) 
  : LS(rhs_LS)
{
	set<Set::Manifold::Point *>::const_iterator pN = N.begin();
	e1 = *pN; pN++; e2 = *pN; 
}

LocalState::LocalState(const LocalState &rhs)
  : e1(rhs.e1), e2(rhs.e2), LS(rhs.LS) {}

LocalState & 
LocalState::operator = (const LocalState &rhs)
{
	if (this == &rhs) return *this; 
	e1 = rhs.e1; e2 = rhs.e2; LS = rhs.LS;
	return *this;
}

void 
LocalState::operator ++ () {
  if ( LS != NULL ) {
    ++(*LS);
  }
}

Potential::Radial::LocalState * 
LocalState::GetLS() {
  return LS;
}

set<Set::Manifold::Point *> 
LocalState::GetNodes() const
{
	set<Set::Manifold::Point *> Nodes;
	Nodes.insert(e1); 
	Nodes.insert(e2);
	return Nodes;
}

set<pair<Set::Manifold::Point *, Set::Manifold::Point *> > 
LocalState::GetNodePairs() const
{
	set<pair<Set::Manifold::Point *, Set::Manifold::Point *> > NodePairs;
	NodePairs.insert(make_pair(e1,e1));
	NodePairs.insert(make_pair(e1,e2));
	NodePairs.insert(make_pair(e2,e1));
	NodePairs.insert(make_pair(e2,e2));
	return NodePairs;
}

void
LocalState::GetNodes(set<Set::Manifold::Point *> & Nodes) const
{
        if ( !Nodes.empty()) Nodes.clear();
	Nodes.insert(e1); 
	Nodes.insert(e2);
	return;
}

void
LocalState::GetNodePairs(set<pair<Set::Manifold::Point *, Set::Manifold::Point *> > &NodePairs) const
{
        if ( !NodePairs.empty()) NodePairs.clear();
	NodePairs.insert(make_pair(e1,e1));
	NodePairs.insert(make_pair(e1,e2));
	NodePairs.insert(make_pair(e2,e1));
	NodePairs.insert(make_pair(e2,e2));
	return;
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
	Potential::Radial::Energy<0> *rhs_V) 
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
Energy<0>::operator () (const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x) const 
{
	typedef Set::VectorSpace::Vector vector_type;
	typedef map<Set::Manifold::Point *, vector_type> map_type;
	map_type::const_iterator px1 = x.begin();
	map_type::const_iterator px2 = px1++;
	Set::VectorSpace::Vector dx = px2->second - px1->second;
	return (*V)(Norm(dx));
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
	Potential::Radial::Energy<1> *rhs_DV) 
	: LS(rhs_LS), DV(rhs_DV) {
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
Energy<1>::operator () (const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x) const
{
	typedef Set::VectorSpace::Vector vector_type;
	typedef map<Set::Manifold::Point *, vector_type> map_type;
	map_type::const_iterator px1 = x.begin();
	map_type::const_iterator px2 = px1++;	
	vector_type dx = px2->second - px1->second;
	double r = Norm(dx);	
	double C = (*DV)(r)/r;
	vector_type f = C*dx; map_type DE;

	//suppose repulsive force is negative
	DE.insert(map_type::value_type(px1->first,-f));
	DE.insert(map_type::value_type(px2->first,f));
	return DE;
}

void
Energy<1>::operator () (const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x,
			map<Set::Manifold::Point *, Set::VectorSpace::Vector> & DE) const
{
        if ( !DE.empty() ) DE.clear();
	typedef Set::VectorSpace::Vector vector_type;
	typedef map<Set::Manifold::Point *, vector_type> map_type;
	map_type::const_iterator px1 = x.begin();
	map_type::const_iterator px2 = px1++;
	vector_type dx = px2->second - px1->second;
	double r = Norm(dx);
	double C = (*DV)(r)/r;
	vector_type f = C*dx;
	
	//suppose repulsive force is negative
	DE.insert(map_type::value_type(px1->first,-f));
	DE.insert(map_type::value_type(px2->first,f));

	// cout << "r=" << r << endl;
	// for ( map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator pDE = DE.begin();
	//       pDE != DE.end(); pDE++ ) {
	//   cout << x.find(pDE->first)->second << "--->\t" << pDE->second << endl;
	// }
	// cout << endl;
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

Energy<2>::Energy(
	LocalState *rhs_LS,
	Potential::Radial::Energy<1> *rhs_DV, 
	Potential::Radial::Energy<2> *rhs_DDV)
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
Energy<2>::operator () (const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x) const
{
	typedef Set::VectorSpace::Vector vector_type;
	typedef map<Set::Manifold::Point *, vector_type> map_type;
	typedef map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, Set::VectorSpace::Hom> dmap_type;
	map_type::const_iterator px1 = x.begin();
	map_type::const_iterator px2 = px1++;
	Set::VectorSpace::Vector dx = px2->second - px1->second;
	double r = Norm(dx); 
	double C = (*DV)(r)/r; 
	double B = (*DDV)(r); 
	double A = B - C;
	unsigned int i, j, n = dx.size();
	Set::VectorSpace::Hom df(n);
	Set::VectorSpace::Vector a = dx/r;
	Set::VectorSpace::Vector b = A*a;
	for (j=0; j<n; j++)
	{
		for (i=0; i<n; i++) df[j][i] = b[j]*a[i];
		df[j][j] += C;
	}
	Set::VectorSpace::Hom dg = -df;
	dmap_type DDE;
	Set::Manifold::Point *e1 = px1->first;
	Set::Manifold::Point *e2 = px2->first;
	DDE.insert(dmap_type::value_type(make_pair(e1,e1),df));
	DDE.insert(dmap_type::value_type(make_pair(e1,e2),dg));
	DDE.insert(dmap_type::value_type(make_pair(e2,e1),dg));
	DDE.insert(dmap_type::value_type(make_pair(e2,e2),df));
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
	Potential::Radial::Jet<0> *rhs_J)
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
Jet<0>::operator () (const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x) const
{
	typedef Set::VectorSpace::Vector vector_type;
	typedef map<Set::Manifold::Point *, vector_type> map_type;
	map_type::const_iterator px1 = x.begin();
	map_type::const_iterator px2 = px1++;
	vector_type dx = px2->second - px1->second;
	double r = Norm(dx); 
	pair<double,double> K = (*J)(r);
	double E = K.first;
	map_type DE; 
	double C = K.second/r;
	vector_type f = C*dx;
	
	//suppose repulsive force is negative
	DE.insert(map_type::value_type(px1->first,-f));
	DE.insert(map_type::value_type(px2->first,f));
	return make_pair(E,DE);
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
	Potential::Radial::Jet<1> *rhs_DJ)
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

pair<map<Set::Manifold::Point *, Set::VectorSpace::Vector>, 
	map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, Set::VectorSpace::Hom> >
Jet<1>::operator () (const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x) const
{
	typedef Set::VectorSpace::Vector vector_type;
	typedef map<Set::Manifold::Point *, vector_type> map_type;
	typedef map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, Set::VectorSpace::Hom> dmap_type;
	map_type::const_iterator px1 = x.begin();
	map_type::const_iterator px2 = px1++;
	vector_type dx = px2->second - px1->second;
	double r = Norm(dx); 
	pair<double,double> DK = (*DJ)(r);
	map_type DE; 
	double C = DK.first/r;
	vector_type f = C*dx;
	
	//suppose repulsive force is negative
	DE.insert(map_type::value_type(px1->first,-f));
	DE.insert(map_type::value_type(px2->first,f));
	dmap_type DDE; 
	double B = DK.second; 
	double A = B - C;
	unsigned int i, j, n = dx.size();
	Set::VectorSpace::Hom df(n);
	Set::VectorSpace::Vector a = dx/r;
	Set::VectorSpace::Vector b = A*a;
	for (j=0; j<n; j++)
	{
		for (i=0; i<n; i++) df[j][i] = b[j]*a[i];
		df[j][j] += C;
	}
	Set::VectorSpace::Hom dg = -df;
	Set::Manifold::Point *e1 = px1->first;
	Set::Manifold::Point *e2 = px2->first;
	DDE.insert(dmap_type::value_type(make_pair(e1,e1),df));
	DDE.insert(dmap_type::value_type(make_pair(e1,e2),dg));
	DDE.insert(dmap_type::value_type(make_pair(e2,e1),dg));
	DDE.insert(dmap_type::value_type(make_pair(e2,e2),df));
	return make_pair(DE,DDE);
}

}

}
