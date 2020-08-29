// EmbeddedAtom.cpp: implementation for the EmbeddedAtom class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include "./EmbeddedAtom.h"

namespace Potential
{
namespace EmbeddedAtom
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

LocalState::LocalState(
	Set::Manifold::Point *rhs_e, const set<Set::Manifold::Point *> &rhs_N) 
	: e(rhs_e), N(rhs_N) {}

LocalState::LocalState(const LocalState &rhs)
	: e(rhs.e), N(rhs.N) {}

LocalState & 
LocalState::operator = (const LocalState &rhs)
{
	if (this == &rhs) return *this; 
	e = rhs.e; N = rhs.N;
	return *this;
}

void 
LocalState::operator ++ () {}

set<Set::Manifold::Point *> 
LocalState::GetNodes() const
{
	set<Set::Manifold::Point *> Nodes = N;
	Nodes.insert(e); 
	return Nodes;
}

set<pair<Set::Manifold::Point *, Set::Manifold::Point *> > 
LocalState::GetNodePairs() const
{
	set<pair<Set::Manifold::Point *, Set::Manifold::Point *> > NodePairs;
	NodePairs.insert(make_pair(e,e));
	set<Set::Manifold::Point *>::const_iterator pN, pN2;
	for (pN=N.begin(); pN!=N.end(); pN++)
	{
		Set::Manifold::Point *d = *pN;
		NodePairs.insert(make_pair(e,d));
		NodePairs.insert(make_pair(d,e));
                for (pN2=N.begin(); pN2!=N.end(); pN2++)
                {
                    Set::Manifold::Point *h = *pN2;
                    NodePairs.insert(make_pair(d,h));
                }
	}
	return NodePairs;
}

Set::VectorSpace::Hom 
LocalState::Dyadic(
	const Set::VectorSpace::Vector &a,
	const Set::VectorSpace::Vector &b)
{
	Set::VectorSpace::Hom ab(a.size(),b.size());
	for (unsigned int j=0; j<b.size(); j++) ab[j] = b[j]*a;
	return ab;
}

//////////////////////////////////////////////////////////////////////
// Class Energy<0>
//////////////////////////////////////////////////////////////////////

Energy<0>::Energy() {}

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

Energy<0>::~Energy() 
{
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<0> *>::iterator pR;
	for (pR=R.begin(); pR!=R.end(); pR++) delete pR->second;
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<0> *>::iterator pV;
	for (pV=V.begin(); pV!=V.end(); pV++) delete pV->second;
	map<Set::Manifold::Point *, Potential::TwoBody::LocalState *>::iterator pLS;
	for (pLS=TBLS.begin(); pLS!=TBLS.end(); pLS++) delete pLS->second;
}

Energy<0>::Energy(
	LocalState *rhs_LS,
	Potential::Radial::Energy<0> *rhs_U,
	Potential::Radial::Energy<0> *R0,
	Potential::Radial::Energy<0> *V0) 
	: LS(rhs_LS), U(rhs_U) 
{
	typedef map<Set::Manifold::Point *, Potential::TwoBody::Energy<0> *> r_type;
	typedef map<Set::Manifold::Point *, Potential::TwoBody::Energy<0> *> v_type;
	typedef map<Set::Manifold::Point *, Potential::TwoBody::LocalState *> ls_type;

	Set::Manifold::Point *e = LS->e;
	set<Set::Manifold::Point *>::iterator pN;
	for (pN=LS->N.begin(); pN!=LS->N.end(); pN++)
	{
		Set::Manifold::Point *d = *pN;
		Potential::TwoBody::LocalState *ed = 
			new Potential::TwoBody::LocalState(e,d);
		TBLS.insert(ls_type::value_type(d,ed));
		R.insert(r_type::value_type(
			d,new Potential::TwoBody::Energy<0>(ed,R0)));
		V.insert(v_type::value_type(
			d,new Potential::TwoBody::Energy<0>(ed,V0)));
	}
}

Energy<0>::Energy(const Energy<0> &rhs) 
	: LS(rhs.LS), U(rhs.U), 
	R(rhs.R), V(rhs.V), TBLS(rhs.TBLS) {}

const Energy<0> & 
Energy<0>::operator = (const Energy<0> &rhs)
{
	if (this == &rhs) return *this; 
	LS = rhs.LS; U = rhs.U; 
	R = rhs.R; V = rhs.V; TBLS = rhs.TBLS;
	return *this;
}

double 
Energy<0>::operator () (const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x0) const
{
	Set::Manifold::Point *e = LS->e; 
	double W=0.0; double rho=0.0; 
	map<Set::Manifold::Point *, Set::VectorSpace::Vector> x=x0;
	Set::VectorSpace::Vector xe = x[e];
	set<Set::Manifold::Point *>::iterator pN;
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<0> *>::const_iterator pR;
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<0> *>::const_iterator pV;
	for (pN=LS->N.begin(), pR=R.begin(), pV=V.begin(); 
		pN!=LS->N.end(); pN++, pR++, pV++)
	{
		Set::Manifold::Point *d = *pN;
		map<Set::Manifold::Point *, Set::VectorSpace::Vector> y;
		y.insert(domain_type::value_type(e,xe));
		y.insert(domain_type::value_type(d,x[d]));
		W += (*pV->second)(y); rho += (*pR->second)(y);
	}
	return W + (*U)(rho);
}

//////////////////////////////////////////////////////////////////////
// Class Energy<1>
//////////////////////////////////////////////////////////////////////

Energy<1>::Energy() {}

Energy<1>::~Energy() 
{
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<0> *>::iterator pR;
	for (pR=R.begin(); pR!=R.end(); pR++) delete pR->second;
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<1> *>::iterator pDR;
	for (pDR=DR.begin(); pDR!=DR.end(); pDR++) delete pDR->second;
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<1> *>::iterator pDV;
	for (pDV=DV.begin(); pDV!=DV.end(); pDV++) delete pDV->second;
	map<Set::Manifold::Point *, Potential::TwoBody::LocalState *>::iterator pLS;
	for (pLS=TBLS.begin(); pLS!=TBLS.end(); pLS++) delete pLS->second;
}

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
	Potential::Radial::Energy<1> *rhs_DU,
	Potential::Radial::Energy<0> *R0,
	Potential::Radial::Energy<1> *DR0,
	Potential::Radial::Energy<1> *DV0) 
	: LS(rhs_LS), DU(rhs_DU) 
{
	typedef map<Set::Manifold::Point *, Potential::TwoBody::Energy<0> *> r_type;
	typedef map<Set::Manifold::Point *, Potential::TwoBody::Energy<1> *> dr_type;
	typedef map<Set::Manifold::Point *, Potential::TwoBody::Energy<1> *> dv_type;
	typedef map<Set::Manifold::Point *, Potential::TwoBody::LocalState *> ls_type;

	Set::Manifold::Point *e = LS->e;
	set<Set::Manifold::Point *>::iterator pN;
	for (pN=LS->N.begin(); pN!=LS->N.end(); pN++)
	{
		Set::Manifold::Point *d = *pN;
		Potential::TwoBody::LocalState *ed = 
			new Potential::TwoBody::LocalState(e,d);
		TBLS.insert(ls_type::value_type(d,ed));
		R.insert(r_type::value_type(
			d,new Potential::TwoBody::Energy<0>(ed,R0)));
		DR.insert(dr_type::value_type(
			d,new Potential::TwoBody::Energy<1>(ed,DR0)));
		DV.insert(dv_type::value_type(
			d,new Potential::TwoBody::Energy<1>(ed,DV0)));
	}
}

Energy<1>::Energy(const Energy<1> &rhs) 
	: LS(rhs.LS), DU(rhs.DU), 
	R(rhs.R), DR(rhs.DR), 
	DV(rhs.DV), TBLS(rhs.TBLS) {}

const Energy<1> & 
Energy<1>::operator = (const Energy<1> &rhs)
{
	if (this == &rhs) return *this; 
	LS = rhs.LS; DU = rhs.DU; 
	DV = rhs.DV; R = rhs.R; 
	DR = rhs.DR; TBLS = rhs.TBLS;
	return *this;
}

map<Set::Manifold::Point *, Set::VectorSpace::Vector>
Energy<1>::operator () (const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x0) const
{
	map<Set::Manifold::Point *, Set::VectorSpace::Vector> x=x0;
	Set::Manifold::Point *e = LS->e; range_type DW;
	Set::VectorSpace::Vector xe = x[e];
	double rho=0.0; set<Set::Manifold::Point *>::iterator pN;
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<0> *>::const_iterator pR;
	for (pN=LS->N.begin(), pR=R.begin(); pN!=LS->N.end(); pN++, pR++)
	{
		Set::Manifold::Point *d = *pN;
		map<Set::Manifold::Point *, Set::VectorSpace::Vector> y;
		y.insert(domain_type::value_type(e,xe));
		y.insert(domain_type::value_type(d,x[d]));
		rho += (*pR->second)(y);
	}
	double du = (*DU)(rho);
	Set::VectorSpace::Vector O(xe.size());
	DW.insert(range_type::value_type(e,O));
	for (pN=LS->N.begin(); pN!=LS->N.end(); pN++)
		DW.insert(range_type::value_type(*pN,O));
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<1> *>::const_iterator pDR;
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<1> *>::const_iterator pDV;
	for (pN=LS->N.begin(), pDR=DR.begin(), pDV=DV.begin(); 
		pN!=LS->N.end(); pN++, pR++, pDR++, pDV++)
	{
		Set::Manifold::Point *d = *pN;
		map<Set::Manifold::Point *, Set::VectorSpace::Vector> y;
		y.insert(domain_type::value_type(e,xe));
		y.insert(domain_type::value_type(d,x[d]));
		map<Set::Manifold::Point *, Set::VectorSpace::Vector> f = (*pDV->second)(y);
		map<Set::Manifold::Point *, Set::VectorSpace::Vector> g = (*pDR->second)(y);
		DW[e] += f[e] + du*g[e];
		DW[d] += f[d] + du*g[d];
	}
	return DW;
}

//////////////////////////////////////////////////////////////////////
// Class Energy<2>
//////////////////////////////////////////////////////////////////////

Energy<2>::Energy() {}

Energy<2>::~Energy() 
{
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<0> *>::iterator pR;
	for (pR=R.begin(); pR!=R.end(); pR++) delete pR->second;
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<1> *>::iterator pDR;
	for (pDR=DR.begin(); pDR!=DR.end(); pDR++) delete pDR->second;
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<2> *>::iterator pDDR;
	for (pDDR=DDR.begin(); pDDR!=DDR.end(); pDDR++) delete pDDR->second;
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<2> *>::iterator pDDV;
	for (pDDV=DDV.begin(); pDDV!=DDV.end(); pDDV++) delete pDDV->second;
	map<Set::Manifold::Point *, Potential::TwoBody::LocalState *>::iterator pLS;
	for (pLS=TBLS.begin(); pLS!=TBLS.end(); pLS++) delete pLS->second;
}

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
	Potential::Radial::Energy<1> *rhs_DU,
	Potential::Radial::Energy<2> *rhs_DDU,
	Potential::Radial::Energy<0> *R0,
	Potential::Radial::Energy<1> *DR0,
	Potential::Radial::Energy<2> *DDR0,
	Potential::Radial::Energy<1> *DV0,
	Potential::Radial::Energy<2> *DDV0) 
	: LS(rhs_LS), DU(rhs_DU), DDU(rhs_DDU) 
{
	typedef map<Set::Manifold::Point *, Potential::TwoBody::Energy<0> *> r_type;
	typedef map<Set::Manifold::Point *, Potential::TwoBody::Energy<1> *> dr_type;
	typedef map<Set::Manifold::Point *, Potential::TwoBody::Energy<2> *> ddr_type;
	typedef map<Set::Manifold::Point *, Potential::TwoBody::Energy<1> *> dv_type;
	typedef map<Set::Manifold::Point *, Potential::TwoBody::Energy<2> *> ddv_type;
	typedef map<Set::Manifold::Point *, Potential::TwoBody::LocalState *> ls_type;

	Set::Manifold::Point *e = LS->e;
	set<Set::Manifold::Point *>::iterator pN;
	for (pN=LS->N.begin(); pN!=LS->N.end(); pN++)
	{
		Set::Manifold::Point *d = *pN;
		Potential::TwoBody::LocalState *ed = 
			new Potential::TwoBody::LocalState(e,d);
		TBLS.insert(ls_type::value_type(d,ed));
		R.insert(r_type::value_type(
			d,new Potential::TwoBody::Energy<0>(ed,R0)));
		DR.insert(dr_type::value_type(
			d,new Potential::TwoBody::Energy<1>(ed,DR0)));
		DDR.insert(ddr_type::value_type(
			d,new Potential::TwoBody::Energy<2>(ed,DR0,DDR0)));
		DDV.insert(ddv_type::value_type(
			d,new Potential::TwoBody::Energy<2>(ed,DV0,DDV0)));
	}
}

Energy<2>::Energy(const Energy<2> &rhs) 
	: LS(rhs.LS), DU(rhs.DU), DDU(rhs.DDU), 
	R(rhs.R), DR(rhs.DR), DDR(rhs.DDR), 
	DDV(rhs.DDV), TBLS(rhs.TBLS) {}

const Energy<2> & 
Energy<2>::operator = (const Energy<2> &rhs)
{
	if (this == &rhs) return *this; 
	LS = rhs.LS; DU = rhs.DU; DDU = rhs.DDU; 
	R = rhs.R; DR = rhs.DR; DDR = rhs.DDR; 
	DDV = rhs.DDV; TBLS = rhs.TBLS; 
	return *this;
}

map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, Set::VectorSpace::Hom>
Energy<2>::operator () (const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x0) const
{
	Set::Manifold::Point *e = LS->e; range_type DDW;
	map<Set::Manifold::Point *, Set::VectorSpace::Vector> x=x0;
	Set::VectorSpace::Vector xe = x[e];
	double rho=0.0;
	set<Set::Manifold::Point *>::iterator pN, pN2;
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<0> *>::const_iterator pR;
	for (pN=LS->N.begin(), pR=R.begin(); pN!=LS->N.end(); pN++, pR++)
	{
		Set::Manifold::Point *d = *pN;
		map<Set::Manifold::Point *, Set::VectorSpace::Vector> y;
		y.insert(domain_type::value_type(e,xe));
		y.insert(domain_type::value_type(d,x[d]));
		rho += (*pR->second)(y);
	}
	double du = (*DU)(rho); 
	double ddu = (*DDU)(rho);
	Set::VectorSpace::Hom O(xe.size());
	DDW.insert(range_type::value_type(make_pair(e,e),O));
	for (pN=LS->N.begin(); pN!=LS->N.end(); pN++)
	{
		Set::Manifold::Point *d = *pN;
		DDW.insert(range_type::value_type(make_pair(e,d),O));
		DDW.insert(range_type::value_type(make_pair(d,e),O));
		for (pN2=LS->N.begin(); pN2!=LS->N.end(); pN2++)
                {
                    Set::Manifold::Point *h = *pN2;
                    DDW.insert(range_type::value_type(make_pair(d,h),O));
                }
	}
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<1> *>::const_iterator pDR, pDR2;
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<2> *>::const_iterator pDDR;
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<2> *>::const_iterator pDDV;
	for (pN=LS->N.begin(), pDR=DR.begin(), 
		pDDR=DDR.begin(), pDDV=DDV.begin(); pN!=LS->N.end(); 
		pN++, pDR++, pDDR++, pDDV++)
	{
		Set::Manifold::Point *d = *pN;
		map<Set::Manifold::Point *, Set::VectorSpace::Vector> y;
		y.insert(domain_type::value_type(e,xe));
		y.insert(domain_type::value_type(d,x[d]));
		Potential::TwoBody::Energy<2>::range_type Df = (*pDDV->second)(y);
		Potential::TwoBody::Energy<1>::range_type g = (*pDR->second)(y);
		Potential::TwoBody::Energy<2>::range_type Dg = (*pDDR->second)(y);
		DDW[make_pair(e,e)] += 
			Df[make_pair(e,e)] + 
			du*Dg[make_pair(e,e)];
		DDW[make_pair(e,d)] += 
			Df[make_pair(e,d)] + 
			du*Dg[make_pair(e,d)];
		DDW[make_pair(d,e)] += 
			Df[make_pair(d,e)] +
			du*Dg[make_pair(d,e)];
		DDW[make_pair(d,d)] += 
			Df[make_pair(d,d)] + 
			du*Dg[make_pair(d,d)];
                
                for (pN2=LS->N.begin(), pDR2=DR.begin(); pN2!=LS->N.end(); 
                     pN2++, pDR2++)
                {
                    Set::Manifold::Point *h = *pN2;
                    map<Set::Manifold::Point *, Set::VectorSpace::Vector> y2;
                    y2.insert(domain_type::value_type(e,xe));
                    y2.insert(domain_type::value_type(h,x[h]));
                    Potential::TwoBody::Energy<1>::range_type g2 = (*pDR2->second)(y2);
                    DDW[make_pair(e,e)] += ddu*LS->Dyadic(g[e],g2[e]);
                    DDW[make_pair(e,h)] += ddu*LS->Dyadic(g[e],g2[h]);
                    DDW[make_pair(d,e)] += ddu*LS->Dyadic(g[d],g2[e]);
                    DDW[make_pair(d,h)] += ddu*LS->Dyadic(g[d],g2[h]);
                }
	}
	return DDW;
}

//////////////////////////////////////////////////////////////////////
// Class Jet<0>
//////////////////////////////////////////////////////////////////////

Jet<0>::Jet() {}

Jet<0>::~Jet() 
{
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<0> *>::iterator pR;
	for (pR=R.begin(); pR!=R.end(); pR++) delete pR->second;
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<1> *>::iterator pDR;
	for (pDR=DR.begin(); pDR!=DR.end(); pDR++) delete pDR->second;
	map<Set::Manifold::Point *, Potential::TwoBody::Jet<0> *>::iterator pVJ;
	for (pVJ=VJ.begin(); pVJ!=VJ.end(); pVJ++) delete pVJ->second;
	map<Set::Manifold::Point *, Potential::TwoBody::LocalState *>::iterator pLS;
	for (pLS=TBLS.begin(); pLS!=TBLS.end(); pLS++) delete pLS->second;
} 

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
	Potential::Radial::Jet<0> *rhs_UJ,
	Potential::Radial::Energy<0> *R0,
	Potential::Radial::Energy<1> *DR0,
	Potential::Radial::Jet<0> *VJ0)
	: LS(rhs_LS), UJ(rhs_UJ) 
{
	typedef map<Set::Manifold::Point *, Potential::TwoBody::Energy<0> *> r_type;
	typedef map<Set::Manifold::Point *, Potential::TwoBody::Energy<1> *> dr_type;
	typedef map<Set::Manifold::Point *, Potential::TwoBody::Jet<0> *> v_type;
	typedef map<Set::Manifold::Point *, Potential::TwoBody::LocalState *> ls_type;

	Set::Manifold::Point *e = LS->e;
	set<Set::Manifold::Point *>::iterator pN;
	for (pN=LS->N.begin(); pN!=LS->N.end(); pN++)
	{
		Set::Manifold::Point *d = *pN;
		Potential::TwoBody::LocalState *ed = 
			new Potential::TwoBody::LocalState(e,d);
		TBLS.insert(ls_type::value_type(d,ed));
		R.insert(r_type::value_type(
			d,new Potential::TwoBody::Energy<0>(ed,R0)));
		DR.insert(dr_type::value_type(
			d,new Potential::TwoBody::Energy<1>(ed,DR0)));
		VJ.insert(v_type::value_type(
			d,new Potential::TwoBody::Jet<0>(ed,VJ0)));
	}
}

Jet<0>::Jet(const Jet<0> &rhs)
	: LS(rhs.LS), UJ(rhs.UJ), 
	R(rhs.R), DR(rhs.DR), 
	VJ(rhs.VJ), TBLS(rhs.TBLS) {}

Jet<0> & 
Jet<0>::operator = (const Jet<0> &rhs)
{
	if (this == &rhs) return *this;
	LS = rhs.LS; UJ = rhs.UJ;
	R = rhs.R; DR = rhs.DR;
	VJ = rhs.VJ; TBLS = rhs.TBLS;
	return *this;
}

pair<double, map<Set::Manifold::Point *, Set::VectorSpace::Vector> > 
Jet<0>::operator () (const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x0) const
{
	map<Set::Manifold::Point *, Set::VectorSpace::Vector> x=x0;
	Set::Manifold::Point *e = LS->e; Set::VectorSpace::Vector xe = x[e];
	double rho=0.0; set<Set::Manifold::Point *>::iterator pN;
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<0> *>::const_iterator pR;
	for (pN=LS->N.begin(), pR=R.begin(); pN!=LS->N.end(); pN++, pR++)
	{
		Set::Manifold::Point *d = *pN;
		map<Set::Manifold::Point *, Set::VectorSpace::Vector> y;
		y.insert(domain_type::value_type(e,xe));
		y.insert(domain_type::value_type(d,x[d]));
		rho += (*pR->second)(y);
	}
	pair<double,double> KU = (*UJ)(rho);
	double W = KU.first;
	double du = KU.second;
	Potential::TwoBody::Energy<1>::range_type DW;
	Set::VectorSpace::Vector O(xe.size());
	DW.insert(Potential::TwoBody::Energy<1>::range_type::value_type(e,O));
	for (pN=LS->N.begin(); pN!=LS->N.end(); pN++)
		DW.insert(Potential::TwoBody::Energy<1>::range_type::value_type(*pN,O));
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<1> *>::const_iterator pDR;
	map<Set::Manifold::Point *, Potential::TwoBody::Jet<0> *>::const_iterator pVJ;
	for (pN=LS->N.begin(), pDR=DR.begin(), pVJ=VJ.begin(); 
		pN!=LS->N.end(); pN++, pR++, pDR++, pVJ++)
	{
		Set::Manifold::Point *d = *pN;
		map<Set::Manifold::Point *, Set::VectorSpace::Vector> y;
		y.insert(domain_type::value_type(e,xe));
		y.insert(domain_type::value_type(d,x[d]));
		Potential::TwoBody::Jet<0>::range_type KV = (*pVJ->second)(y);
		W += KV.first; 
		Potential::TwoBody::Energy<1>::range_type f = KV.second;
		Potential::TwoBody::Energy<1>::range_type g = (*pDR->second)(y);
		DW[e] += f[e] + du*g[e];
		DW[d] += f[d] + du*g[d];
	}
	return make_pair(W,DW);
}

//////////////////////////////////////////////////////////////////////
// Class Jet<1>
//////////////////////////////////////////////////////////////////////

Jet<1>::Jet() {}

Jet<1>::~Jet() 
{
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<0> *>::iterator pR;
	for (pR=R.begin(); pR!=R.end(); pR++) delete pR->second;
	map<Set::Manifold::Point *, Potential::TwoBody::Jet<1> *>::iterator pDRJ;
	for (pDRJ=DRJ.begin(); pDRJ!=DRJ.end(); pDRJ++) delete pDRJ->second;
	map<Set::Manifold::Point *, Potential::TwoBody::Jet<1> *>::iterator pDVJ;
	for (pDVJ=DVJ.begin(); pDVJ!=DVJ.end(); pDVJ++) delete pDVJ->second;
	map<Set::Manifold::Point *, Potential::TwoBody::LocalState *>::iterator pLS;
	for (pLS=TBLS.begin(); pLS!=TBLS.end(); pLS++) delete pLS->second;
}

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
	Potential::Radial::Energy<0> *R0,
	Potential::Radial::Jet<1> *rhs_DUJ,
	Potential::Radial::Jet<1> *DRJ0,
	Potential::Radial::Jet<1> *DVJ0)
	: LS(rhs_LS), DUJ(rhs_DUJ) 
{
	typedef map<Set::Manifold::Point *, Potential::TwoBody::Energy<0> *> r_type;
	typedef map<Set::Manifold::Point *, Potential::TwoBody::Jet<1> *> dr_type;
	typedef map<Set::Manifold::Point *, Potential::TwoBody::Jet<1> *> dv_type;
	typedef map<Set::Manifold::Point *, Potential::TwoBody::LocalState *> ls_type;

	Set::Manifold::Point *e = LS->e;
	set<Set::Manifold::Point *>::iterator pN;
	for (pN=LS->N.begin(); pN!=LS->N.end(); pN++)
	{
		Set::Manifold::Point *d = *pN;
		Potential::TwoBody::LocalState *ed = 
			new Potential::TwoBody::LocalState(e,d);
		TBLS.insert(ls_type::value_type(d,ed));
		R.insert(r_type::value_type(
			d,new Potential::TwoBody::Energy<0>(ed,R0)));
		DRJ.insert(dr_type::value_type(
			d,new Potential::TwoBody::Jet<1>(ed,DRJ0)));
		DVJ.insert(dv_type::value_type(
			d,new Potential::TwoBody::Jet<1>(ed,DVJ0)));
	}
}

Jet<1>::Jet(const Jet<1> &rhs)
	: LS(rhs.LS), DUJ(rhs.DUJ), 
	R(rhs.R), DRJ(rhs.DRJ), 
	DVJ(rhs.DVJ), TBLS(rhs.TBLS) {}

Jet<1> & 
Jet<1>::operator = (const Jet<1> &rhs)
{
	if (this == &rhs) return *this;
	LS = rhs.LS; DUJ = rhs.DUJ;
	R = rhs.R; DRJ = rhs.DRJ;
	DVJ = rhs.DVJ; TBLS = rhs.TBLS;
	return *this;
}

pair<
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>, 
	map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, Set::VectorSpace::Hom> >
Jet<1>::operator () (const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x0) const
{
	map<Set::Manifold::Point *, Set::VectorSpace::Vector> x=x0;
	Set::Manifold::Point *e = LS->e; Set::VectorSpace::Vector xe = x[e];
	double rho=0.0; set<Set::Manifold::Point *>::iterator pN, pN2;
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<0> *>::const_iterator pR;
	for (pN=LS->N.begin(), pR=R.begin(); pN!=LS->N.end(); pN++, pR++)
	{
		Set::Manifold::Point *d = *pN;
		map<Set::Manifold::Point *, Set::VectorSpace::Vector> y;
		y.insert(domain_type::value_type(e,xe));
		y.insert(domain_type::value_type(d,x[d]));
		rho += (*pR->second)(y);
	}
	pair<double,double> DKU = (*DUJ)(rho);
	double du = DKU.first;
	double ddu = DKU.second;
	Potential::TwoBody::Energy<1>::range_type DW;
	Set::VectorSpace::Vector O(xe.size());
	DW.insert(Potential::TwoBody::Energy<1>::range_type::value_type(e,O));
	for (pN=LS->N.begin(); pN!=LS->N.end(); pN++)
		DW.insert(Potential::TwoBody::Energy<1>::range_type::value_type(*pN,O));
	Set::VectorSpace::Hom P(xe.size());
	Potential::TwoBody::Energy<2>::range_type DDW;
	DDW.insert(Potential::TwoBody::Energy<2>::range_type::value_type(make_pair(e,e),P));
	for (pN=LS->N.begin(); pN!=LS->N.end(); pN++)
	{
		Set::Manifold::Point *d = *pN;
		DDW.insert(Potential::TwoBody::Energy<2>::range_type::value_type(make_pair(e,d),P));
		DDW.insert(Potential::TwoBody::Energy<2>::range_type::value_type(make_pair(d,e),P));
		for (pN2=LS->N.begin(); pN2!=LS->N.end(); pN2++)
                {
                    Set::Manifold::Point *h = *pN2;
                    DDW.insert(Potential::TwoBody::Energy<2>::range_type::value_type(make_pair(d,h),P));
                }    

	}
	map<Set::Manifold::Point *, Potential::TwoBody::Jet<1> *>::const_iterator pDRJ, pDRJ2;
	map<Set::Manifold::Point *, Potential::TwoBody::Jet<1> *>::const_iterator pDVJ;
	for (pN=LS->N.begin(), pDRJ=DRJ.begin(), pDVJ=DVJ.begin(); 
		pN!=LS->N.end(); pN++, pDRJ++, pDVJ++)
	{
		Set::Manifold::Point *d = *pN;
		map<Set::Manifold::Point *, Set::VectorSpace::Vector> y;
		y.insert(domain_type::value_type(e,xe));
		y.insert(domain_type::value_type(d,x[d]));
		Potential::TwoBody::Jet<1>::range_type DKV = (*pDVJ->second)(y);
		Potential::TwoBody::Energy<1>::range_type f = DKV.first;
		Potential::TwoBody::Energy<2>::range_type Df = DKV.second;
		Potential::TwoBody::Jet<1>::range_type DKR = (*pDRJ->second)(y);
		Potential::TwoBody::Energy<1>::range_type g = DKR.first;
		Potential::TwoBody::Energy<2>::range_type Dg = DKR.second;
		DW[e] += f[e] + du*g[e];
		DW[d] += f[d] + du*g[d];
		DDW[make_pair(e,e)] += 
			Df[make_pair(e,e)] + 
			du*Dg[make_pair(e,e)];
		DDW[make_pair(e,d)] += 
			Df[make_pair(e,d)] + 
			du*Dg[make_pair(e,d)];
		DDW[make_pair(d,e)] += 
			Df[make_pair(d,e)] + 
			du*Dg[make_pair(d,e)];
		DDW[make_pair(d,d)] += 
			Df[make_pair(d,d)] + 
			du*Dg[make_pair(d,d)];
                
                for (pN2=LS->N.begin(), pDRJ2=DRJ.begin(); pN2!=LS->N.end(); 
                     pN2++, pDRJ2++)
                {
                    Set::Manifold::Point *h = *pN2;
                    map<Set::Manifold::Point *, Set::VectorSpace::Vector> y2;
                    y2.insert(domain_type::value_type(e,xe));
                    y2.insert(domain_type::value_type(h,x[h]));
                    Potential::TwoBody::Jet<1>::range_type DKR2 = (*pDRJ2->second)(y2);
                    Potential::TwoBody::Energy<1>::range_type g2 = DKR2.first;
                    DDW[make_pair(e,e)] += ddu*LS->Dyadic(g[e],g2[e]);
                    DDW[make_pair(e,h)] += ddu*LS->Dyadic(g[e],g2[h]);
                    DDW[make_pair(d,e)] += ddu*LS->Dyadic(g[d],g2[e]);
                    DDW[make_pair(d,h)] += ddu*LS->Dyadic(g[d],g2[h]);
                }
	}
	return make_pair(DW,DDW);
}

}

}
