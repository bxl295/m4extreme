// Isoparametric.cpp: implementation of the Isoparametric class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <cassert>
#include "./Isoparametric.h"
#include "../../Utils/Indexing/Indexing.h"

namespace Element
{
namespace Isoparametric
{
//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

LocalState::LocalState() {}

LocalState::~LocalState() {}

LocalState::LocalState(
	Element::Interpolation::Shape<1> &DS,
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &X,
	const vector<Material::LocalState *> &rhs_MatLS,
	const vector<Set::VectorSpace::Vector> &QP,
	const vector<double> &QW0) 
	: Element::Conforming::LocalState(rhs_MatLS)
{
	assert(QP.size() == QW0.size());
	assert(MatLS.size() == QW0.size());

	unsigned int q;
	unsigned int n = X.begin()->second.size();

	vector<map<Set::Manifold::Point *, Set::VectorSpace::Vector> > DN0(QP.size());
	for (q=0; q<QW0.size(); q++)
	{
		DN0[q] = DS(QP[q]);
		DN.push_back(DN0[q]);
		QW.push_back(QW0[q]);
	}

	double fact = Factorial(n);
	for (q=0; q<QW.size(); q++)
	{
		Set::VectorSpace::Hom J(n);
		map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator pX;
		map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pDN0;
		for (pX=X.begin(), pDN0=DN0[q].begin(); pX!=X.end(); pX++, pDN0++)
			J += Dyadic(pDN0->second,pX->second);
		Set::VectorSpace::Hom Jinv = Inverse(J);
		map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pDN;
		for (pDN0=DN0[q].begin(), pDN=DN[q].begin(); pDN!=DN[q].end(); pDN0++, pDN++)
				pDN->second = Adjoint(Jinv)*pDN0->second;
		QW[q] *= fabs(Jacobian(J))/fact;
	}
}

LocalState::LocalState(const LocalState &rhs)
	: Element::Conforming::LocalState(rhs) {}

LocalState & 
LocalState::operator = (const LocalState &rhs)
{
	if (this == &rhs) return *this; 
	Element::Conforming::LocalState::operator = (rhs);
	return *this;
}

unsigned int 
LocalState::Factorial(unsigned int n)
{
	unsigned int i, fact; fact = 1;
	if (n > 1) for (i=2; i<=n; i++) fact = i*fact;
	return(fact);
}

//////////////////////////////////////////////////////////////////////
// Class Energy<0>
//////////////////////////////////////////////////////////////////////

Energy<0>::Energy() {}

Energy<0>::~Energy() {}

Energy<0>::Energy(
	LocalState *rhs_LS, 
	const vector<Material::Energy<0> *> &rhs_W)
	: Element::Conforming::Energy<0>(rhs_LS, rhs_W) {}

Energy<0>::Energy(const Energy<0> &rhs)
	: Element::Conforming::Energy<0>(rhs) {}

//////////////////////////////////////////////////////////////////////
// Class Energy<1>
//////////////////////////////////////////////////////////////////////

Energy<1>::Energy() {}

Energy<1>::~Energy() {}

Energy<1>::Energy(
	LocalState *rhs_LS, 
	const vector<Material::Energy<1> *> &rhs_DW)
	: Element::Conforming::Energy<1>(rhs_LS, rhs_DW) {}

Energy<1>::Energy(const Energy<1> &rhs)
	: Element::Conforming::Energy<1>(rhs) {}

//////////////////////////////////////////////////////////////////////
// Class Energy<2>
//////////////////////////////////////////////////////////////////////

Energy<2>::Energy() {}

Energy<2>::~Energy() {}

Energy<2>::Energy(
	LocalState *rhs_LS, 
	const vector<Material::Energy<2> *> &rhs_DDW)
	: Element::Conforming::Energy<2>(rhs_LS, rhs_DDW) {}

Energy<2>::Energy(const Energy<2> &rhs)
	: Element::Conforming::Energy<2>(rhs) {}

//////////////////////////////////////////////////////////////////////
// Class Jet<0>
//////////////////////////////////////////////////////////////////////

Jet<0>::Jet() {}

Jet<0>::~Jet() {}

Jet<0>::Jet(
	LocalState *rhs_LS, 
	const vector<Material::Jet<0> *> &rhs_J)
	: Element::Conforming::Jet<0>(rhs_LS,rhs_J) {}

Jet<0>::Jet(const Jet<0> &rhs)
	: Element::Conforming::Jet<0>(rhs) {}

//////////////////////////////////////////////////////////////////////
// Class Jet<1>
//////////////////////////////////////////////////////////////////////

Jet<1>::Jet() {}

Jet<1>::~Jet() {}

Jet<1>::Jet(
	LocalState *rhs_LS, 
	const vector<Material::Jet<1> *> &rhs_DJ)
	: Element::Conforming::Jet<1>(rhs_LS,rhs_DJ) {}

Jet<1>::Jet(const Jet<1> &rhs)
	: Element::Conforming::Jet<1>(rhs) {}

}

}
