// Source.cpp: implementation of the Isoparametric source class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <cassert>
#include "./Source.h"
#include "../../Utils/Indexing/Indexing.h"

namespace Potential
{
namespace Isoparametric
{
//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

LocalState::LocalState() {}

LocalState::~LocalState() {}

LocalState::LocalState(const vector<Material::LocalState *> &rhs_MatLS) 
	: Potential::Conforming::LocalState(rhs_MatLS) {}

LocalState::LocalState(
	Element::Interpolation::Shape<0> &S,
	Element::Interpolation::Shape<1> &DS,
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &X,
	const vector<Material::LocalState *> &rhs_MatLS,
	const vector<Set::VectorSpace::Vector> &QP,
	const vector<double> &QW0) 
	: Potential::Conforming::LocalState(rhs_MatLS)
{
	unsigned int q, m, n;
	assert(QP.size() == QW0.size());
	assert(MatLS.size() == QW0.size());
	vector<map<Set::Manifold::Point *, Set::VectorSpace::Vector> > DN;
	for (q=0; q<QW0.size(); q++)
	{
		N.push_back(S(QP[q]));
		DN.push_back(DS(QP[q]));
		QW.push_back(QW0[q]);
	}

	m = X.begin()->second.size();
	n = DN[0].begin()->second.size();
	double fact = Factorial(n);

	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator pX;
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pDN;
	for (q=0; q<QW.size(); q++)
	{
		Set::VectorSpace::Hom J(m,n);
		for (pDN=DN[q].begin(),pX=X.begin(); pDN!=DN[q].end(); pDN++,pX++)
			J += Dyadic(pDN->second,pX->second); 
		Set::VectorSpace::Hom G = Adjoint(J)*J;
		QW[q] *= sqrt(Jacobian(G))/fact;
	}
}

LocalState::LocalState(const LocalState &rhs)
	: Potential::Conforming::LocalState(rhs) {}

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
	: Potential::Conforming::Energy<0>(rhs_LS, rhs_W) {}

Energy<0>::Energy(const Energy<0> &rhs)
	: Potential::Conforming::Energy<0>(rhs) {}

//////////////////////////////////////////////////////////////////////
// Class Energy<1>
//////////////////////////////////////////////////////////////////////

Energy<1>::Energy() {}

Energy<1>::~Energy() {}

Energy<1>::Energy(
	LocalState *rhs_LS, 
	const vector<Material::Energy<1> *> &rhs_DW)
	: Potential::Conforming::Energy<1>(rhs_LS, rhs_DW) {}

Energy<1>::Energy(const Energy<1> &rhs)
	: Potential::Conforming::Energy<1>(rhs) {}

//////////////////////////////////////////////////////////////////////
// Class Energy<2>
//////////////////////////////////////////////////////////////////////

Energy<2>::Energy() {}

Energy<2>::~Energy() {}

Energy<2>::Energy(
	LocalState *rhs_LS, 
	const vector<Material::Energy<2> *> &rhs_DDW)
	: Potential::Conforming::Energy<2>(rhs_LS, rhs_DDW) {}

Energy<2>::Energy(const Energy<2> &rhs)
	: Potential::Conforming::Energy<2>(rhs) {}

//////////////////////////////////////////////////////////////////////
// Class Jet<0>
//////////////////////////////////////////////////////////////////////

Jet<0>::Jet() {}

Jet<0>::~Jet() {}

Jet<0>::Jet(
	LocalState *rhs_LS, 
	const vector<Material::Jet<0> *> &rhs_J)
	: Potential::Conforming::Jet<0>(rhs_LS,rhs_J) {}

Jet<0>::Jet(const Jet<0> &rhs)
	: Potential::Conforming::Jet<0>(rhs) {}

//////////////////////////////////////////////////////////////////////
// Class Jet<1>
//////////////////////////////////////////////////////////////////////

Jet<1>::Jet() {}

Jet<1>::~Jet() {}

Jet<1>::Jet(
	LocalState *rhs_LS, 
	const vector<Material::Jet<1> *> &rhs_DJ)
	: Potential::Conforming::Jet<1>(rhs_LS,rhs_DJ) {}

Jet<1>::Jet(const Jet<1> &rhs)
	: Potential::Conforming::Jet<1>(rhs) {}

}

}
