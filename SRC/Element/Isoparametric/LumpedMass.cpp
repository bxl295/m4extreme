// LumpedMass.cpp: implementation of the LumpedMass class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <cassert>
#include "./LumpedMass.h"

namespace Element
{
namespace Isoparametric
{
LumpedMass::LumpedMass() {}

LumpedMass::~LumpedMass() {}

LumpedMass::LumpedMass(
	const double &rho_,
	Element::Interpolation::Shape<0> &S,
	Element::Interpolation::Shape<1> &DS,
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &X,
	const vector<Set::VectorSpace::Vector> &QP,
	const vector<double> &QW0) : 
	Element::Conforming::LumpedMass(rho_)
{
	Initialize(S,DS,X,QP,QW0);
}

LumpedMass::LumpedMass(
	const double &rho_,
	Element::Interpolation::Shape<0> &S,
	Element::Interpolation::Shape<1> &DS,
	const Element::Quadrature::Rule &QR,
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &X) : 
	Element::Conforming::LumpedMass(rho_)
{
	vector<Set::VectorSpace::Vector> QP = QR.GetQ();
	vector<double> QW0 = QR.GetW();
	Initialize(S,DS,X,QP,QW0);
}

void
LumpedMass::Initialize(
	Element::Interpolation::Shape<0> &S,
	Element::Interpolation::Shape<1> &DS,
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &X,
	const vector<Set::VectorSpace::Vector> &QP,
	const vector<double> &QW0)
{
	assert(QP.size() == QW0.size());
	unsigned int q, n = X.begin()->second.size();

	vector<map<Set::Manifold::Point *, Set::VectorSpace::Vector> > DN;
	for (q=0; q<QW0.size(); q++)
	{
		N.push_back(S(QP[q]));
		DN.push_back(DS(QP[q]));
		QW.push_back(QW0[q]);
	}

	double fact = Factorial(n);
	for (q=0; q<QW.size(); q++)
	{
		Set::VectorSpace::Hom J(n);
		map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator pX;
		map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pDN;
		for (pX=X.begin(), pDN=DN[q].begin(); pX!=X.end(); pX++, pDN++)
			J += Dyadic(pX->second,pDN->second);
		QW[q] *= fabs(Jacobian(J))/fact;
	}
}

LumpedMass::LumpedMass(const LumpedMass &rhs)
	: Element::Conforming::LumpedMass(rhs) {}

unsigned int 
LumpedMass::Factorial(unsigned int n)
{
	unsigned int i, fact; fact = 1;
	if (n > 1) for (i=2; i<=n; i++) fact = i*fact;
	return(fact);
}

}

}

