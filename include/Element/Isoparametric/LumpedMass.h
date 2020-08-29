// LumpedMass.h: interface for the LumpedMass class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(ELEMENT_ISOPARAMETRIC_LUMPEDMASS_H__INCLUDED_)
#define ELEMENT_ISOPARAMETRIC_LUMPEDMASS_H__INCLUDED_

#pragma once

#include "../Interpolation/Interpolation.h"
#include "../Quadrature/Quadrature.h"
#include "../Conforming/LumpedMass.h"

using namespace std;

namespace Element
{
namespace Isoparametric
{
class LumpedMass : public Element::Conforming::LumpedMass
{
public: 

	LumpedMass();
	virtual ~LumpedMass();
	LumpedMass(
		const double &,
		Element::Interpolation::Shape<0> &,
		Element::Interpolation::Shape<1> &,
		const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &,
		const vector<Set::VectorSpace::Vector> &,
		const vector<double> &);
	LumpedMass(
		const double &,
		Element::Interpolation::Shape<0> &,
		Element::Interpolation::Shape<1> &,
		const Element::Quadrature::Rule &,
		const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &);
	LumpedMass(const LumpedMass &);

private:

	void Initialize(
		Element::Interpolation::Shape<0> &,
		Element::Interpolation::Shape<1> &,
		const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &,
		const vector<Set::VectorSpace::Vector> &,
		const vector<double> &);

private:

	unsigned int Factorial(unsigned int);
};

}

}

#endif // !defined(ELEMENT_ISOPARAMETRIC_LUMPEDMASS_H__INCLUDED_
