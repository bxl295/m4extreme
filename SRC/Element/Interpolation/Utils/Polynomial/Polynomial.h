// Polynomial.h: interface for the Polynomial class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(ELEMENT_INTERPOLATION_UTILS_POLYNOMIAL_H__INCLUDED_)
#define ELEMENT_INTERPOLATION_UTILS_POLYNOMIAL_H__INCLUDED_

#pragma once

#include <set>
#include <map>
#include <cassert>
#include <iostream>
#include "../MultiIndex/MultiIndex.h"
#include "../../../../Set/Algebraic/VectorSpace/Vector/Vector.h"

using namespace std;

namespace Element
{
namespace Interpolation
{
namespace Utils
{
//////////////////////////////////////////////////////////////////////
// Class Polynomial
//////////////////////////////////////////////////////////////////////

class Polynomial
{
public: 

	Polynomial();
	Polynomial(const unsigned int &);
	virtual ~Polynomial();
	Polynomial(const Polynomial &);
	Polynomial & operator = (const Polynomial &);
	const unsigned int & dim() const;
	map<MultiIndex, double> & GetCoefficients();
	const map<MultiIndex, double> & GetCoefficients() const;
	void operator *= (const double &);
	void operator /= (const double &);
	void operator += (const Polynomial &);
	void operator -= (const Polynomial &);
	double operator () (const Set::VectorSpace::Vector &) const;
	Polynomial diff(const unsigned int &) const;

protected:

	unsigned int n;
	map<MultiIndex, double> c;
};

//////////////////////////////////////////////////////////////////////
// Utilities
//////////////////////////////////////////////////////////////////////

set<MultiIndex>
IndexSet(const unsigned int &, const unsigned int &);

set<MultiIndex>
IndexSetInterior(const unsigned int &, const unsigned int &);

unsigned int
IndexSetInteriorSize(const unsigned int &, const unsigned int &);

}

}

}

//////////////////////////////////////////////////////////////////////
// Class Polynomial: Operators
//////////////////////////////////////////////////////////////////////

Element::Interpolation::Utils::Polynomial operator + (
	const Element::Interpolation::Utils::Polynomial &, 
	const Element::Interpolation::Utils::Polynomial &);

Element::Interpolation::Utils::Polynomial operator - (
	const Element::Interpolation::Utils::Polynomial &, 
	const Element::Interpolation::Utils::Polynomial &);

Element::Interpolation::Utils::Polynomial operator * (
	const Element::Interpolation::Utils::Polynomial &, 
	const double &);

Element::Interpolation::Utils::Polynomial operator * (
	const double &, 
	const Element::Interpolation::Utils::Polynomial &);

Element::Interpolation::Utils::Polynomial operator / (
	const Element::Interpolation::Utils::Polynomial &, 
	const double &);

#endif // !defined(ELEMENT_INTERPOLATION_UTILS_POLYNOMIAL_H__INCLUDED_
