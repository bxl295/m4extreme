// Quadrature.h: interface for the Quadrature class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(ELEMENT_QUADRATURE_H__INCLUDED_)
#define ELEMENT_QUADRATURE_H__INCLUDED_

#pragma once

#include <vector>
#include "../../Set/Algebraic/AlgLib.h"

using namespace std;

namespace Element
{
namespace Quadrature
{
//////////////////////////////////////////////////////////////////////
// Class Rule
//////////////////////////////////////////////////////////////////////

class Rule
{
public:

	Rule();
	virtual ~Rule();
	virtual vector<Set::VectorSpace::Vector> GetQ() const=0;
	virtual vector<double> GetW() const=0;
};

}

}

#endif // !defined(ELEMENT_QUADRATURE_H__INCLUDED_)
