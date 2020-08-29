// TensorProduct.h: interface for the TensorProduct class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(ELEMENT_INTERPOLATION_UTILS_TENSORPRODUCT_H__INCLUDED_)
#define ELEMENT_INTERPOLATION_UTILS_TENSORPRODUCT_H__INCLUDED_

#pragma once

#include <map>
#include <vector>
#include <cassert>
#include "../MultiIndex/MultiIndex.h"
#include "../Polynomial/Polynomial.h"
#include "../../../../Set/Algebraic/VectorSpace/Vector/Vector.h"

using namespace std;

namespace Element
{
namespace Interpolation
{
namespace Utils
{
//////////////////////////////////////////////////////////////////////
// Class TensorProduct
//////////////////////////////////////////////////////////////////////

class TensorProduct : public Polynomial
{
public: 

	TensorProduct();
	virtual ~TensorProduct();
	TensorProduct(const Polynomial &);
	TensorProduct(const Polynomial &, const Polynomial &);
	TensorProduct(const vector<Polynomial> &);
	TensorProduct & operator = (const Polynomial &);
	double operator () (const Set::VectorSpace::Vector &) const;
	double operator () (const vector<Set::VectorSpace::Vector> &) const;
};

}

}

}

#endif // !defined(ELEMENT_INTERPOLATION_UTILS_TENSORPRODUCT_H__INCLUDED_
