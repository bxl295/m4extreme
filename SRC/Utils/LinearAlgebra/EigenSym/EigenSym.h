// EigenSym.h: interface for the EigenSym class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(LINEARALGEBRA_EIGENSYM_H__INCLUDED_)
#define LINEARALGEBRA_EIGENSYM_H__INCLUDED_

#pragma once

#include <math.h>
#include "../../../Set/Algebraic/VectorSpace/Vector/Vector.h"
#include "../../../Set/Algebraic/VectorSpace/Category/Category.h"
#include "../../../Set/Algebraic/VectorSpace/Category/Symmetric/Symmetric.h"
#include "../../../Set/Algebraic/VectorSpace/Category/Diagonal/Diagonal.h"
#include "tnt_array1d.h"
#include "tnt_array2d.h"
#include "jama_eig.h"

namespace LinearAlgebra  
{
class EigenSym
{
public:

	EigenSym();
	virtual ~EigenSym();
	void operator () (
		const Set::VectorSpace::Sym &,
		Set::VectorSpace::Diagonal &,
		Set::VectorSpace::Hom &);
};

}

#endif // !defined(LINEARALGEBRA_EIGENSYM_H__INCLUDED_)
