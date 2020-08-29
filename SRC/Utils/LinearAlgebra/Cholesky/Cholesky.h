// Cholesky.h: interface for the Cholesky class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(LINEARALGEBRA_CHOLESKY_H__INCLUDED_)
#define LINEARALGEBRA_CHOLESKY_H__INCLUDED_

#pragma once

#include <math.h>
#include "../../../Set/Algebraic/VectorSpace/Vector/Vector.h"
#include "../../../Set/Algebraic/VectorSpace/Category/Category.h"

namespace LinearAlgebra  
{
class Cholesky
{
public:

	Cholesky();
	Cholesky(const unsigned int &);
	virtual ~Cholesky();

	void Decomposition(Set::VectorSpace::Hom &);
	void Substitution(
		const Set::VectorSpace::Hom &, 
		Set::VectorSpace::Vector &);
	void Substitution(
		const Set::VectorSpace::Hom &, 
		const Set::VectorSpace::Vector &, 
		Set::VectorSpace::Vector &);
	void Solve(
		const Set::VectorSpace::Hom &, 
		Set::VectorSpace::Vector &);
	void Solve(
		const Set::VectorSpace::Hom &, 
		const Set::VectorSpace::Vector &, 
		Set::VectorSpace::Vector &);
	void Invert(
		const Set::VectorSpace::Hom &, 
		Set::VectorSpace::Hom &);
	double Determinant(const Set::VectorSpace::Hom &);

private:

	unsigned int n;
	double *p;
};

}

#endif // !defined(LINEARALGEBRA_CHOLESKY_H__INCLUDED_)
