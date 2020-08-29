// Diagonal.h: interface for the Diagonal class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(SET_VECTORSPACE_CATEGORY_DIAGONAL_H__INCLUDED_)
#define SET_VECTORSPACE_CATEGORY_DIAGONAL_H__INCLUDED_

#pragma once

#include <math.h>
#include <stdlib.h>
#include <cassert>

#include "../Category.h"
#include "../../Vector/Vector.h"

namespace Set
{
namespace VectorSpace
{
//////////////////////////////////////////////////////////////////////
// Class Diagonal
//////////////////////////////////////////////////////////////////////

class Diagonal : public Set::VectorSpace::Vector
{
public:

	Diagonal();
	Diagonal(const unsigned int &);
	Diagonal(const unsigned int &, double * const &);
	virtual ~Diagonal();
	Diagonal(const Diagonal &);
	Diagonal(const Set::VectorSpace::Vector &);
	Diagonal(const Set::VectorSpace::Hom &);
	Diagonal & operator = (const Diagonal &);
	Diagonal & operator = (const Set::VectorSpace::Vector &);
	Set::VectorSpace::Vector 
		operator () (const Set::VectorSpace::Vector &) const;
	double operator () (const Set::VectorSpace::Vector &, 
		const Set::VectorSpace::Vector &) const;
	double operator () (const Diagonal &) const;
	double operator () (const Set::VectorSpace::Hom &) const;
	Set::VectorSpace::Hom Embed() const;
};

}

}

//////////////////////////////////////////////////////////////////////
// Operators
//////////////////////////////////////////////////////////////////////
/*
bool 
operator != (const Set::VectorSpace::Diagonal &, 
			 const Set::VectorSpace::Hom &);

bool 
operator != (const Set::VectorSpace::Hom &,
			 const Set::VectorSpace::Diagonal &);

bool 
operator == (const Set::VectorSpace::Diagonal &,
			 const Set::VectorSpace::Hom &);

bool 
operator == (const Set::VectorSpace::Hom &,
			 const Set::VectorSpace::Diagonal &);
*/
void Random(Set::VectorSpace::Diagonal &);

Set::VectorSpace::Diagonal
operator + (const Set::VectorSpace::Diagonal &, 
			const Set::VectorSpace::Diagonal &);

Set::VectorSpace::Hom
operator + (const Set::VectorSpace::Diagonal &, 
			const Set::VectorSpace::Hom &);

Set::VectorSpace::Hom
operator + (const Set::VectorSpace::Hom &,
			const Set::VectorSpace::Diagonal &);

Set::VectorSpace::Diagonal
operator - (const Set::VectorSpace::Diagonal &, 
			const Set::VectorSpace::Diagonal &);

Set::VectorSpace::Hom
operator - (const Set::VectorSpace::Diagonal &, 
			const Set::VectorSpace::Hom &);

Set::VectorSpace::Hom
operator - (const Set::VectorSpace::Hom &,
			const Set::VectorSpace::Diagonal &); 

Set::VectorSpace::Diagonal 
operator - (const Set::VectorSpace::Diagonal &);

void Null(Set::VectorSpace::Diagonal &);

Set::VectorSpace::Diagonal 
operator * (const Set::VectorSpace::Diagonal &, 
			const double &);

Set::VectorSpace::Diagonal 
operator * (const double &, 
			const Set::VectorSpace::Diagonal &);

Set::VectorSpace::Diagonal 
operator / (const Set::VectorSpace::Diagonal &, 
			const double &);

Set::VectorSpace::Vector
operator * (const Set::VectorSpace::Diagonal &, 
			const Set::VectorSpace::Vector &);

Set::VectorSpace::Diagonal 
operator * (const Set::VectorSpace::Diagonal &, 
			const Set::VectorSpace::Diagonal &); 

Set::VectorSpace::Hom
operator * (const Set::VectorSpace::Diagonal &, 
			const Set::VectorSpace::Hom &);

Set::VectorSpace::Hom
operator * (const Set::VectorSpace::Hom &, 
			const Set::VectorSpace::Diagonal &);

Set::VectorSpace::Diagonal
Adjoint(const Set::VectorSpace::Diagonal &);

double Norm(const Set::VectorSpace::Diagonal &);

Set::VectorSpace::Diagonal
Exp(const Set::VectorSpace::Diagonal &);

Set::VectorSpace::Diagonal
Log(const Set::VectorSpace::Diagonal &);

Set::VectorSpace::Diagonal
Inverse(const Set::VectorSpace::Diagonal &);

void Identity(Set::VectorSpace::Diagonal &);

double Jacobian(const Set::VectorSpace::Diagonal &);
 
double Trace(const Set::VectorSpace::Diagonal &);

#endif // !defined(SET_VECTORSPACE_CATEGORY_DIAGONAL_H__INCLUDED_)
