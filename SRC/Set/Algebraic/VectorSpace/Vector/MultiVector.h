// MultiVector.h: interface for the MultiVector class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(SET_MULTIVECTORSPACE_MULTIVECTOR_H__INCLUDED_)
#define SET_MULTIVECTORSPACE_MULTIVECTOR_H__INCLUDED_

#pragma once

#include<vector>
#include "./Vector.h"

using namespace std;

bool 
operator != (const vector<Set::VectorSpace::Vector> &, 
			 const vector<Set::VectorSpace::Vector> &);
bool 
operator == (const vector<Set::VectorSpace::Vector> &, 
			 const vector<Set::VectorSpace::Vector> &);
vector<Set::VectorSpace::Vector>
operator + (const vector<Set::VectorSpace::Vector> &, 
			const vector<Set::VectorSpace::Vector> &);
vector<Set::VectorSpace::Vector>
operator - (const vector<Set::VectorSpace::Vector> &, 
			const vector<Set::VectorSpace::Vector> &);
void
operator *= (vector<Set::VectorSpace::Vector> &, 
			 const double &);
void
operator /= (vector<Set::VectorSpace::Vector> &, 
			 const double &);
vector<Set::VectorSpace::Vector>
operator * (const vector<Set::VectorSpace::Vector> &, 
			const double &);
vector<Set::VectorSpace::Vector>
operator * (const double &, 
			const vector<Set::VectorSpace::Vector> &);
vector<Set::VectorSpace::Vector>
operator * (const vector<Set::VectorSpace::Vector> &, 
			const double &);
vector<Set::VectorSpace::Vector>
operator / (const vector<Set::VectorSpace::Vector> &, 
			const double &);
double 
operator * (const vector<Set::VectorSpace::Vector> &, 
			const vector<Set::VectorSpace::Vector> &);
void 
Random(vector<Set::VectorSpace::Vector> &);

#endif // !defined(SET_MULTIVECTORSPACE_MULTIVECTOR_H__INCLUDED_)
