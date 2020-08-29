// MultiHom.h: interface for the MultiHom class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(SET_VECTORSPACE_MULTIHOM_H__INCLUDED_)
#define SET_VECTORSPACE_MULTIHOM_H__INCLUDED_

#pragma once

#include<vector>
#include "./Category.h"
#include "../Vector/Vector.h"
#include "../Vector/MultiVector.h"

using namespace std;

bool 
operator != (const vector<vector<Set::VectorSpace::Hom> > &, 
			 const vector<vector<Set::VectorSpace::Hom> > &);
bool 
operator == (const vector<vector<Set::VectorSpace::Hom> > &, 
			 const vector<vector<Set::VectorSpace::Hom> > &);
vector<Set::VectorSpace::Vector>
operator * (const vector<vector<Set::VectorSpace::Hom> > &,
		    const vector<Set::VectorSpace::Vector> &);

#endif // !defined(SET_VECTORSPACE_MULTIHOM_H__INCLUDED_)
