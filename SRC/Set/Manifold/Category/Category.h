// Category.h: Interface for the Category class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(SET_MANIFOLD_CATEGORY_H__INCLUDED_)
#define SET_MANIFOLD_CATEGORY_H__INCLUDED_

#pragma once

#include "../Manifold.h"
#include "../../Algebraic/VectorSpace/Vector/Vector.h"
#include "../../Algebraic/VectorSpace/Category/Category.h"

namespace Set
{
namespace Manifold
{
class Map;
class Field;
class Scalar;
class TMap;

//////////////////////////////////////////////////////////////////////
// Map class
//////////////////////////////////////////////////////////////////////

class Map
{
public: 

	Map();
	virtual ~Map();
	virtual Map *Clone()=0;
	virtual TMap *Diff()=0;
	virtual void Randomize()=0;
	virtual const Point & operator () (const Point &)=0;
	virtual unsigned int size1() const=0;
	virtual unsigned int size2() const=0;
};

//////////////////////////////////////////////////////////////////////
// Field class
//////////////////////////////////////////////////////////////////////

class Field
{
public: 

	Field();
	virtual ~Field();
	virtual Field *Clone()=0;
	virtual TMap *Diff()=0;
	virtual void Randomize()=0;
	virtual const Set::VectorSpace::Vector & operator () (const Point &)=0;
	virtual unsigned int size1() const=0;
	virtual unsigned int size2() const=0;
};

//////////////////////////////////////////////////////////////////////
// Scalar class
//////////////////////////////////////////////////////////////////////

class Scalar
{
public: 

	Scalar();
	virtual ~Scalar();
	virtual Scalar *Clone()=0;
	virtual Field *Diff()=0;
	virtual void Randomize()=0;
	virtual const double & operator () (const Point &)=0;
	virtual unsigned int size1() const=0;
	virtual unsigned int size2() const=0;
};

//////////////////////////////////////////////////////////////////////
// TMap class
//////////////////////////////////////////////////////////////////////

class TMap
{
public: 

	TMap();
	virtual ~TMap();
	virtual TMap *Clone()=0;
	virtual TMap *Diff()=0;
	virtual void Randomize()=0;
	virtual const Set::VectorSpace::Hom & operator () (const Point &)=0;
	virtual unsigned int size1() const=0;
	virtual unsigned int size2() const=0;
};

}

}

#endif // !defined(SET_MANIFOLD_CATEGORY_H__INCLUDED_)
