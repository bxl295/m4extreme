// Factory.h: interface for the Factory class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(SET_MANIFOLD_FACTORY_H__INCLUDED_)
#define SET_MANIFOLD_FACTORY_H__INCLUDED_

#pragma once

#include "./Manifold.h"
#include "./Category/Category.h"

namespace Set
{
namespace Manifold
{	
//////////////////////////////////////////////////////////////////////
// Class Factory
//////////////////////////////////////////////////////////////////////

class Factory
{
public: 

	Factory() {}
	virtual ~Factory() {} 

	virtual Set::Manifold::Map *GetEmb(const Set::Manifold::Point &)=0;
	virtual Set::Manifold::TMap *GetDEmb(const Set::Manifold::Point &)=0;
	virtual Set::Manifold::TMap *GetDDEmb(const Set::Manifold::Point &)=0;
	virtual Set::Manifold::Map *GetSub(const Set::Manifold::Point &)=0;
	virtual Set::Manifold::TMap *GetDSub(const Set::Manifold::Point &)=0;
	virtual Set::Manifold::TMap *GetDDSub(const Set::Manifold::Point &)=0;
	virtual void operator += (const Set::VectorSpace::Vector &)=0;
	virtual void operator -= (const Set::VectorSpace::Vector &)=0;
};

}

}

#endif // !defined(SET_MANIFOLD_FACTORY_H__INCLUDED_
