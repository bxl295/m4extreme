// Factory.h: interface for the Factory class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(MATERIAL_FACTORY_H__INCLUDED_)
#define MATERIAL_FACTORY_H__INCLUDED_

#pragma once

#include "Material.h"

namespace Material
{
//////////////////////////////////////////////////////////////////////
// Class Factory
//////////////////////////////////////////////////////////////////////

class Factory
{
public: 

	Factory() {}
	virtual ~Factory() {}
	virtual Material::LocalState * GetLS()=0;
	virtual Material::Energy<0> * GetW()=0;
	virtual Material::Energy<1> * GetDW()=0;
	virtual Material::Energy<2> * GetDDW()=0;
	virtual Material::Jet<0> * GetJ()=0;
	virtual Material::Jet<1> * GetDJ()=0;
};

//////////////////////////////////////////////////////////////////////
// Class Builder
//////////////////////////////////////////////////////////////////////

class Builder 
{
public: 

	Builder() {}
	virtual ~Builder() {}
	virtual Material::Factory * Build() const=0;
};

//////////////////////////////////////////////////////////////////////
// Class LocalBuilder
//////////////////////////////////////////////////////////////////////

class LocalBuilder 
{
public: 

	LocalBuilder() {}
	virtual ~LocalBuilder() {}
	virtual Material::Factory * Build(const Set::VectorSpace::Vector &)=0;
};

}

#endif // !defined(MATERIAL_FACTORY_H__INCLUDED_
