// Factory.h: interface for the Factory class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(MATERIAL_UNIAXIAL_FACTORY_H__INCLUDED_)
#define MATERIAL_UNIAXIAL_FACTORY_H__INCLUDED_

#pragma once

#include "Uniaxial.h"

namespace Material
{
namespace Uniaxial
{
//////////////////////////////////////////////////////////////////////
// Class Factory
//////////////////////////////////////////////////////////////////////

class Factory
{
public: 

	Factory() {}
	virtual ~Factory() {}
	virtual Material::Uniaxial::LocalState * GetLS()=0;
	virtual Material::Uniaxial::Energy<0> * GetW()=0;
	virtual Material::Uniaxial::Energy<1> * GetDW()=0;
	virtual Material::Uniaxial::Energy<2> * GetDDW()=0;
	virtual Material::Uniaxial::Jet<0> * GetJ()=0;
	virtual Material::Uniaxial::Jet<1> * GetDJ()=0;
};

//////////////////////////////////////////////////////////////////////
// Class Builder
//////////////////////////////////////////////////////////////////////

class Builder 
{
public: 

	Builder() {}
	virtual ~Builder() {}
	virtual Material::Uniaxial::Factory * Build() const=0;
};

}

}

#endif // !defined(MATERIAL_UNIAXIAL_FACTORY_H__INCLUDED_
