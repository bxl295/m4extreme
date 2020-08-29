// Factory.h: interface for the Factory class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(MATERIAL_SHEAR_FACTORY_H__INCLUDED_)
#define MATERIAL_SHEAR_FACTORY_H__INCLUDED_

#pragma once

#include "Shear.h"

namespace Material
{
namespace Shear
{
//////////////////////////////////////////////////////////////////////
// Class Factory
//////////////////////////////////////////////////////////////////////

class Factory
{
public: 

	Factory() {}
	virtual ~Factory() {}
	virtual Material::Shear::LocalState * GetLS()=0;
	virtual Material::Shear::Modulus<0> * GetG()=0;
	virtual Material::Shear::Modulus<1> * GetDG()=0;
	virtual Material::Shear::Modulus<2> * GetDDG()=0;
	virtual Material::Shear::Jet<0> * GetJ()=0;
	virtual Material::Shear::Jet<1> * GetDJ()=0;
	virtual Material::Shear::JetJet<0> * GetJJ()=0;
};

//////////////////////////////////////////////////////////////////////
// Class Builder
//////////////////////////////////////////////////////////////////////

class Builder 
{
public: 

	Builder() {}
	virtual ~Builder() {}
	virtual Material::Shear::Factory * Build() const=0;
};

}

}

#endif // !defined(MATERIAL_SHEAR_FACTORY_H__INCLUDED_
