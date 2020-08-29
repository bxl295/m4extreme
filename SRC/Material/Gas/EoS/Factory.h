// Factory.h: interface for the Factory class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(MATERIAL_GAS_EOS_FACTORY_H__INCLUDED_)
#define MATERIAL_GAS_EOS_FACTORY_H__INCLUDED_

#pragma once

#include "EoS.h"

namespace Material
{
namespace Gas
{
namespace EoS
{
//////////////////////////////////////////////////////////////////////
// Class Factory
//////////////////////////////////////////////////////////////////////

class Factory
{
public: 

	Factory() {}
	virtual ~Factory() {}
	virtual Material::Gas::EoS::LocalState * GetLS()=0;
	virtual Material::Gas::EoS::Energy<0> * GetW()=0;
	virtual Material::Gas::EoS::Energy<1> * GetDW()=0;
	virtual Material::Gas::EoS::Energy<2> * GetDDW()=0;
	virtual Material::Gas::EoS::Jet<0> * GetJ()=0;
	virtual Material::Gas::EoS::Jet<1> * GetDJ()=0;
};

//////////////////////////////////////////////////////////////////////
// Class Builder
//////////////////////////////////////////////////////////////////////

class Builder 
{
public: 

	Builder() {}
	virtual ~Builder() {}
	virtual Material::Gas::EoS::Factory * Build() const=0;
};

}

}

}

#endif // !defined(MATERIAL_GAS_EOS_FACTORY_H__INCLUDED_
