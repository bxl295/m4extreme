// Factory.h: interface for the Factory class.
// Copyright (c) by Bo Li - All rights reserved
//////////////////////////////////////////////////////////////////////

#if !defined(POTENTIAL_RADIAL_FACTORY_H__INCLUDED_)
#define POTENTIAL_RADIAL_FACTORY_H__INCLUDED_

#pragma once

#include "Radial.h"

namespace Potential
{
namespace Radial
{
//////////////////////////////////////////////////////////////////////
// Class Factory
//////////////////////////////////////////////////////////////////////

class Factory
{
public: 

	Factory() {}
	virtual ~Factory() {}
	virtual Potential::Radial::LocalState * GetLS()=0;
	virtual Potential::Radial::Energy<0> * GetE()=0;
	virtual Potential::Radial::Energy<1> * GetDE()=0;
	virtual Potential::Radial::Energy<2> * GetDDE()=0;
	virtual Potential::Radial::Jet<0> * GetJ()=0;
	virtual Potential::Radial::Jet<1> * GetDJ()=0;
};

//////////////////////////////////////////////////////////////////////
// Class Builder
//////////////////////////////////////////////////////////////////////

class Builder 
{
public: 

	Builder() {}
	virtual ~Builder() {}
	virtual Potential::Radial::Factory * Build() const=0;
};

}

}

#endif // !defined(POTENTIAL_RADIAL_FACTORY_H__INCLUDED_
