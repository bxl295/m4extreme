// Factory.h: interface for the Factory class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(POTENTIAL_FIELD_FACTORY_H__INCLUDED_)
#define POTENTIAL_FIELD_FACTORY_H__INCLUDED_

#pragma once

#include "Field.h"

namespace Potential
{
namespace Field
{
//////////////////////////////////////////////////////////////////////
// Class Factory
//////////////////////////////////////////////////////////////////////

class Factory
{
public: 

	Factory() {}
	virtual ~Factory() {}
	virtual Potential::Field::LocalState * GetLS()=0;
	virtual Potential::Field::Energy<0> * GetW()=0;
	virtual Potential::Field::Energy<1> * GetDW()=0;
	virtual Potential::Field::Energy<2> * GetDDW()=0;
	virtual Potential::Field::Jet<0> * GetJ()=0;
	virtual Potential::Field::Jet<1> * GetDJ()=0;
};

//////////////////////////////////////////////////////////////////////
// Class Builder
//////////////////////////////////////////////////////////////////////

class Builder 
{
public: 

	Builder() {}
	virtual ~Builder() {}
	virtual Potential::Field::Factory * Build() const=0;
};

}

}

#endif // !defined(POTENTIAL_FIELD_FACTORY_H__INCLUDED_
