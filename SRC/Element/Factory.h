// Factory.h: interface for the Factory class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(ELEMENT_FACTORY_H__INCLUDED_)
#define ELEMENT_FACTORY_H__INCLUDED_

#pragma once

#include "Element.h"

namespace Element
{
//////////////////////////////////////////////////////////////////////
// Class Factory
//////////////////////////////////////////////////////////////////////

class Factory
{
public: 

	Factory() {}
	virtual ~Factory() {}
	virtual Element::LocalState * GetLS()=0;
	virtual Element::Energy<0> * GetE()=0;
	virtual Element::Energy<1> * GetDE()=0;
	virtual Element::Energy<2> * GetDDE()=0;
	virtual Element::Jet<0> * GetJ()=0;
	virtual Element::Jet<1> * GetDJ()=0;
};

}

#endif // !defined(ELEMENT_FACTORY_H__INCLUDED_
