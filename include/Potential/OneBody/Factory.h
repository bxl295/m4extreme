// Factory.h: interface for the Factory class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(POTENTIAL_ONEBODY_FACTORY_H__INCLUDED_)
#define POTENTIAL_ONEBODY_FACTORY_H__INCLUDED_

#pragma once

#include <set>
#include <map>

#include "../Field/Factory.h"
#include "../../Element/Element.h"
#include "../../Element/Factory.h"
#include "../../Set/Manifold/Manifold.h"

namespace Potential
{
namespace OneBody
{
//////////////////////////////////////////////////////////////////////
// Class Factory
//////////////////////////////////////////////////////////////////////

class Factory : public Element::Factory
{
public: 

	Factory() : PFF(), LS(0), E(0), DE(0), DDE(0), J(0), DJ(0) {}

	virtual ~Factory() 
	{
		if (DJ     != 0) delete DJ;
		if (J      != 0) delete J;
		if (DDE    != 0) delete DDE;
		if (DE     != 0) delete DE;
		if (E      != 0) delete E;
		if (LS     != 0) delete LS;

		delete PFF;
	}

	Factory(Potential::Field::Builder *PFB, Set::Manifold::Point *x) : 
		PFF(PFB->Build()), LS(0), E(0), DE(0), DDE(0), J(0), DJ(0) 
	{
		LS = new Potential::OneBody::LocalState(x);
	}

	Element::LocalState * GetLS()
	{
		return LS;
	}

	Element::Energy<0> * GetE()
	{
		if (E == 0) E = new Potential::OneBody::Energy<0>(LS,PFF->GetW());
		return E;
	}

	Element::Energy<1> * GetDE()
	{
		if (DE == 0) DE = new Potential::OneBody::Energy<1>(LS,PFF->GetDW());
		return DE;
	}

	Element::Energy<2> * GetDDE()
	{
		if (DDE == 0) DDE = new Potential::OneBody::Energy<2>(LS,PFF->GetDDW());
		return DDE;
	}

	Element::Jet<0> * GetJ()
	{
		if (J == 0) J = new Potential::OneBody::Jet<0>(LS,PFF->GetJ());
		return J;
	}

	Element::Jet<1> * GetDJ()
	{
		if (DJ == 0) DJ = new Potential::OneBody::Jet<1>(LS,PFF->GetDJ());
		return DJ;
	}

private:

	Potential::Field::Factory *PFF;
	Potential::OneBody::LocalState *LS;
	Potential::OneBody::Energy<0> *E;
	Potential::OneBody::Energy<1> *DE;
	Potential::OneBody::Energy<2> *DDE;
	Potential::OneBody::Jet<0> *J;
	Potential::OneBody::Jet<1> *DJ;

private:

	Factory(const Factory &);
	Factory & operator = (const Factory &);
};

//////////////////////////////////////////////////////////////////////
// Class Builder
//////////////////////////////////////////////////////////////////////

class Builder
{
public: 

	Builder() {}
	virtual ~Builder() {}

	Builder(Potential::Field::Builder *PFB_) : PFB(PFB_) {}

	Potential::OneBody::Factory * Build(Set::Manifold::Point *x) const
	{
		return new Potential::OneBody::Factory(PFB,x);
	}

private:

	Potential::Field::Builder *PFB;
};

}

}

#endif // !defined(POTENTIAL_ONEBODY_FACTORY_H__INCLUDED_
