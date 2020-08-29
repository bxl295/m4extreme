// Factory.h: Factory for the Coulomb class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(POTENTIAL_FIELD_COULOMB_FACTORY_H__INCLUDED_)
#define POTENTIAL_FIELD_COULOMB_FACTORY_H__INCLUDED_

#pragma once

#include <vector>
#include "./Coulomb.h"
#include "../Factory.h"

namespace Potential
{
namespace Field
{
namespace Coulomb
{
//////////////////////////////////////////////////////////////////////
// Class Factory
//////////////////////////////////////////////////////////////////////

class Factory : public Potential::Field::Factory
{
public: 

	typedef Potential::Field::Coulomb::Data data_type;

	Factory() : LS(0), W(0), DW(0), DDW(0), J(0), DJ(0) {}
	virtual ~Factory() 
	{
		if (LS != 0)  delete LS;
		if (W != 0)   delete W;
		if (DW != 0)  delete DW;
		if (DDW != 0) delete DDW;
		if (J != 0)   delete J;
		if (DJ != 0)  delete DJ;
	}

	Factory(Potential::Field::Coulomb::Data *Dat_) : 
		LS(0), W(0), DW(0), DDW(0), J(0), DJ(0), Dat(Dat_) {}

	Potential::Field::LocalState * GetLS()
	{
		if (LS == 0) LS = new Potential::Field::Coulomb::LocalState(Dat);
		return LS;
	}

	Potential::Field::Energy<0> * GetW()
	{
		if (LS == 0) LS = new Potential::Field::Coulomb::LocalState(Dat);
		if (W == 0) W = new Potential::Field::Coulomb::Energy<0>(LS);
		return W;
	}

	Potential::Field::Energy<1> * GetDW()
	{
		if (LS == 0) LS = new Potential::Field::Coulomb::LocalState(Dat);
		if (DW == 0) DW = new Potential::Field::Coulomb::Energy<1>(LS);
		return DW;
	}

	Potential::Field::Energy<2> * GetDDW()
	{
		if (LS == 0) LS = new Potential::Field::Coulomb::LocalState(Dat);
		if (DDW == 0) DDW = new Potential::Field::Coulomb::Energy<2>(LS);
		return DDW;
	}

	Potential::Field::Jet<0> * GetJ()
	{
		if (LS == 0) LS = new Potential::Field::Coulomb::LocalState(Dat);
		if (J == 0) J = new Potential::Field::Coulomb::Jet<0>(LS);
		return J;
	}

	Potential::Field::Jet<1> * GetDJ()
	{
		if (LS == 0) LS = new Potential::Field::Coulomb::LocalState(Dat);
		if (DJ == 0) DJ = new Potential::Field::Coulomb::Jet<1>(LS);
		return DJ;
	}

private:

	Potential::Field::Coulomb::Data *Dat;
	Potential::Field::Coulomb::LocalState *LS;
	Potential::Field::Coulomb::Energy<0> *W;
	Potential::Field::Coulomb::Energy<1> *DW;
	Potential::Field::Coulomb::Energy<2> *DDW;
	Potential::Field::Coulomb::Jet<0> *J;
	Potential::Field::Coulomb::Jet<1> *DJ;

private:

	Factory(const Factory &);
	Factory & operator = (const Factory &);
};

//////////////////////////////////////////////////////////////////////
// Class Builder
//////////////////////////////////////////////////////////////////////

class Builder : public Potential::Field::Builder
{
public: 

	typedef Potential::Field::Coulomb::Data data_type;

	Builder() {}
	virtual ~Builder() {}

	Builder(Potential::Field::Coulomb::Data *Dat_) : Dat(Dat_) {}

	Potential::Field::Factory * Build() const
	{
		return new Potential::Field::Coulomb::Factory(Dat);
	}

private:

	Potential::Field::Coulomb::Data *Dat;
};

}

}

}

#endif // !defined(POTENTIAL_FIELD_COULOMB_FACTORY_H__INCLUDED_
