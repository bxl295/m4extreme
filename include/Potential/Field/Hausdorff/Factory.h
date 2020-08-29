// Factory.h: Factory for the Hausdorff class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(POTENTIAL_FIELD_HAUSDORFF_FACTORY_H__INCLUDED_)
#define POTENTIAL_FIELD_HAUSDORFF_FACTORY_H__INCLUDED_

#pragma once

#include <vector>
#include "./Hausdorff.h"
#include "Potential/Field/Factory.h"

namespace Potential
{
namespace Field
{
namespace Hausdorff
{
//////////////////////////////////////////////////////////////////////
// Class Factory
//////////////////////////////////////////////////////////////////////

class Factory : public Potential::Field::Factory
{
public: 

	typedef Potential::Field::Hausdorff::Data data_type;

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

	Factory(Potential::Field::Hausdorff::Data *Dat_) : 
		LS(0), W(0), DW(0), DDW(0), J(0), DJ(0), Dat(Dat_) {}

	Potential::Field::LocalState * GetLS()
	{
		if (LS == 0) LS = new Potential::Field::Hausdorff::LocalState(Dat);
		return LS;
	}

	Potential::Field::Energy<0> * GetW()
	{
		if (LS == 0) LS = new Potential::Field::Hausdorff::LocalState(Dat);
		if (W == 0) W = new Potential::Field::Hausdorff::Energy<0>(LS);
		return W;
	}

	Potential::Field::Energy<1> * GetDW()
	{
		if (LS == 0) LS = new Potential::Field::Hausdorff::LocalState(Dat);
		if (DW == 0) DW = new Potential::Field::Hausdorff::Energy<1>(LS);
		return DW;
	}

	Potential::Field::Energy<2> * GetDDW()
	{
		if (LS == 0) LS = new Potential::Field::Hausdorff::LocalState(Dat);
		if (DDW == 0) DDW = new Potential::Field::Hausdorff::Energy<2>(LS);
		return DDW;
	}

	Potential::Field::Jet<0> * GetJ()
	{
		if (LS == 0) LS = new Potential::Field::Hausdorff::LocalState(Dat);
		if (J == 0) J = new Potential::Field::Hausdorff::Jet<0>(LS);
		return J;
	}

	Potential::Field::Jet<1> * GetDJ()
	{
		if (LS == 0) LS = new Potential::Field::Hausdorff::LocalState(Dat);
		if (DJ == 0) DJ = new Potential::Field::Hausdorff::Jet<1>(LS);
		return DJ;
	}

private:

	Potential::Field::Hausdorff::Data *Dat;
	Potential::Field::Hausdorff::LocalState *LS;
	Potential::Field::Hausdorff::Energy<0> *W;
	Potential::Field::Hausdorff::Energy<1> *DW;
	Potential::Field::Hausdorff::Energy<2> *DDW;
	Potential::Field::Hausdorff::Jet<0> *J;
	Potential::Field::Hausdorff::Jet<1> *DJ;

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

	typedef Potential::Field::Hausdorff::Data data_type;

	Builder() {}
	virtual ~Builder() {}

	Builder(Potential::Field::Hausdorff::Data *Dat_) : Dat(Dat_) {}

	Potential::Field::Factory * Build() const
	{
		return new Potential::Field::Hausdorff::Factory(Dat);
	}

private:

	Potential::Field::Hausdorff::Data *Dat;
};

}

}

}

#endif // !defined(POTENTIAL_FIELD_HAUSDORFF_FACTORY_H__INCLUDED_
