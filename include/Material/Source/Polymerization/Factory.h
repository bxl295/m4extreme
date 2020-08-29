// Factory.h: Factory for the Source::Polymerization class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
////////////////////////////////////////////////////////////////////////////

#if !defined(MATERIAL_SOURCE_POLYMERIZATION_FACTORY_H__INCLUDED_)
#define MATERIAL_SOURCE_POLYMERIZATION_FACTORY_H__INCLUDED_

#pragma once

#include <vector>
#include "./Polymerization.h"
#include "Material/Factory.h"

namespace Material
{
namespace Source
{
namespace Polymerization
{
//////////////////////////////////////////////////////////////////////
// Class Factory
//////////////////////////////////////////////////////////////////////

class Factory : public Material::Factory
{
public: 

	typedef Material::Source::Polymerization::Data data_type;

	Factory() : LS(0), W(0), DW(0), DDW(0), J(0), DJ(0) {}
	virtual ~Factory() 
	{
		if (LS != 0)  delete LS;
		if (W != 0)   delete W;
	}

        Factory(double T0_, Clock *Chronos_, Material::Source::Polymerization::Data *Dat_) : 
	T0(T0_), Chronos(Chronos_), LS(0), W(0), DW(0), DDW(0), J(0), DJ(0), Dat(Dat_) {}

	Material::LocalState * GetLS()
	{
	  if (LS == 0) LS = new Material::Source::Polymerization::LocalState(T0, Chronos, Dat);
	  return LS;
	}

	Material::Energy<0> * GetW()
	{
	  if (LS == 0) LS = new Material::Source::Polymerization::LocalState(T0, Chronos, Dat);
	  if (W == 0) W = new Material::Source::Polymerization::Energy<0>(LS);
	  return W;
	}

	Material::Energy<1> * GetDW()
	{
		return NULL;
	}

	Material::Energy<2> * GetDDW()
	{
		return NULL;
	}

	Material::Jet<0> * GetJ()
	{
		return NULL;
	}

	Material::Jet<1> * GetDJ()
	{
		return NULL;
	}

private:
	double T0;
	Clock *Chronos;
	Material::Source::Polymerization::Data *Dat;
	Material::Source::Polymerization::LocalState *LS;
	Material::Source::Polymerization::Energy<0> *W;
	Material::Source::Polymerization::Energy<1> *DW;
	Material::Source::Polymerization::Energy<2> *DDW;
	Material::Source::Polymerization::Jet<0> *J;
	Material::Source::Polymerization::Jet<1> *DJ;

private:

	Factory(const Factory &);
	Factory & operator = (const Factory &);
};

//////////////////////////////////////////////////////////////////////
// Class Builder
//////////////////////////////////////////////////////////////////////

class Builder : public Material::Builder
{
public: 

	typedef Material::Source::Polymerization::Data data_type;
	typedef Material::Source::Polymerization::Factory fact_type;

	Builder() {}
	virtual ~Builder() {}

        Builder(double T0_, Clock *Chronos_, data_type *Dat_) : 
	Dat(Dat_), Chronos(Chronos_), T0(T0_) {}

	Material::Factory * Build() const
	{
	  return new fact_type(T0, Chronos, Dat);
	}

private:
	double T0;
	Clock *Chronos;
	data_type *Dat;
};

}

}

}

#endif // !defined(MATERIAL_SOURCE_POLYMERIZATION_FACTORY_H__INCLUDED_
