// Factory.h: Factory for the LennardJones class.
// Copyright (c) by Bo Li - All rights reserved
//////////////////////////////////////////////////////////////////////

#if !defined(MATERIAL_RADIAL_LENNARDJONES_FACTORY_H__INCLUDED_)
#define MATERIAL_RADIAL_LENNARDJONES_FACTORY_H__INCLUDED_

#pragma once

#include <vector>
#include "./LennardJones.h"
#include "../Factory.h"

namespace Potential
{
namespace Radial
{
namespace LennardJones
{
//////////////////////////////////////////////////////////////////////
// Class Factory
//////////////////////////////////////////////////////////////////////

class Factory : public Potential::Radial::Factory
{
public: 

	typedef Potential::Radial::LennardJones::Data data_type;

	Factory() : LS(0), E(0), DE(0), DDE(0), J(0), DJ(0) {}
	virtual ~Factory() 
	{
		if (LS != 0)  delete LS;
		if (E != 0)   delete E;
		if (DE != 0)  delete DE;
		if (DDE != 0) delete DDE;
		if (J != 0)   delete J;
		if (DJ != 0)  delete DJ;
	}

	Factory(Potential::Radial::LennardJones::Data *Dat_) : 
		LS(0), E(0), DE(0), DDE(0), J(0), DJ(0), Dat(Dat_) {}

	Potential::Radial::LocalState * GetLS()
	{
		if (LS == 0) LS = new Potential::Radial::LennardJones::LocalState(Dat);
		return LS;
	}

	Potential::Radial::Energy<0> * GetE()
	{
		if (LS == 0) LS = new Potential::Radial::LennardJones::LocalState(Dat);
		if (E == 0) E = new Potential::Radial::LennardJones::Energy<0>(LS);
		return E;
	}

	Potential::Radial::Energy<1> * GetDE()
	{
		if (LS == 0) LS = new Potential::Radial::LennardJones::LocalState(Dat);
		if (DE == 0) DE = new Potential::Radial::LennardJones::Energy<1>(LS);
		return DE;
	}

	Potential::Radial::Energy<2> * GetDDE()
	{
		if (LS == 0) LS = new Potential::Radial::LennardJones::LocalState(Dat);
		if (DDE == 0) DDE = new Potential::Radial::LennardJones::Energy<2>(LS);
		return DDE;
	}

	Potential::Radial::Jet<0> * GetJ()
	{
		if (LS == 0) LS = new Potential::Radial::LennardJones::LocalState(Dat);
		if (J == 0) J = new Potential::Radial::LennardJones::Jet<0>(LS);
		return J;
	}

	Potential::Radial::Jet<1> * GetDJ()
	{
		if (LS == 0) LS = new Potential::Radial::LennardJones::LocalState(Dat);
		if (DJ == 0) DJ = new Potential::Radial::LennardJones::Jet<1>(LS);
		return DJ;
	}

private:

	Potential::Radial::LennardJones::Data *Dat;
	Potential::Radial::LennardJones::LocalState *LS;
	Potential::Radial::LennardJones::Energy<0> *E;
	Potential::Radial::LennardJones::Energy<1> *DE;
	Potential::Radial::LennardJones::Energy<2> *DDE;
	Potential::Radial::LennardJones::Jet<0> *J;
	Potential::Radial::LennardJones::Jet<1> *DJ;

private:

	Factory(const Factory &);
	Factory & operator = (const Factory &);
};

//////////////////////////////////////////////////////////////////////
// Class Builder
//////////////////////////////////////////////////////////////////////

class Builder : public Potential::Radial::Builder
{
public: 

	typedef Potential::Radial::LennardJones::Data data_type;

	Builder() {}
	virtual ~Builder() {}

	Builder(Potential::Radial::LennardJones::Data *Dat_) : Dat(Dat_) {}

	Potential::Radial::Factory * Build() const
	{
		return new Potential::Radial::LennardJones::Factory(Dat);
	}

private:

	Potential::Radial::LennardJones::Data *Dat;
};

}

}

}

#endif // !defined(MATERIAL_RADIAL_LENNARDJONES_FACTORY_H__INCLUDED_
