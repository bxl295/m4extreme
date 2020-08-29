// Factory.h: Factory for the ShearViscosity class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(MATERIAL_UNIAXIAL_SHEARVISCOSITY_FACTORY_H__INCLUDED_)
#define MATERIAL_UNIAXIAL_SHEARVISCOSITY_FACTORY_H__INCLUDED_

#pragma once

#include <vector>
#include "./ShearViscosity.h"
#include "Material/Uniaxial/Factory.h"

namespace Material
{
namespace Uniaxial
{
namespace ShearViscosity
{
//////////////////////////////////////////////////////////////////////
// Class Factory
//////////////////////////////////////////////////////////////////////

class Factory : public Material::Uniaxial::Factory
{
public: 

	typedef Material::Uniaxial::ShearViscosity::Data data_type;

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

	Factory(Material::Uniaxial::ShearViscosity::Data *Dat_) :
		LS(0), W(0), DW(0), DDW(0), J(0), DJ(0), Dat(Dat_) {}

	Material::Uniaxial::LocalState * GetLS()
	{
		if (LS == 0) LS = new Material::Uniaxial::ShearViscosity::LocalState(Dat);
		return LS;
	}

	Material::Uniaxial::Energy<0> * GetW()
	{
		if (LS == 0) LS = new Material::Uniaxial::ShearViscosity::LocalState(Dat);
		if (W == 0) W = new Material::Uniaxial::ShearViscosity::Energy<0>(LS);
		return W;
	}

	Material::Uniaxial::Energy<1> * GetDW()
	{
		if (LS == 0) LS = new Material::Uniaxial::ShearViscosity::LocalState(Dat);
		if (DW == 0) DW = new Material::Uniaxial::ShearViscosity::Energy<1>(LS);
		return DW;
	}

	Material::Uniaxial::Energy<2> * GetDDW()
	{
		if (LS == 0) LS = new Material::Uniaxial::ShearViscosity::LocalState(Dat);
		if (DDW == 0) DDW = new Material::Uniaxial::ShearViscosity::Energy<2>(LS);
		return DDW;
	}

	Material::Uniaxial::Jet<0> * GetJ()
	{
		if (LS == 0) LS = new Material::Uniaxial::ShearViscosity::LocalState(Dat);
		if (J == 0) J = new Material::Uniaxial::ShearViscosity::Jet<0>(LS);
		return J;
	}

	Material::Uniaxial::Jet<1> * GetDJ()
	{
		if (LS == 0) LS = new Material::Uniaxial::ShearViscosity::LocalState(Dat);
		if (DJ == 0) DJ = new Material::Uniaxial::ShearViscosity::Jet<1>(LS);
		return DJ;
	}

private:

	Material::Uniaxial::ShearViscosity::Data *Dat;
	Material::Uniaxial::ShearViscosity::LocalState *LS;
	Material::Uniaxial::ShearViscosity::Energy<0> *W;
	Material::Uniaxial::ShearViscosity::Energy<1> *DW;
	Material::Uniaxial::ShearViscosity::Energy<2> *DDW;
	Material::Uniaxial::ShearViscosity::Jet<0> *J;
	Material::Uniaxial::ShearViscosity::Jet<1> *DJ;

private:

	Factory(const Factory &);
	Factory & operator = (const Factory &);
};

//////////////////////////////////////////////////////////////////////
// Class Builder
//////////////////////////////////////////////////////////////////////

class Builder : public Material::Uniaxial::Builder
{
public: 

	typedef Material::Uniaxial::ShearViscosity::Data data_type;

	Builder() {}
	virtual ~Builder() {}

	Builder(Material::Uniaxial::ShearViscosity::Data *Dat_) : Dat(Dat_) {}

	Material::Uniaxial::Factory * Build() const
	{
		return new Material::Uniaxial::ShearViscosity::Factory(Dat);
	}

private:

	Material::Uniaxial::ShearViscosity::Data *Dat;
};

}

}

}

#endif // !defined(MATERIAL_UNIAXIAL_SHEARVISCOSITY_FACTORY_H__INCLUDED_
