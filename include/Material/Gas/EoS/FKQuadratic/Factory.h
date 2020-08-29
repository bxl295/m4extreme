// Factory.h: Factory for the FKQuadratic class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(MATERIAL_GAS_EOS_FKQUADRATIC_FACTORY_H__INCLUDED_)
#define MATERIAL_GAS_EOS_FKQUADRATIC_FACTORY_H__INCLUDED_

#pragma once

#include <vector>
#include "./FKQuadratic.h"
#include "Material/Gas/EoS/Factory.h"

namespace Material
{
namespace Gas
{
namespace EoS
{
namespace FKQuadratic
{
//////////////////////////////////////////////////////////////////////
// Class Factory
//////////////////////////////////////////////////////////////////////

class Factory : public Material::Gas::EoS::Factory
{
public: 

	typedef Material::Gas::EoS::FKQuadratic::Data data_type;

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

	Factory(Material::Gas::EoS::FKQuadratic::Data *Dat_) : 
		LS(0), W(0), DW(0), DDW(0), J(0), DJ(0), Dat(Dat_) {}

	Material::Gas::EoS::LocalState * GetLS()
	{
		if (LS == 0) LS = new Material::Gas::EoS::FKQuadratic::LocalState(Dat);
		return LS;
	}

	Material::Gas::EoS::Energy<0> * GetW()
	{
		if (LS == 0) LS = new Material::Gas::EoS::FKQuadratic::LocalState(Dat);
		if (W == 0) W = new Material::Gas::EoS::FKQuadratic::Energy<0>(LS);
		return W;
	}

	Material::Gas::EoS::Energy<1> * GetDW()
	{
		if (LS == 0) LS = new Material::Gas::EoS::FKQuadratic::LocalState(Dat);
		if (DW == 0) DW = new Material::Gas::EoS::FKQuadratic::Energy<1>(LS);
		return DW;
	}

	Material::Gas::EoS::Energy<2> * GetDDW()
	{
		if (LS == 0) LS = new Material::Gas::EoS::FKQuadratic::LocalState(Dat);
		if (DDW == 0) DDW = new Material::Gas::EoS::FKQuadratic::Energy<2>(LS);
		return DDW;
	}

	Material::Gas::EoS::Jet<0> * GetJ()
	{
		if (LS == 0) LS = new Material::Gas::EoS::FKQuadratic::LocalState(Dat);
		if (J == 0) J = new Material::Gas::EoS::FKQuadratic::Jet<0>(LS);
		return J;
	}

	Material::Gas::EoS::Jet<1> * GetDJ()
	{
		if (LS == 0) LS = new Material::Gas::EoS::FKQuadratic::LocalState(Dat);
		if (DJ == 0) DJ = new Material::Gas::EoS::FKQuadratic::Jet<1>(LS);
		return DJ;
	}

private:

	Material::Gas::EoS::FKQuadratic::Data *Dat;
	Material::Gas::EoS::FKQuadratic::LocalState *LS;
	Material::Gas::EoS::FKQuadratic::Energy<0> *W;
	Material::Gas::EoS::FKQuadratic::Energy<1> *DW;
	Material::Gas::EoS::FKQuadratic::Energy<2> *DDW;
	Material::Gas::EoS::FKQuadratic::Jet<0> *J;
	Material::Gas::EoS::FKQuadratic::Jet<1> *DJ;

private:

	Factory(const Factory &);
	Factory & operator = (const Factory &);
};

//////////////////////////////////////////////////////////////////////
// Class Builder
//////////////////////////////////////////////////////////////////////

class Builder : public Material::Gas::EoS::Builder
{
public: 

	typedef Material::Gas::EoS::FKQuadratic::Data data_type;

	Builder() {}
	virtual ~Builder() {}

	Builder(Material::Gas::EoS::FKQuadratic::Data *Dat_) : Dat(Dat_) {}

	Material::Gas::EoS::Factory * Build() const
	{
		return new Material::Gas::EoS::FKQuadratic::Factory(Dat);
	}

private:

	Material::Gas::EoS::FKQuadratic::Data *Dat;
};

}

}

}

}

#endif // !defined(MATERIAL_GAS_EOS_FKQUADRATIC_FACTORY_H__INCLUDED_
