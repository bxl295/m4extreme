// Factory.h: Factory for the Incremental class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(MATERIAL_INCREMENTAL_FACTORY_H__INCLUDED_)
#define MATERIAL_INCREMENTAL_FACTORY_H__INCLUDED_

#pragma once

#include <vector>
#include "Incremental.h"
#include "../Factory.h"
#include "../Gas/EoS/Factory.h"

namespace Material
{
namespace Incremental
{
//////////////////////////////////////////////////////////////////////
// Class Factory
//////////////////////////////////////////////////////////////////////

class Factory : public Material::Factory
{
public: 

	Factory() : LS(0), W(0), DW(0), DDW(0), J(0), DJ(0) {}

	virtual ~Factory() 
	{
		if (DJ     != 0) delete DJ;
		if (J      != 0) delete J;
		if (DDW    != 0) delete DDW;
		if (DW     != 0) delete DW;
		if (W      != 0) delete W;
		if (LS     != 0) delete LS;
	}

	Factory(Material::Factory *Fact_) : Fact(Fact_),
		LS(0), W(0), DW(0), DDW(0), J(0), DJ(0) {}

	void NewLS()
	{
		LS = new Material::Incremental::LocalState(Fact->GetLS());
	}

	Material::LocalState * GetLS()
	{
		if (LS == 0) NewLS();
		return LS;
	}

	Material::Energy<0> * GetW()
	{
		if (LS == 0) NewLS();
		if (W == 0) W = new Material::Incremental::Energy<0>(LS,Fact->GetW());
		return W;
	}

	Material::Energy<1> * GetDW()
	{
		if (LS == 0) NewLS();
		if (DW == 0) DW = new Material::Incremental::Energy<1>(LS,Fact->GetDW());
		return DW;
	}

	Material::Energy<2> * GetDDW()
	{
		if (LS == 0) NewLS();
		if (DDW == 0) DDW = new Material::Incremental::Energy<2>(LS,Fact->GetDDW());
		return DDW;
	}

	Material::Jet<0> * GetJ()
	{
		if (LS == 0) NewLS();
		if (J == 0) J = new Material::Incremental::Jet<0>(LS,Fact->GetJ());
		return J;
	}

	Material::Jet<1> * GetDJ()
	{
		if (LS == 0) NewLS();
		if (DJ == 0) DJ = new Material::Incremental::Jet<1>(LS,Fact->GetDJ());
		return DJ;
	}

private:

	Material::Factory *Fact;
	Material::Incremental::LocalState *LS;
	Material::Incremental::Energy<0> *W;
	Material::Incremental::Energy<1> *DW;
	Material::Incremental::Energy<2> *DDW;
	Material::Incremental::Jet<0> *J;
	Material::Incremental::Jet<1> *DJ;

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

	Builder() {}
	virtual ~Builder() {}

	Builder(Material::Factory *Fact_) : Fact(Fact_) {}

	Material::Factory * Build() const
	{
		return new Material::Incremental::Factory(Fact);
	}

private:

	Material::Factory *Fact;
};

}

}

#endif // !defined(MATERIAL_INCREMENTAL_FACTORY_H__INCLUDED_
