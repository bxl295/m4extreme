// Factory.h: Factory for the Gas class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(MATERIAL_GAS_FACTORY_H__INCLUDED_)
#define MATERIAL_GAS_FACTORY_H__INCLUDED_

#pragma once

#include <vector>
#include "Gas.h"
#include "../Factory.h"
#include "./EoS/Factory.h"

namespace Material
{
namespace Gas
{
//////////////////////////////////////////////////////////////////////
// Class Factory
//////////////////////////////////////////////////////////////////////

class Factory : public Material::Factory
{
public:
        typedef Material::Gas::EoS::Builder eos_build_type;
        typedef Material::Gas::EoS::Factory eos_fact_type;

	Factory() : LS(0), W(0), DW(0), DDW(0), J(0), DJ(0), EoSFact(0) {}
	virtual ~Factory() 
	{
	        if (EoSFact != 0) delete EoSFact;
		if (LS != 0)  delete LS;
		if (W != 0)   delete W;
		if (DW != 0)  delete DW;
		if (DDW != 0) delete DDW;
		if (J != 0)   delete J;
		if (DJ != 0)  delete DJ;                
	}

	Factory(eos_build_type *EoSB_) :
		LS(0), W(0), DW(0), DDW(0), J(0), DJ(0), EoSFact(EoSB_->Build()) {}

	Factory(const Factory &rhs) : 
		LS(0), W(0), DW(0), DDW(0), J(0), DJ(0), EoSFact(rhs.EoSFact) {}

	Material::Factory *Clone() const{return new Factory(*this);}

	Material::LocalState * GetLS()
	{
		if (LS == 0) LS = new Material::Gas::LocalState(EoSFact->GetLS());
		return LS;
	}

	Material::Energy<0> * GetW()
	{
		if (LS == 0) LS = new Material::Gas::LocalState(EoSFact->GetLS());
		if (W == 0) W = new Material::Gas::Energy<0>(LS, EoSFact->GetW());
		return W;
	}

	Material::Energy<1> * GetDW()
	{
		if (LS == 0) LS = new Material::Gas::LocalState(EoSFact->GetLS());
		if (DW == 0) DW = new Material::Gas::Energy<1>(LS, EoSFact->GetDW());
		return DW;
	}

	Material::Energy<2> * GetDDW()
	{
		if (LS == 0) LS = new Material::Gas::LocalState(EoSFact->GetLS());
		if (DDW == 0) DDW = new Material::Gas::Energy<2>(LS, EoSFact->GetDW(), EoSFact->GetDDW());
		return DDW;
	}

	Material::Jet<0> * GetJ()
	{
		if (LS == 0) LS = new Material::Gas::LocalState(EoSFact->GetLS());
		if (J == 0) J = new Material::Gas::Jet<0>(LS, EoSFact->GetJ());
		return J;
	}

	Material::Jet<1> * GetDJ()
	{
		if (LS == 0) LS = new Material::Gas::LocalState(EoSFact->GetLS());
		if (DJ == 0) DJ = new Material::Gas::Jet<1>(LS, EoSFact->GetDJ());
		return DJ;
	}

private:

	eos_fact_type *EoSFact;
	Material::Gas::LocalState *LS;
	Material::Gas::Energy<0> *W;
	Material::Gas::Energy<1> *DW;
	Material::Gas::Energy<2> *DDW;
	Material::Gas::Jet<0> *J;
	Material::Gas::Jet<1> *DJ;

private:

	Factory & operator = (const Factory &);
};

//////////////////////////////////////////////////////////////////////
// Class Builder
//////////////////////////////////////////////////////////////////////

class Builder : public Material::Builder
{
public: 

	typedef Material::Gas::EoS::Builder eos_builder_type;

	Builder() {}
	virtual ~Builder() {}

	Builder(eos_builder_type *EoSB_) : EoSB(EoSB_) {}

	Material::Factory * Build() const
	{
		return new Material::Gas::Factory(EoSB);
	}

private:
        eos_builder_type *EoSB;
};

}

}

#endif // !defined(MATERIAL_GAS_FACTORY_H__INCLUDED_
