// Factory.h: Factory for the Hadamard class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(MATERIAL_HADAMARD_FACTORY_H__INCLUDED_)
#define MATERIAL_HADAMARD_FACTORY_H__INCLUDED_

#pragma once

#include <vector>
#include "Hadamard.h"
#include "../Factory.h"
#include "../Gas/EoS/Factory.h"

namespace Material
{
namespace Hadamard
{
//////////////////////////////////////////////////////////////////////
// Class Factory
//////////////////////////////////////////////////////////////////////

class Factory : public Material::Factory
{
public: 

	Factory() : 
		GasLS(0), GasW(0), GasDW(0), GasDDW(0), GasJ(0), GasDJ(0),
		LS(0), W(0), DW(0), DDW(0), J(0), DJ(0) {}

	virtual ~Factory() 
	{
		if (DJ     != 0) delete DJ;
		if (J      != 0) delete J;
		if (DDW    != 0) delete DDW;
		if (DW     != 0) delete DW;
		if (W      != 0) delete W;
		if (LS     != 0) delete LS;
		if (GasDJ  != 0) delete GasDJ;
		if (GasJ   != 0) delete GasJ;
		if (GasDDW != 0) delete GasDDW;
		if (GasDW  != 0) delete GasDW;
		if (GasW   != 0) delete GasW;
		if (GasLS  != 0) delete GasLS;

		delete EoSFact;
	}

	Factory(Material::Hadamard::Data *HadDat_,
		Material::Gas::EoS::Builder *EoSB) : EoSFact(EoSB->Build()),
		GasLS(0), GasW(0), GasDW(0), GasDDW(0), GasJ(0), GasDJ(0),
		LS(0), W(0), DW(0), DDW(0), J(0), DJ(0), HadDat(HadDat_) {}

	void NewLS()
	{
		GasLS = new Material::Gas::LocalState(EoSFact->GetLS());
		LS = new Material::Hadamard::LocalState(HadDat,GasLS);
	}

	Material::LocalState * GetLS()
	{
		if (LS == 0) NewLS();
		return LS;
	}

	Material::Energy<0> * GetW()
	{
		if (LS == 0) NewLS();
		if (GasW == 0) GasW  = new Material::Gas::Energy<0>(GasLS,EoSFact->GetW()); 
		if (W    == 0) W     = new Material::Hadamard::Energy<0>(LS,GasW);
		return W;
	}

	Material::Energy<1> * GetDW()
	{
		if (LS == 0) NewLS();
		if (GasDW == 0) GasDW  = new Material::Gas::Energy<1>(GasLS,EoSFact->GetDW()); 
		if (DW    == 0) DW     = new Material::Hadamard::Energy<1>(LS,GasDW);
		return DW;
	}

	Material::Energy<2> * GetDDW()
	{
		if (LS == 0) NewLS();
		if (GasDDW == 0) GasDDW = new Material::Gas::Energy<2>(GasLS,EoSFact->GetDW(),EoSFact->GetDDW()); 
		if (DDW    == 0) DDW =   new Material::Hadamard::Energy<2>(LS,GasDDW);
		return DDW;
	}

	Material::Jet<0> * GetJ()
	{
		if (LS == 0) NewLS();
		if (GasJ == 0) GasJ  = new Material::Gas::Jet<0>(GasLS,EoSFact->GetJ()); 
		if (J    == 0) J     = new Material::Hadamard::Jet<0>(LS,GasJ);
		return J;
	}

	Material::Jet<1> * GetDJ()
	{
		if (LS == 0) NewLS();
		if (GasDJ == 0) GasDJ  = new Material::Gas::Jet<1>(GasLS,EoSFact->GetDJ()); 
		if (DJ    == 0) DJ     = new Material::Hadamard::Jet<1>(LS,GasDJ);
		return DJ;
	}

private:

	Material::Hadamard::Data *HadDat;
	Material::Gas::EoS::Factory *EoSFact;
	Material::Gas::LocalState *GasLS;
	Material::Gas::Energy<0> *GasW;
	Material::Gas::Energy<1> *GasDW;
	Material::Gas::Energy<2> *GasDDW;
	Material::Gas::Jet<0> *GasJ;
	Material::Gas::Jet<1> *GasDJ;
	Material::Hadamard::LocalState *LS;
	Material::Hadamard::Energy<0> *W;
	Material::Hadamard::Energy<1> *DW;
	Material::Hadamard::Energy<2> *DDW;
	Material::Hadamard::Jet<0> *J;
	Material::Hadamard::Jet<1> *DJ;

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

	typedef Material::Hadamard::Data data_type;

	Builder() {}
	virtual ~Builder() {}

	Builder(Material::Hadamard::Data *Dat_,
		Material::Gas::EoS::Builder *EoSB_) : 
		Dat(Dat_), EoSB(EoSB_) {}

	Material::Factory * Build() const
	{
		return new Material::Hadamard::Factory(Dat,EoSB);
	}

private:

	Material::Hadamard::Data *Dat;
	Material::Gas::EoS::Builder *EoSB;
};

}

}

#endif // !defined(MATERIAL_HADAMARD_FACTORY_H__INCLUDED_
