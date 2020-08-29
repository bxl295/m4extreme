// Factory.h: Factory for the LogLaw class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(MATERIAL_UNIAXIAL_LOGLAW_FACTORY_H__INCLUDED_)
#define MATERIAL_UNIAXIAL_LOGLAW_FACTORY_H__INCLUDED_

#pragma once

#include <vector>
#include "LogLaw.h"
#include "Material/Uniaxial/Factory.h"

namespace Material
{
namespace Uniaxial
{
namespace LogLaw
{
//////////////////////////////////////////////////////////////////////
// Class Factory
//////////////////////////////////////////////////////////////////////

class Factory : public Material::Uniaxial::Factory
{
public: 

	typedef Material::Uniaxial::LogLaw::Data data_type;

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

        Factory(double T0_, Material::Uniaxial::LogLaw::Data *Dat_) :
	T0(T0_), LS(0), W(0), DW(0), DDW(0), J(0), DJ(0), Dat(Dat_) {}

	void NewLS() {
	  if (LS == 0) 
	    LS = new Material::Uniaxial::LogLaw::LocalState(T0, Dat);
	}

	Material::Uniaxial::LocalState * GetLS()
	{
	  NewLS();
		return LS;
	}

	Material::Uniaxial::Energy<0> * GetW()
	{
	  NewLS();
		if (W == 0) W = new Material::Uniaxial::LogLaw::Energy<0>(LS);
		return W;
	}

	Material::Uniaxial::Energy<1> * GetDW()
	{
	  NewLS();
		if (DW == 0) DW = new Material::Uniaxial::LogLaw::Energy<1>(LS);
		return DW;
	}

	Material::Uniaxial::Energy<2> * GetDDW()
	{
	  NewLS();
		if (DDW == 0) DDW = new Material::Uniaxial::LogLaw::Energy<2>(LS);
		return DDW;
	}

	Material::Uniaxial::Jet<0> * GetJ()
	{
	  NewLS();
		if (J == 0) J = new Material::Uniaxial::LogLaw::Jet<0>(LS);
		return J;
	}

	Material::Uniaxial::Jet<1> * GetDJ()
	{
	  NewLS();
		if (DJ == 0) DJ = new Material::Uniaxial::LogLaw::Jet<1>(LS);
		return DJ;
	}

private:
	double T0;
	Material::Uniaxial::LogLaw::Data *Dat;
	Material::Uniaxial::LogLaw::LocalState *LS;
	Material::Uniaxial::LogLaw::Energy<0> *W;
	Material::Uniaxial::LogLaw::Energy<1> *DW;
	Material::Uniaxial::LogLaw::Energy<2> *DDW;
	Material::Uniaxial::LogLaw::Jet<0> *J;
	Material::Uniaxial::LogLaw::Jet<1> *DJ;

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

	typedef Material::Uniaxial::LogLaw::Data data_type;

	Builder() {}
	virtual ~Builder() {}

        Builder(double T0_, Material::Uniaxial::LogLaw::Data *Dat_) : Dat(Dat_), T0(T0_) {}

	Material::Uniaxial::Factory * Build() const
	{
	  return new Material::Uniaxial::LogLaw::Factory(T0, Dat);
	}

private:
	double T0;
	Material::Uniaxial::LogLaw::Data *Dat;
};

}

}

}

#endif // !defined(MATERIAL_UNIAXIAL_LOGLAW_FACTORY_H__INCLUDED_
