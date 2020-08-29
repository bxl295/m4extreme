// Factory.h: Factory for the Exponential class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(MATERIAL_SHEAR_EXPONENTIAL_FACTORY_H__INCLUDED_)
#define MATERIAL_SHEAR_EXPONENTIAL_FACTORY_H__INCLUDED_

#pragma once

#include <vector>
#include "./Exponential.h"
#include "../Factory.h"

namespace Material
{
namespace Shear
{
namespace Exponential
{
//////////////////////////////////////////////////////////////////////
// Class Factory
//////////////////////////////////////////////////////////////////////

class Factory : public Material::Shear::Factory
{
public: 

	typedef Material::Shear::Exponential::Data data_type;

	Factory() : LS(0), G(0), DG(0), DDG(0), J(0), DJ(0), JJ(0) {}
	virtual ~Factory() 
	{
		if (LS  != 0) delete LS;
		if (G   != 0) delete G;
		if (DG  != 0) delete DG;
		if (DDG != 0) delete DDG;
		if (J   != 0) delete J;
		if (DJ  != 0) delete DJ;
		if (JJ  != 0) delete JJ;
	}

	Factory(Material::Shear::Exponential::Data *Dat_) : 
		LS(0), G(0), DG(0), DDG(0), J(0), DJ(0), JJ(0), Dat(Dat_) {}

	Material::Shear::LocalState * GetLS()
	{
		if (LS == 0) LS = new Material::Shear::Exponential::LocalState(Dat);
		return LS;
	}

	Material::Shear::Modulus<0> * GetG()
	{
		if (LS == 0) LS = new Material::Shear::Exponential::LocalState(Dat);
		if (G  == 0) G  = new Material::Shear::Exponential::Modulus<0>(LS);
		return G;
	}

	Material::Shear::Modulus<1> * GetDG()
	{
		if (LS == 0) LS = new Material::Shear::Exponential::LocalState(Dat);
		if (DG == 0) DG = new Material::Shear::Exponential::Modulus<1>(LS);
		return DG;
	}

	Material::Shear::Modulus<2> * GetDDG()
	{
		if (LS  == 0) LS  = new Material::Shear::Exponential::LocalState(Dat);
		if (DDG == 0) DDG = new Material::Shear::Exponential::Modulus<2>(LS);
		return DDG;
	}

	Material::Shear::Jet<0> * GetJ()
	{
		if (LS == 0) LS = new Material::Shear::Exponential::LocalState(Dat);
		if (J  == 0) J  = new Material::Shear::Exponential::Jet<0>(LS);
		return J;
	}

	Material::Shear::Jet<1> * GetDJ()
	{
		if (LS == 0) LS = new Material::Shear::Exponential::LocalState(Dat);
		if (DJ == 0) DJ = new Material::Shear::Exponential::Jet<1>(LS);
		return DJ;
	}

	Material::Shear::JetJet<0> * GetJJ()
	{
		if (LS == 0) LS = new Material::Shear::Exponential::LocalState(Dat);
		if (JJ == 0) JJ = new Material::Shear::Exponential::JetJet<0>(LS);
		return JJ;
	}

private:

	Material::Shear::Exponential::Data *Dat;
	Material::Shear::Exponential::LocalState *LS;
	Material::Shear::Exponential::Modulus<0> *G;
	Material::Shear::Exponential::Modulus<1> *DG;
	Material::Shear::Exponential::Modulus<2> *DDG;
	Material::Shear::Exponential::Jet<0> *J;
	Material::Shear::Exponential::Jet<1> *DJ;
	Material::Shear::Exponential::JetJet<0> *JJ;

private:

	Factory(const Factory &);
	Factory & operator = (const Factory &);
};

//////////////////////////////////////////////////////////////////////
// Class Builder
//////////////////////////////////////////////////////////////////////

class Builder : public Material::Shear::Builder
{
public: 

	typedef Material::Shear::Exponential::Data data_type;

	Builder() {}
	virtual ~Builder() {}

	Builder(Material::Shear::Exponential::Data *Dat_) : Dat(Dat_) {}

	Material::Shear::Factory * Build() const
	{
		return new Material::Shear::Exponential::Factory(Dat);
	}

private:

	Material::Shear::Exponential::Data *Dat;
};

}

}

}

#endif // !defined(MATERIAL_SHEAR_EXPONENTIAL_FACTORY_H__INCLUDED_
