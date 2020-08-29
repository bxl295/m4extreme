// Factory.h: Factory for the Constant class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(MATERIAL_SHEAR_CONSTANT_FACTORY_H__INCLUDED_)
#define MATERIAL_SHEAR_CONSTANT_FACTORY_H__INCLUDED_

#pragma once

#include <vector>
#include "./Constant.h"
#include "../Factory.h"

namespace Material
{
namespace Shear
{
namespace Constant
{
//////////////////////////////////////////////////////////////////////
// Class Factory
//////////////////////////////////////////////////////////////////////

class Factory : public Material::Shear::Factory
{
public: 

	typedef Material::Shear::Constant::Data data_type;

	Factory() : LS(0), G(0), DG(0), DDG(0), J(0), DJ(0), JJ(0) {}
	virtual ~Factory() 
	{
		if (LS != 0)  delete LS;
		if (G != 0)   delete G;
		if (DG != 0)  delete DG;
		if (DDG != 0) delete DDG;
		if (J != 0)   delete J;
		if (DJ != 0)  delete DJ;
		if (JJ != 0)  delete JJ;
	}

	Factory(Material::Shear::Constant::Data *Dat_) : 
		LS(0), G(0), DG(0), DDG(0), J(0), DJ(0), JJ(0), Dat(Dat_) {}

	Material::Shear::LocalState * GetLS()
	{
		if (LS == 0) LS = new Material::Shear::Constant::LocalState(Dat);
		return LS;
	}

	Material::Shear::Modulus<0> * GetG()
	{
		if (LS == 0) LS = new Material::Shear::Constant::LocalState(Dat);
		if (G  == 0) G  = new Material::Shear::Constant::Modulus<0>(LS);
		return G;
	}

	Material::Shear::Modulus<1> * GetDG()
	{
		if (LS == 0) LS = new Material::Shear::Constant::LocalState(Dat);
		if (DG == 0) DG = new Material::Shear::Constant::Modulus<1>(LS);
		return DG;
	}

	Material::Shear::Modulus<2> * GetDDG()
	{
		if (LS  == 0) LS  = new Material::Shear::Constant::LocalState(Dat);
		if (DDG == 0) DDG = new Material::Shear::Constant::Modulus<2>(LS);
		return DDG;
	}

	Material::Shear::Jet<0> * GetJ()
	{
		if (LS == 0) LS = new Material::Shear::Constant::LocalState(Dat);
		if (J  == 0) J  = new Material::Shear::Constant::Jet<0>(LS);
		return J;
	}

	Material::Shear::Jet<1> * GetDJ()
	{
		if (LS == 0) LS = new Material::Shear::Constant::LocalState(Dat);
		if (DJ == 0) DJ = new Material::Shear::Constant::Jet<1>(LS);
		return DJ;
	}

	Material::Shear::JetJet<0> * GetJJ()
	{
		if (LS == 0) LS = new Material::Shear::Constant::LocalState(Dat);
		if (JJ == 0) JJ = new Material::Shear::Constant::JetJet<0>(LS);
		return JJ;
	}

private:

	Material::Shear::Constant::Data *Dat;
	Material::Shear::Constant::LocalState *LS;
	Material::Shear::Constant::Modulus<0> *G;
	Material::Shear::Constant::Modulus<1> *DG;
	Material::Shear::Constant::Modulus<2> *DDG;
	Material::Shear::Constant::Jet<0> *J;
	Material::Shear::Constant::Jet<1> *DJ;
	Material::Shear::Constant::JetJet<0> *JJ;

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

	typedef Material::Shear::Constant::Data data_type;

	Builder() {}
	virtual ~Builder() {}

	Builder(Material::Shear::Constant::Data *Dat_) : Dat(Dat_) {}

	Material::Shear::Factory * Build() const
	{
		return new Material::Shear::Constant::Factory(Dat);
	}

private:

	Material::Shear::Constant::Data *Dat;
};

}

}

}

#endif // !defined(MATERIAL_SHEAR_CONSTANT_FACTORY_H__INCLUDED_
