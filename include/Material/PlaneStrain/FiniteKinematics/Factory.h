// Factory.h: Factory for the PlaneStrain::FiniteKinematics class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(MATERIAL_PLANESTRAIN_FINITEKINEMATICS_FACTORY_H__INCLUDED_)
#define MATERIAL_PLANESTRAIN_FINITEKINEMATICS_FACTORY_H__INCLUDED_

#pragma once

#include <vector>
#include "./FiniteKinematics.h"
#include "../../Material.h"
#include "../../Factory.h"

namespace Material
{
namespace PlaneStrain
{
namespace FiniteKinematics
{
//////////////////////////////////////////////////////////////////////
// Class Factory
//////////////////////////////////////////////////////////////////////

class Factory : public Material::Factory
{
public: 

	Factory() : W(0), DW(0), DDW(0), J(0), DJ(0) {}

	virtual ~Factory() 
	{
		if (DJ     != 0) delete DJ;
		if (J      != 0) delete J;
		if (DDW    != 0) delete DDW;
		if (DW     != 0) delete DW;
		if (W      != 0) delete W;

		delete MF;
	}

	Factory(Material::Builder *MB) : MF(MB->Build()),
		W(0), DW(0), DDW(0), J(0), DJ(0) {}


	Material::Energy<0> * GetW()
	{
		if (W == 0) W = 
			new Material::PlaneStrain::FiniteKinematics::Energy<0>(MF->GetW());
		return W;
	}

	Material::Energy<1> * GetDW()
	{
		if (DW == 0) DW = 
			new Material::PlaneStrain::FiniteKinematics::Energy<1>(MF->GetDW());
		return DW;
	}

	Material::Energy<2> * GetDDW()
	{
		if (DDW == 0) DDW = 
			new Material::PlaneStrain::FiniteKinematics::Energy<2>(MF->GetDDW());
		return DDW;
	}

	Material::Jet<0> * GetJ()
	{
		if (J == 0) J = 
			new Material::PlaneStrain::FiniteKinematics::Jet<0>(MF->GetJ());
		return J;
	}

	Material::Jet<1> * GetDJ()
	{
		if (DJ == 0) DJ = 
			new Material::PlaneStrain::FiniteKinematics::Jet<1>(MF->GetDJ());
		return DJ;
	}

private:

	Material::Factory *MF;
	Material::PlaneStrain::FiniteKinematics::Energy<0> *W;
	Material::PlaneStrain::FiniteKinematics::Energy<1> *DW;
	Material::PlaneStrain::FiniteKinematics::Energy<2> *DDW;
	Material::PlaneStrain::FiniteKinematics::Jet<0> *J;
	Material::PlaneStrain::FiniteKinematics::Jet<1> *DJ;

private:

	Material::LocalState * GetLS(){return 0;}
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

	Builder(Material::Builder *B_) : B(B_) {}

	Material::Factory * Build() const
	{
		return new Material::PlaneStrain::FiniteKinematics::Factory(B);
	}

private:

	Material::Builder *B;
};

}

}

}

#endif // !defined(MATERIAL_PLANESTRAIN_FINITEKINEMATICS_FACTORY_H__INCLUDED_
