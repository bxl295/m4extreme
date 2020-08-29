// Factory.h: Factory for the FiniteKinematics class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(MATERIAL_DEVIATORIC_FINITEKINEMATICS_FACTORY_H__INCLUDED_)
#define MATERIAL_DEVIATORIC_FINITEKINEMATICS_FACTORY_H__INCLUDED_

#pragma once

#include <vector>
#include "FiniteKinematics.h"
#include "Material/Factory.h"

namespace Material
{
namespace Deviatoric
{
namespace FiniteKinematics
{
//////////////////////////////////////////////////////////////////////
// Class Factory
//////////////////////////////////////////////////////////////////////

class Factory : public Material::Factory
{
public:
        typedef Material::Builder mat_build_type;
        typedef Material::Factory mat_fact_type;

	Factory() : LS(0), W(0), DW(0), DDW(0), J(0), DJ(0), MatFact(0) {}
	virtual ~Factory() 
	{
	        if (MatFact != 0) delete MatFact;
		if (LS != 0)  delete LS;
		if (W != 0)   delete W;
		if (DW != 0)  delete DW;
		if (DDW != 0) delete DDW;
		if (J != 0)   delete J;
		if (DJ != 0)  delete DJ;                
	}

	Factory(mat_build_type *MatB_) :
		LS(0), W(0), DW(0), DDW(0), J(0), DJ(0), MatFact(MatB_->Build()) {}

	Factory(const Factory &rhs) : 
		LS(0), W(0), DW(0), DDW(0), J(0), DJ(0), MatFact(rhs.MatFact) {}

	Material::Factory *Clone() const{return new Factory(*this);}

	Material::LocalState * GetLS()
	{
		if (LS == 0) LS = new Material::Deviatoric::FiniteKinematics::LocalState(MatFact->GetLS());
		return LS;
	}

	Material::Energy<0> * GetW()
	{
		if (LS == 0) LS = new Material::Deviatoric::FiniteKinematics::LocalState(MatFact->GetLS());
		if (W == 0) W = new Material::Deviatoric::FiniteKinematics::Energy<0>(LS, MatFact->GetW());
		return W;
	}

	Material::Energy<1> * GetDW()
	{
		if (LS == 0) LS = new Material::Deviatoric::FiniteKinematics::LocalState(MatFact->GetLS());
		if (DW == 0) DW = new Material::Deviatoric::FiniteKinematics::Energy<1>(LS, MatFact->GetDW());
		return DW;
	}

	Material::Energy<2> * GetDDW()
	{
		if (LS == 0) LS = new Material::Deviatoric::FiniteKinematics::LocalState(MatFact->GetLS());
		if (DDW == 0) DDW = new Material::Deviatoric::FiniteKinematics::Energy<2>(LS, MatFact->GetDW(), MatFact->GetDDW());
		return DDW;
	}

	Material::Jet<0> * GetJ()
	{
		if (LS == 0) LS = new Material::Deviatoric::FiniteKinematics::LocalState(MatFact->GetLS());
		if (J == 0) J = new Material::Deviatoric::FiniteKinematics::Jet<0>(LS, MatFact->GetJ());
		return J;
	}

	Material::Jet<1> * GetDJ()
	{
		if (LS == 0) LS = new Material::Deviatoric::FiniteKinematics::LocalState(MatFact->GetLS());
		if (DJ == 0) DJ = new Material::Deviatoric::FiniteKinematics::Jet<1>(LS, MatFact->GetDJ());
		return DJ;
	}

private:

	mat_fact_type *MatFact;
	Material::Deviatoric::FiniteKinematics::LocalState *LS;
	Material::Deviatoric::FiniteKinematics::Energy<0> *W;
	Material::Deviatoric::FiniteKinematics::Energy<1> *DW;
	Material::Deviatoric::FiniteKinematics::Energy<2> *DDW;
	Material::Deviatoric::FiniteKinematics::Jet<0> *J;
	Material::Deviatoric::FiniteKinematics::Jet<1> *DJ;

private:

	Factory & operator = (const Factory &);
};

//////////////////////////////////////////////////////////////////////
// Class Builder
//////////////////////////////////////////////////////////////////////

class Builder : public Material::Builder
{
public: 

	typedef Material::Builder mat_builder_type;

	Builder() {}
	virtual ~Builder() {}

	Builder(mat_builder_type *MatB_) : MatB(MatB_) {}

	Material::Factory * Build() const
	{
		return new Material::Deviatoric::FiniteKinematics::Factory(MatB);
	}

private:
        mat_builder_type *MatB;
};

}

}

}

#endif // !defined(MATERIAL_DEVIATORIC_FINITEKINEMATICS_FACTORY_H__INCLUDED_
