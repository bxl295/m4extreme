// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
////////////////////////////////////////////////////////////////////////////

#if !defined(MATERIAL_TWOSTAGES_FACTORY_H_INCLUDED_)
#define MATERIAL_TWOSTAGES_FACTORY_H_INCLUDED_

#pragma once

#include <vector>
#include "./TwoStages.h"
#include "Material/Factory.h"

namespace Material
{
namespace TwoStages
{
//////////////////////////////////////////////////////////////////////
// Class Factory
//////////////////////////////////////////////////////////////////////

class Factory : public Material::Factory
{
public: 

        Factory() : Fcomp(0), Ften(0), LS(0), W(0), DW(0), DDW(0), J(0), DJ(0) {}
	virtual ~Factory() 
	{
	        if (Fcomp != 0 ) delete Fcomp;
		if (Ften != 0 ) delete Ften;
		if (LS != 0)  delete LS;
		if (W != 0)   delete W;
		if (DW != 0)  delete DW;
		if (DDW != 0) delete DDW;
		if (J != 0)   delete J;
		if (DJ != 0)  delete DJ;
	}

	Factory(
		Material::Builder * Bcomp_, 
		Material::Builder *Bten_) : 
	  Fcomp(Bcomp_->Build()), Ften(Bten_->Build()),
	  LS(0), W(0), DW(0), DDW(0), J(0), DJ(0) {}

	Material::LocalState * GetLS()
	{
	    if (LS == 0) LS = new Material::TwoStages::LocalState(Fcomp->GetLS(), Ften->GetLS());
	    return LS;
	}

	Material::Energy<0> * GetW()
	{
	    if (LS == 0) LS = new Material::TwoStages::LocalState(Fcomp->GetLS(), Ften->GetLS());
	    if (W == 0) W = new Material::TwoStages::Energy<0>(LS, Fcomp->GetW(), Ften->GetW());
	    return W;
	}

	Material::Energy<1> * GetDW()
	{
	    if (LS == 0) LS = new Material::TwoStages::LocalState(Fcomp->GetLS(), Ften->GetLS());
	    if (DW == 0) DW = new Material::TwoStages::Energy<1>(LS, Fcomp->GetDW(),Ften->GetDW());		
	    return DW;
	}

	Material::Energy<2> * GetDDW()
	{
	    if (LS == 0) LS = new Material::TwoStages::LocalState(Fcomp->GetLS(), Ften->GetLS());
	    if (DDW == 0) DDW = new Material::TwoStages::Energy<2>(LS, Fcomp->GetDDW(), Ften->GetDDW());
	    return DDW;
	}

	Material::Jet<0> * GetJ()
	{
	    if (LS == 0) LS = new Material::TwoStages::LocalState(Fcomp->GetLS(), Ften->GetLS());
	    if (J == 0) J = new Material::TwoStages::Jet<0>(LS, Fcomp->GetJ(), Ften->GetJ());
	    return J;
	}

	Material::Jet<1> * GetDJ()
	{
	    if (LS == 0) LS = new Material::TwoStages::LocalState(Fcomp->GetLS(), Ften->GetLS());
	    if (DJ == 0) DJ = new Material::TwoStages::Jet<1>(LS, Fcomp->GetDJ(),Ften->GetDJ());
	    return DJ;
	}

private:
	Material::Factory *Fcomp;
	Material::Factory *Ften;
	Material::TwoStages::LocalState *LS;
	Material::TwoStages::Energy<0> *W;
	Material::TwoStages::Energy<1> *DW;
	Material::TwoStages::Energy<2> *DDW;
	Material::TwoStages::Jet<0> *J;
	Material::TwoStages::Jet<1> *DJ;

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

	Builder(
		Material::Builder *Bcomp_, 
		Material::Builder *Bten_) : Bcomp(Bcomp_), Bten(Bten_) {}

	Material::Factory * Build() const
	{
	  return new Material::TwoStages::Factory(Bcomp, Bten);
	}

private:
	Material::Builder *Bcomp; 
	Material::Builder *Bten;
};

}

}

#endif // !defined(MATERIAL_TWOSTAGES_FACTORY_H_INCLUDED_)
