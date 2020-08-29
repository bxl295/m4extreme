// Factory.h: interface for the Factory class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(ELEMENT_ISOPARAMETRIC_FACTORY_H__INCLUDED_)
#define ELEMENT_ISOPARAMETRIC_FACTORY_H__INCLUDED_

#pragma once

#include <set>
#include <map>

#include "./Isoparametric.h"	
#include "./Source.h"	
#include "../Factory.h"	
#include "../Element.h"
#include "../LumpedMass.h"
#include "../Interpolation/Interpolation.h"
#include "../Quadrature/Quadrature.h"
#include "../../Material/Material.h"
#include "../../Material/Factory.h"
#include "../../Set/Manifold/Manifold.h"
#include "../../Set/Algebraic/VectorSpace/Vector/Vector.h"

namespace Element
{
namespace Isoparametric
{
//////////////////////////////////////////////////////////////////////
// Class Factory
//////////////////////////////////////////////////////////////////////

class Factory : public Element::Factory
{
public: 

	Factory() : LS(0), E(0), DE(0), DDE(0), J(0), DJ(0) {}

	virtual ~Factory() 
	{
		if (DJ     != 0) delete DJ;
		if (J      != 0) delete J;
		if (DDE    != 0) delete DDE;
		if (DE     != 0) delete DE;
		if (E      != 0) delete E;
		if (LS     != 0) delete LS;
	}

	Factory(Material::Factory *MF_,
		Element::Interpolation::Shape<1> &DN,
		const Element::Quadrature::Rule &QR,
		const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &Xloc) : 
		MF(vector<Material::Factory *>(QR.GetW().size(),MF_)),
		LS(0), E(0), DE(0), DDE(0), J(0), DJ(0) 
	{
		vector<Material::LocalState *> MatLS(QR.GetW().size(),MF_->GetLS());
		LS = new Element::Isoparametric::LocalState(DN,Xloc,MatLS,QR.GetQ(),QR.GetW());
	}

	Factory(const vector<Material::Factory *> &MF_,
		Element::Interpolation::Shape<1> &DN,
		const Element::Quadrature::Rule &QR,
		const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &Xloc) : 
		MF(MF_), LS(0), E(0), DE(0), DDE(0), J(0), DJ(0) 
	{
		vector<Material::LocalState *> MatLS;
		vector<Material::Factory *>::iterator pMF;
		for (pMF=MF.begin(); pMF!=MF.end(); pMF++) MatLS.push_back((*pMF)->GetLS());
		LS = new Element::Isoparametric::LocalState(DN,Xloc,MatLS,QR.GetQ(),QR.GetW());
	}

	Element::LocalState * GetLS()
	{
		return LS;
	}

	Element::Energy<0> * GetE()
	{
		if (E == 0) 
		{
			vector<Material::Energy<0> *> MatW;
			vector<Material::Factory *>::iterator pMF;
			for (pMF=MF.begin(); pMF!=MF.end(); pMF++) MatW.push_back((*pMF)->GetW());
			E = new Element::Isoparametric::Energy<0>(LS,MatW);
		}
		return E;
	}

	Element::Energy<1> * GetDE()
	{
		if (DE == 0) 
		{
			vector<Material::Energy<1> *> MatDW;
			vector<Material::Factory *>::iterator pMF;
			for (pMF=MF.begin(); pMF!=MF.end(); pMF++) MatDW.push_back((*pMF)->GetDW());
			DE = new Element::Isoparametric::Energy<1>(LS,MatDW);
		}
		return DE;
	}

	Element::Energy<2> * GetDDE()
	{
		if (DDE == 0) 
		{
			vector<Material::Energy<2> *> MatDDW;
			vector<Material::Factory *>::iterator pMF;
			for (pMF=MF.begin(); pMF!=MF.end(); pMF++) MatDDW.push_back((*pMF)->GetDDW());
			DDE = new Element::Isoparametric::Energy<2>(LS,MatDDW);
		}
		return DDE;
	}

	Element::Jet<0> * GetJ()
	{
		if (J == 0) 
		{
			vector<Material::Jet<0> *> MatJ;
			vector<Material::Factory *>::iterator pMF;
			for (pMF=MF.begin(); pMF!=MF.end(); pMF++) MatJ.push_back((*pMF)->GetJ());
			J = new Element::Isoparametric::Jet<0>(LS,MatJ);
		}
		return J;
	}

	Element::Jet<1> * GetDJ()
	{
		if (DJ == 0) 
		{
			vector<Material::Jet<1> *> MatDJ;
			vector<Material::Factory *>::iterator pMF;
			for (pMF=MF.begin(); pMF!=MF.end(); pMF++) MatDJ.push_back((*pMF)->GetDJ());
			DJ = new Element::Isoparametric::Jet<1>(LS,MatDJ);
		}
		return DJ;
	}

private:

	vector<Material::Factory *> MF;
	Element::Isoparametric::LocalState *LS;
	Element::Isoparametric::Energy<0> *E;
	Element::Isoparametric::Energy<1> *DE;
	Element::Isoparametric::Energy<2> *DDE;
	Element::Isoparametric::Jet<0> *J;
	Element::Isoparametric::Jet<1> *DJ;

private:

	Factory(const Factory &);
	Factory & operator = (const Factory &);
};

class LocalBuilder
{
public: 

	LocalBuilder() {}
	virtual ~LocalBuilder() {}

	LocalBuilder(Material::Factory *MF_,
		Element::Interpolation::Shape<1> *DN_, 
		Element::Quadrature::Rule *QR_) : 
		MF(QR_->GetW().size(),MF_), DN(DN_), QR(QR_) {}

	LocalBuilder(const vector<Material::Factory *> &MF_,
		Element::Interpolation::Shape<1> *DN_, 
		Element::Quadrature::Rule *QR_) : MF(MF_), DN(DN_), QR(QR_) {}

	Element::Factory * Build(
		const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &Xloc)
	{
		return new Element::Isoparametric::Factory(MF,*DN,*QR,Xloc);
	}

private:

	vector<Material::Factory *> MF;
	Element::Interpolation::Shape<1> *DN;
	Element::Quadrature::Rule *QR;
};

}

}

namespace Potential
{
namespace Isoparametric
{
//////////////////////////////////////////////////////////////////////
// Class Factory
//////////////////////////////////////////////////////////////////////

class Factory : public Element::Factory
{
public: 

	Factory() : LS(0), E(0), DE(0), DDE(0), J(0), DJ(0) {}

	virtual ~Factory() 
	{
		if (DJ     != 0) delete DJ;
		if (J      != 0) delete J;
		if (DDE    != 0) delete DDE;
		if (DE     != 0) delete DE;
		if (E      != 0) delete E;
		if (LS     != 0) delete LS;
	}

	Factory(Material::Factory *MF_,
		Element::Interpolation::Shape<0> &N,
		Element::Interpolation::Shape<1> &DN,
		const Element::Quadrature::Rule &QR,
		const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &Xloc) : 
		MF(vector<Material::Factory *>(QR.GetW().size(),MF_)),
		LS(0), E(0), DE(0), DDE(0), J(0), DJ(0) 
	{
		vector<Material::LocalState *> MatLS(QR.GetW().size(),MF_->GetLS());
		LS = new Potential::Isoparametric::LocalState(N,DN,Xloc,MatLS,QR.GetQ(),QR.GetW());
	}

	Factory(const vector<Material::Factory *> &MF_,
		Element::Interpolation::Shape<0> &N,
		Element::Interpolation::Shape<1> &DN,
		const Element::Quadrature::Rule &QR,
		const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &Xloc) : 
		MF(MF_), LS(0), E(0), DE(0), DDE(0), J(0), DJ(0) 
	{
		vector<Material::LocalState *> MatLS;
		vector<Material::Factory *>::iterator pMF;
		for (pMF=MF.begin(); pMF!=MF.end(); pMF++) MatLS.push_back((*pMF)->GetLS());
		LS = new Potential::Isoparametric::LocalState(N,DN,Xloc,MatLS,QR.GetQ(),QR.GetW());
	}

	Element::LocalState * GetLS()
	{
		return LS;
	}

	Element::Energy<0> * GetE()
	{
		if (E == 0) 
		{
			vector<Material::Energy<0> *> MatW;
			vector<Material::Factory *>::iterator pMF;
			for (pMF=MF.begin(); pMF!=MF.end(); pMF++) MatW.push_back((*pMF)->GetW());
			E = new Potential::Isoparametric::Energy<0>(LS,MatW);
		}
		return E;
	}

	Element::Energy<1> * GetDE()
	{
		if (DE == 0) 
		{
			vector<Material::Energy<1> *> MatDW;
			vector<Material::Factory *>::iterator pMF;
			for (pMF=MF.begin(); pMF!=MF.end(); pMF++) MatDW.push_back((*pMF)->GetDW());
			DE = new Potential::Isoparametric::Energy<1>(LS,MatDW);
		}
		return DE;
	}

	Element::Energy<2> * GetDDE()
	{
		if (DDE == 0) 
		{
			vector<Material::Energy<2> *> MatDDW;
			vector<Material::Factory *>::iterator pMF;
			for (pMF=MF.begin(); pMF!=MF.end(); pMF++) MatDDW.push_back((*pMF)->GetDDW());
			DDE = new Potential::Isoparametric::Energy<2>(LS,MatDDW);
		}
		return DDE;
	}

	Element::Jet<0> * GetJ()
	{
		if (J == 0) 
		{
			vector<Material::Jet<0> *> MatJ;
			vector<Material::Factory *>::iterator pMF;
			for (pMF=MF.begin(); pMF!=MF.end(); pMF++) MatJ.push_back((*pMF)->GetJ());
			J = new Potential::Isoparametric::Jet<0>(LS,MatJ);
		}
		return J;
	}

	Element::Jet<1> * GetDJ()
	{
		if (DJ == 0) 
		{
			vector<Material::Jet<1> *> MatDJ;
			vector<Material::Factory *>::iterator pMF;
			for (pMF=MF.begin(); pMF!=MF.end(); pMF++) MatDJ.push_back((*pMF)->GetDJ());
			DJ = new Potential::Isoparametric::Jet<1>(LS,MatDJ);
		}
		return DJ;
	}

private:

	vector<Material::Factory *> MF;
	Potential::Isoparametric::LocalState *LS;
	Potential::Isoparametric::Energy<0> *E;
	Potential::Isoparametric::Energy<1> *DE;
	Potential::Isoparametric::Energy<2> *DDE;
	Potential::Isoparametric::Jet<0> *J;
	Potential::Isoparametric::Jet<1> *DJ;

private:

	Factory(const Factory &);
	Factory & operator = (const Factory &);
};

}

}

#endif // !defined(ELEMENT_ISOPARAMETRIC_FACTORY_H__INCLUDED_
