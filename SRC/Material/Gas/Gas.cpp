// Gas.cpp: implementation of the Gas class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include "./Gas.h"
#include "../../Utils/Indexing/Indexing.h"

namespace Material
{
namespace Gas
{
//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

LocalState::~LocalState() {}

Material::LocalState *
LocalState::Clone() const 
{
	return new LocalState(*this);
}

LocalState::LocalState(Material::Gas::EoS::LocalState *EoSLS_) :
	EoSLS(EoSLS_) {}

LocalState::LocalState(const LocalState &rhs) : 
	EoSLS(rhs.EoSLS) {}

void 
LocalState::operator ++ () { return EoSLS->operator++(); }

//////////////////////////////////////////////////////////////////////
// Class Energy<0>
//////////////////////////////////////////////////////////////////////

Energy<0>::~Energy() {}

Material::Energy<0> *
Energy<0>::Clone() const
{
	return new Energy<0>(*this);
}

Energy<0>::Energy(LocalState *LS_,
	Material::Gas::EoS::Energy<0> *W_) : 
	LS(LS_), W(W_) {}

Energy<0>::Energy(const Energy<0> &rhs) : 
	LS(rhs.LS), W(rhs.W){}

double 
Energy<0>::operator () (const Set::VectorSpace::Vector &Dy) const
{
	assert (Dy.size() == 9);
	Set::VectorSpace::Hom F(3,3,Dy.begin());
	return (*W)(Jacobian(F));
}

double 
Energy<0>::operator () (const Set::VectorSpace::Vector &Dy, double T) const
{
	assert (Dy.size() == 9);
	Set::VectorSpace::Hom F(3,3,Dy.begin());
	return (*W)(Jacobian(F), T);
}

//////////////////////////////////////////////////////////////////////
// Class Energy<1>
//////////////////////////////////////////////////////////////////////

Energy<1>::~Energy() {}

Material::Energy<1> *
Energy<1>::Clone() const
{
	return new Energy<1>(*this);
}

Energy<1>::Energy(LocalState *LS_,
	Material::Gas::EoS::Energy<1> *DW_) :
	LS(LS_), DW(DW_) {}

Energy<1>::Energy(const Energy<1> &rhs) : 
	LS(rhs.LS), DW(rhs.DW) {}

Set::VectorSpace::Vector 
Energy<1>::operator () (const Set::VectorSpace::Vector &Dy) const
{
	assert (Dy.size() == 9);
	Set::VectorSpace::Hom F(3,3,Dy.begin());
	double J = Jacobian(F);
	return (*DW)(J)*J*Adjoint(Inverse(F));
}

double
Energy<1>::operator () (const Set::VectorSpace::Vector &Dy, double T) const
{
	assert (Dy.size() == 9);
	return 0.0;
	//Set::VectorSpace::Hom F(3,3,Dy.begin());
	//return (*DW)(Jacobian(F), T);
}

//////////////////////////////////////////////////////////////////////
// Class Energy<2>
//////////////////////////////////////////////////////////////////////

Energy<2>::~Energy() {}

Material::Energy<2> *
Energy<2>::Clone() const
{
	return new Energy<2>(*this);
}

Energy<2>::Energy(LocalState *LS_,
	Material::Gas::EoS::Energy<1> *DW_, 
	Material::Gas::EoS::Energy<2> *DDW_)
	: LS(LS_), DW(DW_), DDW(DDW_) {}

Energy<2>::Energy(const Energy<2> &rhs) 
	: LS(rhs.LS), DW(rhs.DW), DDW(rhs.DDW){}

Set::VectorSpace::Hom 
Energy<2>::operator () (const Set::VectorSpace::Vector &Dy) const
{
	assert (Dy.size() == 9);
	Set::VectorSpace::Hom F(3,3,Dy.begin());
	Set::VectorSpace::Hom FInv = Inverse(F);
	double J = Jacobian(F);
	double dw = (*DW)(J), ddw = (*DDW)(J);
	double Lambda = (ddw*J + dw)*J, Mu = - dw*J;

	Set::VectorSpace::Hom DDE(9);
	double ****C = Indexing::New(DDE.begin(),3,3,3,3);

	unsigned int i, j, k, l;
	double factor;

	for (i=0; i<3; i++)
		for (j=0; j<3; j++) 
		{
			factor = Lambda*FInv[i][j];
			for (k=0; k<3; k++)
				for (l=0; l<3; l++) 
					C[l][k][j][i] = factor*FInv[k][l];
		}

	for (k=0; k<3; k++) 
		for (j=0; j<3; j++)
		{
			factor = Mu*FInv[k][j];
			for (i=0; i<3; i++) 
				for (l=0; l<3; l++)
					C[l][k][j][i] += factor*FInv[i][l];
		}

	Indexing::Delete(C,3,3,3,3);

	return DDE;
}

double
Energy<2>::operator () (const Set::VectorSpace::Vector &Dy, double T) const
{
	assert (Dy.size() == 9);
	return 0.0;
	//Set::VectorSpace::Hom F(3,3,Dy.begin());
	//return (*DDW)(Jacobian(F), T);
}

//////////////////////////////////////////////////////////////////////
// Class Jet<0>
//////////////////////////////////////////////////////////////////////

Jet<0>::~Jet() {}

Material::Jet<0> *
Jet<0>::Clone() const
{
	return new Jet<0>(*this);
}

Jet<0>::Jet(LocalState *LS_,
	Material::Gas::EoS::Jet<0> *V_) : 
	LS(LS_), V(V_) {}

Jet<0>::Jet(const Jet<0> &rhs) : LS(rhs.LS), V(rhs.V) {}

pair<double,Set::VectorSpace::Vector>
Jet<0>::operator () (const Set::VectorSpace::Vector &Dy) const
{
	assert (Dy.size() == 9);
	Set::VectorSpace::Hom F(3,3,Dy.begin());
	double J = Jacobian(F);
	pair<double,double> K = (*V)(J);
	return make_pair(K.first,K.second*J*Adjoint(Inverse(F)));
}

//////////////////////////////////////////////////////////////////////
// Class Jet<1>
//////////////////////////////////////////////////////////////////////

Jet<1>::~Jet() {}

Material::Jet<1> *
Jet<1>::Clone() const
{
	return new Jet<1>(*this);
}

Jet<1>::Jet(LocalState *LS_,
	Material::Gas::EoS::Jet<1> *DV_) : 
	LS(LS_), DV(DV_) {}

Jet<1>::Jet(const Jet<1> &rhs) : LS(rhs.LS), DV(rhs.DV) {}

pair<Set::VectorSpace::Vector,Set::VectorSpace::Hom> 
Jet<1>::operator () (const Set::VectorSpace::Vector &Dy) const
{
	assert (Dy.size() == 9);
	Set::VectorSpace::Hom F(3,3,Dy.begin());
	double J = Jacobian(F);
	pair<double,double> DK = (*DV)(J);
	double dw = DK.first, ddw = DK.second;
	double Lambda = (ddw*J + dw)*J, Mu = - dw*J;
	Set::VectorSpace::Hom FInv = Inverse(F);
	Set::VectorSpace::Vector DE = dw*J*Adjoint(FInv);

	Set::VectorSpace::Hom DDE(9);
	double ****C = Indexing::New(DDE.begin(),3,3,3,3);

	unsigned int i, j, k, l;
	double factor;

	for (i=0; i<3; i++)
		for (j=0; j<3; j++) 
		{
			factor = Lambda*FInv[i][j];
			for (k=0; k<3; k++)
				for (l=0; l<3; l++) 
					C[l][k][j][i] = factor*FInv[k][l];
		}

	for (k=0; k<3; k++) 
		for (j=0; j<3; j++)
		{
			factor = Mu*FInv[k][j];
			for (i=0; i<3; i++) 
				for (l=0; l<3; l++)
					C[l][k][j][i] += factor*FInv[i][l];
		}

	Indexing::Delete(C,3,3,3,3);

	return make_pair(DE,DDE);
}

}

}
