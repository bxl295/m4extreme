// FiniteKinematics.cpp: implementation of the FiniteKinematics class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include "./FiniteKinematics.h"
#include "Utils/Indexing/Indexing.h"

namespace Material
{
namespace Deviatoric
{
namespace FiniteKinematics
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

LocalState::LocalState(Material::LocalState *MLS_) :
	MLS(MLS_) {}

LocalState::LocalState(const LocalState &rhs) : 
	MLS(rhs.MLS) {}

void 
LocalState::operator ++ () { MLS->operator++();}

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
	Material::Energy<0> *W_) : 
	LS(LS_), W(W_) {}

Energy<0>::Energy(const Energy<0> &rhs) : 
	LS(rhs.LS), W(rhs.W){}

double 
Energy<0>::operator () (const Set::VectorSpace::Vector &Dy) const
{
	assert (Dy.size() == 9);
	Set::VectorSpace::Hom F(3,3,Dy.begin());
	double alpha = 1.0 / pow(Jacobian(F), 1.0/3.0);
	return (*W)(alpha * F);
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
	Material::Energy<1> *DW_) :
	LS(LS_), DW(DW_) {}

Energy<1>::Energy(const Energy<1> &rhs) : 
	LS(rhs.LS), DW(rhs.DW) {}

Set::VectorSpace::Vector 
Energy<1>::operator () (const Set::VectorSpace::Vector &Dy) const
{
	assert (Dy.size() == 9);
	Set::VectorSpace::Hom F(3,3,Dy.begin());
	double alpha = 1.0 / pow(Jacobian(F), 1.0/3.0);
	Set::VectorSpace::Hom Fdev = alpha * F;
	Set::VectorSpace::Hom Pdev(3, 3, (*DW)(Fdev).begin()); 
	Set::VectorSpace::Hom tau = Pdev * Adjoint(Fdev);
	double p = Trace(tau) / 3.0;
	double * thead = tau.begin();
	*(thead) -= p;
	*(thead + 4) -= p;
	*(thead + 8) -= p;
	return tau * Adjoint(Inverse(F));
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
	Material::Energy<1> *DW_, 
	Material::Energy<2> *DDW_)
	: LS(LS_), DW(DW_), DDW(DDW_) {}

Energy<2>::Energy(const Energy<2> &rhs) 
	: LS(rhs.LS), DW(rhs.DW), DDW(rhs.DDW){}

Set::VectorSpace::Hom 
Energy<2>::operator () (const Set::VectorSpace::Vector &Dy) const
{
	Set::VectorSpace::Hom DDE(9);
	cerr << "Material::Deviatoric::FiniteKinematics::Energy<2> underconstruction" << endl;
	assert(false);
	return DDE;
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
	Material::Jet<0> *V_) : 
	LS(LS_), V(V_) {}

Jet<0>::Jet(const Jet<0> &rhs) : LS(rhs.LS), V(rhs.V) {}

pair<double,Set::VectorSpace::Vector>
Jet<0>::operator () (const Set::VectorSpace::Vector &Dy) const
{
	assert (Dy.size() == 9);
	Set::VectorSpace::Hom F(3,3,Dy.begin());
	double alpha = 1.0 / pow(Jacobian(F), 1.0/3.0);
	Set::VectorSpace::Hom Fdev = alpha * F;
	pair<double, Set::VectorSpace::Vector> K = (*V)(Fdev);
	Set::VectorSpace::Hom Pdev(3, 3, K.second.begin());
	Set::VectorSpace::Hom tau = Pdev * Adjoint(Fdev);
	double p = Trace(tau) / 3.0;
	double * thead = tau.begin();
	*(thead) -= p;
	*(thead + 4) -= p;
	*(thead + 8) -= p;
	return make_pair(K.first, tau * Adjoint(Inverse(F)));
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
	Material::Jet<1> *DV_) : 
	LS(LS_), DV(DV_) {}

Jet<1>::Jet(const Jet<1> &rhs) : LS(rhs.LS), DV(rhs.DV) {}

pair<Set::VectorSpace::Vector,Set::VectorSpace::Hom> 
Jet<1>::operator () (const Set::VectorSpace::Vector &Dy) const
{
        Set::VectorSpace::Vector DE(9);
	Set::VectorSpace::Hom DDE(9);
	cerr << "Material::Deviatoric::FiniteKinematics::Jet<1> underconstruction" << endl;
	assert(false);
	return make_pair(DE,DDE);
}

}

}

}
