// Isotropic.cpp: implementation of the Isotropic class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <cassert>
#include "./Isotropic.h"
#include "../../../Utils/Indexing/Indexing.h"

namespace Material
{
namespace Hookean
{
namespace Isotropic
{
//////////////////////////////////////////////////////////////////////
// Class Data
//////////////////////////////////////////////////////////////////////

Data::Data() : Lambda(0.0), Mu(0.0) {}

Data::~Data(){}

Data::Data(const double * const rhs) 
	: Lambda(rhs[0]), Mu(rhs[1]) {}

Data::Data(const double &rhs_Lambda, const double &rhs_Mu) 
: Lambda(rhs_Lambda), Mu(rhs_Mu) {}

Data::Data(const Data &rhs) 
: Lambda(rhs.Lambda), Mu(rhs.Mu) {}

Data &
Data::operator = (const Data &rhs)
{
	if (this == &rhs) return *this;
	Lambda = rhs.Lambda; Mu = rhs.Mu;
	return *this;
}

const double &
Data::GetLambda() const
{
	return Lambda;
}

double & 
Data::GetLambda()
{
	return Lambda;
}

const double &
Data::GetMu() const
{
	return Mu;
}

double & 
Data::GetMu()
{
	return Mu;
}

void 
Data::Randomize()
{
	Lambda = (double)rand()/(double)RAND_MAX;
	Mu     = (double)rand()/(double)RAND_MAX;
}

//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

LocalState::~LocalState() {}

Material::LocalState *
LocalState::Clone() const 
{
	return new LocalState(*this);
}

LocalState::LocalState(Data *Properties_) : 
	Properties(Properties_) {}

LocalState::LocalState(const LocalState &rhs) :
	Properties(rhs.Properties) {}

void 
LocalState::operator ++ () {}

//////////////////////////////////////////////////////////////////////
// Class Energy<0>
//////////////////////////////////////////////////////////////////////

Energy<0>::~Energy() {}

Material::Energy<0> *
Energy<0>::Clone() const
{
	return new Energy<0>(*this);
}

Energy<0>::Energy(LocalState *LS_) : LS(LS_) {}

Energy<0>::Energy(const Energy<0> &rhs) : LS(rhs.LS) {}

double 
Energy<0>::operator () (const Set::VectorSpace::Vector &Du) const
{
	assert (Du.size() == 9);
	const double &Lambda = LS->Properties->Lambda;
	const double &Mu = LS->Properties->Mu;
	Set::VectorSpace::Hom B(3,3,Du.begin());
	Set::VectorSpace::Hom Eps = 0.5*(B + Adjoint(B));
	double Theta = Trace(Eps);
	double W = 0.5*Lambda*Theta*Theta + Mu*Eps(Eps);
	return W;
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

Energy<1>::Energy(LocalState *LS_) : LS(LS_) {}

Energy<1>::Energy(const Energy<1> &rhs) : LS(rhs.LS) {}

Set::VectorSpace::Vector 
Energy<1>::operator () (const Set::VectorSpace::Vector &Du) const
{
	assert (Du.size() == 9);
	const double &Lambda = LS->Properties->Lambda;
	const double &Mu = LS->Properties->Mu;
	Set::VectorSpace::Hom B(3,3,Du.begin());
	Set::VectorSpace::Hom Eps = 0.5*(B + Adjoint(B));
	Set::VectorSpace::Hom Sig(3,3);
	double Theta = Trace(Eps);
	Sig = 2.0*Mu*Eps; 
	Sig[0][0] += Lambda*Theta;
	Sig[1][1] += Lambda*Theta;
	Sig[2][2] += Lambda*Theta;
	return Sig;
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

Energy<2>::Energy(LocalState *LS_) : LS(LS_) {}

Energy<2>::Energy(const Energy<2> &rhs) : LS(rhs.LS) {}

Set::VectorSpace::Hom 
Energy<2>::operator () (const Set::VectorSpace::Vector &Du) const
{
	assert (Du.size() == 9);
	unsigned int i, j; 
	const double &Lambda = LS->Properties->Lambda;
	const double &Mu = LS->Properties->Mu;
	Set::VectorSpace::Hom DDW(9);
	double **** C = Indexing::New(DDW.begin(),3,3,3,3);
	for (j=0; j<3; j++)
	{
		for (i=0; i<3; i++)
		{
			C[j][j][i][i] += Lambda;
			C[j][i][j][i] += Mu;
			C[i][j][j][i] += Mu;
		}
	}
	Indexing::Delete(C,3,3,3,3);
	return DDW;
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

Jet<0>::Jet(LocalState *LS_) : LS(LS_) {}

Jet<0>::Jet(const Jet<0> &rhs) : LS(rhs.LS) {}

pair<double,Set::VectorSpace::Vector>
Jet<0>::operator () (const Set::VectorSpace::Vector &Du) const
{
	assert (Du.size() == 9);
	const double &Lambda = LS->Properties->Lambda;
	const double &Mu = LS->Properties->Mu;
	Set::VectorSpace::Hom B(3,3,Du.begin());
	Set::VectorSpace::Hom Eps = 0.5*(B + Adjoint(B));
	double Theta = Trace(Eps);
	double W = 0.5*Lambda*Theta*Theta + Mu*Eps(Eps);
	Set::VectorSpace::Vector DW(9);
	Set::VectorSpace::Hom Sig(3,3,DW.begin());
	Sig = 2.0*Mu*Eps; 
	Sig[0][0] += Lambda*Theta;
	Sig[1][1] += Lambda*Theta;
	Sig[2][2] += Lambda*Theta;
	return make_pair(W,DW);
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

Jet<1>::Jet(LocalState *LS_) : LS(LS_) {}

Jet<1>::Jet(const Jet<1> &rhs) : LS(rhs.LS) {}

pair<Set::VectorSpace::Vector,Set::VectorSpace::Hom> 
Jet<1>::operator () (const Set::VectorSpace::Vector &Du) const
{
	assert (Du.size() == 9);
	const double &Lambda = LS->Properties->Lambda;
	const double &Mu = LS->Properties->Mu;
	Set::VectorSpace::Hom B(3,3,Du.begin());
	Set::VectorSpace::Hom Eps = 0.5*(B + Adjoint(B));
	double Theta = Trace(Eps);
	Set::VectorSpace::Vector DW(9);
	Set::VectorSpace::Hom Sig(3,3,DW.begin());
	Sig = 2.0*Mu*Eps; 
	Sig[0][0] += Lambda*Theta;
	Sig[1][1] += Lambda*Theta;
	Sig[2][2] += Lambda*Theta;
	unsigned int i, j; 
	Set::VectorSpace::Hom DDW(9);
	double **** C = Indexing::New(DDW.begin(),3,3,3,3);
	for (j=0; j<3; j++)
	{
		for (i=0; i<3; i++)
		{
			C[j][j][i][i] += Lambda;
			C[j][i][j][i] += Mu;
			C[i][j][j][i] += Mu;
		}
	}
	Indexing::Delete(C,3,3,3,3);
	return make_pair(DW,DDW);
}

}

}

}
