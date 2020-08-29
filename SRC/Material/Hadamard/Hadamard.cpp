// Hadamard.cpp: implementation of the Hadamard class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <cassert>
#include "./Hadamard.h"
#include "../../Utils/Indexing/Indexing.h"

namespace Material
{
namespace Hadamard
{
//////////////////////////////////////////////////////////////////////
// Class Data
//////////////////////////////////////////////////////////////////////

Data::Data() : Mu(0.0) {}

Data::Data(const double &rhs_Mu) : Mu(rhs_Mu) {}

Data::~Data(){}

Data::Data(const Data &rhs) : Mu(rhs.Mu) {}

Data &
Data::operator = (const Data &rhs)
{
	if (this == &rhs) return *this;
	Mu = rhs.Mu; return *this;
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

LocalState::LocalState(Data *Properties_, Material::Gas::LocalState *LS_)
	: Properties(Properties_), LS(LS_) {}

LocalState::LocalState(const LocalState &rhs) : 
	Properties(rhs.Properties), LS(rhs.LS) {}

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

Energy<0>::Energy(LocalState *LS_, 
	Material::Gas::Energy<0> *f_) : LS(LS_), f(f_) {}

Energy<0>::Energy(const Energy<0> &rhs) : LS(rhs.LS), f(rhs.f) {}

double 
Energy<0>::operator () (const Set::VectorSpace::Vector &Dy) const
{
	assert (Dy.size() == 9);
	const double &Mu = LS->Properties->Mu;
	Set::VectorSpace::Hom F(3,3,Dy.begin());
	double J = Jacobian(F);
	Set::VectorSpace::Hom Fdev = F/pow(J,1.0/3.0);
	Set::VectorSpace::Hom Cdev = Adjoint(Fdev)*Fdev;
	double Wdev = 0.5*Mu*(Trace(Cdev) - 3.0);
	if (Wdev < 0.0) Wdev = 0.0; 
	return (*f)(Dy) + Wdev;
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
	Material::Gas::Energy<1> *Df_) : LS(LS_), Df(Df_) {}

Energy<1>::Energy(const Energy<1> &rhs) : LS(rhs.LS), Df(rhs.Df) {}

Set::VectorSpace::Vector 
Energy<1>::operator () (const Set::VectorSpace::Vector &Dy) const
{
	assert (Dy.size() == 9);
	const double &Mu = LS->Properties->Mu;
	Set::VectorSpace::Hom F(3,3,Dy.begin());
	double J = Jacobian(F);
	double fact = Mu/pow(J,2.0/3.0);
	Set::VectorSpace::Hom C = Adjoint(F)*F;
	Set::VectorSpace::Hom Finv = Inverse(F);
	return (*Df)(Dy) + fact*(F - (Trace(C)/3.0)*Adjoint(Finv));
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
	Material::Gas::Energy<2> *DDf_) : LS(LS_), DDf(DDf_) {}

Energy<2>::Energy(const Energy<2> &rhs) : LS(rhs.LS), DDf(rhs.DDf) {}

Set::VectorSpace::Hom 
Energy<2>::operator () (const Set::VectorSpace::Vector &Dy) const
{
	assert (Dy.size() == 9);
	const double &Mu = LS->Properties->Mu;
	Set::VectorSpace::Hom F(3,3,Dy.begin());
	double J = Jacobian(F);
	Set::VectorSpace::Hom Finv = Inverse(F);
	Set::VectorSpace::Hom C = Adjoint(F)*F;

	Set::VectorSpace::Hom DDW = (*DDf)(Dy);
	double **** D = Indexing::New(DDW.begin(),3,3,3,3);

	unsigned int i, j, k, l;
	double fact = Mu/pow(J,2.0/3.0);
	double t = Trace(C)/3.0;

	{
		double fact2 = fact*t, fact3;

		for (k=0; k<3; k++) 
			for (j=0; j<3; j++)
			{
				fact3 = fact2*Finv[k][j];
				for (i=0; i<3; i++) 
					for (l=0; l<3; l++)
						D[l][k][j][i] += fact3*Finv[i][l];
			}
	}

	{
		double fact2 = fact*(2.0/3.0)*t, fact3;

		for (l=0; l<3; l++) 
			for (k=0; k<3; k++)
			{
				fact3 = fact2*Finv[k][l];
				for (j=0; j<3; j++) 
					for (i=0; i<3; i++)
						D[l][k][j][i] += fact3*Finv[i][j];
			}
	}

	{
		double fact2 = (2.0/3.0)*fact, fact3;

		for (l=0; l<3; l++) 
			for (k=0; k<3; k++)
			{
				fact3 = fact2*Finv[k][l];
				for (j=0; j<3; j++) 
					for (i=0; i<3; i++)
						D[l][k][j][i] -= fact3*F[j][i];
			}

		for (l=0; l<3; l++) 
			for (k=0; k<3; k++)
			{
				fact3 = fact2*F[l][k];
				for (j=0; j<3; j++) 
					for (i=0; i<3; i++)
						D[l][k][j][i] -= fact3*Finv[i][j];
			}
	}

	for (k=0; k<3; k++)
		for (i=0; i<3; i++) 
			D[k][i][k][i] += fact;

	Indexing::Delete(D,3,3,3,3);

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

Jet<0>::Jet(LocalState *LS_, 
	Material::Gas::Jet<0> *g_) : LS(LS_), g(g_) {}

Jet<0>::Jet(const Jet<0> &rhs) : LS(rhs.LS), g(rhs.g) {}

pair<double,Set::VectorSpace::Vector>
Jet<0>::operator () (const Set::VectorSpace::Vector &Dy) const
{
	assert (Dy.size() == 9);
	const double &Mu = LS->Properties->Mu;
	Set::VectorSpace::Hom F(3,3,Dy.begin());
	double J = Jacobian(F); 
	double fact = 1.0/pow(J,1.0/3.0);
	Set::VectorSpace::Hom Fdev = fact*F;
	Set::VectorSpace::Hom Cdev = Adjoint(Fdev)*Fdev;
	pair<double, Set::VectorSpace::Vector> K = (*g)(Dy);
	double Wdev = 0.5*Mu*(Trace(Cdev) - 3.0);
	if (Wdev < 0.0) Wdev = 0.0; 
	double W = K.first + Wdev;
	Set::VectorSpace::Hom Finv = Inverse(Fdev);
	Set::VectorSpace::Vector DW = K.second + 
		fact*Mu*(Fdev - (Trace(Cdev)/3.0)*Adjoint(Finv));
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

Jet<1>::Jet(LocalState *LS_, 
	Material::Gas::Jet<1> *Dg_) : LS(LS_), Dg(Dg_) {}

Jet<1>::Jet(const Jet<1> &rhs) : LS(rhs.LS), Dg(rhs.Dg) {}

pair<Set::VectorSpace::Vector,Set::VectorSpace::Hom> 
Jet<1>::operator () (const Set::VectorSpace::Vector &Dy) const
{
	assert (Dy.size() == 9);
	const double &Mu = LS->Properties->Mu;
	Set::VectorSpace::Hom F(3,3,Dy.begin());
	double J = Jacobian(F);
	Set::VectorSpace::Hom C = Adjoint(F)*F;
	Set::VectorSpace::Hom Finv = Inverse(F);
	double fact = Mu/pow(J,2.0/3.0);

	pair<Set::VectorSpace::Vector, Set::VectorSpace::Hom> DK = (*Dg)(Dy);

	Set::VectorSpace::Vector DW = DK.first +
		fact*(F - (Trace(C)/3.0)*Adjoint(Finv));

	Set::VectorSpace::Hom DDW = DK.second;
	double **** D = Indexing::New(DDW.begin(),3,3,3,3);

	unsigned int i, j, k, l;
	double t = Trace(C)/3.0;

	{
		double fact2 = fact*t, fact3;

		for (k=0; k<3; k++) 
			for (j=0; j<3; j++)
			{
				fact3 = fact2*Finv[k][j];
				for (i=0; i<3; i++) 
					for (l=0; l<3; l++)
						D[l][k][j][i] += fact3*Finv[i][l];
			}
	}

	{
		double fact2 = fact*(2.0/3.0)*t, fact3;

		for (l=0; l<3; l++) 
			for (k=0; k<3; k++)
			{
				fact3 = fact2*Finv[k][l];
				for (j=0; j<3; j++) 
					for (i=0; i<3; i++)
						D[l][k][j][i] += fact3*Finv[i][j];
			}
	}

	{
		double fact2 = (2.0/3.0)*fact, fact3;

		for (l=0; l<3; l++) 
			for (k=0; k<3; k++)
			{
				fact3 = fact2*Finv[k][l];
				for (j=0; j<3; j++) 
					for (i=0; i<3; i++)
						D[l][k][j][i] -= fact3*F[j][i];
			}

		for (l=0; l<3; l++) 
			for (k=0; k<3; k++)
			{
				fact3 = fact2*F[l][k];
				for (j=0; j<3; j++) 
					for (i=0; i<3; i++)
						D[l][k][j][i] -= fact3*Finv[i][j];
			}
	}

	for (k=0; k<3; k++)
		for (i=0; i<3; i++) 
			D[k][i][k][i] += fact;

	Indexing::Delete(D,3,3,3,3);

	return make_pair(DW,DDW);
}

}

}
