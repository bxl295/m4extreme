// Coulomb.cpp: implementation of the Coulomb class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <cassert>
#include "./Coulomb.h"

namespace Potential
{
namespace Field
{
namespace Coulomb
{
//////////////////////////////////////////////////////////////////////
// Class Data
//////////////////////////////////////////////////////////////////////

Data::Data() {}

Data::Data(const map<Set::Euclidean::Orthonormal::Point *,double> &C_) : C(C_)
{
	map<Set::Euclidean::Orthonormal::Point *,double>::iterator pC;
	for (pC=C.begin(); pC!=C.end(); pC++) assert (pC->first->size() == 3);
}

Data::~Data(){}

Data::Data(const Data &rhs) : C(rhs.C) {}

Data &
Data::operator = (const Data &rhs)
{
	if (this == &rhs) return *this;
	C = rhs.C; return *this;
}

const map<Set::Euclidean::Orthonormal::Point *,double> & 
Data::GetC() const
{
	return C;
}
	
map<Set::Euclidean::Orthonormal::Point *,double> & 
Data::GetC()
{
	return C;
}

void 
Data::Randomize()
{
	map<Set::Euclidean::Orthonormal::Point *,double>::iterator pC;
	for (pC=C.begin(); pC!=C.end(); pC++) 
		pC->second = (double)rand()/(double)RAND_MAX;
}

//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

LocalState::~LocalState() {}

Potential::Field::LocalState *
LocalState::Clone() const 
{
	return new LocalState(*this);
}

LocalState::LocalState(Data *Prop_) : Prop(Prop_) {}

LocalState::LocalState(const LocalState &rhs) : Prop(rhs.Prop) {}

void 
LocalState::operator ++ () {}

//////////////////////////////////////////////////////////////////////
// Class Energy<0>
//////////////////////////////////////////////////////////////////////

Energy<0>::~Energy() {}

Potential::Field::Energy<0> *
Energy<0>::Clone() const
{
	return new Energy<0>(*this);
}

Energy<0>::Energy(LocalState *LS_) : LS(LS_) {}

Energy<0>::Energy(const Energy<0> &rhs) : LS(rhs.LS) {}

double 
Energy<0>::operator () (const Set::VectorSpace::Vector &x) const
{
	double E=0.0;
	const map<Set::Euclidean::Orthonormal::Point *,double> &C = LS->Prop->C;
	map<Set::Euclidean::Orthonormal::Point *,double>::const_iterator pC;
	for (pC=C.begin(); pC!=C.end(); pC++) 
		E += pC->second/Norm(x - *pC->first);
	return E;
}

//////////////////////////////////////////////////////////////////////
// Class Energy<1>
//////////////////////////////////////////////////////////////////////

Energy<1>::~Energy() {}

Potential::Field::Energy<1> *
Energy<1>::Clone() const
{
	return new Energy<1>(*this);
}

Energy<1>::Energy(LocalState *LS_) : LS(LS_) {}

Energy<1>::Energy(const Energy<1> &rhs) : LS(rhs.LS) {}

Set::VectorSpace::Vector 
Energy<1>::operator () (const Set::VectorSpace::Vector &x) const
{
	Set::VectorSpace::Vector DE(3);
	const map<Set::Euclidean::Orthonormal::Point *,double> &C = LS->Prop->C;
	map<Set::Euclidean::Orthonormal::Point *,double>::const_iterator pC;
	for (pC=C.begin(); pC!=C.end(); pC++)
	{
		Set::VectorSpace::Vector dx = x - *pC->first;
		DE += -pC->second*dx/pow(Norm(dx),3);
	}
	return DE;
}

//////////////////////////////////////////////////////////////////////
// Class Energy<2>
//////////////////////////////////////////////////////////////////////

Energy<2>::~Energy() {}

Potential::Field::Energy<2> *
Energy<2>::Clone() const
{
	return new Energy<2>(*this);
}

Energy<2>::Energy(LocalState *LS_) : LS(LS_) {}

Energy<2>::Energy(const Energy<2> &rhs) : LS(rhs.LS) {}

Set::VectorSpace::Hom 
Energy<2>::operator () (const Set::VectorSpace::Vector &x) const
{
	Set::VectorSpace::Hom DDE(3);
	const map<Set::Euclidean::Orthonormal::Point *,double> &C = LS->Prop->C;
	map<Set::Euclidean::Orthonormal::Point *,double>::const_iterator pC;
	for (pC=C.begin(); pC!=C.end(); pC++)
	{
		Set::VectorSpace::Vector dx = x - *pC->first;
		double dxdx = dx(dx);
		DDE += (pC->second/pow(dxdx,1.5))*
			(Dyadic(dx/(dxdx/3.0),dx)-Set::VectorSpace::HomId(3));
	}
	return DDE;
}

//////////////////////////////////////////////////////////////////////
// Class Jet<0>
//////////////////////////////////////////////////////////////////////

Jet<0>::~Jet() {}

Potential::Field::Jet<0> *
Jet<0>::Clone() const
{
	return new Jet<0>(*this);
}

Jet<0>::Jet(LocalState *LS_) : LS(LS_) {}

Jet<0>::Jet(const Jet<0> &rhs) : LS(rhs.LS) {}

pair<double,Set::VectorSpace::Vector>
Jet<0>::operator () (const Set::VectorSpace::Vector &x) const
{
	double E=0.0;
	Set::VectorSpace::Vector DE(3);
	const map<Set::Euclidean::Orthonormal::Point *,double> &C = LS->Prop->C;
	map<Set::Euclidean::Orthonormal::Point *,double>::const_iterator pC;
	for (pC=C.begin(); pC!=C.end(); pC++) 
	{
		Set::VectorSpace::Vector dx = x - *pC->first;
		double norm = Norm(dx);
		E += pC->second/norm;
		DE += -pC->second*dx/pow(norm,3);
	}
	return make_pair(E,DE);
}

//////////////////////////////////////////////////////////////////////
// Class Jet<1>
//////////////////////////////////////////////////////////////////////

Jet<1>::~Jet() {}

Potential::Field::Jet<1> *
Jet<1>::Clone() const
{
	return new Jet<1>(*this);
}

Jet<1>::Jet(LocalState *LS_) : LS(LS_) {}

Jet<1>::Jet(const Jet<1> &rhs) : LS(rhs.LS) {}

pair<Set::VectorSpace::Vector,Set::VectorSpace::Hom> 
Jet<1>::operator () (const Set::VectorSpace::Vector &x) const
{
	Set::VectorSpace::Vector DE(3);
	Set::VectorSpace::Hom DDE(3);
	const map<Set::Euclidean::Orthonormal::Point *,double> &C = LS->Prop->C;
	map<Set::Euclidean::Orthonormal::Point *,double>::const_iterator pC;
	for (pC=C.begin(); pC!=C.end(); pC++) 
	{
		Set::VectorSpace::Vector dx = x - *pC->first;
		double dxdx = dx(dx), den = pow(dxdx,1.5);
		DE += -pC->second*dx/den;
		DDE += (pC->second/den)*
			(Dyadic(dx/(dxdx/3.0),dx)-Set::VectorSpace::HomId(3));
	}
	return make_pair(DE,DDE);
}

}

}

}
