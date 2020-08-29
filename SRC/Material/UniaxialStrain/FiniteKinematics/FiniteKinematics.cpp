// FiniteKinematics.cpp: implementation of the FiniteKinematics class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <cassert>
#include "./FiniteKinematics.h"

namespace Material
{
namespace UniaxialStrain
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

LocalState::LocalState(Material::LocalState *LS_) : LS(LS_) {}

LocalState::LocalState(const LocalState &rhs) : LS(rhs.LS) {}

void 
LocalState::operator ++ () 
{
	++(*LS);
}

void 
LocalState::Reset(const Set::VectorSpace::Hom &Du)
{
	assert (Du.size() == 1);
	Set::VectorSpace::Hom Dv(3);
	Dv[0][0] = Du[0][0];
	Dv[1][1] = 1.0;
	Dv[2][2] = 1.0;
	LS->Reset(Dv);
}

//////////////////////////////////////////////////////////////////////
// Class Energy<0>
//////////////////////////////////////////////////////////////////////

Energy<0>::Energy(Material::Energy<0> *W_) : W(W_) {}

Energy<0>::~Energy() {}

Material::Energy<0> *
Energy<0>::Clone() const
{
	return new Energy<0>(*this);
}

Energy<0>::Energy(const Energy<0> &rhs) : W(rhs.W) {}

double 
Energy<0>::operator () (const Set::VectorSpace::Vector &Du) const
{
	assert (Du.size() == 1);
	Set::VectorSpace::Vector Dv(9);
	Dv[0] = Du[0];
	Dv[4] = 1.0;
	Dv[8] = 1.0;
	return (*W)(Dv);
}

//////////////////////////////////////////////////////////////////////
// Class Energy<1>
//////////////////////////////////////////////////////////////////////

Energy<1>::Energy(Material::Energy<1> *DW_) : DW(DW_) {}

Energy<1>::~Energy() {}

Material::Energy<1> *
Energy<1>::Clone() const
{
	return new Energy<1>(*this);
}

Energy<1>::Energy(const Energy<1> &rhs) : DW(rhs.DW) {}

Set::VectorSpace::Vector 
Energy<1>::operator () (const Set::VectorSpace::Vector &Du) const
{
	assert (Du.size() == 1);
	Set::VectorSpace::Vector Dv(9);
	Dv[0] = Du[0];
	Dv[4] = 1.0;
	Dv[8] = 1.0;
	Set::VectorSpace::Vector P = (*DW)(Dv);
	Set::VectorSpace::Vector Q(1);
	Q[0] = P[0]; 
	return Q;
}

//////////////////////////////////////////////////////////////////////
// Class Energy<2>
//////////////////////////////////////////////////////////////////////

Energy<2>::Energy(Material::Energy<2> *DDW_) : DDW(DDW_) {}

Energy<2>::~Energy() {}

Material::Energy<2> *
Energy<2>::Clone() const
{
	return new Energy<2>(*this);
}

Energy<2>::Energy(const Energy<2> &rhs) : DDW(rhs.DDW) {}

Set::VectorSpace::Hom 
Energy<2>::operator () (const Set::VectorSpace::Vector &Du) const
{
	assert (Du.size() == 1);
	Set::VectorSpace::Vector Dv(9);
	Dv[0] = Du[0];
	Dv[4] = 1.0;
	Dv[8] = 1.0;
	Set::VectorSpace::Hom DP = (*DDW)(Dv);
	Set::VectorSpace::Hom DQ(1);
	DQ[0][0] = DP[0][0]; 
	return DQ;
}

//////////////////////////////////////////////////////////////////////
// Class Jet<0>
//////////////////////////////////////////////////////////////////////

Jet<0>::~Jet() {}

Jet<0>::Jet(Material::Jet<0> *J_) : J(J_) {}

Material::Jet<0> *
Jet<0>::Clone() const
{
	return new Jet<0>(*this);
}

Jet<0>::Jet(const Jet<0> &rhs) : J(rhs.J) {}

pair<double,Set::VectorSpace::Vector>
Jet<0>::operator () (const Set::VectorSpace::Vector &Du) const
{
	assert (Du.size() == 1);
	Set::VectorSpace::Vector Dv(9);
	Dv[0] = Du[0];
	Dv[4] = 1.0;
	Dv[8] = 1.0;
	pair<double,Set::VectorSpace::Vector> K = (*J)(Dv);
	Set::VectorSpace::Vector P = K.second;
	Set::VectorSpace::Vector Q(1);
	Q[0] = P[0]; 
	return make_pair(K.first,Q);
}

//////////////////////////////////////////////////////////////////////
// Class Jet<1>
//////////////////////////////////////////////////////////////////////

Jet<1>::~Jet() {}

Jet<1>::Jet(Material::Jet<1> *DJ_) : DJ(DJ_) {}

Material::Jet<1> *
Jet<1>::Clone() const
{
	return new Jet<1>(*this);
}

Jet<1>::Jet(const Jet<1> &rhs) : DJ(rhs.DJ) {}

pair<Set::VectorSpace::Vector,Set::VectorSpace::Hom> 
Jet<1>::operator () (const Set::VectorSpace::Vector &Du) const
{
	assert (Du.size() == 1);
	Set::VectorSpace::Vector Dv(9);
	Dv[0] = Du[0];
	Dv[4] = 1.0;
	Dv[8] = 1.0;
	pair<
		Set::VectorSpace::Vector, 
		Set::VectorSpace::Hom> DK = (*DJ)(Dv);
	Set::VectorSpace::Vector P = DK.first;
	Set::VectorSpace::Vector Q(1);
	Q[0] = P[0]; 
	Set::VectorSpace::Hom DP = DK.second;
	Set::VectorSpace::Hom DQ(1);
	DQ[0][0] = DP[0][0]; 
	return make_pair(Q,DQ);
}

}

}

}
