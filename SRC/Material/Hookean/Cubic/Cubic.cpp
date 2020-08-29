// Cubic.cpp: implementation of the Cubic class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <cassert>
#include "./Cubic.h"
#include "../../../Utils/Indexing/Indexing.h"

namespace Material
{
namespace Hookean
{
namespace Cubic
{
//////////////////////////////////////////////////////////////////////
// Class Data
//////////////////////////////////////////////////////////////////////

Data::Data() : C11(0.0), C12(0.0), C44(0.0) {}

Data::~Data(){}

Data::Data(const double * const rhs) 
	: C11(rhs[0]), C12(rhs[1]), C44(rhs[2]) {}

Data::Data(
	const double &rhs_C11, 
	const double &rhs_C12, 
	const double &rhs_C44) 
: C11(rhs_C11), C12(rhs_C12), C44(rhs_C44) {}

Data::Data(const Data &rhs) 
: C11(rhs.C11), C12(rhs.C12), C44(rhs.C44) {}

Data &
Data::operator = (const Data &rhs)
{
	if (this == &rhs) return *this;
	C11 = rhs.C11; C12 = rhs.C12; C44 = rhs.C44;
	return *this;
}

const double &
Data::GetC11() const
{
	return C11;
}

double & 
Data::GetC11()
{
	return C11;
}

const double &
Data::GetC12() const
{
	return C12;
}

double & 
Data::GetC12()
{
	return C12;
}

const double &
Data::GetC44() const
{
	return C44;
}

double & 
Data::GetC44()
{
	return C44;
}

void 
Data::Randomize()
{
	C12 = (double)rand()/(double)RAND_MAX;
	C11 = C12 + 2.0*(double)rand()/(double)RAND_MAX;
	C44 = (double)rand()/(double)RAND_MAX;
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
	const double &C11 = LS->Properties->C11;
	const double &C12 = LS->Properties->C12;
	const double &C44 = LS->Properties->C44;
	Set::VectorSpace::Vector Eps(6);
	Eps[0] = Du[0];
	Eps[1] = Du[4];
	Eps[2] = Du[8];
	Eps[3] = Du[1] + Du[3];
	Eps[4] = Du[2] + Du[6];
	Eps[5] = Du[5] + Du[7];
	Set::VectorSpace::Vector Sig(6);
	Sig[0] = C11*Eps[0] + C12*(Eps[1] + Eps[2]);
	Sig[1] = C11*Eps[1] + C12*(Eps[0] + Eps[2]);
	Sig[2] = C11*Eps[2] + C12*(Eps[0] + Eps[1]);
	Sig[3] = C44*Eps[3];
	Sig[4] = C44*Eps[4];
	Sig[5] = C44*Eps[5];
	return 0.5*Sig(Eps);
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
	const double &C11 = LS->Properties->C11;
	const double &C12 = LS->Properties->C12;
	const double &C44 = LS->Properties->C44;
	Set::VectorSpace::Vector Eps(6);
	Eps[0] = Du[0];
	Eps[1] = Du[4];
	Eps[2] = Du[8];
	Eps[3] = Du[1] + Du[3];
	Eps[4] = Du[2] + Du[6];
	Eps[5] = Du[5] + Du[7];
	Set::VectorSpace::Vector Sig(6);
	Sig[0] = C11*Eps[0] + C12*(Eps[1] + Eps[2]);
	Sig[1] = C11*Eps[1] + C12*(Eps[0] + Eps[2]);
	Sig[2] = C11*Eps[2] + C12*(Eps[0] + Eps[1]);
	Sig[3] = C44*Eps[3];
	Sig[4] = C44*Eps[4];
	Sig[5] = C44*Eps[5];
	Set::VectorSpace::Vector DW(9);
	DW[0] = Sig[0];
	DW[4] = Sig[1];
	DW[8] = Sig[2];
	DW[1] = DW[3] = Sig[3];
	DW[2] = DW[6] = Sig[4];
	DW[5] = DW[7] = Sig[5];
	return DW;
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
	const double &C11 = LS->Properties->C11;
	const double &C12 = LS->Properties->C12;
	const double &C44 = LS->Properties->C44;
	Set::VectorSpace::Hom DDW(9);
	DDW[0][0] = DDW[4][4] = DDW[8][8] = C11;
	DDW[0][4] = DDW[4][0] = C12;
	DDW[0][8] = DDW[8][0] = C12;
	DDW[4][8] = DDW[8][4] = C12;
	DDW[1][1] = DDW[1][3] = DDW[3][1] = DDW[3][3] = C44;
	DDW[2][2] = DDW[2][6] = DDW[6][2] = DDW[6][6] = C44;
	DDW[5][5] = DDW[5][7] = DDW[7][5] = DDW[7][7] = C44;
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
	const double &C11 = LS->Properties->C11;
	const double &C12 = LS->Properties->C12;
	const double &C44 = LS->Properties->C44;
	Set::VectorSpace::Vector Eps(6);
	Eps[0] = Du[0];
	Eps[1] = Du[4];
	Eps[2] = Du[8];
	Eps[3] = Du[1] + Du[3];
	Eps[4] = Du[2] + Du[6];
	Eps[5] = Du[5] + Du[7];
	Set::VectorSpace::Vector Sig(6);
	Sig[0] = C11*Eps[0] + C12*(Eps[1] + Eps[2]);
	Sig[1] = C11*Eps[1] + C12*(Eps[0] + Eps[2]);
	Sig[2] = C11*Eps[2] + C12*(Eps[0] + Eps[1]);
	Sig[3] = C44*Eps[3];
	Sig[4] = C44*Eps[4];
	Sig[5] = C44*Eps[5];
	Set::VectorSpace::Vector DW(9);
	DW[0] = Sig[0];
	DW[4] = Sig[1];
	DW[8] = Sig[2];
	DW[1] = DW[3] = Sig[3];
	DW[2] = DW[6] = Sig[4];
	DW[5] = DW[7] = Sig[5];
	return make_pair(0.5*Sig(Eps),DW);
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
	const double &C11 = LS->Properties->C11;
	const double &C12 = LS->Properties->C12;
	const double &C44 = LS->Properties->C44;
	Set::VectorSpace::Vector Eps(6);
	Eps[0] = Du[0];
	Eps[1] = Du[4];
	Eps[2] = Du[8];
	Eps[3] = Du[1] + Du[3];
	Eps[4] = Du[2] + Du[6];
	Eps[5] = Du[5] + Du[7];
	Set::VectorSpace::Vector Sig(6);
	Sig[0] = C11*Eps[0] + C12*(Eps[1] + Eps[2]);
	Sig[1] = C11*Eps[1] + C12*(Eps[0] + Eps[2]);
	Sig[2] = C11*Eps[2] + C12*(Eps[0] + Eps[1]);
	Sig[3] = C44*Eps[3];
	Sig[4] = C44*Eps[4];
	Sig[5] = C44*Eps[5];
	Set::VectorSpace::Vector DW(9);
	DW[0] = Sig[0];
	DW[4] = Sig[1];
	DW[8] = Sig[2];
	DW[1] = DW[3] = Sig[3];
	DW[2] = DW[6] = Sig[4];
	DW[5] = DW[7] = Sig[5];
	Set::VectorSpace::Hom DDW(9);
	DDW[0][0] = DDW[4][4] = DDW[8][8] = C11;
	DDW[0][4] = DDW[4][0] = C12;
	DDW[0][8] = DDW[8][0] = C12;
	DDW[4][8] = DDW[8][4] = C12;
	DDW[1][1] = DDW[1][3] = DDW[3][1] = DDW[3][3] = C44;
	DDW[2][2] = DDW[2][6] = DDW[6][2] = DDW[6][6] = C44;
	DDW[5][5] = DDW[5][7] = DDW[7][5] = DDW[7][7] = C44;
	return make_pair(DW,DDW);
}

}

}

}
