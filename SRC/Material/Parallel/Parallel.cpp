// Parallel.cpp: implementation of the Parallel class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <cassert>
#include "./Parallel.h"
#include "../../Utils/Indexing/Indexing.h"

namespace Material
{
namespace Parallel
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

LocalState::LocalState(const vector<Material::LocalState *> &LS_) : LS(LS_) {}

LocalState::LocalState(const LocalState &rhs) : LS(rhs.LS) {}

void 
LocalState::operator ++ () 
{
	vector<Material::LocalState *>::iterator pLS;
	for (pLS=LS.begin(); pLS!=LS.end(); pLS++) ++(*(*pLS));
}

const vector<Material::LocalState *> & 
LocalState::GetLS() const
{
	return LS;
}

vector<Material::LocalState *> & 
LocalState::GetLS()
{
	return LS;
}

//////////////////////////////////////////////////////////////////////
// Class Energy<0>
//////////////////////////////////////////////////////////////////////

Energy<0>::~Energy() {}

Material::Energy<0> *
Energy<0>::Clone() const
{
	return new Energy<0>(*this);
}

Energy<0>::Energy(const vector<Material::Energy<0> *> &f_) : f(f_) {}

Energy<0>::Energy(const Energy<0> &rhs) : f(rhs.f) {}

double 
Energy<0>::operator () (const Set::VectorSpace::Vector &Dy) const
{
	assert (Dy.size() == 9); double W = 0.0;
	vector<Material::Energy<0> *>::const_iterator pf;
	for (pf=f.begin(); pf!=f.end(); pf++) W += (*(*pf))(Dy);
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

Energy<1>::Energy(const vector<Material::Energy<1> *> &Df_) : Df(Df_) {}

Energy<1>::Energy(const Energy<1> &rhs) : Df(rhs.Df) {}

Set::VectorSpace::Vector 
Energy<1>::operator () (const Set::VectorSpace::Vector &Dy) const
{
	assert (Dy.size() == 9);
	Set::VectorSpace::Vector DW(9);
	vector<Material::Energy<1> *>::const_iterator pDf;
	for (pDf=Df.begin(); pDf!=Df.end(); pDf++) DW += (*(*pDf))(Dy);
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

Energy<2>::Energy(const vector<Material::Energy<2> *> &DDf_) : DDf(DDf_) {}

Energy<2>::Energy(const Energy<2> &rhs) : DDf(rhs.DDf) {}

Set::VectorSpace::Hom 
Energy<2>::operator () (const Set::VectorSpace::Vector &Dy) const
{
	assert (Dy.size() == 9);
	Set::VectorSpace::Hom DDW(9);
	vector<Material::Energy<2> *>::const_iterator pDDf;
	for (pDDf=DDf.begin(); pDDf!=DDf.end(); pDDf++) DDW += (*(*pDDf))(Dy);
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

Jet<0>::Jet(const vector<Material::Jet<0> *> &g_) : g(g_) {}

Jet<0>::Jet(const Jet<0> &rhs) : g(rhs.g) {}

pair<double,Set::VectorSpace::Vector>
Jet<0>::operator () (const Set::VectorSpace::Vector &Dy) const
{
	assert (Dy.size() == 9);
	double W = 0.0; Set::VectorSpace::Vector DW(9);
	vector<Material::Jet<0> *>::const_iterator pg;
	for (pg=g.begin(); pg!=g.end(); pg++) 
	{
		pair<double,Set::VectorSpace::Vector> J = (*(*pg))(Dy);
		W += J.first; DW += J.second;
	}
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

Jet<1>::Jet(const vector<Material::Jet<1> *> &Dg_) : Dg(Dg_) {}

Jet<1>::Jet(const Jet<1> &rhs) : Dg(rhs.Dg) {}

pair<Set::VectorSpace::Vector,Set::VectorSpace::Hom> 
Jet<1>::operator () (const Set::VectorSpace::Vector &Dy) const
{
	assert (Dy.size() == 9);
	Set::VectorSpace::Vector DW(9);
	Set::VectorSpace::Hom DDW(9);
	vector<Material::Jet<1> *>::const_iterator pDg;
	for (pDg=Dg.begin(); pDg!=Dg.end(); pDg++) 
	{
		pair<Set::VectorSpace::Vector, Set::VectorSpace::Hom> DJ = (*(*pDg))(Dy);
		DW += DJ.first; DDW += DJ.second;
	}
	return make_pair(DW,DDW);
}

}

}
