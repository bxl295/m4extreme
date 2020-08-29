// EmbeddedAtom.h: interface for the EmbeddedAtom class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(POTENTIAL_EMBEDDEDATOM_H__INCLUDED_)
#define POTENTIAL_EMBEDDEDATOM_H__INCLUDED_

#pragma once

#include "../Radial/Radial.h"
#include "../TwoBody/TwoBody.h"
#include "../../Element/Element.h"

using namespace std;

namespace Potential
{
namespace EmbeddedAtom
{
class LocalState;
template <unsigned int> class Energy;
template <unsigned int> class Jet;

//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

class LocalState : public Element::LocalState
{
friend class Energy<0>;
friend class Energy<1>;
friend class Energy<2>;
friend class Jet<0>;
friend class Jet<1>;

public: 

	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> domain_type;
	typedef vector<Set::VectorSpace::Vector> range_type;

	LocalState();
	virtual ~LocalState();
	Element::LocalState *Clone() const;
	LocalState(Set::Manifold::Point *, const set<Set::Manifold::Point *> &);
	LocalState(const LocalState &);
	LocalState & operator = (const LocalState &);
	void operator ++ ();
	set<Set::Manifold::Point *> GetNodes() const;
	set<pair<Set::Manifold::Point *, Set::Manifold::Point *> > GetNodePairs() const;

protected:

	Set::Manifold::Point *e;
	set<Set::Manifold::Point *> N;

	Set::VectorSpace::Hom Dyadic(
		const Set::VectorSpace::Vector &,
		const Set::VectorSpace::Vector &);
};

//////////////////////////////////////////////////////////////////////
// Class Energy<0>
//////////////////////////////////////////////////////////////////////

template<> 
class Energy<0> : public Element::Energy<0>
{
public: 

	typedef Potential::EmbeddedAtom::Energy<1> tangent_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> domain_type;
	typedef double range_type;

	Energy();
	virtual ~Energy();
	Element::Energy<0> *Clone() const;
	Element::LocalState *GetLocalState() const;
	Energy(
		LocalState *,
		Potential::Radial::Energy<0> *,
		Potential::Radial::Energy<0> *,
		Potential::Radial::Energy<0> *);
	Energy(const Energy<0> &);
	const Energy<0> & operator = (const Energy<0> &);
	range_type operator () (const domain_type &) const;

private:

	LocalState *LS;
	Potential::Radial::Energy<0> *U;				//	Embedding function
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<0> *> R;	//	Electron density function
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<0> *> V;	//	Pair potential
	map<Set::Manifold::Point *, Potential::TwoBody::LocalState *> TBLS;
};

//////////////////////////////////////////////////////////////////////
// Class Energy<1>
//////////////////////////////////////////////////////////////////////

template<> 
class Energy<1> : public Element::Energy<1>
{
public: 

	typedef Potential::EmbeddedAtom::Energy<2> tangent_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> domain_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> range_type;

	Energy();
	virtual ~Energy();
	Element::Energy<1> *Clone() const;
	Element::LocalState *GetLocalState() const;
	Energy(
		LocalState *,
		Potential::Radial::Energy<1> *,
		Potential::Radial::Energy<0> *,
		Potential::Radial::Energy<1> *,
		Potential::Radial::Energy<1> *);
	Energy(const Energy<1> &);
	const Energy<1> & operator = (const Energy<1> &);
	range_type operator () (const domain_type &) const;

private:

	LocalState *LS;
	Potential::Radial::Energy<1> *DU;
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<0> *> R;
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<1> *> DR;
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<1> *> DV;
	map<Set::Manifold::Point *, Potential::TwoBody::LocalState *> TBLS;
};

//////////////////////////////////////////////////////////////////////
// Class Energy<2>
//////////////////////////////////////////////////////////////////////

template<> 
class Energy<2> : public Element::Energy<2>
{
public: 

	typedef Potential::EmbeddedAtom::Energy<2> tangent_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> domain_type;
	typedef map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, 
		Set::VectorSpace::Hom> range_type;

	Energy();
	virtual ~Energy();
	Element::Energy<2> *Clone() const;
	Element::LocalState *GetLocalState() const;
	Energy(
		LocalState *,
		Potential::Radial::Energy<1> *,
		Potential::Radial::Energy<2> *,
		Potential::Radial::Energy<0> *,
		Potential::Radial::Energy<1> *,
		Potential::Radial::Energy<2> *,
		Potential::Radial::Energy<1> *,
		Potential::Radial::Energy<2> *);
	Energy(const Energy<2> &);
	const Energy<2> & operator = (const Energy<2> &);
	range_type operator () (const domain_type &) const;

private:

	LocalState *LS;
	Potential::Radial::Energy<1> *DU;
	Potential::Radial::Energy<2> *DDU;
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<0> *> R;
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<1> *> DR;
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<2> *> DDR;
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<2> *> DDV;
	map<Set::Manifold::Point *, Potential::TwoBody::LocalState *> TBLS;
};

//////////////////////////////////////////////////////////////////////
// Class Jet<p>
//////////////////////////////////////////////////////////////////////

template<unsigned int p> class Jet;

//////////////////////////////////////////////////////////////////////
// Class Jet<0>
//////////////////////////////////////////////////////////////////////

template<> 
class Jet<0> : public Element::Jet<0>
{
public: 

	typedef Element::Jet<1> tangent_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> domain_type;
	typedef pair<double, map<Set::Manifold::Point *, Set::VectorSpace::Vector> > range_type;

	Jet();
	virtual ~Jet();
	Element::Jet<0> *Clone() const;
	Element::LocalState *GetLocalState() const;
	Jet(
		LocalState *,
		Potential::Radial::Jet<0> *rhs_UJ,
		Potential::Radial::Energy<0> *R0,
		Potential::Radial::Energy<1> *DR0,
		Potential::Radial::Jet<0> *VJ0);
	Jet(const Jet<0> &);
	Jet<0> & operator = (const Jet<0> &);
	range_type operator () (const domain_type &) const;

protected:

	LocalState *LS;
	Potential::Radial::Jet<0> *UJ;
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<0> *> R;
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<1> *> DR;
	map<Set::Manifold::Point *, Potential::TwoBody::Jet<0> *> VJ;
	map<Set::Manifold::Point *, Potential::TwoBody::LocalState *> TBLS;
};

//////////////////////////////////////////////////////////////////////
// Class Jet<1>
//////////////////////////////////////////////////////////////////////

template <> 
class Jet<1> : public Element::Jet<1>
{
public: 

	typedef Element::Jet<2> tangent_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> domain_type;
	typedef pair<
		map<Set::Manifold::Point *, Set::VectorSpace::Vector>, 
		map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, 
		Set::VectorSpace::Hom> > range_type;

	Jet();
	virtual ~Jet();
	Element::Jet<1> *Clone() const;
	Element::LocalState *GetLocalState() const;
	Jet(LocalState *,
		Potential::Radial::Energy<0> *,
		Potential::Radial::Jet<1> *,
		Potential::Radial::Jet<1> *,
		Potential::Radial::Jet<1> *);
	Jet(const Jet<1> &);
	Jet<1> & operator = (const Jet<1> &);
	range_type operator () (const domain_type &) const;

protected:

	LocalState *LS;
	Potential::Radial::Jet<1> *DUJ;
	map<Set::Manifold::Point *, Potential::TwoBody::Energy<0> *> R;
	map<Set::Manifold::Point *, Potential::TwoBody::Jet<1> *> DRJ;
	map<Set::Manifold::Point *, Potential::TwoBody::Jet<1> *> DVJ;
	map<Set::Manifold::Point *, Potential::TwoBody::LocalState *> TBLS;
};

}

}

#endif // !defined(POTENTIAL_EMBEDDEDATOM_H__INCLUDED_)

