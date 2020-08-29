// TwoBody.h: interface for the TwoBody class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(POTENTIAL_TWOBODY_H__INCLUDED_)
#define POTENTIAL_TWOBODY_H__INCLUDED_

#pragma once

#include <vector>
#include <utility>
#include "../Radial/Radial.h"
#include "../../Element/Element.h"

using namespace std;

namespace Potential
{
namespace TwoBody
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
	LocalState(Set::Manifold::Point *, Set::Manifold::Point *,
		   Potential::Radial::LocalState *rhs_LS=NULL);
	LocalState(const set<Set::Manifold::Point *> &,
		   Potential::Radial::LocalState *rhs_LS=NULL);
	LocalState(const LocalState &);
	LocalState & operator = (const LocalState &);
	void operator ++ ();
	Potential::Radial::LocalState * GetLS();
	set<Set::Manifold::Point *> GetNodes() const;
	set<pair<Set::Manifold::Point *, Set::Manifold::Point *> > GetNodePairs() const;

	void GetNodes(set<Set::Manifold::Point *> & Nodes) const;
	void GetNodePairs(set<pair<Set::Manifold::Point *, Set::Manifold::Point *> > &NodePairs) const;

 protected:

	Set::Manifold::Point *e1;
	Set::Manifold::Point *e2;
	Potential::Radial::LocalState *LS;
};

//////////////////////////////////////////////////////////////////////
// Class Energy<0>
//////////////////////////////////////////////////////////////////////

template<> 
class Energy<0> : public Element::Energy<0>
{
public: 

	typedef Potential::TwoBody::Energy<1> tangent_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> domain_type;
	typedef double range_type;

	Energy();
	virtual ~Energy();
	Element::Energy<0> *Clone() const;
	Element::LocalState *GetLocalState() const;
	Energy(
		LocalState *,
		Potential::Radial::Energy<0> *);
	Energy(const Energy<0> &);
	const Energy<0> & operator = (const Energy<0> &);
	range_type operator () (const domain_type &) const;

private:

	LocalState *LS;
	Potential::Radial::Energy<0> *V;
};

//////////////////////////////////////////////////////////////////////
// Class Energy<1>
//////////////////////////////////////////////////////////////////////

template<> 
class Energy<1> : public Element::Energy<1>
{
public: 

	typedef Potential::TwoBody::Energy<2> tangent_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> domain_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> range_type;

	Energy();
	virtual ~Energy();
	Element::Energy<1> *Clone() const;
	Element::LocalState *GetLocalState() const;
	Energy(
		LocalState *,
		Potential::Radial::Energy<1> *);
	Energy(const Energy<1> &);
	const Energy<1> & operator = (const Energy<1> &);
	range_type operator () (const domain_type &) const;

	void operator () (const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &x,
			  map<Set::Manifold::Point *, Set::VectorSpace::Vector> & DE) const;
	
private:

	LocalState *LS;
	Potential::Radial::Energy<1> *DV;
};

//////////////////////////////////////////////////////////////////////
// Class Energy<2>
//////////////////////////////////////////////////////////////////////

template<> 
class Energy<2> : public Element::Energy<2>
{
public: 

	typedef Potential::TwoBody::Energy<2> tangent_type;
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
		Potential::Radial::Energy<2> *);
	Energy(const Energy<2> &);
	const Energy<2> & operator = (const Energy<2> &);
	range_type operator () (const domain_type &) const;

private:

	LocalState *LS;
	Potential::Radial::Energy<1> *DV;
	Potential::Radial::Energy<2> *DDV;
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
	Jet(LocalState *, Potential::Radial::Jet<0> *rhs_J);
	Jet(const Jet<0> &);
	Jet<0> & operator = (const Jet<0> &);
	range_type operator () (const domain_type &) const;

protected:

	LocalState *LS;
	Potential::Radial::Jet<0> *J;
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
	Jet(LocalState *, Potential::Radial::Jet<1> *);
	Jet(const Jet<1> &);
	Jet<1> & operator = (const Jet<1> &);
	range_type operator () (const domain_type &) const;

protected:

	LocalState *LS;
	Potential::Radial::Jet<1> *DJ;
};

}

}

#endif // !defined(POTENTIAL_TWOBODY_H__INCLUDED_)

