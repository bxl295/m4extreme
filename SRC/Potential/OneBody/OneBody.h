// OneBody.h: interface for the OneBody class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(POTENTIAL_ONEBODY_H__INCLUDED_)
#define POTENTIAL_ONEBODY_H__INCLUDED_

#pragma once

#include <vector>
#include <utility>
#include "../Field/Field.h"
#include "../../Element/Element.h"

using namespace std;

namespace Potential
{
namespace OneBody
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
	LocalState(Set::Manifold::Point *);
	LocalState(const LocalState &);
	LocalState & operator = (const LocalState &);
	void operator ++ ();
	set<Set::Manifold::Point *> GetNodes() const;
	set<pair<Set::Manifold::Point *, Set::Manifold::Point *> > GetNodePairs() const;

	void GetNodes(set<Set::Manifold::Point *>  &) const;

protected:

	Set::Manifold::Point *e0;
};

//////////////////////////////////////////////////////////////////////
// Class Energy<0>
//////////////////////////////////////////////////////////////////////

template<> 
class Energy<0> : public Element::Energy<0>
{
public: 

	typedef Potential::OneBody::Energy<1> tangent_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> domain_type;
	typedef double range_type;

	Energy();
	virtual ~Energy();
	Element::Energy<0> *Clone() const;
	Element::LocalState *GetLocalState() const;
	Energy(LocalState *, Potential::Field::Energy<0> *);
	Energy(const Energy<0> &);
	const Energy<0> & operator = (const Energy<0> &);
	range_type operator () (const domain_type &) const;

private:

	LocalState *LS;
	Potential::Field::Energy<0> *V;
};

//////////////////////////////////////////////////////////////////////
// Class Energy<1>
//////////////////////////////////////////////////////////////////////

template<> 
class Energy<1> : public Element::Energy<1>
{
public: 

	typedef Potential::OneBody::Energy<2> tangent_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> domain_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> range_type;

	Energy();
	virtual ~Energy();
	Element::Energy<1> *Clone() const;
	Element::LocalState *GetLocalState() const;
	Energy(LocalState *, Potential::Field::Energy<1> *);
	Energy(const Energy<1> &);
	const Energy<1> & operator = (const Energy<1> &);
	range_type operator () (const domain_type &) const;

	void operator () (const domain_type &, range_type &) const;

private:

	LocalState *LS;
	Potential::Field::Energy<1> *DV;
};

//////////////////////////////////////////////////////////////////////
// Class Energy<2>
//////////////////////////////////////////////////////////////////////

template<> 
class Energy<2> : public Element::Energy<2>
{
public: 

	typedef Potential::OneBody::Energy<2> tangent_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> domain_type;
	typedef map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, 
		Set::VectorSpace::Hom> range_type;

	Energy();
	virtual ~Energy();
	Element::Energy<2> *Clone() const;
	Element::LocalState *GetLocalState() const;
	Energy(LocalState *, Potential::Field::Energy<2> *);
	Energy(const Energy<2> &);
	const Energy<2> & operator = (const Energy<2> &);
	range_type operator () (const domain_type &) const;

private:

	LocalState *LS;
	Potential::Field::Energy<2> *DDV;
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
	Jet(LocalState *, Potential::Field::Jet<0> *rhs_J);
	Jet(const Jet<0> &);
	Jet<0> & operator = (const Jet<0> &);
	range_type operator () (const domain_type &) const;

protected:

	LocalState *LS;
	Potential::Field::Jet<0> *J;
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
	Jet(LocalState *, Potential::Field::Jet<1> *);
	Jet(const Jet<1> &);
	Jet<1> & operator = (const Jet<1> &);
	range_type operator () (const domain_type &) const;

protected:

	LocalState *LS;
	Potential::Field::Jet<1> *DJ;
};

}

}

#endif // !defined(POTENTIAL_ONEBODY_H__INCLUDED_)

