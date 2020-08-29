// Source.h: interface for the Source class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(POTENTIAL_CONFORMING_H__INCLUDED_)
#define POTENTIAL_CONFORMING_H__INCLUDED_

#pragma once

#include <set>
#include <map>
#include <vector>
#include "../Element.h"
#include "../Quadrature/Quadrature.h"
#include "../Interpolation/Interpolation.h"
#include "../../Material/Material.h"

using namespace std;

namespace Potential
{
namespace Conforming
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
	virtual Element::LocalState *Clone() const;
	LocalState(const vector<Material::LocalState *> &);
	LocalState(
		const vector<Material::LocalState *> &,
		const vector<double> &);
	LocalState(
		const vector<Material::LocalState *> &,
		const vector<map<Set::Manifold::Point *, double> > &,
		const vector<double> &);
	LocalState(const LocalState &);
	LocalState & operator = (const LocalState &);
	vector<map<Set::Manifold::Point *, double> > & GetN();
	vector<double> & GetQW();
	range_type operator () (const domain_type &) const;
	void operator ++ ();
	set<Set::Manifold::Point *> GetNodes() const;
	set<pair<Set::Manifold::Point *, Set::Manifold::Point *> > GetNodePairs() const;

        void operator () (const domain_type &, range_type &) const;
        
protected:

	vector<Material::LocalState *> MatLS;
	vector<map<Set::Manifold::Point *, double> > N;
	vector<double> QW;
};

//////////////////////////////////////////////////////////////////////
// Class Energy<p>
//////////////////////////////////////////////////////////////////////

template<unsigned int p> class Energy;

//////////////////////////////////////////////////////////////////////
// Class Energy<0>
//////////////////////////////////////////////////////////////////////

template<>
class Energy<0> : public Element::Energy<0>
{
public: 

	typedef Potential::Conforming::Energy<1> tangent_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> domain_type;
	typedef double range_type;

	Energy();
	virtual ~Energy();
	virtual Element::Energy<0> *Clone() const;
	virtual Element::LocalState *GetLocalState() const;
	Energy(LocalState *, const vector<Material::Energy<0> *> &);
	Energy(const Energy<0> &);
	Energy<0> & operator = (const Energy<0> &);
	range_type operator () (const domain_type &) const;

protected:

	LocalState *LS;
	vector<Material::Energy<0> *> W;
};

//////////////////////////////////////////////////////////////////////
// Class Energy<1>
//////////////////////////////////////////////////////////////////////

template<>
class Energy<1> : public Element::Energy<1>
{
public: 

	typedef Potential::Conforming::Energy<2> tangent_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> domain_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> range_type;

	Energy();
	virtual ~Energy();
	virtual Element::Energy<1> *Clone() const;
	virtual Element::LocalState *GetLocalState() const;
	Energy(LocalState *, const vector<Material::Energy<1> *> &);
	Energy(const Energy<1> &);
	Energy<1> & operator = (const Energy<1> &);
	range_type operator () (const domain_type &) const;
        
        void operator()(const domain_type &, range_type &) const;

protected:

	LocalState *LS;
	vector<Material::Energy<1> *> DW;
};

//////////////////////////////////////////////////////////////////////
// Class Energy<2>
//////////////////////////////////////////////////////////////////////

template<>
class Energy<2> : public Element::Energy<2>
{
public: 

	typedef Potential::Conforming::Energy<3> tangent_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> domain_type;
	typedef map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, 
		Set::VectorSpace::Hom> range_type;

	Energy();
	virtual ~Energy();
	virtual Element::Energy<2> *Clone() const;
	virtual Element::LocalState *GetLocalState() const;
	Energy(LocalState *, const vector<Material::Energy<2> *> &);
	Energy(const Energy<2> &);
	Energy<2> & operator = (const Energy<2> &);
	range_type operator () (const domain_type &) const;

protected:

	LocalState *LS;
	vector<Material::Energy<2> *> DDW;
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

	typedef Potential::Conforming::Jet<1> tangent_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> domain_type;
	typedef pair<double, map<Set::Manifold::Point *, Set::VectorSpace::Vector> > range_type;

	Jet();
	virtual ~Jet();
	virtual Element::Jet<0> *Clone() const;
	virtual Element::LocalState *GetLocalState() const;
	Jet(LocalState *, const vector<Material::Jet<0> *> &);
	Jet(const Jet<0> &);
	Jet<0> & operator = (const Jet<0> &);
	range_type operator () (const domain_type &) const;

protected:

	LocalState *LS;
	vector<Material::Jet<0> *> J;
};

//////////////////////////////////////////////////////////////////////
// Class Jet<1>
//////////////////////////////////////////////////////////////////////

template <> 
class Jet<1> : public Element::Jet<1>
{
public: 

	typedef Potential::Conforming::Jet<2> tangent_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> domain_type;
	typedef pair<
		map<Set::Manifold::Point *, Set::VectorSpace::Vector>, 
		map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, 
		Set::VectorSpace::Hom> > range_type;

	Jet();
	virtual ~Jet();
	virtual Element::Jet<1> *Clone() const;
	virtual Element::LocalState *GetLocalState() const;
	Jet(LocalState *, const vector<Material::Jet<1> *> &);
	Jet(const Jet<1> &);
	Jet<1> & operator = (const Jet<1> &);
	range_type operator () (const domain_type &) const;

protected:

	LocalState *LS;
	vector<Material::Jet<1> *> DJ;
};

}

}

#endif // !defined(POTENTIAL_CONFORMING_H__INCLUDED_
