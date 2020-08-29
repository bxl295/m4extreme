// Conforming.h: interface for the Conforming class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(ELEMENT_CONFORMING_H__INCLUDED_)
#define ELEMENT_CONFORMING_H__INCLUDED_

#pragma once

#include <set>
#include <map>
#include <vector>
#include "../Element.h"
#include "../Quadrature/Quadrature.h"
#include "../Interpolation/Interpolation.h"
#include "../../Material/Material.h"

using namespace std;

namespace Element
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
		const vector<map<Set::Manifold::Point *, Set::VectorSpace::Vector> > &,
		const vector<double> &);
	LocalState(const LocalState &);
	virtual LocalState & operator = (const LocalState &);
	virtual const vector<map<Set::Manifold::Point *, Set::VectorSpace::Vector> > & GetDN() const;
	virtual vector<map<Set::Manifold::Point *, Set::VectorSpace::Vector> > & GetDN();
	virtual const vector<double> & GetQW() const;
	virtual vector<double> & GetQW();
	virtual range_type operator () (const domain_type &) const;
	virtual void operator ++ ();
	virtual set<Set::Manifold::Point *> GetNodes() const;
	virtual set<pair<Set::Manifold::Point *, 
		Set::Manifold::Point *> > GetNodePairs() const;

	virtual void GetNodes(set<Set::Manifold::Point *> &) const;


protected:

	vector<Material::LocalState *> MatLS;
	vector<map<Set::Manifold::Point *, Set::VectorSpace::Vector> > DN;
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

	typedef Element::Conforming::Energy<1> tangent_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> domain_type;
	typedef double range_type;

	Energy();
	virtual ~Energy();
	virtual Element::Energy<0> *Clone() const;
	virtual Element::LocalState *GetLocalState() const;
	Energy(LocalState *, const vector<Material::Energy<0> *> &);
	Energy(const Energy<0> &);
	virtual Energy<0> & operator = (const Energy<0> &);
	virtual range_type operator () (const domain_type &) const;

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

	typedef Element::Conforming::Energy<2> tangent_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> domain_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> range_type;

	Energy();
	virtual ~Energy();
	virtual Element::Energy<1> *Clone() const;
	virtual Element::LocalState *GetLocalState() const;
	Energy(LocalState *, const vector<Material::Energy<1> *> &);
	Energy(const Energy<1> &);
	virtual Energy<1> & operator = (const Energy<1> &);
	virtual range_type operator () (const domain_type &) const;
	virtual void operator () (const domain_type &, range_type &) const;

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

	typedef Element::Conforming::Energy<3> tangent_type;
	typedef map<Set::Manifold::Point *, 
		Set::VectorSpace::Vector> domain_type;
	typedef map<pair<Set::Manifold::Point *, 
		Set::Manifold::Point *>, 
		Set::VectorSpace::Hom> range_type;

	Energy();
	virtual ~Energy();
	virtual Element::Energy<2> *Clone() const;
	virtual Element::LocalState *GetLocalState() const;
	Energy(LocalState *, const vector<Material::Energy<2> *> &);
	Energy(const Energy<2> &);
	virtual Energy<2> & operator = (const Energy<2> &);
	virtual range_type operator () (const domain_type &) const;

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

	typedef Element::Jet<1> tangent_type;
	typedef map<Set::Manifold::Point *, 
		Set::VectorSpace::Vector> domain_type;
	typedef pair<double, map<Set::Manifold::Point *, 
		Set::VectorSpace::Vector> > range_type;

	Jet();
	virtual ~Jet();
	virtual Element::Jet<0> *Clone() const;
	virtual Element::LocalState *GetLocalState() const;
	Jet(LocalState *, const vector<Material::Jet<0> *> &);
	Jet(const Jet<0> &);
	virtual Jet<0> & operator = (const Jet<0> &);
	virtual range_type operator () (const domain_type &) const;

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

	typedef Element::Jet<2> tangent_type;
	typedef map<Set::Manifold::Point *, 
		Set::VectorSpace::Vector> domain_type;
	typedef pair<
		map<Set::Manifold::Point *, 
		Set::VectorSpace::Vector>, 
		map<pair<Set::Manifold::Point *, 
		Set::Manifold::Point *>, 
		Set::VectorSpace::Hom> > range_type;

	Jet();
	virtual ~Jet();
	virtual Element::Jet<1> *Clone() const;
	virtual Element::LocalState *GetLocalState() const;
	Jet(LocalState *, const vector<Material::Jet<1> *> &);
	Jet(const Jet<1> &);
	virtual Jet<1> & operator = (const Jet<1> &);
	virtual range_type operator () (const domain_type &) const;

protected:

	LocalState *LS;
	vector<Material::Jet<1> *> DJ;
};

}

}

#endif // !defined(ELEMENT_CONFORMING_H__INCLUDED_
