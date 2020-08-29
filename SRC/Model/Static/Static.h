// Static.h: interface for the Static class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(MODEL_STATIC_H__INCLUED_)
#define MODEL_STATIC_H__INCLUED_

#pragma once

#include "../Model.h"
#include "../../Clock/Clock.h"
#include "../../Element/Element.h"
#include "../../Set/Manifold/Category/Category.h"
#include "../../Set/Manifold/Euclidean/Orthonormal/Orthonormal.h"
#include "../../Solver/Linear/Linear.h"

namespace Model
{
namespace Static
{
class LocalState;
template <unsigned int> class Energy;
template <unsigned int> class Jet;

//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

class LocalState : public Model::LocalState
{
friend class Energy<0>;
friend class Energy<1>;
friend class Energy<2>;
friend class Jet<0>;
friend class Jet<1>;

public: 

	typedef set<Set::Manifold::Point *> domain_type;        
        typedef Set::VectorSpace::Vector vector_type;
        
	LocalState();
	virtual ~LocalState();
	Model::LocalState *Clone() const;
	LocalState(Clock *, set<Element::LocalState *> *,
		map<Set::Manifold::Point *, Set::Manifold::Map *> *,
		map<Set::Manifold::Point *, Set::Manifold::TMap *> *);
	LocalState(const LocalState &rhs);
	LocalState & operator = (const LocalState &);
        virtual void Embed(const domain_type &,
			   map <Set::Manifold::Point *, Set::Euclidean::Orthonormal::Point> &);
	virtual void Embed(const domain_type &,
                           map <Set::Manifold::Point *, vector_type> &);
	virtual map <Set::Manifold::Point *, vector_type> Embed(const domain_type &);
	virtual void Embed(const map<Set::Manifold::Point *, vector_type> &,
			   map<Set::Manifold::Point *, vector_type> &);
	virtual void Submerge(const map<Set::Manifold::Point *, vector_type> &,
			      map<Set::Manifold::Point *, vector_type> &);
	virtual void Submerge( const map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, Set::VectorSpace::Hom> &,
			       map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, Set::VectorSpace::Hom> &);
	virtual void operator ++ ();
	set<DOF> GetDOFs() const; 
	set<pair<DOF,DOF> > GetDOFPairs() const;
	void operator () (Solver::Linear::System<DOF> &, domain_type &) const;
	void operator () (Solver::Linear::System<DOF> &, 
		map<Set::Manifold::Point *, vector_type> &) const;
	set<Element::LocalState *> *GetELS() const;

	const std::map<Element::LocalState*, int> & GetComputationalCost() const {
	      return _computational_cost;
	}


protected:

	Clock *T; set<Element::LocalState *> *ELS;
	map<Set::Manifold::Point *, Set::Manifold::Map *> *Emb;
	map<Set::Manifold::Point *, Set::Manifold::TMap *> *DEmb;

	std::map<Element::LocalState*, int> _computational_cost;
};

//////////////////////////////////////////////////////////////////////
// Class Energy<p>
//////////////////////////////////////////////////////////////////////

template<unsigned int p> class Energy;

//////////////////////////////////////////////////////////////////////
// Class Energy<0>
//////////////////////////////////////////////////////////////////////

template<> class Energy<0> : public Model::Energy<0>
{
public: 

	typedef Model::Static::Energy<1> tangent_type;
	typedef set<Set::Manifold::Point *> domain_type;
	typedef double range_type;

	Energy();
	virtual ~Energy();
	Model::Energy<0> *Clone() const;
	Energy(LocalState *, 
		set<Element::Energy<0> *> *);
	Energy(const Energy<0> &);
	Energy<0> & operator = (const Energy<0> &);
	virtual range_type operator () (const domain_type &) const;

protected:

	LocalState *LS;
	set<Element::Energy<0> *> *E;
};

//////////////////////////////////////////////////////////////////////
// Class Energy<1>
//////////////////////////////////////////////////////////////////////

template<> class Energy<1> : public Model::Energy<1>
{
friend class Jet<1>;

public: 

	typedef Model::Static::Energy<2> tangent_type;
	typedef set<Set::Manifold::Point *> domain_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> range_type;

	Energy();
	virtual ~Energy();
	Model::Energy<1> *Clone() const;
	Energy(LocalState *,
		set<Element::Energy<1> *> *);
	Energy(const Energy<1> &);
	Energy<1> & operator = (const Energy<1> &);

	virtual double operator () (const domain_type &, 
		const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &) const;
	virtual void operator () (const domain_type &, Solver::Linear::System<DOF> &) const;

	virtual range_type operator () (const domain_type &) const;
        virtual void operator () (const domain_type &, range_type &);

	const std::map<Element::LocalState*, int> & GetComputationalCost() const{
	      return _computational_cost;
	}

protected:

	LocalState *LS;
	set<Element::Energy<1> *> *DE;

	std::map<Element::LocalState*, int> _computational_cost;
};

//////////////////////////////////////////////////////////////////////
// Class Energy<2>
//////////////////////////////////////////////////////////////////////

template<> class Energy<2> : public Model::Energy<2>
{
public: 

	typedef Model::Static::Energy<3> tangent_type;
	typedef set<Set::Manifold::Point *> domain_type;
	typedef map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, 
		Set::VectorSpace::Hom> range_type;

	Energy();
	virtual ~Energy();
	Model::Energy<2> *Clone() const;
	Energy(LocalState *,
		set<Element::Energy<2> *> *);
	Energy(const Energy<2> &);
	Energy<2> & operator = (const Energy<2> &);
	virtual range_type operator () (const domain_type &) const;
	virtual double operator () (
		const domain_type &,
		const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &) const;
	virtual void operator () (const domain_type &, Solver::Linear::System<DOF> &) const;

protected:

	LocalState *LS;
	set<Element::Energy<2> *> *DDE;
};

//////////////////////////////////////////////////////////////////////
// Class Jet<p>
//////////////////////////////////////////////////////////////////////

template<unsigned int p> class Jet;

//////////////////////////////////////////////////////////////////////
// Class Jet<0>
//////////////////////////////////////////////////////////////////////

template<> class Jet<0> : public Model::Jet<0>
{
public: 

	typedef Model::Static::Jet<0> tangent_type;
	typedef set<Set::Manifold::Point *> domain_type;
	typedef pair<double, map<Set::Manifold::Point *, 
		Set::VectorSpace::Vector> > range_type;

	Jet();
	virtual ~Jet();
	Model::Jet<0> *Clone() const;
	Jet(LocalState *, set<Element::Jet<0> *> *);
	Jet(const Jet<0> &);
	Jet<0> & operator = (const Jet<0> &);
	range_type operator () (const domain_type &) const;
	pair<double,double> operator () (const domain_type &, 
		const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &) const;

protected:

	LocalState *LS;
	set<Element::Jet<0> *> *J;
};

//////////////////////////////////////////////////////////////////////
// Class Jet<1>
//////////////////////////////////////////////////////////////////////

template<> class Jet<1> : public Model::Jet<1>
{
public: 

	typedef Model::Static::Jet<2> tangent_type;
	typedef set<Set::Manifold::Point *> domain_type;
	typedef pair<map<Set::Manifold::Point *, Set::VectorSpace::Vector>,
		map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, 
		Set::VectorSpace::Hom> > range_type;

	Jet();
	virtual ~Jet();
	Model::Jet<1> *Clone() const;
	Jet(LocalState *, set<Element::Jet<1> *> *);
	Jet(const Jet<1> &);
	Jet<1> & operator = (const Jet<1> &);
	range_type operator () (const domain_type &) const;
	pair<double,double> operator () (const domain_type &,
		const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &) const;
	void operator () (const domain_type &, Solver::Linear::System<DOF> &) const;
	void operator () (Solver::Linear::System<DOF> &, domain_type &) const;

protected:

	LocalState *LS;
	set<Element::Jet<1> *> *DJ;
};

}

}

#endif // !defined(MODEL_STATIC_H__INCLUED_
