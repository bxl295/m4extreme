// Model.h: interface for the Model class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(MODEL_H__INCLUDED_)
#define MODEL_H__INCLUDED_

#pragma once

#include "../Set/Manifold/Manifold.h"
#include "../Set/Manifold/Euclidean/Orthonormal/Orthonormal.h"
#include "../Set/Algebraic/VectorSpace/VectorSpace.h"
#include "../Set/Algebraic/VectorSpace/Category/Category.h"
#include "../Solver/Linear/Linear.h"

typedef pair<Set::Manifold::Point *, unsigned int> DOF;

namespace Model
{
//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

class LocalState
{
public: 

	typedef set<Set::Manifold::Point *> domain_type;

	LocalState(){}
	virtual ~LocalState(){}
	virtual LocalState *Clone() const=0;
	virtual void operator ++ ()=0;
	virtual set<DOF> GetDOFs() const=0;
	virtual set<pair<DOF,DOF> > GetDOFPairs() const=0;
	virtual void operator () (Solver::Linear::System<DOF> &, domain_type &) const=0;
	virtual void operator () (Solver::Linear::System<DOF> &, 
				  map<Set::Manifold::Point *, Set::VectorSpace::Vector> &) const=0;
	virtual void Embed(const domain_type &,
			   map <Set::Manifold::Point *, Set::Euclidean::Orthonormal::Point> &) = 0;
	virtual void Embed(const domain_type &,
			   map <Set::Manifold::Point *, Set::VectorSpace::Vector> &) = 0;
	virtual void Embed(const map <Set::Manifold::Point *, Set::VectorSpace::Vector> &,
			   map <Set::Manifold::Point *, Set::VectorSpace::Vector> &) = 0;
	virtual void Submerge(const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &,
			      map<Set::Manifold::Point *, Set::VectorSpace::Vector> &) = 0;
        virtual void Reset(const set<Set::Manifold::Point *> &) {return;}
	virtual void Embed(const domain_type &,
			   map <Set::Manifold::Point *, double> &) {return;}
};

//////////////////////////////////////////////////////////////////////
// Class Energy<p>
//////////////////////////////////////////////////////////////////////

template<unsigned int p> class Energy;

//////////////////////////////////////////////////////////////////////
// Class Energy<0>
//////////////////////////////////////////////////////////////////////

template<>
class Energy<0>
{
public: 

	typedef Model::Energy<1> tangent_type;
	typedef set<Set::Manifold::Point *> domain_type;
	typedef double range_type;

	Energy(){}
	virtual ~Energy(){}
	virtual Energy<0> *Clone() const=0;
	virtual range_type operator () (const domain_type &) const=0;
};

//////////////////////////////////////////////////////////////////////
// Class Energy<1>
//////////////////////////////////////////////////////////////////////

template<>
class Energy<1>
{
public: 

	typedef Model::Energy<2> tangent_type;
	typedef set<Set::Manifold::Point *> domain_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> range_type;

	Energy(){}
	virtual ~Energy(){}
	virtual Energy<1> *Clone() const=0;
	virtual range_type operator () (const domain_type &) const=0;
        
        virtual void operator () (const domain_type &, range_type &) {
            return;
        };
        
	virtual double operator () (const domain_type &, 
		const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &) const=0;
	virtual void operator () (const domain_type &, Solver::Linear::System<DOF> &) const=0;
};

//////////////////////////////////////////////////////////////////////
// Class Energy<2>
//////////////////////////////////////////////////////////////////////

template<>
class Energy<2>
{
public: 

	typedef Model::Energy<3> tangent_type;
	typedef set<Set::Manifold::Point *> domain_type;
	typedef map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, 
		Set::VectorSpace::Hom> range_type;

	Energy(){}
	virtual ~Energy(){}
	virtual Energy<2> *Clone() const=0;
	virtual range_type operator () (const domain_type &) const=0;
	virtual double operator () (const domain_type &,
		const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &) const=0;
	virtual void operator () (const domain_type &, Solver::Linear::System<DOF> &) const=0;
};

//////////////////////////////////////////////////////////////////////
// Class Jet<p>
//////////////////////////////////////////////////////////////////////

template<unsigned int p> class Jet;

//////////////////////////////////////////////////////////////////////
// Class Jet<0>
//////////////////////////////////////////////////////////////////////

template<> 
class Jet<0>
{
public: 

	typedef Model::Jet<0> tangent_type;
	typedef set<Set::Manifold::Point *> domain_type;
	typedef pair<double, map<Set::Manifold::Point *, 
		Set::VectorSpace::Vector> > range_type;

	Jet(){}
	virtual ~Jet(){}
	virtual Jet<0> *Clone() const=0;
	virtual range_type operator () (const domain_type &) const=0;
	virtual pair<double,double> operator () (const domain_type &,
		const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &) const=0;
};

//////////////////////////////////////////////////////////////////////
// Class Jet<1>
//////////////////////////////////////////////////////////////////////

template<> 
class Jet<1>
{
public: 

	typedef Model::Jet<2> tangent_type;
	typedef set<Set::Manifold::Point *> domain_type;
	typedef pair<map<Set::Manifold::Point *, Set::VectorSpace::Vector>,
		map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, 
		Set::VectorSpace::Hom> > range_type;

	Jet(){}
	virtual ~Jet(){}
	virtual Jet<1> *Clone() const=0;
	virtual range_type operator () (const domain_type &) const=0;
	virtual void operator () (const domain_type &, Solver::Linear::System<DOF> &) const=0;
	virtual void operator () (Solver::Linear::System<DOF> &, domain_type &) const=0;
	virtual pair<double,double> operator () (const domain_type &,
		const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &) const=0;
};

}

#endif // !defined(MODEL_H__INCLUDED_)
