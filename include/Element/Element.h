// Element.h: interface for the Element class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(ELEMENT_H__INCLUDED_)
#define ELEMENT_H__INCLUDED_

#pragma once

#include <set>
#include <map>
#include <cassert>
#include "../Set/Algebraic/AlgLib.h"
#include "../Set/Manifold/Manifold.h"

using namespace std;

namespace Element
{
//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

class LocalState
{
public: 
       typedef Set::Manifold::Point dof_type;
       typedef Set::VectorSpace::Vector vector_type;
       typedef map<dof_type*, vector_type> vectorset_type;
       typedef map<dof_type*, double> scalarset_type;
	
	LocalState(){
            _isActivated = true;
            _isStretched   = false;
        }
	virtual ~LocalState(){}
	virtual LocalState *Clone() const=0;
	virtual void operator ++ ()=0;
	virtual set<Set::Manifold::Point *> GetNodes() const=0;

        virtual void GetNodes(set<Set::Manifold::Point *> &) const {
            assert(false);
        }

	virtual void GetNodePairs(set<pair<Set::Manifold::Point *, 
				  Set::Manifold::Point *> > &) const {
	  assert(false);
	}

	virtual set<pair<Set::Manifold::Point *, 
		Set::Manifold::Point *> > GetNodePairs() const=0;
	virtual void Reset(const vectorset_type &) {}
	virtual void Reset(const scalarset_type &) {}

        virtual bool isActivated() const { return _isActivated; }
        virtual bool & isActivated() { return _isActivated; }
        bool isStretched() const { return _isStretched; }
        bool & isStretched() { return _isStretched; }

	virtual void write(ostream & os, const std::map<dof_type*, int> & idmap) const {
	  os.write((char*)&_isActivated, sizeof(bool));
	  os.write((char*)&_isStretched, sizeof(bool));
	}

	virtual void read(istream & is, const std::map<int, dof_type*> & idmap, 
			  const std::map<dof_type*, vector_type> & x) {
	  is.read((char*)&_isActivated, sizeof(bool));
	  is.read((char*)&_isStretched, sizeof(bool));
	}
       
protected:
        bool _isActivated;
        bool _isStretched;
       
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

	typedef Element::Energy<1> tangent_type;
	typedef map<Set::Manifold::Point *, 
		Set::VectorSpace::Vector> domain_type;
	typedef double range_type;

	Energy(){}
	virtual ~Energy(){}
	virtual Energy<0> *Clone() const=0;
	virtual LocalState *GetLocalState() const=0;
	virtual range_type operator () (const domain_type &) const=0;
};

//////////////////////////////////////////////////////////////////////
// Class Energy<1>
//////////////////////////////////////////////////////////////////////

template<>
class Energy<1>
{
public: 

	typedef Element::Energy<2> tangent_type;
	typedef map<Set::Manifold::Point *, 
		Set::VectorSpace::Vector> domain_type;
	typedef map<Set::Manifold::Point *, 
		Set::VectorSpace::Vector> range_type;

	typedef map<Set::Manifold::Point *, double> scalar_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> vector_type;
	typedef map<Set::Manifold::Point *, double> scalarset_type;

	Energy():_ELS(NULL){}
	virtual ~Energy(){}
	virtual Energy<1> *Clone() const=0;
	virtual LocalState *GetLocalState() const=0;
	virtual range_type operator () (const domain_type &) const {
	  assert(false);
	}
        
        virtual void operator () (const domain_type &, range_type &) const {
            assert(false);
        };

        virtual void operator () (const vector_type &, const scalar_type &, scalarset_type &) const {
	  return;
	}

public:
        LocalState * _ELS;
};

//////////////////////////////////////////////////////////////////////
// Class Energy<2>
//////////////////////////////////////////////////////////////////////

template<>
class Energy<2>
{
public: 

	typedef Element::Energy<3> tangent_type;
	typedef map<Set::Manifold::Point *, 
		Set::VectorSpace::Vector> domain_type;
	typedef map<pair<Set::Manifold::Point *, 
		Set::Manifold::Point *>, 
		Set::VectorSpace::Hom> range_type;

	typedef map<Set::Manifold::Point *, double> scalar_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> vector_type;
	typedef map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, double> scalarset_type;

	Energy(){}
	virtual ~Energy(){}
	virtual Energy<2> *Clone() const=0;
	virtual LocalState *GetLocalState() const=0;
	virtual range_type operator () (const domain_type &) const {
	  assert(false);
	}

        virtual void operator () (const domain_type &, range_type &) const {
            assert(false);
        };

        virtual void operator () (const vector_type &, const scalar_type &, scalarset_type &) const {
	  return;
	}

};

//////////////////////////////////////////////////////////////////////
// Class Jet<p>
//////////////////////////////////////////////////////////////////////

template<unsigned int p> class Jet;

//////////////////////////////////////////////////////////////////////
// Class Jet<0>
//////////////////////////////////////////////////////////////////////

template<> class Jet<0>
{
public: 

	typedef Element::Jet<1> tangent_type;
	typedef map<Set::Manifold::Point *, 
		Set::VectorSpace::Vector> domain_type;
	typedef pair<double, map<Set::Manifold::Point *, 
		Set::VectorSpace::Vector> > range_type;

	Jet(){}
	virtual ~Jet(){}
	virtual Jet<0> *Clone() const=0;
	virtual LocalState *GetLocalState() const=0;
	virtual range_type operator () (const domain_type &) const=0;
};

//////////////////////////////////////////////////////////////////////
// Class Jet<1>
//////////////////////////////////////////////////////////////////////

template <> class Jet<1>
{
public: 

	typedef Element::Jet<2> tangent_type;
	typedef map<Set::Manifold::Point *, 
		Set::VectorSpace::Vector> domain_type;
	typedef pair<
		map<Set::Manifold::Point *, Set::VectorSpace::Vector>, 
		map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, 
		Set::VectorSpace::Hom> > range_type;

	Jet(){}
	virtual ~Jet(){}
	virtual Jet<1> *Clone() const=0;
	virtual LocalState *GetLocalState() const=0;
	virtual range_type operator () (const domain_type &) const=0;

	virtual void operator()(const domain_type &, range_type &) const {
	  assert(false);
	}
};

}

#endif // !defined(ELEMENT_H__INCLUDED_
