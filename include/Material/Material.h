// Material.h: interface for the Material class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(SET_MATERIAL_H__INCLUDED_)
#define SET_MATERIAL_H__INCLUDED_

#pragma once

#include <iostream>
#include "../Set/Algebraic/VectorSpace/Vector/Vector.h"
#include "../Set/Algebraic/VectorSpace/Category/Category.h"

using namespace std;

namespace Material
{
//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////
    
class Data {
public:
      Data(){}
      virtual ~Data(){}        
};

class LocalState
{
public: 

	typedef Set::VectorSpace::Vector domain_type;

	LocalState(){}
	virtual ~LocalState(){}
	virtual LocalState *Clone() const=0;

	// update local state
	virtual void operator ++ ()=0;

	// reset current deformation
	virtual void Reset(const Set::VectorSpace::Hom &) { return; }

	// reset current temperature
	virtual void Reset(double) { return; }

	// get current temperature
        virtual double GetTemperature() const { return 0.0; }

	// serialization
        virtual void write(ostream & os) const {}
	virtual void read(istream & is) {}
};

//////////////////////////////////////////////////////////////////////
// Class Energy<p>
//////////////////////////////////////////////////////////////////////

template<unsigned int p> class Energy;

//////////////////////////////////////////////////////////////////////
// Class Energy<0>
//////////////////////////////////////////////////////////////////////

template<> class Energy<0>
{
public: 

	typedef Material::Energy<1> tangent_type;
	typedef Set::VectorSpace::Vector domain_type;
	typedef double range_type;

	Energy(){}
	virtual ~Energy(){}
	virtual Energy<0> *Clone() const=0;

	// W(F)
	virtual range_type operator () (const domain_type &) const=0;

	// W(F, T)
	virtual range_type operator () (const domain_type &, double) const {return 0.0;}
};

//////////////////////////////////////////////////////////////////////
// Class Energy<1>
//////////////////////////////////////////////////////////////////////

template <> class Energy<1>
{
public: 

	typedef Material::Energy<2> tangent_type;
	typedef Set::VectorSpace::Vector domain_type;
	typedef Set::VectorSpace::Vector range_type;

	Energy(){}
	virtual ~Energy(){}
	virtual Energy<1> *Clone() const=0;

	// P = dW/dF
	virtual range_type operator () (const domain_type &) const=0;

	// N = dW/dT
	virtual double operator() (const domain_type &, double) const {return 0.0;}
};

//////////////////////////////////////////////////////////////////////
// Class Energy<2>
//////////////////////////////////////////////////////////////////////

template <> class Energy<2>
{
public: 

	typedef Material::Energy<3> tangent_type;
	typedef Set::VectorSpace::Vector domain_type;
	typedef Set::VectorSpace::Hom range_type;

	Energy(){}
	virtual ~Energy(){}
	virtual Energy<2> *Clone() const=0;

	// ddW/dFdF
	virtual range_type operator () (const domain_type &) const=0;

	// ddW/dTdT
	virtual double operator () (const domain_type &, double) const {return 0.0;}
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

	typedef Material::Jet<1> tangent_type;
	typedef Set::VectorSpace::Vector domain_type;
	typedef pair<double,Set::VectorSpace::Vector> range_type;

	Jet(){}
	virtual ~Jet(){}
	virtual Jet<0> *Clone() const=0;
	virtual range_type operator () (const domain_type &) const=0;
};

//////////////////////////////////////////////////////////////////////
// Class Jet<1>
//////////////////////////////////////////////////////////////////////

template <> class Jet<1>
{
public: 

	typedef Material::Jet<2> tangent_type;
	typedef Set::VectorSpace::Vector domain_type;
	typedef pair<Set::VectorSpace::Vector,Set::VectorSpace::Hom> range_type;

	Jet(){}
	virtual ~Jet(){}
	virtual Jet<1> *Clone() const=0;
	virtual range_type operator () (const domain_type &) const=0;
};

}

#endif // !defined(SET_MATERIAL_H__INCLUDED_)
