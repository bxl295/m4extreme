// Angular.h: interface for the Angular class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(POTENTIAL_ANGULAR__INCLUDED_)
#define POTENTIAL_ANGULAR__INCLUDED_

#pragma once

#include <utility>

using namespace std;

namespace Potential
{
namespace Angular
{
//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

class LocalState
{
public: 

	LocalState();
	virtual ~LocalState();
	virtual LocalState *Clone() const=0;
	virtual void operator ++ ()=0;
};

//////////////////////////////////////////////////////////////////////
// Class Energy<p>
//////////////////////////////////////////////////////////////////////

template<unsigned int> class Energy;

//////////////////////////////////////////////////////////////////////
// Class Energy<0>
//////////////////////////////////////////////////////////////////////

template<> class Energy<0>
{
public: 

	typedef Potential::Angular::Energy<1> tangent_type;
	typedef double domain_type;
	typedef double range_type;

	Energy();
	virtual ~Energy();
	virtual Energy<0> *Clone() const=0;
	virtual const double operator () (const double &) const=0;
};

//////////////////////////////////////////////////////////////////////
// Class Energy<1>
//////////////////////////////////////////////////////////////////////

template<> class Energy<1>
{
public: 

	typedef Potential::Angular::Energy<2> tangent_type;
	typedef double domain_type;
	typedef double range_type;

	Energy();
	virtual ~Energy();
	virtual Energy<1> *Clone() const=0;
	virtual const double operator () (const double &) const=0;
};

//////////////////////////////////////////////////////////////////////
// Class Energy<2>
//////////////////////////////////////////////////////////////////////

template<> class Energy<2>
{
public: 

	typedef Potential::Angular::Energy<2> tangent_type;
	typedef double domain_type;
	typedef double range_type;

	Energy();
	virtual ~Energy();
	virtual Energy<2> *Clone() const=0;
	virtual const double operator () (const double &) const=0;
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

	typedef Potential::Angular::Jet<1> tangent_type;
	typedef double domain_type;
	typedef pair<double, double> range_type;

	Jet();
	virtual ~Jet();
	virtual Jet<0> *Clone() const=0;
	virtual const range_type operator () (const domain_type &)=0;
};

//////////////////////////////////////////////////////////////////////
// Class Jet<1>
//////////////////////////////////////////////////////////////////////

template <> class Jet<1>
{
public: 

	typedef Potential::Angular::Jet<2> tangent_type;
	typedef double domain_type;
	typedef pair<double,double> range_type;

	Jet();
	virtual ~Jet();
	virtual Jet<1> *Clone() const=0;
	virtual const range_type operator () (const domain_type &)=0;
};

}

}

#endif // !defined(POTENTIAL_ANGULAR__INCLUDED_)
