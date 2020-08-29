// Spring.h: interface for the Spring class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(POTENTIAL_ANGULAR_SPRING__INCLUDED_)
#define POTENTIAL_ANGULAR_SPRING__INCLUDED_

#pragma once

#include "../Angular.h"

namespace Potential
{
namespace Angular
{
namespace Spring
{
class Data;
class LocalState;
template<unsigned int> class Energy;
template<unsigned int> class Jet;

//////////////////////////////////////////////////////////////////////
// Class Data
//////////////////////////////////////////////////////////////////////

class Data
{
friend class LocalState;
friend class Energy<0>;
friend class Energy<1>;
friend class Energy<2>;
friend class Jet<0>;
friend class Jet<1>;

public:

	Data();
	virtual ~Data();
	Data(const double &, const double &);
	Data(const Data &);
	Data & operator = (const Data &);
	void Randomize();

private:

	double C0;	// Spring constant at c0
	double c0;	// Cosine of equilibrium angle
};

//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

class LocalState : public Potential::Angular::LocalState
{
friend class Energy<0>;
friend class Energy<1>;
friend class Energy<2>;
friend class Jet<0>;
friend class Jet<1>;

public: 

	LocalState();
	virtual ~LocalState();
	Potential::Angular::LocalState *Clone() const;
	LocalState(const Data *);
	LocalState(const LocalState &);
	LocalState & operator = (const LocalState &);
	void operator ++ ();

private:

	const Data *Prop;
};

//////////////////////////////////////////////////////////////////////
// Class Energy<p>
//////////////////////////////////////////////////////////////////////

template<unsigned int p> class Energy;

//////////////////////////////////////////////////////////////////////
// Class Energy<0>
//////////////////////////////////////////////////////////////////////

template<> 
class Energy<0> : public Potential::Angular::Energy<0>
{
public: 

	typedef Potential::Angular::Spring::Energy<1> tangent_type;
	typedef double domain_type;
	typedef double range_type;

	Energy();
	virtual ~Energy();
	Potential::Angular::Energy<0> *Clone() const;
	Energy(LocalState *);
	Energy(const Energy<0> &);
	const Energy<0> & operator = (const Energy<0> &);
	const range_type operator () (const domain_type &) const;

private:

	LocalState *LS;
};

//////////////////////////////////////////////////////////////////////
// Class Energy<1>
//////////////////////////////////////////////////////////////////////

template<> 
class Energy<1> : public Potential::Angular::Energy<1>
{
public: 

	typedef Potential::Angular::Spring::Energy<2> tangent_type;
	typedef double domain_type;
	typedef double range_type;

	Energy();
	virtual ~Energy();
	Potential::Angular::Energy<1> *Clone() const;
	Energy(LocalState *);
	Energy(const Energy<1> &);
	const Energy<1> & operator = (const Energy<1> &);
	const range_type operator () (const domain_type &) const;

private:

	LocalState *LS;
};

//////////////////////////////////////////////////////////////////////
// Class Energy<2>
//////////////////////////////////////////////////////////////////////

template<> 
class Energy<2> : public Potential::Angular::Energy<2>
{
public: 

	typedef Potential::Angular::Spring::Energy<3> tangent_type;
	typedef double domain_type;
	typedef double range_type;

	Energy();
	virtual ~Energy();
	Potential::Angular::Energy<2> *Clone() const;
	Energy(LocalState *);
	Energy(const Energy<2> &);
	const Energy<2> & operator = (const Energy<2> &);
	const range_type operator () (const domain_type &) const;

private:

	LocalState *LS;
};

//////////////////////////////////////////////////////////////////////
// Class Jet<p>
//////////////////////////////////////////////////////////////////////

template<unsigned int p> class Jet;

//////////////////////////////////////////////////////////////////////
// Class Jet<0>
//////////////////////////////////////////////////////////////////////

template<>
class Jet<0> : public Potential::Angular::Jet<0>
{
public: 

	typedef Potential::Angular::Jet<1> tangent_type;
	typedef double domain_type;
	typedef pair<double,double> range_type;

	Jet();
	virtual ~Jet();
	Potential::Angular::Jet<0> *Clone() const;
	Jet(LocalState *);
	Jet(const Jet<0> &);
	Jet<0> & operator = (const Jet<0> &);
	const range_type operator () (const domain_type &);

private:

	LocalState *LS;
};

//////////////////////////////////////////////////////////////////////
// Class Jet<1>
//////////////////////////////////////////////////////////////////////

template <>
class Jet<1> : public Potential::Angular::Jet<1>
{
public: 

	typedef Potential::Angular::Jet<2> tangent_type;
	typedef double domain_type;
	typedef pair<double,double> range_type;

	Jet();
	virtual ~Jet();
	Potential::Angular::Jet<1> *Clone() const;
	Jet(LocalState *);
	Jet(const Jet<1> &);
	Jet<1> & operator = (const Jet<1> &);
	const range_type operator () (const domain_type &);

private:

	LocalState *LS;
};

}

}

}

#endif // !defined(POTENTIAL_ANGULAR_SPRING__INCLUDED_)
