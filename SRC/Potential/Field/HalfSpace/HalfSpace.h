// HalfSpace.h: interface for the HalfSpace class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(POTENTIAL_FIELD_HALFSPACE_H__INCLUDED_)
#define POTENTIAL_FIELD_HALFSPACE_H__INCLUDED_

#pragma once

#include <utility>
#include "../Field.h"
#include "../../../Set/Algebraic/AlgLib.h"

namespace Potential
{
namespace Field
{
namespace HalfSpace	
{
class Data;
class LocalState;
template <unsigned int> class Energy;
template <unsigned int> class Jet;

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
	Data(const Set::VectorSpace::Vector &, 
	     const Set::VectorSpace::Vector &, const double &);
	Data(const Set::VectorSpace::Vector &,
	     const Set::VectorSpace::Vector &, 
	     const Set::VectorSpace::Vector &, const double &);
	Data(const vector<Set::VectorSpace::Vector *> &,	     
	     const Set::VectorSpace::Vector &, const double &);
	Data(const Data &);
	Data & operator = (const Data &);
	void Randomize();

        Set::VectorSpace::Vector & GetY() { return y; }
        Set::VectorSpace::Vector & GetN() { return N; }
        double & GetC() { return C; }

private:

	Set::VectorSpace::Vector y;
	Set::VectorSpace::Vector N;
	double C;

	bool Finite;
	double y0_max, y0_min, y1_max, y1_min;
	
	Set::VectorSpace::Vector anchor;
	vector<Set::VectorSpace::Vector> vertices;
	vector<double> theta;	
};

//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

class LocalState : public Potential::Field::LocalState
{
friend class Energy<0>;
friend class Energy<1>;
friend class Energy<2>;
friend class Jet<0>;
friend class Jet<1>;

public: 

	typedef Potential::Field::HalfSpace::Data data_type;
	typedef Potential::Field::HalfSpace::Energy<0> energy_type;
	typedef Potential::Field::HalfSpace::Jet<0> jet_type;
	typedef Set::VectorSpace::Vector domain_type;

	virtual ~LocalState();
	Potential::Field::LocalState *Clone() const;
	LocalState(Data *);
	LocalState(const LocalState &);
	void operator ++ ();

	bool IsOutofBound(double xn, const Set::VectorSpace::Vector & dx) const;

private:

	Data *Prop;

private:

	LocalState & operator = (const LocalState &);
};

//////////////////////////////////////////////////////////////////////
// Class Energy<p>
//////////////////////////////////////////////////////////////////////

template<unsigned int p> class Energy;

//////////////////////////////////////////////////////////////////////
// Class Energy<0>
//////////////////////////////////////////////////////////////////////

template<>
class Energy<0> : public Potential::Field::Energy<0>
{
public: 

	typedef Potential::Field::HalfSpace::Energy<1> tangent_type;
	typedef Set::VectorSpace::Vector domain_type;
	typedef double range_type;

	virtual ~Energy();
	Potential::Field::Energy<0> *Clone() const;
	Energy(LocalState *);
	Energy(const Energy<0> &);
	range_type operator () (const domain_type &) const;

private:

	LocalState *LS;

private:

	Energy<0> & operator = (const Energy<0> &);
};

//////////////////////////////////////////////////////////////////////
// Class Energy<1>
//////////////////////////////////////////////////////////////////////

template <>
class Energy<1> : public Potential::Field::Energy<1>
{
public: 

	typedef Potential::Field::HalfSpace::Energy<2> tangent_type;
	typedef Set::VectorSpace::Vector domain_type;
	typedef Set::VectorSpace::Vector range_type;

	virtual ~Energy();
	Potential::Field::Energy<1> *Clone() const;
	Energy(LocalState *);
	Energy(const Energy<1> &);
	range_type operator () (const domain_type &) const;

private:

	LocalState *LS;

private:

	Energy<1> & operator = (const Energy<1> &);
};

//////////////////////////////////////////////////////////////////////
// Class Energy<2>
//////////////////////////////////////////////////////////////////////

template <>
class Energy<2> : public Potential::Field::Energy<2>
{
public: 

	typedef Potential::Field::HalfSpace::Energy<3> tangent_type;
	typedef Set::VectorSpace::Vector domain_type;
	typedef Set::VectorSpace::Hom range_type;

	virtual ~Energy();
	Potential::Field::Energy<2> *Clone() const;
	Energy(LocalState *);
	Energy(const Energy<2> &);
	range_type operator () (const domain_type &) const;

private:

	LocalState *LS;

private:

	Energy<2> & operator = (const Energy<2> &);
};

//////////////////////////////////////////////////////////////////////
// Class Jet<p>
//////////////////////////////////////////////////////////////////////

template<unsigned int p> class Jet;

//////////////////////////////////////////////////////////////////////
// Class Jet<0>
//////////////////////////////////////////////////////////////////////

template<>
class Jet<0> : public Potential::Field::Jet<0>
{
public: 

	typedef Potential::Field::HalfSpace::Jet<1> tangent_type;
	typedef Set::VectorSpace::Vector domain_type;
	typedef pair<double,Set::VectorSpace::Vector> range_type;

	virtual ~Jet();
	Potential::Field::Jet<0> *Clone() const;
	Jet(LocalState *);
	Jet(const Jet<0> &);
	range_type operator () (const domain_type &) const;

private:

	LocalState *LS;

private:

	Jet<0> & operator = (const Jet<0> &);
};

//////////////////////////////////////////////////////////////////////
// Class Jet<1>
//////////////////////////////////////////////////////////////////////

template <>
class Jet<1> : public Potential::Field::Jet<1>
{
public: 

	typedef Potential::Field::HalfSpace::Jet<2> tangent_type;
	typedef Set::VectorSpace::Vector domain_type;
	typedef pair<Set::VectorSpace::Vector,Set::VectorSpace::Hom> range_type;

	virtual ~Jet();
	Potential::Field::Jet<1> *Clone() const;
	Jet(LocalState *);
	Jet(const Jet<1> &);
	range_type operator () (const domain_type &) const;

private:

	LocalState *LS;

private:

	Jet<1> & operator = (const Jet<1> &);
};

}

}

}

#endif // !defined(POTENTIAL_FIELD_HALFSPACE_H__INCLUDED_
