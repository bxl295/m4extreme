// TimeDependent.h: interface for the TimeDependent class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(POTENTIAL_FIELD_TIMEDEPENDENT_H__INCLUDED_)
#define POTENTIAL_FIELD_TIMEDEPENDENT_H__INCLUDED_

#pragma once

#include <utility>
#include "Clock/Clock.h"
#include "Potential/Field/Field.h"
#include "Set/Algebraic/AlgLib.h"

namespace Potential
{
namespace Field
{
namespace TimeDependent	
{
class TimeFunction;
class Data;
class LocalState;
template <unsigned int> class Energy;
template <unsigned int> class Jet;

//////////////////////////////////////////////////////////////////////
// Class Data
//////////////////////////////////////////////////////////////////////
 class TimeFunction {
 public:
 TimeFunction(int Dim_, Clock * Chronos_):
   Chronos(Chronos_), B(Dim_) {}
   
   virtual ~TimeFunction() {};
   
   const Set::VectorSpace::Vector & GetDensity() {
     return B;
   }
   
 protected:
   Clock *Chronos;
   Set::VectorSpace::Vector B;
 };

 class LinearLoading : public TimeFunction {
 public:
 LinearLoading(int Dim_,
	       Clock * Chronos_,
	       double Dt_,
	       const Set::VectorSpace::Vector & B0_,
	       const Set::VectorSpace::Vector & Bf_):
   dB(Dim_), tf( Chronos->Time() + Dt_),
   TimeFunction(Dim_, Chronos_) {
     if (Dt_ > 0.) {
       dB = (Bf_ - B0_) / Dt_;
     }
     
     B = B0_;
   }
   
   virtual ~LinearLoading() {}
   
   void Update() {
     if ( Chronos->Time() < tf ) {
       B += dB * Chronos->DTime();
     }
   }
   
   void Reset( double Dt_,
	       const Set::VectorSpace::Vector & B0_,
	       const Set::VectorSpace::Vector & Bf_ ) {
     if (Dt_ > 0. ) {
       dB = (Bf_ - B0_) / Dt_;
     }
     tf = Chronos->Time() + Dt_;
   }
   
 private:
   double tf;
   Set::VectorSpace::Vector dB;
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
	typedef Potential::Field::TimeDependent::Energy<0> energy_type;
	typedef Set::VectorSpace::Vector domain_type;

	virtual ~LocalState();
	Potential::Field::LocalState *Clone() const;
	LocalState(double, TimeFunction *);
	LocalState(const LocalState &);
	void operator ++ ();
	const double & GetQW() const { return QW; }
	double & GetQW() { return QW; }
	const TimeFunction *GetProp() const { return Prop; }
	TimeFunction* GetProp() { return Prop; }
private:
	double QW;
	TimeFunction *Prop;

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

	typedef Potential::Field::TimeDependent::Energy<1> tangent_type;
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

	typedef Potential::Field::TimeDependent::Energy<2> tangent_type;
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

	typedef Potential::Field::TimeDependent::Energy<3> tangent_type;
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

	typedef Potential::Field::TimeDependent::Jet<1> tangent_type;
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

	typedef Potential::Field::TimeDependent::Jet<2> tangent_type;
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

#endif // !defined(POTENTIAL_FIELD_TIMEDEPENDENT_H__INCLUDED_
