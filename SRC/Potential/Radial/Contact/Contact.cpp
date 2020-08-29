// Contact.cpp: implementation for the Contact class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <cassert>
#include <iostream>
#include "./Contact.h"

namespace Potential
{
  namespace Radial
  {
    namespace Contact
    {
      //////////////////////////////////////////////////////////////////////
      // Class Data
      //////////////////////////////////////////////////////////////////////

      Data::Data(){}

      Data::Data(Clock *rhs_Chronos, 
		 const double &rhs_kp, const double &rhs_RMin, 
		 const double &rhs_m,  const double &rhs_C) :
	Chronos(rhs_Chronos), kp(rhs_kp), RMin(rhs_RMin), 
	m0(rhs_m), C(rhs_C) {
	assert(m0>0.0 && kp>0.0 && RMin>0.0);
      }
	
      Data::Data(const Data &rhs) :
	kp(rhs.kp), RMin(rhs.RMin), m0(rhs.m0), C(rhs.C), Chronos(rhs.Chronos){}

      Data::~Data(){}

      Data & 
      Data::operator = (const Data &rhs)
      {
	if (this == &rhs) return *this;
	kp = rhs.kp; RMin = rhs.RMin; C = rhs.C; 
	m0 = rhs.m0; Chronos = rhs.Chronos;
	return *this;
      }

      void 
      Data::Randomize()
      {
        kp = (double)rand()/(double)RAND_MAX;
	RMin = (double)rand()/(double)RAND_MAX;
	m0 =  (double)rand()/(double)RAND_MAX;
	C =  (double)rand()/(double)RAND_MAX;
      }

      //////////////////////////////////////////////////////////////////////
      // Class LocalState
      //////////////////////////////////////////////////////////////////////

      LocalState::LocalState() {Prop = 0; mass = 0.0;}

      LocalState::~LocalState() {}

      Potential::Radial::LocalState *
      LocalState::Clone() const 
      {
	return new LocalState(*this);
      }

      LocalState::LocalState(const Data *rhs_Prop)
	: Prop(rhs_Prop), mass(rhs_Prop->m0) {}

      LocalState::LocalState(const LocalState &) {}

      LocalState & 
      LocalState::operator = (const LocalState &rhs)
      {
	if (this == &rhs) return *this; 
	Prop = rhs.Prop; 
	mass = rhs.mass;
	return *this;
      }

      void 
      LocalState::operator ++ () {}
      
      void 
      LocalState::SetMass(double m_) { 
	mass = m_; 
      }

      const double & 
      LocalState::GetMass() const { 
	return mass; 
      }

      //////////////////////////////////////////////////////////////////////
      // Class Energy<0>
      //////////////////////////////////////////////////////////////////////

      Energy<0>::Energy() {LS=0;}

      Energy<0>::~Energy(){}

      Potential::Radial::Energy<0> *
      Energy<0>::Clone() const
      {
	return new Energy<0>(*this);
      }

      Energy<0>::Energy(LocalState *rhs_LS) : LS(rhs_LS){}

      Energy<0>::Energy(const Energy<0> &rhs) : LS(rhs.LS){}

      const Energy<0> & 
      Energy<0>::operator = (const Energy<0> &rhs)
      {
	if (this == &rhs) return *this; 
	LS = rhs.LS; return *this;
      }

      double 
      Energy<0>::operator () (const double &r) const
      {
	assert(false);
	return 0.0;
      }

      //////////////////////////////////////////////////////////////////////
      // Class Energy<1>
      //////////////////////////////////////////////////////////////////////

      Energy<1>::Energy() {LS=0;}

      Energy<1>::~Energy(){}

      Potential::Radial::Energy<1> *
      Energy<1>::Clone() const
      {
	return new Energy<1>(*this);
      }

      Energy<1>::Energy(LocalState *rhs_LS) : LS(rhs_LS){}

      Energy<1>::Energy(const Energy<1> &rhs) : LS(rhs.LS){}

      const Energy<1> & 
      Energy<1>::operator = (const Energy<1> &rhs)
      {
	if (this == &rhs) return *this; 
	LS = rhs.LS; return *this;
      }

      double 
      Energy<1>::operator () (const double &r) const
      {
        const Data *Prop = LS->Prop;
	double x = Prop->RMin - r;
	const double & dt = Prop->Chronos->DTime();
	if ( x <= 0.0 ) {
	  return 0.0;
	}
	else {
	  double F1 = x / dt / dt;
	  double F2 = 2.0 * Prop->C * x / dt / r;
	  return -Prop->kp * LS->mass * (F1>F2?F1:F2);
	}
      }

      //////////////////////////////////////////////////////////////////////
      // Class Energy<2>
      //////////////////////////////////////////////////////////////////////

      Energy<2>::Energy() {LS=0;}

      Energy<2>::~Energy(){}

      Potential::Radial::Energy<2> *
      Energy<2>::Clone() const
      {
	return new Energy<2>(*this);
      }

      Energy<2>::Energy(LocalState *rhs_LS) : LS(rhs_LS){}

      Energy<2>::Energy(const Energy<2> &rhs) : LS(rhs.LS){}

      const Energy<2> & 
      Energy<2>::operator = (const Energy<2> &rhs)
      {
	if (this == &rhs) return *this; 
	LS = rhs.LS; return *this;
      }

      double 
      Energy<2>::operator () (const double &r) const
      {
        const Data *Prop = LS->Prop;
	double x = Prop->RMin - r;
	const double & dt = Prop->Chronos->DTime();
	if ( x <= 0.0 ) {
	  return 0.0;
	}
	else {
	  double dFdr1 = 1.0 / dt / dt;
	  double dFdr2 = 2.0 * Prop->C * Prop->RMin / dt / r / r;
	  return Prop->kp * LS->mass * (dFdr1>dFdr2?dFdr1:dFdr2);
	}
      }

      //////////////////////////////////////////////////////////////////////
      // Class Jet<0>
      //////////////////////////////////////////////////////////////////////

      Jet<0>::Jet() {LS=0;}

      Jet<0>::~Jet() {}

      Potential::Radial::Jet<0> *
      Jet<0>::Clone() const
      {
	return new Jet<0>(*this);
      }

      Jet<0>::Jet(LocalState *rhs_LS) : LS(rhs_LS) {}

      Jet<0>::Jet(const Jet<0> &rhs) : LS(rhs.LS) {}

      Jet<0> & 
      Jet<0>::operator = (const Jet<0> &rhs)
      {
	if (&rhs == this) return *this;
	LS = rhs.LS; return *this;
      }

      pair<double,double>
      Jet<0>::operator () (const double &r) const
      {
	assert(false);
	return make_pair(0.0, 0.0);
      }

      //////////////////////////////////////////////////////////////////////
      // Class Jet<1>
      //////////////////////////////////////////////////////////////////////

      Jet<1>::Jet() {LS=0;}

      Jet<1>::~Jet() {}

      Potential::Radial::Jet<1> *
      Jet<1>::Clone() const
      {
	return new Jet<1>(*this);
      }

      Jet<1>::Jet(LocalState *rhs_LS) : LS(rhs_LS) {}

      Jet<1>::Jet(const Jet<1> &rhs) : LS(rhs.LS) {}

      Jet<1> & 
      Jet<1>::operator = (const Jet<1> &rhs)
      {
	if (&rhs == this) return *this;
	LS = rhs.LS; return *this;
      }

      pair<double,double> 
      Jet<1>::operator () (const double &r) const
      {
        const Data *Prop = LS->Prop;
	double x = Prop->RMin - r;
	const double & dt = Prop->Chronos->DTime();
	if ( x <= 0.0 ) {
	  return make_pair(0.0, 0.0);
	}
	else {
	  double dt2 = dt * dt;
	  double F1 = x / dt2;	  
	  double dFdr1 = 1.0 / dt2;

	  double F2 = 2.0 * Prop->C * x / dt / r;
	  double dFdr2 = 2.0 * Prop->C * Prop->RMin / dt / r / r;

	  double km = Prop->kp * LS->mass;
	  return make_pair(-km * (F1>F2?F1:F2), 
			   km * (dFdr1>dFdr2?dFdr1:dFdr2));
	}
      }

    }

  }

}
