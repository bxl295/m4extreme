// Polymerization.cpp: implementation of the Polymerization class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <cassert>
#include "./Polymerization.h"

namespace Material
{
namespace Source
{
namespace Polymerization
{
//////////////////////////////////////////////////////////////////////
// Class Data
//////////////////////////////////////////////////////////////////////
	Data::Data() : 
	  c0(0.0), c1(0.0), c2(0.0), 
	  Tc(0.0),
	  a1(0.0), b1(0.0), a2(0.0), b2(0.0),
	  m(0.0), n(0.0), R(0.0)
	  {}

	Data::~Data(){}

	Data::Data(const double & rhs_c0,
		   const double & rhs_c1,
		   const double & rhs_c2,
		   const double & rhs_Tc, 
		   const double & rhs_a1,
		   const double & rhs_b1, 
		   const double & rhs_a2, 
		   const double & rhs_b2, 
		   const double & rhs_m, 
		   const double & rhs_n, 
		   const double & rhs_R):
	  c0(rhs_c0), c1(rhs_c1), c2(rhs_c2),
	  Tc(rhs_Tc),
	  a1(rhs_a1), b1(rhs_b1), a2(rhs_a2), b2(rhs_b2),
	  m(rhs_m), n(rhs_n), R(rhs_R) {
	}

	Data::Data(const Data &rhs) :
	  c0(rhs.c0), c1(rhs.c1), c2(rhs.c2),
	  Tc(rhs.Tc),
	  a1(rhs.a1), b1(rhs.b1), a2(rhs.a2), b2(rhs.b2),
	  m(rhs.m), n(rhs.n), R(rhs.R)
	{}

	Data &
	Data::operator = (const Data &rhs)
	{
	  if (this == &rhs) return *this;
	  c0 = rhs.c0; c1 = rhs.c1; c2 = rhs.c2; 
	  Tc = rhs.Tc;
	  a1 = rhs.a1; b1 = rhs.b1; a2 = rhs.a2; b2 = rhs.b2;
	  m = rhs.m; n = rhs.n; R = rhs.R;
	  return *this;
	}

//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

LocalState::~LocalState() {}

Material::LocalState *
LocalState::Clone() const 
{
  return new LocalState(*this);
}

LocalState::LocalState(double T0, Clock *Chronos_, Data *Properties_) : 
  Properties(Properties_), Chronos(Chronos_),T(T0), 
  CureDegree(0.0), CureDegreeOld(0.0), CureRate(0.0){}

LocalState::LocalState(const LocalState &rhs) : 
  Properties(rhs.Properties), Chronos(rhs.Chronos), T(rhs.T),
  CureDegree(rhs.CureDegree), CureDegreeOld(rhs.CureDegreeOld) {}

void 
LocalState::operator ++ () {
  CureDegreeOld = CureDegree;
}

void
LocalState::Reset(double T_) {
  T = T_;
}

void 
LocalState::ComputeCureRate() {
  if ( T >= Properties->Tc ) {
    const double & a1 = Properties->a1;
    const double & b1 = Properties->b1;
    const double & a2 = Properties->a2;
    const double & b2 = Properties->b2;
    const double & m  = Properties->m;
    const double & n  = Properties->n;
    const double & R = Properties->R;
    CureRate =  (a1*exp(-b1/R/T) + a2*exp(-b2/R/T)*pow(CureDegreeOld, m)) * pow(1.0-CureDegreeOld, n);
    if ( CureRate < 0.0 ) { 
      CureRate = 0.0;
    }
    
    CureDegree = CureDegreeOld + Chronos->DTime() * CureRate / 60.0;
    if ( CureDegree > 1.0 ) {
      CureDegree = 1.0;
    }
  }
  
  return;
}

//////////////////////////////////////////////////////////////////////
// Class Energy<0>
//////////////////////////////////////////////////////////////////////

Energy<0>::~Energy() {}

Material::Energy<0> *
Energy<0>::Clone() const
{
  return new Energy<0>(*this);
}

Energy<0>::Energy(LocalState *LS_) : LS(LS_) {}

Energy<0>::Energy(const Energy<0> &rhs) : LS(rhs.LS) {}

double 
Energy<0>::operator () (const Set::VectorSpace::Vector &qp) const
{
  Data * D = LS->Properties;
  LS->ComputeCureRate();
  const double & T = LS->T;
  return -(D->c0 + D->c1 * T + D->c2 * T * T) * LS->CureRate;
}

}

}

}
