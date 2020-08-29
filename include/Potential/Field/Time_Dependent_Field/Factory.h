// Factory.h: Factory for the TimeDependent class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(POTENTIAL_FIELD_TIMEDEPENDENT_FACTORY_H__INCLUDED_)
#define POTENTIAL_FIELD_TIMEDEPENDENT_FACTORY_H__INCLUDED_

#pragma once

#include <vector>
#include "./TimeDependent.h"
#include "Potential/Field/Factory.h"

namespace Potential
{
  namespace Field
  {
    namespace TimeDependent
    {
      //////////////////////////////////////////////////////////////////////
      // Class Factory
      //////////////////////////////////////////////////////////////////////

      class Factory : public Potential::Field::Factory
      {
      public: 

	typedef Potential::Field::TimeDependent::TimeFunction data_type;

      Factory() : LS(0), W(0), DW(0), DDW(0), J(0), DJ(0) {}
	virtual ~Factory() 
	  {
	    if (LS != 0)  delete LS;
	    if (W != 0)   delete W;
	    if (DW != 0)  delete DW;
	    if (DDW != 0) delete DDW;
	    if (J != 0)   delete J;
	    if (DJ != 0)  delete DJ;
	  }

      Factory(double QW_, Potential::Field::TimeDependent::TimeFunction *Dat_) : 
	LS(0), W(0), DW(0), DDW(0), J(0), DJ(0), Dat(Dat_), QW(QW_) {}

	Potential::Field::LocalState * GetLS()
	  {
	    if (LS == 0) LS = new Potential::Field::TimeDependent::LocalState(QW, Dat);
	    return LS;
	  }

	Potential::Field::Energy<0> * GetW()
	  {
	    if (LS == 0) GetLS();
	    if (W == 0) W = new Potential::Field::TimeDependent::Energy<0>(LS);
	    return W;
	  }

	Potential::Field::Energy<1> * GetDW()
	  {
	    if (LS == 0) GetLS();
	    if (DW == 0) DW = new Potential::Field::TimeDependent::Energy<1>(LS);
	    return DW;
	  }

	Potential::Field::Energy<2> * GetDDW()
	  {
	    if (LS == 0) GetLS();
	    if (DDW == 0) DDW = new Potential::Field::TimeDependent::Energy<2>(LS);
	    return DDW;
	  }

	Potential::Field::Jet<0> * GetJ()
	  {
	    if (LS == 0) GetLS();
	    if (J == 0) J = new Potential::Field::TimeDependent::Jet<0>(LS);
	    return J;
	  }

	Potential::Field::Jet<1> * GetDJ()
	  {
	    if (LS == 0) GetLS();
	    if (DJ == 0) DJ = new Potential::Field::TimeDependent::Jet<1>(LS);
	    return DJ;
	  }

      private:
	double QW;
	Potential::Field::TimeDependent::TimeFunction *Dat;
	Potential::Field::TimeDependent::LocalState *LS;
	Potential::Field::TimeDependent::Energy<0> *W;
	Potential::Field::TimeDependent::Energy<1> *DW;
	Potential::Field::TimeDependent::Energy<2> *DDW;
	Potential::Field::TimeDependent::Jet<0> *J;
	Potential::Field::TimeDependent::Jet<1> *DJ;

      private:

	Factory(const Factory &);
	Factory & operator = (const Factory &);
      };

      //////////////////////////////////////////////////////////////////////
      // Class Builder
      //////////////////////////////////////////////////////////////////////

      class Builder : public Potential::Field::Builder
      {
      public: 

	typedef Potential::Field::TimeDependent::TimeFunction data_type;

	Builder() {}
	virtual ~Builder() {}

      Builder(double QW_, Potential::Field::TimeDependent::TimeFunction *Dat_) : Dat(Dat_), QW(QW_) {}

	Potential::Field::Factory * Build() const
	  {
	    return new Potential::Field::TimeDependent::Factory(QW, Dat);
	  }

      private:
	double QW;
	Potential::Field::TimeDependent::TimeFunction *Dat;
      };

    }

  }

}

#endif // !defined(POTENTIAL_FIELD_TIMEDEPENDENT_FACTORY_H__INCLUDED_
