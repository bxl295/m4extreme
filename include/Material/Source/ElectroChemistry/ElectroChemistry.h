// ElectroChemistry.h: interface for the ElectroChemistry class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(MATERIAL_SOURCE_ELECTROCHEMISTRY_H__INCLUDED_)
#define MATERIAL_SOURCE_ELECTROCHEMISTRY_H__INCLUDED_

#pragma once

#include <utility>
#include "Clock/Clock.h"
#include "Material/Material.h"
#include "Set/Algebraic/AlgLib.h"

namespace Material
{
  namespace Source
  {
    namespace ElectroChemistry
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
	Data(const double & Ik_,
	     const double & tmax_,
	     const double & tcycle_,
	     const double & sigma0_,
	     const double & Phif0_,
	     const double & alpha_,
	     const double & a0_,
	     const double & a1_,
	     const double & a2_,
	     const double & a3_,
	     const double & a4_,
	     const double & a5_,
	     const double & a6_,
	     const double & n_,
	     const double & F_);
	virtual ~Data();
	Data(const Data &);
	Data & operator = (const Data &);
	
      private:
	double Ik;
	double tmax;
	double tcycle;
	double sigma0;
	double Phif0;
	double alpha;
	double a0;
	double a1;
	double a2;
	double a3;
	double a4;
	double a5;
	double a6;
	double n;
	double F;
      };

      //////////////////////////////////////////////////////////////////////
      // Class LocalState
      //////////////////////////////////////////////////////////////////////

      class LocalState : public Material::LocalState
      {
	friend class Energy<0>;
	friend class Energy<1>;
	friend class Energy<2>;
	friend class Jet<0>;
	friend class Jet<1>;

      public: 

	typedef Material::Source::ElectroChemistry::Data data_type;
	typedef Material::Source::ElectroChemistry::Energy<0> energy_type;
	typedef Set::VectorSpace::Vector domain_type;

	virtual ~LocalState();
	Material::LocalState *Clone() const;
	LocalState(double, Clock *, Data *);
	LocalState(const LocalState &);
	void operator ++ ();
	void Reset(double);

	void SetPorosity(double);
	double GetEntropy();
	double GetCurrent();
	double GetConductivity();
	double GetVoltage();

      private:
	Data *Properties;
	Clock *Chronos;
	double T;
	double SOC;
	double Phif;
	
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
	class Energy<0> : public Material::Energy<0>
      {
      public: 

	typedef Material::Source::ElectroChemistry::Energy<1> tangent_type;
	typedef Set::VectorSpace::Vector domain_type;
	typedef double range_type;

	virtual ~Energy();
	Material::Energy<0> *Clone() const;
	Energy(LocalState *);
	Energy(const Energy<0> &);
	range_type operator () (const domain_type &) const;

      private:

	LocalState *LS;

      private:

	Energy<0> & operator = (const Energy<0> &);
      };


      //////////////////////////////////////////////////////////////////////
      // Class Jet<p>
      //////////////////////////////////////////////////////////////////////

      template<unsigned int p> class Jet;


    }

  }

}

#endif // !defined(MATERIAL_SOURCE_ELECTROCHEMISTRY_H__INCLUDED_
