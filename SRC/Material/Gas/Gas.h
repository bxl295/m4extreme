// Gas.h: interface for the Gas class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(MATERIAL_LAW_GAS_H__INCLUDED_)
#define MATERIAL_LAW_GAS_H__INCLUDED_

#pragma once

#include <cassert>
#include "./EoS/EoS.h"
#include "../Material.h"
#include "../../Set/Algebraic/AlgLib.h"

namespace Material
{
namespace Gas
{
class LocalState;
template <unsigned int> class Energy;
template <unsigned int> class Jet;

//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

class LocalState : public Material::LocalState
{
public: 

	typedef Material::Gas::Energy<0> energy_type;
	typedef Set::VectorSpace::Vector domain_type;

	virtual ~LocalState();
	Material::LocalState *Clone() const;
	LocalState(Material::Gas::EoS::LocalState *);
	LocalState(const LocalState &);
	void operator ++ ();
	void Reset(double T) { return EoSLS->Reset(T); }
	Material::Gas::EoS::LocalState * GetMLS() { return EoSLS; }
	double GetTemperature() const { return EoSLS->GetTemperature(); }
        double GetPressure() const { return EoSLS->GetPressure(); }
private:

	Material::Gas::EoS::LocalState *EoSLS;

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

	typedef Material::Gas::Energy<1> tangent_type;
	typedef Set::VectorSpace::Vector domain_type;
	typedef double range_type;

	virtual ~Energy();
	Material::Energy<0> *Clone() const;
	Energy(LocalState *, Material::Gas::EoS::Energy<0> *);
	Energy(const Energy<0> &);
	range_type operator () (const domain_type &) const;
	range_type operator () (const domain_type &, double) const;
private:

	LocalState *LS;
	Material::Gas::EoS::Energy<0> *W;

private:

	const Energy<0> & operator = (const Energy<0> &);
};

//////////////////////////////////////////////////////////////////////
// Class Energy<1>
//////////////////////////////////////////////////////////////////////

template <>
class Energy<1> : public Material::Energy<1>
{
public: 

	typedef Material::Gas::Energy<2> tangent_type;
	typedef Set::VectorSpace::Vector domain_type;
	typedef Set::VectorSpace::Vector range_type;

	virtual ~Energy();
	Material::Energy<1> *Clone() const;
	Energy(LocalState *, Material::Gas::EoS::Energy<1> *);
	Energy(const Energy<1> &);
	range_type operator () (const domain_type &) const;
        double operator () (const domain_type &, double) const;
private:

	LocalState *LS;
	Material::Gas::EoS::Energy<1> *DW;

private:

	const Energy<1> & operator = (const Energy<1> &);
};

//////////////////////////////////////////////////////////////////////
// Class Energy<2>
//////////////////////////////////////////////////////////////////////

template <>
class Energy<2> : public Material::Energy<2>
{
public: 

	typedef Material::Gas::Energy<3> tangent_type;
	typedef Set::VectorSpace::Vector domain_type;
	typedef Set::VectorSpace::Hom range_type;

	virtual ~Energy();
	Material::Energy<2> *Clone() const;
	Energy(LocalState *,
		Material::Gas::EoS::Energy<1> *,
		Material::Gas::EoS::Energy<2> *);
	Energy(const Energy<2> &);
	range_type operator () (const domain_type &) const;
	double operator () (const domain_type &, double) const;
private:

	LocalState *LS;
	Material::Gas::EoS::Energy<1> *DW;
	Material::Gas::EoS::Energy<2> *DDW;

private:

	const Energy<2> & operator = (const Energy<2> &);
};

//////////////////////////////////////////////////////////////////////
// Class Jet<p>
//////////////////////////////////////////////////////////////////////

template<unsigned int p> class Jet;

//////////////////////////////////////////////////////////////////////
// Class Jet<0>
//////////////////////////////////////////////////////////////////////

template<>
class Jet<0> : public Material::Jet<0>
{
public: 

	typedef Material::Gas::Jet<1> tangent_type;
	typedef Set::VectorSpace::Vector domain_type;
	typedef pair<double,Set::VectorSpace::Vector> range_type;

	virtual ~Jet();
	Material::Jet<0> *Clone() const;
	Jet(LocalState *, Material::Gas::EoS::Jet<0> *);
	Jet(const Jet<0> &);
	range_type operator () (const domain_type &) const;

private:

	LocalState *LS;
	Material::Gas::EoS::Jet<0> *V;

private:

	Jet<0> & operator = (const Jet<0> &);
};

//////////////////////////////////////////////////////////////////////
// Class Jet<1>
//////////////////////////////////////////////////////////////////////

template <>
class Jet<1> : public Material::Jet<1>
{
public: 

	typedef Material::Gas::Jet<2> tangent_type;
	typedef Set::VectorSpace::Vector domain_type;
	typedef pair<Set::VectorSpace::Vector,Set::VectorSpace::Hom> range_type;

	virtual ~Jet();
	Material::Jet<1> *Clone() const;
	Jet(LocalState *, Material::Gas::EoS::Jet<1> *);
	Jet(const Jet<1> &);
	range_type operator () (const domain_type &) const;

private:

	LocalState *LS;
	Material::Gas::EoS::Jet<1> *DV;

private:

	Jet<1> & operator = (const Jet<1> &);
};

}

}

#endif // !defined(MATERIAL_LAW_GAS_H__INCLUDED_
