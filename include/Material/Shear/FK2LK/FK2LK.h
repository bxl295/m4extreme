// FK2LK.h: interface for the FK2LK class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(MATERIAL_SHEAR_FK2LK__INCLUDED_)
#define MATERIAL_SHEAR_FK2LK__INCLUDED_

#pragma once

#include <utility>
#include "../Shear.h"

namespace Material
{
namespace Shear
{
namespace FK2LK
{
class LocalState;
template<unsigned int> class Modulus;
template<unsigned int> class Jet;
template<unsigned int> class JetJet;

//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

class LocalState : public Material::Shear::LocalState
{
friend class Modulus<0>;
friend class Modulus<1>;
friend class Modulus<2>;
friend class Jet<0>;
friend class Jet<1>;
friend class JetJet<0>;

public: 

	typedef Material::Shear::FK2LK::Modulus<0> energy_type;
	typedef Material::Shear::FK2LK::Jet<0> jet_type;
	typedef Material::Shear::FK2LK::JetJet<0> jetjet_type;
	typedef double domain_type;

	virtual ~LocalState();
	Material::Shear::LocalState *Clone() const;
	LocalState(Material::Shear::LocalState *);
	LocalState(const LocalState &);
	void operator ++ ();

private:

	Material::Shear::LocalState *LS;

private:

	const LocalState & operator = (const LocalState &);
};

//////////////////////////////////////////////////////////////////////
// Class Modulus<0>
//////////////////////////////////////////////////////////////////////

template<> 
class Modulus<0> : public Material::Shear::Modulus<0>
{
public: 

	typedef Material::Shear::FK2LK::Modulus<1> tangent_type;
	typedef double domain_type;
	typedef double range_type;

	virtual ~Modulus();
	Material::Shear::Modulus<0> *Clone() const;
	Modulus(Material::Shear::Modulus<0> *);
	Modulus(const Modulus<0> &);
	range_type operator () (const domain_type &) const;

private:

	Material::Shear::Modulus<0> *f;

private:

	const Modulus<0> & operator = (const Modulus<0> &);
};

//////////////////////////////////////////////////////////////////////
// Class Modulus<1>
//////////////////////////////////////////////////////////////////////

template<> 
class Modulus<1> : public Material::Shear::Modulus<1>
{
public: 

	typedef Material::Shear::FK2LK::Modulus<2> tangent_type;
	typedef double domain_type;
	typedef double range_type;

	virtual ~Modulus();
	Material::Shear::Modulus<1> *Clone() const;
	Modulus(Material::Shear::Modulus<1> *);
	Modulus(const Modulus<1> &);
	range_type operator () (const domain_type &) const;

private:

	Material::Shear::Modulus<1> *Df;

private:

	const Modulus<1> & operator = (const Modulus<1> &);
};

//////////////////////////////////////////////////////////////////////
// Class Modulus<2>
//////////////////////////////////////////////////////////////////////

template<> 
class Modulus<2> : public Material::Shear::Modulus<2>
{
public: 

	typedef Material::Shear::FK2LK::Modulus<3> tangent_type;
	typedef double domain_type;
	typedef double range_type;

	virtual ~Modulus();
	Material::Shear::Modulus<2> *Clone() const;
	Modulus(Material::Shear::Jet<1> *);
	Modulus(const Modulus<2> &);
	range_type operator () (const domain_type &) const;

private:

	Material::Shear::Jet<1> *Dg;

private:

	const Modulus<2> & operator = (const Modulus<2> &);
};

//////////////////////////////////////////////////////////////////////
// Class Jet<p>
//////////////////////////////////////////////////////////////////////

template<unsigned int p> class Jet;

//////////////////////////////////////////////////////////////////////
// Class Jet<0>
//////////////////////////////////////////////////////////////////////

template<>
class Jet<0> : public Material::Shear::Jet<0>
{
public: 

	typedef Material::Shear::Jet<1> tangent_type;
	typedef double domain_type;
	typedef pair<double,double> range_type;

	virtual ~Jet();
	Material::Shear::Jet<0> *Clone() const;
	Jet(Material::Shear::Jet<0> *);
	Jet(const Jet<0> &);
	range_type operator () (const domain_type &) const;

private:

	Material::Shear::Jet<0> *g;

private:

	const Jet<0> & operator = (const Jet<0> &);
};

//////////////////////////////////////////////////////////////////////
// Class Jet<1>
//////////////////////////////////////////////////////////////////////

template <>
class Jet<1> : public Material::Shear::Jet<1>
{
public: 

	typedef Material::Shear::Jet<2> tangent_type;
	typedef double domain_type;
	typedef pair<double,double> range_type;

	virtual ~Jet();
	Material::Shear::Jet<1> *Clone() const;
	Jet(Material::Shear::Jet<1> *);
	Jet(const Jet<1> &);
	range_type operator () (const domain_type &) const;

private:

	Material::Shear::Jet<1> *Dg;

private:

	const Jet<1> & operator = (const Jet<1> &);
};

//////////////////////////////////////////////////////////////////////
// Class JetJet<p>
//////////////////////////////////////////////////////////////////////

template<unsigned int p> class JetJet;

//////////////////////////////////////////////////////////////////////
// Class JetJet<0>
//////////////////////////////////////////////////////////////////////

template<>
class JetJet<0> : public Material::Shear::JetJet<0>
{
public: 

	typedef Material::Shear::JetJet<1> tangent_type;
	typedef double domain_type;
	typedef triplet<double,double,double> range_type;

	virtual ~JetJet();
	Material::Shear::JetJet<0> *Clone() const;
	JetJet(Material::Shear::JetJet<0> *);
	JetJet(const JetJet<0> &);
	range_type operator () (const domain_type &) const;

private:

	Material::Shear::JetJet<0> *h;

private:

	const JetJet<0> & operator = (const JetJet<0> &);
};

}

}

}

#endif // !defined(MATERIAL_SHEAR_FK2LK__INCLUDED_)
