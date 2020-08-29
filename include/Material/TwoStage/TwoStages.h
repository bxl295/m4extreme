// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
////////////////////////////////////////////////////////////////////////////

#if !defined(MATERIAL_TWOSTAGES_INCLUDED_)
#define MATERIAL_TWOSTAGES_INCLUDED_

#pragma once

#include <map>
#include <utility>

#include "Material/Material.h"

namespace Material
{
namespace TwoStages
{

class LocalState;
template<unsigned int> class Energy;
template<unsigned int> class Jet;


// Class LocalState

class LocalState : public Material::LocalState
{
friend class Energy<0>;
friend class Energy<1>;
friend class Energy<2>;
friend class Jet<0>;
friend class Jet<1>;

public: 

	typedef Material::TwoStages::Energy<0> energy_type;
	typedef Set::VectorSpace::Vector domain_type;

	LocalState();
	LocalState(Material::LocalState *,
		   Material::LocalState *);
	virtual ~LocalState();
	Material::LocalState *Clone() const;
	LocalState(const LocalState &);
	void operator ++ ();

private:

	Material::LocalState *LScomp;
	Material::LocalState *LSten;

private:

	LocalState & operator = (const LocalState &);
};


// Class Energy<0>

template<>
class Energy<0> : public Material::Energy<0>
{
public:

        typedef Material::TwoStages::Energy<1> tangent_type;
	typedef Set::VectorSpace::Vector domain_type;
        typedef double range_type;

	Energy();
        virtual ~Energy();
        Material::Energy<0> *Clone() const;
        Energy(LocalState *, 
	       Material::Energy<0> *, 
	       Material::Energy<0> *);
        Energy(const Energy<0> &);
        range_type operator () (const domain_type &) const;

private:

	LocalState * LS;
        Material::Energy<0> * Wcomp;
	Material::Energy<0> * Wten;

private:

        const Energy<0> & operator = (const Energy<0> &);
};

// Class Energy<1>

template<>
class Energy<1> : public Material::Energy<1>
{
public:

        typedef Material::TwoStages::Energy<2> tangent_type;
	typedef Set::VectorSpace::Vector domain_type;
	typedef Set::VectorSpace::Vector range_type;

	Energy();
        virtual ~Energy();
        Material::Energy<1> *Clone() const;
        Energy(LocalState *,
	       Material::Energy<1> *, 
	       Material::Energy<1> *);
        Energy(const Energy<1> &);
        range_type operator () (const domain_type &) const;

private:

	LocalState * LS;
        Material::Energy<1> * DWcomp;
	Material::Energy<1> * DWten;

private:

        const Energy<1> & operator = (const Energy<1> &);
};

// Class Energy<2>

template<>
class Energy<2> : public Material::Energy<2>
{
public:

        typedef Material::TwoStages::Energy<3> tangent_type;
	typedef Set::VectorSpace::Vector domain_type;
	typedef Set::VectorSpace::Hom range_type;

	Energy();
        virtual ~Energy();
        Material::Energy<2> *Clone() const;
        Energy(LocalState *,	       
	       Material::Energy<2> * ,
	       Material::Energy<2> * );
        Energy(const Energy<2> &);
        range_type operator () (const domain_type &) const;

private:

	LocalState * LS;	
	Material::Energy<2> * DDWcomp;
	Material::Energy<2> * DDWten;

private:

        const Energy<2> & operator = (const Energy<2> &);
};

// Class Jet<p>

template<unsigned int p> class Jet;

// Class Jet<0>

template<>
class Jet<0> : public Material::Jet<0>
{
public:

        typedef Material::TwoStages::Jet<1> tangent_type;
	typedef Set::VectorSpace::Vector domain_type;
	typedef pair<double, Set::VectorSpace::Vector> range_type;

	Jet();
        virtual ~Jet();
        Material::Jet<0> *Clone() const;
        Jet(LocalState *,
	    Material::Jet<0> *,	
	    Material::Jet<0> *);
        Jet(const Jet<0> &);
        range_type operator () (const domain_type &) const;

private:

	LocalState * LS;
        Material::Jet<0> * Jcomp;
	Material::Jet<0> * Jten;

private:

        const Jet<0> & operator = (const Jet<0> &);
};

// Class Jet<1>
template <>
class Jet<1> : public Material::Jet<1>
{
public:

        typedef Material::TwoStages::Jet<2> tangent_type;
	typedef Set::VectorSpace::Vector domain_type;
	typedef pair<Set::VectorSpace::Vector, 
		Set::VectorSpace::Hom> range_type;

	Jet();
        virtual ~Jet();
        Material::Jet<1> *Clone() const;
        Jet(LocalState *,
	    Material::Jet<1> *,
	    Material::Jet<1> *);
        Jet(const Jet<1> &);
        range_type operator () (const domain_type &) const;

private:

	LocalState * LS;
        Material::Jet<1> * DJcomp;
	Material::Jet<1> * DJten;

private:

        const Jet<1> & operator = (const Jet<1> &);
};

}

}

#endif // !defined(MATERIAL_TWOSTAGES_INCLUDED_)
