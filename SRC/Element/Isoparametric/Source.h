// Source.h: interface for the Isoparametric source class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(POTENTIAL_ISOPARAMETRIC_H__INCLUDED_)
#define POTENTIAL_ISOPARAMETRIC_H__INCLUDED_

#pragma once

#include "../Conforming/Source.h"

using namespace std;

namespace Potential
{
namespace Isoparametric
{
class LocalState;
template <unsigned int> class Energy;
template <unsigned int> class Jet;

//////////////////////////////////////////////////////////////////////
// Class LocalState
//////////////////////////////////////////////////////////////////////

class LocalState : public Potential::Conforming::LocalState
{
friend class Energy<0>;
friend class Energy<1>;
friend class Energy<2>;
friend class Jet<0>;
friend class Jet<1>;

public: 

	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> domain_type;
	typedef vector<Set::VectorSpace::Vector> range_type;

	LocalState();
	virtual ~LocalState();
	LocalState(const vector<Material::LocalState *> &);
	LocalState(
		Element::Interpolation::Shape<0> &,
		Element::Interpolation::Shape<1> &,
		const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &,
		const vector<Material::LocalState *> &,
		const vector<Set::VectorSpace::Vector> &,
		const vector<double> &); 
	LocalState(const LocalState &);

private:

	unsigned int Factorial(unsigned int);
};

//////////////////////////////////////////////////////////////////////
// Class Energy<p>
//////////////////////////////////////////////////////////////////////

template<unsigned int p> class Energy;

//////////////////////////////////////////////////////////////////////
// Class Energy<0>
//////////////////////////////////////////////////////////////////////

template<>
class Energy<0> : public Potential::Conforming::Energy<0>
{
public: 

	typedef Potential::Isoparametric::Energy<1> tangent_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> domain_type;
	typedef double range_type;

	Energy();
	virtual ~Energy();
	Energy(LocalState *, const vector<Material::Energy<0> *> &);
	Energy(const Energy<0> &);
};

//////////////////////////////////////////////////////////////////////
// Class Energy<1>
//////////////////////////////////////////////////////////////////////

template<>
class Energy<1> : public Potential::Conforming::Energy<1>
{
public: 

	typedef Potential::Isoparametric::Energy<2> tangent_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> domain_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> range_type;

	Energy();
	virtual ~Energy();
	Energy(LocalState *, const vector<Material::Energy<1> *> &);
	Energy(const Energy<1> &);
};

//////////////////////////////////////////////////////////////////////
// Class Energy<2>
//////////////////////////////////////////////////////////////////////

template<>
class Energy<2> : public Potential::Conforming::Energy<2>
{
public: 

	typedef Potential::Isoparametric::Energy<3> tangent_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> domain_type;
	typedef map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, 
		Set::VectorSpace::Hom> range_type;

	Energy();
	virtual ~Energy();
	Energy(LocalState *, const vector<Material::Energy<2> *> &);
	Energy(const Energy<2> &);
};

//////////////////////////////////////////////////////////////////////
// Class Jet<p>
//////////////////////////////////////////////////////////////////////

template<unsigned int p> class Jet;

//////////////////////////////////////////////////////////////////////
// Class Jet<0>
//////////////////////////////////////////////////////////////////////

template<> 
class Jet<0> : public Potential::Conforming::Jet<0>
{
public: 

	typedef Potential::Isoparametric::Jet<1> tangent_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> domain_type;
	typedef pair<double, map<Set::Manifold::Point *, Set::VectorSpace::Vector> > range_type;

	Jet();
	virtual ~Jet();
	Jet(LocalState *, const vector<Material::Jet<0> *> &);
	Jet(const Jet<0> &);
};

//////////////////////////////////////////////////////////////////////
// Class Jet<1>
//////////////////////////////////////////////////////////////////////

template <> 
class Jet<1> : public Potential::Conforming::Jet<1>
{
public: 

	typedef Potential::Isoparametric::Jet<2> tangent_type;
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> domain_type;
	typedef pair<
		map<Set::Manifold::Point *, Set::VectorSpace::Vector>, 
		map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, 
		Set::VectorSpace::Hom> > range_type;

	Jet();
	virtual ~Jet();
	Jet(LocalState *, const vector<Material::Jet<1> *> &);
	Jet(const Jet<1> &);
};

}

}

#endif // !defined(POTENTIAL_ISOPARAMETRIC_H__INCLUDED_
