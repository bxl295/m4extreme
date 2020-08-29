// Hausdorff.h: interface for the Hausdorff class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(POTENTIAL_FIELD_HAUSDORFF_H__INCLUDED_)
#define POTENTIAL_FIELD_HAUSDORFF_H__INCLUDED_

#pragma once

#include <set>
#include <vector>
#include <utility>
#include "Potential/Field/Field.h"
#include "Set/SetLib.h"

namespace Potential
{
namespace Field
{
namespace Hausdorff	
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
	Data(vector<set<Set::Manifold::Point*> > *,	     
	     map<Set::Manifold::Point*, Set::VectorSpace::Vector> *,
	     vector<Set::VectorSpace::Vector> *,
	     const double &,
	     bool = false);
	Data(const Data &);
	Data & operator = (const Data &);

	void Initialize();
	void ResetBounds();
       
	bool & IsDeformable();
	bool IsDeformable() const;
	bool & IsUpdated();	
	bool IsUpdated() const;

private:
	void _calculateBounds();

private:
	vector<set<Set::Manifold::Point*> > *BC;
	map<Set::Manifold::Point*, Set::VectorSpace::Vector> *yemb;
	vector<Set::VectorSpace::Vector> *N;
	bool isDeformable;
	bool isUpdated;
	double C;
	vector<Set::VectorSpace::Vector*> y;
	vector<double> xmin, xmax, ymin, ymax;
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

	typedef Potential::Field::Hausdorff::Data data_type;
	typedef Potential::Field::Hausdorff::Energy<0> energy_type;
	typedef Potential::Field::Hausdorff::Jet<0> jet_type;
	typedef Set::VectorSpace::Vector domain_type;

	virtual ~LocalState();
	Potential::Field::LocalState *Clone() const;
	LocalState(Data *);
	LocalState(const LocalState &);
	void operator ++ ();
	bool IsOutofBound(int index, double xn, const Set::VectorSpace::Vector & dx) const;
	bool FindNearestSurface(const Set::VectorSpace::Vector &x,
				int & index, double & xn, 
				Set::VectorSpace::Vector & dx) const;
private:
	double _det(const Set::VectorSpace::Vector & u, 
		    const Set::VectorSpace::Vector & v) const;
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

	typedef Potential::Field::Hausdorff::Energy<1> tangent_type;
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

	typedef Potential::Field::Hausdorff::Energy<2> tangent_type;
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

	typedef Potential::Field::Hausdorff::Energy<3> tangent_type;
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

	typedef Potential::Field::Hausdorff::Jet<1> tangent_type;
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

	typedef Potential::Field::Hausdorff::Jet<2> tangent_type;
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

#endif // !defined(POTENTIAL_FIELD_HAUSDORFF_H__INCLUDED_
