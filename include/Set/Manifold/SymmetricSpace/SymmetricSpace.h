// SymmetricSpace.h: interface for the SymmetricSpace class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(MANIFOLD_SYMMETRICSPACE_H__INCLUDED_)
#define MANIFOLD_SYMMETRICSPACE_H__INCLUDED_

#pragma once

#include <utility>
#include "../Manifold.h"
#include "../Category/Category.h"
#include "../../Algebraic/AlgLib.h"
#include "../../Indexed/Table/Table.h"
#include "../../../Utils/LinearAlgebra/EigenSym/EigenSym.h"
#include "../../../Utils/Triplet/Triplet.h"

using namespace std;

namespace Set
{
namespace SymmetricSpace
{
class Point;
class Vector;
class VectorZero;
class Covector;
class CovectorZero;
template<unsigned int> class ExpMap;
template<unsigned int> class ExpJet;
template<unsigned int> class ExpJetJet;
template<unsigned int> class LogMap;
template<unsigned int> class LogJet;
template<unsigned int> class LogJetJet;

//////////////////////////////////////////////////////////////////////
// Class Point
//////////////////////////////////////////////////////////////////////

class Point : public Set::Manifold::Point, public Set::VectorSpace::Sym
{
public:

	typedef Set::SymmetricSpace::Vector vector_type;
	typedef Set::SymmetricSpace::Covector covector_type;

	Point();
	Point(const unsigned int &);
	Point(const unsigned int &, double * const &u);
	Set::Manifold::Point *Clone() const;
	virtual ~Point();
	Point(const Point &);
	Point(const Set::VectorSpace::Sym &);
	Point(const Set::VectorSpace::Hom &);
	Point & operator = (const Point &);
	Point & operator = (const Set::VectorSpace::Sym &);
	Point & operator = (const Set::VectorSpace::Hom &);
	Set::Manifold::Point & operator = (const Set::Manifold::Point &);
	void operator += (const Set::VectorSpace::Vector &);
	void operator -= (const Set::VectorSpace::Vector &);
	bool operator != (const Set::Manifold::Point &) const;
	bool operator == (const Set::Manifold::Point &) const;
	void Randomize();
	unsigned int size() const;
	void print(ostream *);
};

//////////////////////////////////////////////////////////////////////
// Class Vector
//////////////////////////////////////////////////////////////////////

class Vector : public Set::VectorSpace::Sym, public Set::Manifold::Point
{
public:

	typedef Set::SymmetricSpace::Covector dual_type;

	Vector();
	Vector(const unsigned int &);
	Vector(const unsigned int &, double * const &);
	Set::Manifold::Point *Clone() const;
	virtual ~Vector();
	Vector(const Vector &);
	Vector(const Set::VectorSpace::Sym &);
	Vector & operator = (const Vector &);
	Vector & operator = (const Set::VectorSpace::Sym &);
	Set::Manifold::Point & operator = (const Set::Manifold::Point &);
	void operator += (const Set::VectorSpace::Vector &);
	void operator -= (const Set::VectorSpace::Vector &);
	void operator += (const Set::VectorSpace::Sym &);
	void operator -= (const Set::VectorSpace::Sym &);
	bool operator != (const Set::Manifold::Point &) const;
	bool operator == (const Set::Manifold::Point &) const;
	void Randomize();
	unsigned int size() const;
	void print(ostream *);
};

//////////////////////////////////////////////////////////////////////
// Class VectorZero
//////////////////////////////////////////////////////////////////////

class VectorZero : public Set::VectorSpace::SymZero
{
public:

	typedef Set::SymmetricSpace::CovectorZero dual_type;

	VectorZero();
	VectorZero(const unsigned int &);
	VectorZero(const unsigned int &, double * const &);
	virtual ~VectorZero();
	VectorZero(const VectorZero &);
	VectorZero(const Set::VectorSpace::SymZero &);
	VectorZero(const CovectorZero &);
	VectorZero & operator = (const Set::VectorSpace::SymZero &);
	VectorZero & operator = (const VectorZero &);
};

//////////////////////////////////////////////////////////////////////
// Class Covector
//////////////////////////////////////////////////////////////////////

class Covector : public Set::VectorSpace::Sym
{
public:

	typedef Set::SymmetricSpace::Vector dual_type;

	Covector();
	Covector(const unsigned int &);
	Covector(const unsigned int &, double * const &);
	virtual ~Covector();
	Covector(const Set::VectorSpace::Sym &);
	Covector(const Covector &);
	Covector(const Set::SymmetricSpace::Vector &);
	Covector & operator = (const Set::VectorSpace::Sym &);
	Covector & operator = (const Covector &);
};

//////////////////////////////////////////////////////////////////////
// Class CovectorZero
//////////////////////////////////////////////////////////////////////

class CovectorZero : public Set::VectorSpace::SymZero
{
public:

	typedef Set::SymmetricSpace::VectorZero dual_type;

	CovectorZero();
	CovectorZero(const unsigned int &);
	CovectorZero(const unsigned int &, double * const &);
	virtual ~CovectorZero();
	CovectorZero(const Set::VectorSpace::SymZero &);
	CovectorZero(const CovectorZero &);
	CovectorZero(const Set::SymmetricSpace::VectorZero &);
	CovectorZero & operator = (const Set::VectorSpace::SymZero &);
	CovectorZero & operator = (const CovectorZero &);
};

//////////////////////////////////////////////////////////////////////
// Class ExpMap<0>
//////////////////////////////////////////////////////////////////////

template<> class ExpMap<0>
{
public: 

	typedef Set::SymmetricSpace::ExpMap<1> tangent_type;
	typedef Set::SymmetricSpace::Vector domain_type;
	typedef Set::SymmetricSpace::Point range_type;

	ExpMap();
	virtual ~ExpMap();
	range_type operator () (const domain_type &);
	range_type operator () (
		const Set::VectorSpace::Diagonal &, 
		const Set::VectorSpace::Hom &);
	range_type operator () (const Set::VectorSpace::Vector &);
};

//////////////////////////////////////////////////////////////////////
// Class ExpMap<1>
//////////////////////////////////////////////////////////////////////

template<> class ExpMap<1>
{
public: 

	typedef Set::SymmetricSpace::ExpMap<2> tangent_type;
	typedef Set::SymmetricSpace::Vector domain_type;
	typedef Set::VectorSpace::Hom range_type;

	ExpMap();
	virtual ~ExpMap();
	ExpMap(const ExpMap<0> &);
	range_type operator () (const domain_type &);
	range_type operator () (
		const Set::VectorSpace::Diagonal &, 
		const Set::VectorSpace::Hom &);
	range_type operator () (const Set::VectorSpace::Vector &);
};

//////////////////////////////////////////////////////////////////////
// Class ExpMap<2>
//////////////////////////////////////////////////////////////////////

template<> class ExpMap<2>
{
public: 

	typedef Set::SymmetricSpace::ExpMap<3> tangent_type;
	typedef Set::SymmetricSpace::Vector domain_type;
	typedef Set::VectorSpace::Hom range_type;

	ExpMap();
	virtual ~ExpMap();
	ExpMap(const ExpMap<1> &);
	range_type operator () (const domain_type &);
	range_type operator () (
		const Set::VectorSpace::Diagonal &, 
		const Set::VectorSpace::Hom &);
	range_type operator () (const Set::VectorSpace::Vector &);
};

//////////////////////////////////////////////////////////////////////
// Class ExpJet<0>
//////////////////////////////////////////////////////////////////////

template<> class ExpJet<0>
{
public: 

	typedef Set::SymmetricSpace::ExpJet<1> tangent_type;
	typedef Set::SymmetricSpace::Vector domain_type;
	typedef pair<Set::SymmetricSpace::Point,
		Set::VectorSpace::Hom> range_type;

	ExpJet();
	virtual ~ExpJet();
	range_type operator () (const domain_type &);
	range_type operator () (
		const Set::VectorSpace::Diagonal &, 
		const Set::VectorSpace::Hom &);
	range_type operator () (const Set::VectorSpace::Vector &);
};

//////////////////////////////////////////////////////////////////////
// Class ExpJet<1>
//////////////////////////////////////////////////////////////////////

template<> class ExpJet<1>
{
public: 

	typedef Set::SymmetricSpace::ExpJet<2> tangent_type;
	typedef Set::SymmetricSpace::Vector domain_type;
	typedef pair<Set::VectorSpace::Hom,
		Set::VectorSpace::Hom> range_type;

	ExpJet();
	virtual ~ExpJet();
	ExpJet(const ExpJet<0> &);
	range_type operator () (const domain_type &);
	range_type operator () (
		const Set::VectorSpace::Diagonal &, 
		const Set::VectorSpace::Hom &);
	range_type operator () (const Set::VectorSpace::Vector &);
};

//////////////////////////////////////////////////////////////////////
// Class ExpJetJet<0>
//////////////////////////////////////////////////////////////////////

template<> class ExpJetJet<0>
{
public: 

	typedef Set::SymmetricSpace::ExpJetJet<1> tangent_type;
	typedef Set::SymmetricSpace::Vector domain_type;
	typedef triplet<Set::SymmetricSpace::Point,
		Set::VectorSpace::Hom,
		Set::VectorSpace::Hom> range_type;

	ExpJetJet();
	virtual ~ExpJetJet();
	range_type operator () (const domain_type &);
	range_type operator () (
		const Set::VectorSpace::Diagonal &, 
		const Set::VectorSpace::Hom &);
	range_type operator () (const Set::VectorSpace::Vector &);
};

//////////////////////////////////////////////////////////////////////
// Class LogMap<0>
//////////////////////////////////////////////////////////////////////

template<> class LogMap<0>
{
public: 

	typedef Set::SymmetricSpace::LogMap<1> tangent_type;
	typedef Set::SymmetricSpace::Point domain_type;
	typedef Set::VectorSpace::Sym range_type;

	LogMap();
	virtual ~LogMap();
	range_type operator () (const domain_type &);
	range_type operator () (
		const Set::VectorSpace::Diagonal &, 
		const Set::VectorSpace::Hom &);
	range_type operator () (const Set::VectorSpace::Vector &);
};

//////////////////////////////////////////////////////////////////////
// Class LogMap<1>
//////////////////////////////////////////////////////////////////////

template<> class LogMap<1>
{
public: 

	typedef Set::SymmetricSpace::LogMap<2> tangent_type;
	typedef Set::SymmetricSpace::Point domain_type;
	typedef Set::VectorSpace::Hom range_type;

	LogMap();
	virtual ~LogMap();
	LogMap(const LogMap<0> &);
	range_type operator () (const domain_type &);
	range_type operator () (
		const Set::VectorSpace::Diagonal &, 
		const Set::VectorSpace::Hom &);
	range_type operator () (const Set::VectorSpace::Vector &);
 private:
	range_type operatorGeneral(
		const Set::VectorSpace::Diagonal &, 
		const Set::VectorSpace::Hom &);
	range_type operator3D(
		const Set::VectorSpace::Diagonal &, 
		const Set::VectorSpace::Hom &);
};

//////////////////////////////////////////////////////////////////////
// Class LogMap<2>
//////////////////////////////////////////////////////////////////////

template<> class LogMap<2>
{
public: 

	typedef Set::SymmetricSpace::LogMap<3> tangent_type;
	typedef Set::SymmetricSpace::Point domain_type;
	typedef Set::VectorSpace::Hom range_type;

	LogMap();
	virtual ~LogMap();
	LogMap(const LogMap<1> &);
	range_type operator () (const domain_type &);
	range_type operator () (
		const Set::VectorSpace::Diagonal &, 
		const Set::VectorSpace::Hom &);
	range_type operator () (const Set::VectorSpace::Vector &);
 private:
	range_type operatorGeneral(
		const Set::VectorSpace::Diagonal &, 
		const Set::VectorSpace::Hom &);
	range_type operator3D(
		const Set::VectorSpace::Diagonal &, 
		const Set::VectorSpace::Hom &);
};

//////////////////////////////////////////////////////////////////////
// Class LogJet<0>
//////////////////////////////////////////////////////////////////////

template<> class LogJet<0>
{
public: 

	typedef Set::SymmetricSpace::LogJet<1> tangent_type;
	typedef Set::SymmetricSpace::Point domain_type;
	typedef pair<Set::VectorSpace::Sym,
		Set::VectorSpace::Hom> range_type;

	LogJet();
	virtual ~LogJet();
	range_type operator () (const domain_type &);
	range_type operator () (
		const Set::VectorSpace::Diagonal &, 
		const Set::VectorSpace::Hom &);
	range_type operator () (const Set::VectorSpace::Vector &);
};

//////////////////////////////////////////////////////////////////////
// Class LogJet<1>
//////////////////////////////////////////////////////////////////////

template<> class LogJet<1>
{
public: 

	typedef Set::SymmetricSpace::LogJet<2> tangent_type;
	typedef Set::SymmetricSpace::Point domain_type;
	typedef pair<Set::VectorSpace::Hom,
		Set::VectorSpace::Hom> range_type;

	LogJet();
	virtual ~LogJet();
	LogJet(const LogJet<0> &);
	range_type operator () (const domain_type &);
	range_type operator () (
		const Set::VectorSpace::Diagonal &, 
		const Set::VectorSpace::Hom &);
	range_type operator () (const Set::VectorSpace::Vector &);
};

//////////////////////////////////////////////////////////////////////
// Class LogJetJet<0>
//////////////////////////////////////////////////////////////////////

template<> class LogJetJet<0>
{
public: 

	typedef Set::SymmetricSpace::LogJetJet<1> tangent_type;
	typedef Set::SymmetricSpace::Point domain_type;
	typedef triplet<Set::VectorSpace::Sym,
		Set::VectorSpace::Hom,
		Set::VectorSpace::Hom> range_type;

	LogJetJet();
	virtual ~LogJetJet();
	range_type operator () (const domain_type &);
	range_type operator () (
		const Set::VectorSpace::Diagonal &, 
		const Set::VectorSpace::Hom &);
	range_type operator () (const Set::VectorSpace::Vector &);
};

}

}

Set::SymmetricSpace::Point
operator + (const Set::SymmetricSpace::Point &,
			const Set::SymmetricSpace::Vector &);
Set::SymmetricSpace::Point
operator - (const Set::SymmetricSpace::Point &,
			const Set::SymmetricSpace::Vector &);
Set::SymmetricSpace::Vector
operator - (const Set::SymmetricSpace::Point &,
			const Set::SymmetricSpace::Point &);
Set::SymmetricSpace::Point
operator + (const Set::SymmetricSpace::Point &,
			const Set::VectorSpace::Vector &);
Set::SymmetricSpace::Point
operator - (const Set::SymmetricSpace::Point &,
			const Set::VectorSpace::Vector &);
void Random(Set::SymmetricSpace::Point &);	
//
Set::SymmetricSpace::Vector
operator + (const Set::SymmetricSpace::Vector &,
			const Set::VectorSpace::Vector &);
Set::SymmetricSpace::Vector
operator - (const Set::SymmetricSpace::Vector &,
			const Set::VectorSpace::Vector &);
//
Set::SymmetricSpace::Vector
operator + (const Set::SymmetricSpace::Vector &,
			const Set::VectorSpace::Sym &);
Set::SymmetricSpace::Vector
operator - (const Set::SymmetricSpace::Vector &,
			const Set::VectorSpace::Sym &);
//
Set::SymmetricSpace::Vector
operator + (const Set::SymmetricSpace::Vector &,
			const Set::VectorSpace::SymZero &);
Set::SymmetricSpace::Vector
operator - (const Set::SymmetricSpace::Vector &,
			const Set::VectorSpace::SymZero &);
//
Set::SymmetricSpace::Vector
operator + (const Set::SymmetricSpace::Vector &,
			const Set::SymmetricSpace::Vector &);
Set::SymmetricSpace::Vector
operator - (const Set::SymmetricSpace::Vector &,
			const Set::SymmetricSpace::Vector &);

#endif // !defined(MANIFOLD_SYMMETRICSPACE_H__INCLUDED_)
