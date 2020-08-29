// LinearMapping.h: interface for the LinearMapping class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(MANIFOLD_LINEARMAPPING_H__INCLUDED_)
#define MANIFOLD_LINEARMAPPING_H__INCLUDED_

#pragma once

#include <cassert>
#include "../Manifold.h"
#include "../../Algebraic/AlgLib.h"

namespace Set
{
namespace LinearMapping
{
class Point;
class Vector;
class VectorZero;
class Covector;
class CovectorZero;

//////////////////////////////////////////////////////////////////////
// Class Point
//////////////////////////////////////////////////////////////////////

class Point : public Set::Manifold::Point, public Set::VectorSpace::Hom
{
public:

	typedef Set::LinearMapping::Vector vector_type;
	typedef Set::LinearMapping::Covector covector_type;

	Point();
	Point(const unsigned int &, const unsigned int &);
	Point(const unsigned int &, const unsigned int &, double * const &u);
	Point(const unsigned int &);
	Point(const unsigned int &, double * const &u);
	Set::Manifold::Point *Clone() const;
	virtual ~Point();
	Point(const Point &);
	Point & operator = (const Point &);
	Point & operator = (const Set::VectorSpace::Hom &);
	void operator += (const Set::LinearMapping::Vector &);
	void operator -= (const Set::LinearMapping::Vector &);
	void Randomize();
	unsigned int size() const;
	void print(ostream *);
	Set::Manifold::Point & operator = (const Set::Manifold::Point &);
	bool operator != (const Set::Manifold::Point &) const;
	bool operator == (const Set::Manifold::Point &) const;
	void operator += (const Set::VectorSpace::Vector &);
	void operator -= (const Set::VectorSpace::Vector &);
	void operator += (const Set::VectorSpace::Hom &);
	void operator -= (const Set::VectorSpace::Hom &);
};

//////////////////////////////////////////////////////////////////////
// Class Vector
//////////////////////////////////////////////////////////////////////

class Vector : public Set::VectorSpace::Hom, public Set::Manifold::Point
{
public:

	typedef Set::LinearMapping::Covector dual_type;

	Vector();
	Vector(const unsigned int &, const unsigned int &);
	Vector(const unsigned int &, const unsigned int &, double * const &);
	Vector(const unsigned int &);
	Vector(const unsigned int &, double * const &);
	Set::Manifold::Point *Clone() const;
	virtual ~Vector();
	Vector(const Vector &);
	Vector(const Set::VectorSpace::Hom &);
	Vector(const Covector &);
	Vector(const Set::LinearMapping::Point &);
	Vector & operator = (const Set::VectorSpace::Vector &);
	Vector & operator = (const Set::VectorSpace::Hom &);
	Vector & operator = (const Vector &);
	Vector & operator = (const Set::LinearMapping::Point &);
	Set::Manifold::Point & operator = (const Set::Manifold::Point &);
	bool operator != (const Set::Manifold::Point &) const;
	bool operator == (const Set::Manifold::Point &) const;
	void Randomize();
	void operator += (const Set::VectorSpace::Vector &);
	void operator -= (const Set::VectorSpace::Vector &);
	void operator += (const Set::VectorSpace::Hom &);
	void operator -= (const Set::VectorSpace::Hom &);
	unsigned int size() const;
	void print(ostream *);
};

//////////////////////////////////////////////////////////////////////
// Class VectorZero
//////////////////////////////////////////////////////////////////////

class VectorZero : public Set::VectorSpace::HomZero
{
public:

	typedef Set::LinearMapping::CovectorZero dual_type;

	VectorZero();
	VectorZero(const unsigned int &, const unsigned int &);
	VectorZero(const unsigned int &, const unsigned int &, double * const &);
	VectorZero(const unsigned int &);
	VectorZero(const unsigned int &, double * const &);
	virtual ~VectorZero();
	VectorZero(const VectorZero &);
	VectorZero(const Set::VectorSpace::HomZero &);
	VectorZero(const CovectorZero &);
	VectorZero & operator = (const Set::VectorSpace::HomZero &);
	VectorZero & operator = (const VectorZero &);
};

//////////////////////////////////////////////////////////////////////
// Class Covector
//////////////////////////////////////////////////////////////////////

class Covector : public Set::VectorSpace::Hom
{
public:

	typedef Set::LinearMapping::Vector dual_type;

	Covector();
	Covector(const unsigned int &, const unsigned int &);
	Covector(const unsigned int &, const unsigned int &, double * const &);
	Covector(const unsigned int &);
	Covector(const unsigned int &, double * const &);
	virtual ~Covector();
	Covector(const Set::VectorSpace::Hom &);
	Covector(const Covector &);
	Covector(const Set::LinearMapping::Vector &);
	Covector & operator = (const Set::VectorSpace::Vector &);
	Covector & operator = (const Set::VectorSpace::Hom &);
	Covector & operator = (const Covector &);
};

//////////////////////////////////////////////////////////////////////
// Class CovectorZero
//////////////////////////////////////////////////////////////////////

class CovectorZero : public Set::VectorSpace::HomZero
{
public:

	typedef Set::LinearMapping::VectorZero dual_type;

	CovectorZero();
	CovectorZero(const unsigned int &, const unsigned int &);
	CovectorZero(const unsigned int &, const unsigned int &, double * const &);
	CovectorZero(const unsigned int &);
	CovectorZero(const unsigned int &, double * const &);
	virtual ~CovectorZero();
	CovectorZero(const Set::VectorSpace::HomZero &);
	CovectorZero(const CovectorZero &);
	CovectorZero(const Set::LinearMapping::VectorZero &);
	CovectorZero & operator = (const Set::VectorSpace::HomZero &);
	CovectorZero & operator = (const CovectorZero &);
};

}

}

Set::LinearMapping::Point
operator + (const Set::LinearMapping::Point &,
			const Set::LinearMapping::Vector &);
Set::LinearMapping::Point
operator - (const Set::LinearMapping::Point &,
			const Set::LinearMapping::Vector &);
Set::LinearMapping::Vector
operator - (const Set::LinearMapping::Point &,
			const Set::LinearMapping::Point &);
Set::LinearMapping::Point
operator + (const Set::LinearMapping::Point &,
			const Set::VectorSpace::Vector &);
Set::LinearMapping::Point
operator - (const Set::LinearMapping::Point &,
			const Set::VectorSpace::Vector &);
Set::LinearMapping::Point
operator + (const Set::LinearMapping::Point &,
			const Set::VectorSpace::Hom &);
Set::LinearMapping::Point
operator - (const Set::LinearMapping::Point &,
			const Set::VectorSpace::Hom &);
void Random(Set::LinearMapping::Point &);

#endif // !defined(MANIFOLD_LINEARMAPPING_H__INCLUDED_)
