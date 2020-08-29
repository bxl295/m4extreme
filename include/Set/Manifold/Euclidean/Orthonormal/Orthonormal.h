// Orthonormal.h: interface for the Orthonormal class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(MANIFOLD_EUCLIDEAN_ORTHONORMAL_H__INCLUDED_)
#define MANIFOLD_EUCLIDEAN_ORTHONORMAL_H__INCLUDED_

#pragma once

#include "../../Manifold.h"
#include "../../../Algebraic/AlgLib.h"

namespace Set
{
namespace Euclidean
{	
namespace Orthonormal
{
class Point;
class Vector;
class VectorZero;
class Covector;
class CovectorZero;

//////////////////////////////////////////////////////////////////////
// Class Point
//////////////////////////////////////////////////////////////////////

class Point : public Set::Manifold::Point, public Set::Array
{
public:

	typedef Set::Euclidean::Orthonormal::Vector vector_type;
	typedef Set::Euclidean::Orthonormal::Covector covector_type;

	Point();
	Point(const unsigned int &);
	Point(const unsigned int &, double * const &u);
	Set::Manifold::Point *Clone() const;
	virtual ~Point();
	Point(const Point &);
	Point & operator = (const Point &);
	void operator += (const Vector &);
	void operator -= (const Vector &);
	void operator *= (const double &);
	void operator /= (const double &);
	void Randomize();
	unsigned int size() const;
	void print(ostream *);
	Set::Manifold::Point & operator = (const Set::Manifold::Point &);
	bool operator != (const Set::Manifold::Point &) const;
	bool operator == (const Set::Manifold::Point &) const;
	void operator += (const Set::VectorSpace::Vector &);
	void operator -= (const Set::VectorSpace::Vector &);
};

//////////////////////////////////////////////////////////////////////
// Class Vector
//////////////////////////////////////////////////////////////////////

class Vector : public Set::VectorSpace::Vector, public Set::Manifold::Point
{
public:

	typedef Set::Euclidean::Orthonormal::Covector dual_type;

	Vector();
	Vector(const unsigned int &);
	Vector(const unsigned int &, double * const &);
	Set::Manifold::Point *Clone() const;
	virtual ~Vector();
	Vector(const Vector &);
	Vector(const Set::VectorSpace::Vector &);
	Vector(const Covector &);
	Vector(const Set::Euclidean::Orthonormal::Point &);
	Vector & operator = (const Set::VectorSpace::Vector &);
	Vector & operator = (const Vector &);
	Vector & operator = (const Set::Euclidean::Orthonormal::Point &);
	Set::Manifold::Point & operator = (const Set::Manifold::Point &);
	bool operator != (const Set::Manifold::Point &) const;
	bool operator == (const Set::Manifold::Point &) const;
	void Randomize();
	void operator += (const Set::VectorSpace::Vector &);
	void operator -= (const Set::VectorSpace::Vector &);
	unsigned int size() const;
	void print(ostream *);
};

//////////////////////////////////////////////////////////////////////
// Class VectorZero
//////////////////////////////////////////////////////////////////////

class VectorZero : public Set::VectorSpace::VectorZero
{
public:

	typedef Set::Euclidean::Orthonormal::CovectorZero dual_type;

	VectorZero();
	VectorZero(const unsigned int &);
	VectorZero(const unsigned int &, double * const &);
	virtual ~VectorZero();
	VectorZero(const VectorZero &);
	VectorZero(const Set::VectorSpace::VectorZero &);
	VectorZero(const CovectorZero &);
	VectorZero & operator = (const Set::VectorSpace::VectorZero &);
	VectorZero & operator = (const VectorZero &);
};

//////////////////////////////////////////////////////////////////////
// Class Covector
//////////////////////////////////////////////////////////////////////

class Covector : public Set::VectorSpace::Vector
{
public:

	typedef Set::Euclidean::Orthonormal::Vector dual_type;

	Covector();
	Covector(const unsigned int &);
	Covector(const unsigned int &, double * const &);
	virtual ~Covector();
	Covector(const Set::VectorSpace::Vector &);
	Covector(const Covector &);
	Covector(const Set::Euclidean::Orthonormal::Vector &);
	Covector & operator = (const Set::VectorSpace::Vector &);
	Covector & operator = (const Covector &);
};

//////////////////////////////////////////////////////////////////////
// Class CovectorZero
//////////////////////////////////////////////////////////////////////

class CovectorZero : public Set::VectorSpace::VectorZero
{
public:

	typedef Set::Euclidean::Orthonormal::VectorZero dual_type;

	CovectorZero();
	CovectorZero(const unsigned int &);
	CovectorZero(const unsigned int &, double * const &);
	virtual ~CovectorZero();
	CovectorZero(const Set::VectorSpace::VectorZero &);
	CovectorZero(const CovectorZero &);
	CovectorZero(const Set::Euclidean::Orthonormal::VectorZero &);
	CovectorZero & operator = (const Set::VectorSpace::VectorZero &);
	CovectorZero & operator = (const CovectorZero &);
};

}

}

}

Set::Euclidean::Orthonormal::Point
operator + (const Set::Euclidean::Orthonormal::Point &,
			const Set::Euclidean::Orthonormal::Vector &);
Set::Euclidean::Orthonormal::Point
operator - (const Set::Euclidean::Orthonormal::Point &,
			const Set::Euclidean::Orthonormal::Vector &);
Set::Euclidean::Orthonormal::Vector
operator - (const Set::Euclidean::Orthonormal::Point &,
			const Set::Euclidean::Orthonormal::Point &);
Set::Euclidean::Orthonormal::Point
operator + (const Set::Euclidean::Orthonormal::Point &,
			const Set::VectorSpace::Vector &);
Set::Euclidean::Orthonormal::Point
operator - (const Set::Euclidean::Orthonormal::Point &,
			const Set::VectorSpace::Vector &);
void Random(Set::Euclidean::Orthonormal::Point &);

#endif // !defined(MANIFOLD_EUCLIDEAN_ORTHONORMAL_H__INCLUDED_)
