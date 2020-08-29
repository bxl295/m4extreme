// Polar.h: interface for the Polar class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(MANIFOLD_EUCLIDEAN_POLAR_H__INCLUDED_)
#define MANIFOLD_EUCLIDEAN_POLAR_H__INCLUDED_

#pragma once

#include <cassert>
#include "../Orthonormal/Orthonormal.h"
#include "../../Manifold.h"
#include "../../Category/Category.h"
#include "../../../Algebraic/AlgLib.h"

namespace Set
{
namespace Euclidean
{	
namespace Polar
{
#define PI2 6.2831853071795862

class Point;
class Vector;
class VectorZero;
class Covector;
class CovectorZero;
template<unsigned int> class Embedding;
template<unsigned int> class Submersion;

//////////////////////////////////////////////////////////////////////
// Class Point
//////////////////////////////////////////////////////////////////////

class Point : public Set::Manifold::Point, public Set::Array
{
public:

	typedef Set::Euclidean::Polar::Vector vector_type;
	typedef Set::Euclidean::Polar::Covector covector_type;

	Point();
	Point(double * const &u);
	Set::Manifold::Point *Clone() const;
	virtual ~Point();
	Point(const Point &);
	Point & operator = (const Point &);
	void operator += (const Vector &);
	void operator -= (const Vector &);
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

	typedef Set::Euclidean::Polar::Covector dual_type;

	Vector();
	Vector(double * const &);
	Set::Manifold::Point *Clone() const;
	virtual ~Vector();
	Vector(const Vector &);
	Vector(const Set::VectorSpace::Vector &);
	Vector(const Covector &);
	Vector(const Set::Euclidean::Polar::Point &);
	Vector & operator = (const Set::VectorSpace::Vector &);
	Vector & operator = (const Vector &);
	Vector & operator = (const Set::Euclidean::Polar::Point &);
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

	typedef Set::Euclidean::Polar::CovectorZero dual_type;

	VectorZero();
	VectorZero(double * const &);
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

	typedef Set::Euclidean::Polar::Vector dual_type;

	Covector();
	Covector(double * const &);
	virtual ~Covector();
	Covector(const Set::VectorSpace::Vector &);
	Covector(const Covector &);
	Covector(const Set::Euclidean::Polar::Vector &);
	Covector & operator = (const Set::VectorSpace::Vector &);
	Covector & operator = (const Covector &);
};

//////////////////////////////////////////////////////////////////////
// Class CovectorZero
//////////////////////////////////////////////////////////////////////

class CovectorZero : public Set::VectorSpace::VectorZero
{
public:

	typedef Set::Euclidean::Polar::VectorZero dual_type;

	CovectorZero();
	CovectorZero(double * const &);
	virtual ~CovectorZero();
	CovectorZero(const Set::VectorSpace::VectorZero &);
	CovectorZero(const CovectorZero &);
	CovectorZero(const Set::Euclidean::Polar::VectorZero &);
	CovectorZero & operator = (const Set::VectorSpace::VectorZero &);
	CovectorZero & operator = (const CovectorZero &);
};

//////////////////////////////////////////////////////////////////////
// Class Embedding<0>
//////////////////////////////////////////////////////////////////////

template<>
class Embedding<0> : public Set::Manifold::Map
{
public: 

	typedef Set::Euclidean::Polar::Embedding<1> tangent_type;
	typedef Set::Euclidean::Polar::Point domain_type;
	typedef Set::Euclidean::Orthonormal::Point range_type;

	Embedding();
	Embedding(const Embedding<0> &);
	virtual ~Embedding();
	Embedding<0> & operator = (const Embedding<0> &);
	Set::Manifold::Map *Clone();
	Set::Manifold::TMap *Diff();
	void Randomize();
	const range_type & operator () (const domain_type &);
	unsigned int size1() const;
	unsigned int size2() const;
	const Set::Manifold::Point & 
		operator () (const Set::Manifold::Point &);

private:

	Set::Euclidean::Orthonormal::Point y; 
};

//////////////////////////////////////////////////////////////////////
// Class Embedding<1>
//////////////////////////////////////////////////////////////////////

template <>
class Embedding<1> : public Set::Manifold::TMap
{
public: 

	typedef Set::Euclidean::Polar::Embedding<2> tangent_type;
	typedef Set::Euclidean::Polar::Point domain_type;
	typedef Set::VectorSpace::Hom range_type;

	Embedding();
	Embedding(const Embedding<0> &);
	Embedding(const Embedding<1> &);
	virtual ~Embedding();
	Embedding<1> & operator = (const Embedding<1> &);
	Set::Manifold::TMap *Clone();
	Set::Manifold::TMap *Diff();
	void Randomize();
	const range_type & operator () (const domain_type &);
	unsigned int size1() const;
	unsigned int size2() const;
	const Set::VectorSpace::Hom & 
		operator () (const Set::Manifold::Point &);

private:

	Set::VectorSpace::Hom A;
};

//////////////////////////////////////////////////////////////////////
// Class Embedding<2>
//////////////////////////////////////////////////////////////////////

template <>
class Embedding<2> : public Set::Manifold::TMap
{
public: 

	typedef Set::Euclidean::Polar::Embedding<2> tangent_type;
	typedef Set::Euclidean::Polar::Point domain_type;
	typedef Set::VectorSpace::Hom range_type;

	Embedding();
	Embedding(const Embedding<1> &);
	Embedding(const Embedding<2> &);
	virtual ~Embedding();
	Embedding<2> & operator = (const Embedding<2> &);
	Set::Manifold::TMap *Clone();
	void Randomize();
	const range_type & operator () (const domain_type &);
	unsigned int size1() const;
	unsigned int size2() const;
	const Set::VectorSpace::Hom & 
		operator () (const Set::Manifold::Point &);

private:

	Set::VectorSpace::Hom A;

private:

	Set::Manifold::TMap *Diff(){return 0;}
};

//////////////////////////////////////////////////////////////////////
// Class Submersion<0>
//////////////////////////////////////////////////////////////////////

template<>
class Submersion<0> : public Set::Manifold::Map
{
public: 

	typedef Set::Euclidean::Polar::Submersion<1> tangent_type;
	typedef Set::Euclidean::Orthonormal::Point domain_type;
	typedef Set::Euclidean::Polar::Point range_type;

	Submersion();
	Submersion(const Submersion<0> &);
	virtual ~Submersion();
	Submersion<0> & operator = (const Submersion<0> &);
	Set::Manifold::Map *Clone();
	Set::Manifold::TMap *Diff();
	void Randomize();
	const range_type & operator () (const domain_type &);
	unsigned int size1() const;
	unsigned int size2() const;
	const Set::Manifold::Point & 
		operator () (const Set::Manifold::Point &);

private:

	Set::Euclidean::Polar::Point x; 
};

//////////////////////////////////////////////////////////////////////
// Class Submersion<1>
//////////////////////////////////////////////////////////////////////

template <>
class Submersion<1> : public Set::Manifold::TMap
{
public: 

	typedef Set::Euclidean::Polar::Submersion<2> tangent_type;
	typedef Set::Euclidean::Orthonormal::Point domain_type;
	typedef Set::VectorSpace::Hom range_type;

	Submersion();
	Submersion(const Submersion<0> &);
	Submersion(const Submersion<1> &);
	virtual ~Submersion();
	Submersion<1> & operator = (const Submersion<1> &);
	Set::Manifold::TMap *Clone();
	Set::Manifold::TMap *Diff();
	void Randomize();
	const range_type & operator () (const domain_type &);
	unsigned int size1() const;
	unsigned int size2() const;
	const Set::VectorSpace::Hom & 
		operator () (const Set::Manifold::Point &);

private:

	Set::VectorSpace::Hom A;
};

//////////////////////////////////////////////////////////////////////
// Class Submersion<2>
//////////////////////////////////////////////////////////////////////

template <>
class Submersion<2> : public Set::Manifold::TMap
{
public: 

	typedef Set::Euclidean::Polar::Submersion<2> tangent_type;
	typedef Set::Euclidean::Orthonormal::Point domain_type;
	typedef Set::VectorSpace::Hom range_type;

	Submersion();
	Submersion(const Submersion<1> &);
	Submersion(const Submersion<2> &);
	virtual ~Submersion();
	Submersion<2> & operator = (const Submersion<2> &);
	Set::Manifold::TMap *Clone();
	void Randomize();
	const range_type & operator () (const domain_type &);
	unsigned int size1() const;
	unsigned int size2() const;
	const Set::VectorSpace::Hom & 
		operator () (const Set::Manifold::Point &);

private:

	Set::VectorSpace::Hom A;

private:

	Set::Manifold::TMap *Diff(){return 0;}
};

}

}

}

Set::Euclidean::Polar::Point
operator + (const Set::Euclidean::Polar::Point &,
			const Set::Euclidean::Polar::Vector &);
Set::Euclidean::Polar::Point
operator - (const Set::Euclidean::Polar::Point &,
			const Set::Euclidean::Polar::Vector &);
Set::Euclidean::Polar::Vector
operator - (const Set::Euclidean::Polar::Point &,
			const Set::Euclidean::Polar::Point &);
Set::Euclidean::Polar::Point
operator + (const Set::Euclidean::Polar::Point &,
			const Set::VectorSpace::Vector &);
Set::Euclidean::Polar::Point
operator - (const Set::Euclidean::Polar::Point &,
			const Set::VectorSpace::Vector &);
void Random(Set::Euclidean::Polar::Point &);

#endif // !defined(MANIFOLD_EUCLIDEAN_POLAR_H__INCLUDED_)
