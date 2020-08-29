// Axisymmetric.h: interface for the Axisymmetric class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(MANIFOLD_EUCLIDEAN_AXISYMMETRIC_H__INCLUDED_)
#define MANIFOLD_EUCLIDEAN_AXISYMMETRIC_H__INCLUDED_

#pragma once

#include <cassert>
#include "./Cylindrical.h"
#include "../Orthonormal/Orthonormal.h"
#include "../../Manifold.h"
#include "../../Category/Category.h"
#include "../../../Algebraic/AlgLib.h"

namespace Set	
{
namespace Euclidean
{	
namespace Axisymmetric
{
#define PI2 6.2831853071795862

template<unsigned int> class Embedding;
template<unsigned int> class Submersion;

//////////////////////////////////////////////////////////////////////
// Class Embedding<0>
//////////////////////////////////////////////////////////////////////

template<>
class Embedding<0> : public Set::Manifold::Map
{
public: 

	typedef Set::Euclidean::Axisymmetric::Embedding<1> tangent_type;
	typedef Set::Euclidean::Orthonormal::Point domain_type;
	typedef Set::Euclidean::Cylindrical::Point range_type;

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

	Set::Euclidean::Cylindrical::Point y; 
};

//////////////////////////////////////////////////////////////////////
// Class Embedding<1>
//////////////////////////////////////////////////////////////////////

template <>
class Embedding<1> : public Set::Manifold::TMap
{
public: 

	typedef Set::Euclidean::Axisymmetric::Embedding<2> tangent_type;
	typedef Set::Euclidean::Orthonormal::Point domain_type;
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

	typedef Set::Euclidean::Axisymmetric::Embedding<2> tangent_type;
	typedef Set::Euclidean::Orthonormal::Point domain_type;
	typedef Set::VectorSpace::HomZero range_type;

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

	Set::VectorSpace::HomZero A;

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

	typedef Set::Euclidean::Axisymmetric::Submersion<1> tangent_type;
	typedef Set::Euclidean::Cylindrical::Point domain_type;
	typedef Set::Euclidean::Orthonormal::Point range_type;

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

	Set::Euclidean::Orthonormal::Point x; 
};

//////////////////////////////////////////////////////////////////////
// Class Submersion<1>
//////////////////////////////////////////////////////////////////////

template <>
class Submersion<1> : public Set::Manifold::TMap
{
public: 

	typedef Set::Euclidean::Axisymmetric::Submersion<2> tangent_type;
	typedef Set::Euclidean::Cylindrical::Point domain_type;
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

	typedef Set::Euclidean::Axisymmetric::Submersion<2> tangent_type;
	typedef Set::Euclidean::Cylindrical::Point domain_type;
	typedef Set::VectorSpace::HomZero range_type;

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

	Set::VectorSpace::HomZero A;

private:

	Set::Manifold::TMap *Diff(){return 0;}
};

}

}

}

#endif // !defined(MANIFOLD_EUCLIDEAN_AXISYMMETRIC_H__INCLUDED_)
