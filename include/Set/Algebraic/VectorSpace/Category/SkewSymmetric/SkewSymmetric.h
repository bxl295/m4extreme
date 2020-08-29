// SkewSymmetric.h: interface for the Skw class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(SET_VECTORSPACE_CATEGORY_SKEWSYMMETRIC_H__INCLUDED_)
#define SET_VECTORSPACE_CATEGORY_SKEWSYMMETRIC_H__INCLUDED_

#pragma once

#include "../Category.h"
#include "../../Vector/Vector.h"

namespace Set
{
namespace VectorSpace
{
class Skw;
class SkwZero;

//////////////////////////////////////////////////////////////////////
// Class Skw
//////////////////////////////////////////////////////////////////////

class Skw : public Vector
{
public:

	typedef Set::VectorSpace::SkwZero zero_type;

	Skw();
	Skw(const unsigned int &);
	Skw(const unsigned int &, double * const &);
	virtual ~Skw();
	Skw(const Vector &);
	Skw(const Skw &);
	Skw(const Hom &);
	Skw & operator = (const Vector &);
	Skw & operator = (const Skw &);
	Skw & operator = (const Hom &);
	double operator () (const Vector &) const;
	double operator () (const Skw &) const;
	double operator () (const Hom &) const;
	Hom Embed() const;
	unsigned int size1() const;
	unsigned int size2() const;	

private:

	unsigned int n1;
};

//////////////////////////////////////////////////////////////////////
// Class SkwZero
//////////////////////////////////////////////////////////////////////

class SkwZero : public VectorZero
{
public:

	SkwZero();
	SkwZero(const unsigned int &);
	SkwZero(const unsigned int &, double * const &);
	virtual ~SkwZero();
	SkwZero(const VectorZero &);
	SkwZero(const SkwZero &);
	SkwZero(const HomZero &);
	SkwZero & operator = (const VectorZero &);
	SkwZero & operator = (const SkwZero &);
	SkwZero & operator = (const HomZero &);
	double operator () (const Vector &) const;
	double operator () (const Skw &) const;
	double operator () (const Hom &) const;
	HomZero Embed() const;
	unsigned int size1() const;
	unsigned int size2() const;	

private:

	unsigned int n1;
};

namespace SkewSymmetric
{
//////////////////////////////////////////////////////////////////////
// Class Exp
//////////////////////////////////////////////////////////////////////

class Exp
{
public: 

	typedef Set::VectorSpace::Skw domain_type;
	typedef Set::VectorSpace::Hom range_type;

	Exp();
	virtual ~Exp();
	Exp(const Exp &);
	Exp & operator = (const Exp &);
	range_type operator () (const domain_type &);
};

//////////////////////////////////////////////////////////////////////
// Class Log
//////////////////////////////////////////////////////////////////////

class Log
{
public: 

	typedef Set::VectorSpace::Hom domain_type;
	typedef Set::VectorSpace::Skw range_type;

	Log();
	virtual ~Log();
	Log(const Log &);
	Log & operator = (const Log &);
	range_type operator () (const domain_type &);
};

}

}

}

#endif // !defined(SET_VECTORSPACE_CATEGORY_SKEWSYMMETRIC_H__INCLUDED_)
