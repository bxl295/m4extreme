// Symmetric.h: interface for the Sym class. 
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(SET_VECTORSPACE_CATEGORY_SYMMETRIC_H__INCLUDED_)
#define SET_VECTORSPACE_CATEGORY_SYMMETRIC_H__INCLUDED_

#pragma once

#include "../Category.h"
#include "../../Vector/Vector.h"
#include "../../../../../Utils/Indexing/Indexing.h"

namespace Set
{
namespace VectorSpace
{
class Sym;
class SymDual;
class SymZero;
class SymId;
template<unsigned int> class SymSub;
template<unsigned int> class SymEmb;

//////////////////////////////////////////////////////////////////////
// Class Sym
//////////////////////////////////////////////////////////////////////

class Sym : public Vector
{
public:

	typedef Set::VectorSpace::SymZero zero_type;
	typedef Set::VectorSpace::SymId id_type;

	inline
	Sym():Vector(), n1(0){}

	inline
	Sym(unsigned int m1):Vector(m1*(m1+1)/2), n1(m1){}

	inline
	Sym(unsigned int m1, double * u):Vector(m1*(m1+1)/2,u), n1(m1){}

	inline
	virtual ~Sym(){}

	Sym(const Vector &);
	Sym(const Hom &);

	inline
	Sym(const Sym &A) : Vector(A), n1(A.n1){}

	Sym(const SymDual &);
	Sym & operator = (const Vector &);
	Sym & operator = (const Sym &);
	Sym & operator = (const Hom &);
	double operator () (const Vector &) const;
	double operator () (const Sym &) const;
	double operator () (const SymDual &) const;
	double operator () (const Hom &) const;
	double operator () (const Vector &, const Vector &) const;
	Hom Embed() const;

	inline
	unsigned int size1() const {return n1;}

	inline
	unsigned int size2() const {return n1;}

	virtual void write(ostream & os) const;
	virtual void read(istream & is);

protected:

	unsigned int n1;
};

//////////////////////////////////////////////////////////////////////
// Class SymDual
//////////////////////////////////////////////////////////////////////

class SymDual : public Vector
{
public:

	typedef Set::VectorSpace::SymZero zero_type;
	typedef Set::VectorSpace::SymId id_type;

	SymDual();
	SymDual(const unsigned int &);
	SymDual(const unsigned int &, double * const &);
	virtual ~SymDual();
	SymDual(const Vector &);
	SymDual(const SymDual &);
	SymDual(const Hom &);
	SymDual(const Sym &);
	SymDual & operator = (const Vector &);
	SymDual & operator = (const SymDual &);
	SymDual & operator = (const Hom &);
	double operator () (const Vector &) const;
	double operator () (const Sym &) const;
	double operator () (const SymDual &) const;
	double operator () (const Hom &) const;
	double operator () (const Vector &, const Vector &) const;
	Hom Embed() const;
	unsigned int size1() const;
	unsigned int size2() const;	

	virtual void write(ostream & os) const;
	virtual void read(istream & is);
private:

	unsigned int n1;
};

//////////////////////////////////////////////////////////////////////
// Class SymZero
//////////////////////////////////////////////////////////////////////

class SymZero : public VectorZero
{
public:

	SymZero();
	SymZero(const unsigned int &);
	SymZero(const unsigned int &, double * const &);
	virtual ~SymZero();
	SymZero(const VectorZero &);
	SymZero(const SymZero &);
	SymZero(const HomZero &);
	SymZero & operator = (const VectorZero &);
	SymZero & operator = (const SymZero &);
	SymZero & operator = (const HomZero &);
	double operator () (const Vector &) const;
	double operator () (const Sym &) const;
	double operator () (const SymDual &) const;
	double operator () (const Hom &) const;
	double operator () (const Vector &, const Vector &) const;
	HomZero Embed() const;
	unsigned int size1() const;
	unsigned int size2() const;	

private:

	unsigned int n1;
};

//////////////////////////////////////////////////////////////////////
// Class SymId
//////////////////////////////////////////////////////////////////////

class SymId : public Vector
{
public:

	SymId();
	SymId(const unsigned int &);
	SymId(const unsigned int &, double * const &);
	virtual ~SymId();
	SymId(const SymId &);
	SymId(const HomId &);
	SymId & operator = (const SymId &);
	SymId & operator = (const HomId &);
	double operator () (const Vector &) const;
	double operator () (const Sym &) const;
	double operator () (const SymDual &) const;
	double operator () (const Hom &) const;
	double operator () (const Vector &, const Vector &) const;
	HomId Embed() const;
	unsigned int size1() const;
	unsigned int size2() const;	

private:

	unsigned int n1;
};

//////////////////////////////////////////////////////////////////////
// Class SymSub<0>
//////////////////////////////////////////////////////////////////////

template<>
class SymSub<0>
{
public: 

	SymSub();
	SymSub(const unsigned int &);
	virtual ~SymSub();
	Set::VectorSpace::Sym 
		operator () (const Set::VectorSpace::Hom &);

private:

	unsigned int d;
};

//////////////////////////////////////////////////////////////////////
// Class SymSub<1>
//////////////////////////////////////////////////////////////////////

template<>
class SymSub<1>
{
public: 

	SymSub();
	SymSub(const unsigned int &);
	virtual ~SymSub();
	Set::VectorSpace::Hom 
	  operator () (const Set::VectorSpace::Hom &);
	void operator3D(double*, Set::VectorSpace::Hom &);

private:

	unsigned int d;
};

//////////////////////////////////////////////////////////////////////
// Class SymSub<2>
//////////////////////////////////////////////////////////////////////

template<>
class SymSub<2>
{
public: 

	SymSub();
	SymSub(const unsigned int &);
	virtual ~SymSub();
	Set::VectorSpace::Hom 
		operator () (const Set::VectorSpace::Hom &);

private:

	unsigned int d;
};

//////////////////////////////////////////////////////////////////////
// Class SymEmb<0>
//////////////////////////////////////////////////////////////////////

template<>
class SymEmb<0>
{
public: 

	SymEmb();
	SymEmb(const unsigned int &);
	virtual ~SymEmb();
	Set::VectorSpace::Hom 
		operator () (const Set::VectorSpace::Sym &);

private:

	unsigned int d;
};

//////////////////////////////////////////////////////////////////////
// Class SymEmb<1>
//////////////////////////////////////////////////////////////////////

template<>
class SymEmb<1>
{
public: 

	SymEmb();
	SymEmb(const unsigned int &);
	virtual ~SymEmb();
	Set::VectorSpace::Hom 
		operator () (const Set::VectorSpace::Hom &);

private:

	unsigned int d;
};

//////////////////////////////////////////////////////////////////////
// Class SymEmb<2>
//////////////////////////////////////////////////////////////////////

template<>
class SymEmb<2>
{
public: 

	SymEmb();
	SymEmb(const unsigned int &);
	virtual ~SymEmb();
	Set::VectorSpace::Hom 
		operator () (const Set::VectorSpace::Hom &);

private:

	unsigned int d;
};

}

}

//////////////////////////////////////////////////////////////////////
// Sym Operators
//////////////////////////////////////////////////////////////////////

Set::VectorSpace::Sym 
operator + (const Set::VectorSpace::Sym &, 
			const Set::VectorSpace::Sym &);
Set::VectorSpace::Sym 
operator + (const Set::VectorSpace::Sym &, 
			const Set::VectorSpace::SymId &);
Set::VectorSpace::Sym 
operator + (const Set::VectorSpace::SymId &, 
			const Set::VectorSpace::Sym &);
Set::VectorSpace::Sym 
operator + (const Set::VectorSpace::Sym &, 
			const Set::VectorSpace::SymZero &);
Set::VectorSpace::Sym 
operator + (const Set::VectorSpace::SymZero &, 
			const Set::VectorSpace::Sym &);
Set::VectorSpace::Sym 
operator - (const Set::VectorSpace::Sym &, 
			const Set::VectorSpace::Sym &);
Set::VectorSpace::Sym 
operator - (const Set::VectorSpace::Sym &, 
			const Set::VectorSpace::SymId &);
Set::VectorSpace::Sym 
operator - (const Set::VectorSpace::SymId &, 
			const Set::VectorSpace::Sym &);
Set::VectorSpace::Sym 
operator - (const Set::VectorSpace::Sym &, 
			const Set::VectorSpace::SymZero &);
Set::VectorSpace::Sym 
operator - (const Set::VectorSpace::SymZero &, 
			const Set::VectorSpace::Sym &);
Set::VectorSpace::Sym 
operator - (const Set::VectorSpace::Sym &);
void Null(Set::VectorSpace::Sym &);
Set::VectorSpace::Sym 
operator * (const Set::VectorSpace::Sym &, 
			const double &);
Set::VectorSpace::Sym 
operator * (const double &, 
			const Set::VectorSpace::Sym &);
Set::VectorSpace::Sym 
operator / (const Set::VectorSpace::Sym &, 
			const double &);
Set::VectorSpace::Vector 
operator * (const Set::VectorSpace::Sym &, 
			const Set::VectorSpace::Vector &);
Set::VectorSpace::Sym
operator * (const Set::VectorSpace::Sym &, 
			const Set::VectorSpace::Sym &);
void Identity(Set::VectorSpace::Sym &);
Set::VectorSpace::SymDual Inverse(const Set::VectorSpace::Sym &);
double Norm(const Set::VectorSpace::Sym &);
Set::VectorSpace::Sym Adjoint(const Set::VectorSpace::Sym &);
double Jacobian(const Set::VectorSpace::Sym &);
double Trace(const Set::VectorSpace::Sym &);
Set::VectorSpace::Sym 
Deviator(const Set::VectorSpace::Sym &);

//////////////////////////////////////////////////////////////////////
// SymDual Operators
//////////////////////////////////////////////////////////////////////

Set::VectorSpace::SymDual 
operator + (const Set::VectorSpace::SymDual &, 
			const Set::VectorSpace::SymDual &);
Set::VectorSpace::SymDual 
operator + (const Set::VectorSpace::SymDual &, 
			const Set::VectorSpace::SymId &);
Set::VectorSpace::SymDual 
operator + (const Set::VectorSpace::SymId &, 
			const Set::VectorSpace::SymDual &);
Set::VectorSpace::SymDual 
operator + (const Set::VectorSpace::SymDual &, 
			const Set::VectorSpace::SymZero &);
Set::VectorSpace::SymDual 
operator + (const Set::VectorSpace::SymZero &, 
			const Set::VectorSpace::SymDual &);
Set::VectorSpace::SymDual 
operator - (const Set::VectorSpace::SymDual &, 
			const Set::VectorSpace::SymDual &);
Set::VectorSpace::SymDual 
operator - (const Set::VectorSpace::SymId &, 
			const Set::VectorSpace::SymDual &);
Set::VectorSpace::SymDual 
operator - (const Set::VectorSpace::SymDual &, 
			const Set::VectorSpace::SymId &);
Set::VectorSpace::SymDual 
operator - (const Set::VectorSpace::SymZero &, 
			const Set::VectorSpace::SymDual &);
Set::VectorSpace::SymDual 
operator - (const Set::VectorSpace::SymDual &, 
			const Set::VectorSpace::SymZero &);
Set::VectorSpace::SymDual 
operator - (const Set::VectorSpace::SymDual &);
void Null(Set::VectorSpace::SymDual &);
Set::VectorSpace::SymDual 
operator * (const Set::VectorSpace::SymDual &, 
			const double &);
Set::VectorSpace::SymDual 
operator * (const double &, 
			const Set::VectorSpace::SymDual &);
Set::VectorSpace::SymDual 
operator / (const Set::VectorSpace::SymDual &, 
			const double &);
Set::VectorSpace::Vector 
operator * (const Set::VectorSpace::SymDual &, 
			const Set::VectorSpace::Vector &);
Set::VectorSpace::SymDual
operator * (const Set::VectorSpace::SymDual &, 
			const Set::VectorSpace::SymDual &);
void Identity(Set::VectorSpace::SymDual &);
Set::VectorSpace::Sym Inverse(const Set::VectorSpace::SymDual &);
double Norm(const Set::VectorSpace::SymDual &);
Set::VectorSpace::SymDual Adjoint(const Set::VectorSpace::SymDual &);
double Jacobian(const Set::VectorSpace::SymDual &);
double Trace(const Set::VectorSpace::SymDual &);
Set::VectorSpace::SymDual
Deviator(const Set::VectorSpace::SymDual &);

//////////////////////////////////////////////////////////////////////
// SymZero Operators
//////////////////////////////////////////////////////////////////////

Set::VectorSpace::SymZero 
operator + (const Set::VectorSpace::SymZero &, 
			const Set::VectorSpace::SymZero &);

Set::VectorSpace::SymZero 
operator - (const Set::VectorSpace::SymZero &, 
			const Set::VectorSpace::SymZero &);

#endif // !defined(SET_VECTORSPACE_CATEGORY_SYMMETRIC_H__INCLUDED_)
