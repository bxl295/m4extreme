// Vector.h: interface for the Vector class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(SET_VECTORSPACE_VECTOR_H__INCLUDED_)
#define SET_VECTORSPACE_VECTOR_H__INCLUDED_

#pragma once

#include "../../../Indexed/Array/Array.h"

using namespace std;

namespace Set
{
namespace VectorSpace
{
class Vector;
class VectorZero;

//////////////////////////////////////////////////////////////////////
// Class Vector
//////////////////////////////////////////////////////////////////////

class Vector : public Set::Array
{
public:

	typedef Set::VectorSpace::VectorZero zero_type;

	inline
	Vector(){}

	inline
	Vector(unsigned int n):Array(n){}

	inline
	Vector(unsigned int n, double u):Array(n,u){}

	inline
	Vector(unsigned int n, double * u):Array(n,u){}

	inline
	virtual ~Vector(){}

	inline
	Vector(const Vector & A):Array(A){}

	inline
	Vector(const Array & A):Array(A) {}
	
	Vector(const VectorZero &);
	
	inline
	Vector & operator = (const Vector &A){
	  if (this == &A) return *this;
	  Array::operator = (A);
	  return *this;
	}
	
	inline
	Vector & operator = (const Array & A){
	  Array::operator = (A);
	  return *this;
	}

	Vector & operator = (const VectorZero &);

	inline
	void operator += (const Vector & A){
	  assert(n == A.n); double *p, *q; 
	  for (p=this->head, q=A.head; p<this->tail; p++, q++) *p += *q;
	}

	inline
	void operator -= (const Vector & A){
	  assert(n == A.n); double *p, *q; 
	  for (p=this->head, q=A.head; p<this->tail; p++, q++) *p -= *q;
	}
	
	inline
	void operator *= (double a){
	  double *p; 
	  for (p=this->head; p<this->tail; p++) *p *= a;
	}

	inline
	void operator /= (double a){
	  double *p; 
	  for (p=this->head; p<this->tail; p++) *p /= a;
	}

	inline
	double operator () (const Vector & A) const{
	  assert(n == A.n); double a = 0.0; double *p, *q;
	  for (p=this->head, q=A.head; p<this->tail; p++, q++) a += (*q)*(*p);
	  return a;
	}

	void operator += (const VectorZero &);
	void operator -= (const VectorZero &);
	double operator () (const VectorZero &) const;
};

class VectorZero : public Vector
{
public:

	VectorZero();
	VectorZero(const unsigned int &);
	VectorZero(const unsigned int &, double * const &);
	virtual ~VectorZero();
	VectorZero(const Vector &);
	VectorZero(const VectorZero &);
	VectorZero & operator = (const VectorZero &);
	void operator += (const VectorZero &);
	void operator -= (const VectorZero &);
	void operator *= (const double &);
	void operator /= (const double &);
	double operator () (const Vector &) const;

private:

	void operator += (const Vector &);
	void operator -= (const Vector &);
};

}

}

inline
Set::VectorSpace::Vector 
operator + (const Set::VectorSpace::Vector & A, 
	    const Set::VectorSpace::Vector & B){
	Set::VectorSpace::Vector C=A; C += B; return C;
}

inline
Set::VectorSpace::Vector 
operator - (const Set::VectorSpace::Vector & A, 
	    const Set::VectorSpace::Vector & B){
        Set::VectorSpace::Vector C=A; C -= B; return C;
}

inline
Set::VectorSpace::Vector 
operator - (const Set::VectorSpace::Vector & A){
	Set::VectorSpace::Vector B=A; double *p; 
	for (p=B.begin(); p!=B.end(); p++) *p = -(*p);
	return B;
}

inline
void Null(Set::VectorSpace::Vector & A){
	double *p; 
	for (p=A.begin(); p!=A.end(); p++) *p = 0.0;
}

Set::VectorSpace::VectorZero 
operator + (const Set::VectorSpace::VectorZero &, 
			const Set::VectorSpace::VectorZero &);
Set::VectorSpace::Vector 
operator + (const Set::VectorSpace::Vector &, 
			const Set::VectorSpace::VectorZero &);
Set::VectorSpace::Vector 
operator + (const Set::VectorSpace::VectorZero &, 
			const Set::VectorSpace::Vector &);
Set::VectorSpace::VectorZero 
operator - (const Set::VectorSpace::VectorZero &, 
			const Set::VectorSpace::VectorZero &);
Set::VectorSpace::Vector 
operator - (const Set::VectorSpace::Vector &, 
			const Set::VectorSpace::VectorZero &);
Set::VectorSpace::Vector 
operator - (const Set::VectorSpace::VectorZero &, 
			const Set::VectorSpace::Vector &);
Set::VectorSpace::VectorZero 
operator - (const Set::VectorSpace::VectorZero &);

inline
Set::VectorSpace::Vector 
operator * (const Set::VectorSpace::Vector & A, 
			double a){
	Set::VectorSpace::Vector B=A; B *= a; return B;
}

inline
Set::VectorSpace::Vector 
operator * (double a, 
			const Set::VectorSpace::Vector & A){
	Set::VectorSpace::Vector B=A; B *= a; return B;
}

inline
Set::VectorSpace::Vector 
operator / (const Set::VectorSpace::Vector & A, 
			double a){
	Set::VectorSpace::Vector B=A; B /= a; return B;
}

Set::VectorSpace::VectorZero 
operator * (const Set::VectorSpace::VectorZero &, 
			const double &);
Set::VectorSpace::VectorZero 
operator * (const double &, 
			const Set::VectorSpace::VectorZero &);
Set::VectorSpace::VectorZero 
operator / (const Set::VectorSpace::VectorZero &, 
			const double &);

inline
double Norm(const Set::VectorSpace::Vector & A){
	double *p; double a=0.0;
	for (p=A.begin(); p!=A.end(); p++) a += (*p) * (*p);
	return sqrt(a);
}

double Norm(const Set::VectorSpace::VectorZero &);

#endif // !defined(SET_VECTORSPACE_VECTOR_H__INCLUDED_)
