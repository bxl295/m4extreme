// Category.h: interface for the Category class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(SET_VECTORSPACE_CATEGORY_H__INCLUDED_)
#define SET_VECTORSPACE_CATEGORY_H__INCLUDED_

#pragma once

#include <iostream>
#include "../Vector/Vector.h"
#include "../../../Indexed/Array/Array.h"

using namespace std;

namespace Set
{
namespace VectorSpace
{
class Hom;
class HomZero;
class HomId;

//////////////////////////////////////////////////////////////////////
// Class Hom
//////////////////////////////////////////////////////////////////////

class Hom : public Vector
{
public:

	typedef Set::VectorSpace::HomZero zero_type;
	typedef Set::VectorSpace::HomId id_type;

	inline
	Hom() : n1(0), n2(0), L(0){}

	inline
	Hom( unsigned int m1,  unsigned int m2)
	  : Vector(m1*m2), n1(m1), n2(m2), L(new Vector * [m2])
	  {
	    unsigned int i; double *p;
	    for (i=0, p=this->head; i<n2; i++, p+=n1) {L[i] = new Vector(n1,p);}
	  }

	inline  
	Hom( unsigned int m1,  unsigned int m2, double * u)
	  : Vector(m1*m2,u), n1(m1), n2(m2), L(new Vector * [m2])
	  {
	    unsigned int i; double *p;
	    for (i=0, p=this->head; i<n2; i++, p+=n1) {L[i] = new Vector(n1,p);}
	  }

	inline
	Hom( unsigned int m)
	  : Vector(m*m), n1(m), n2(m), L(new Vector * [n2])
	  {
	    unsigned int i; double *p;
	    for (i=0, p=this->head; i<n2; i++, p+=n1) {L[i] = new Vector(n1,p);}
	  }

	inline  
	Hom( unsigned int m, double * u)
	  : Vector(m*m,u), n1(m), n2(m), L(new Vector * [n2])
	  {
	    unsigned int i; double *p;
	    for (i=0, p=this->head; i<n2; i++, p+=n1) {L[i] = new Vector(n1,p);}
	  }

	inline  
	virtual ~Hom()
	  {
	    for (unsigned int i=0; i<n2; i++) delete L[i];
	    delete [] L; 
	  }

	inline
	Hom(const Hom &A)
	  : Vector(A), n1(A.n1), n2(A.n2), L(new Vector * [n2])
	  {
	    unsigned int i; double *p;
	    for (i=0, p=this->head; i<n2; i++, p+=n1) {L[i] = new Vector(n1,p);}
	  }

	inline
	Hom & operator = (const Hom &A)
	  {
	    assert(n1 == A.n1); assert(n2 == A.n2); 
	    if (this == &A) return *this;
	    Vector::operator = (A);
	    return *this;
	  }

	inline
	Hom & operator = (const Vector &A)
	  {
	    assert(Vector::n == A.size()); 
	    Vector::operator = (A);
	    return *this;
	  }

	inline
	Hom & operator = (const Set::Array &);

	inline
	Vector operator () (const Vector &A) const
	  {
	    assert (A.size() == n2); 
	    unsigned int i, j; 	
	    const double *p;
	    const double *r;
	    double *q;
	    Vector B(n1); 
	    for (i=0, r=Vector::head, p=A.begin(); i<n2; i++, p++)
	      for (j=0, q=B.begin(); j<n1; j++, q++, r++)
		*q += (*p)*(*r);
	    return B;
	  }

	// this[i][j] * A[j] * B[i]
	inline
	double operator () (const Vector &A, const Vector &B) const
	  {
	    assert (A.size() == n2); 
	    assert (B.size() == n1); 
	    unsigned int i, j; 
	    double *q;
	    const double *p;
	    const double *r;
	    Vector C(n1);
	    for (i=0, r=Vector::head, p=A.begin(); i<n2; i++, p++)
	      for (j=0, q=C.begin(); j<n1; j++, q++, r++)
		*q += (*p)*(*r);
	    return C(B);
	  }

	// u[i] = this[i][j] * v[j]
	inline
	  Vector & operator () (const Vector & v, Vector & u) const {
	  assert(v.size() == n2 && u.size() == n1);
	  unsigned int i, j;
	  double * uh = u.begin();
	  const double * vh = v.begin();

	  for ( i=0; i < n1; ++i ) *(uh+i) = 0.0;

	  for ( i=0; i < n2; ++i ) {
	    for ( j=0; j < n1; ++j ) {
	      *(uh+j) += *(head + i * n1 + j) * *(vh+i);
	    }
	  }

	  return u;
	}

	// u[i] = a * this[i][j] * v[j]
	inline
	  Vector & operator () (double a, const Vector & v, Vector & u) const {
	  assert(v.size() == n2 && u.size() == n1);
	  this->operator()(v, u);
	  double * uh = u.begin();
	  for ( int i=0; i < n1; ++i ) *(uh+i) *= a;

	  return u;
	}

	/*
	inline
	  Hom & Adjoint() {
	  if ( n1 == n2 ) {
	    if ( n1 == 3 ) {
	      swap(*(head+1), *(head+3));
	      swap(*(head+2), *(head+6));
	      swap(*(head+5), *(head+7));
	    }
	    else if ( n1 == 2 ) {
	      swap(*(head+1), *(head+2));
	    }
	    else{
	      *this = Adjoint(*this);
	    }
	  }
	  else {
	    *this = Adjoint(*this);
	  }
	  return *this;
	}
	*/

	inline
	void operator += (const Vector &A)
	  {
	    assert(n1*n2 == A.size()); 
	    Vector::operator +=(A);
	  }

	inline
	void operator -= (const Vector &A)
	  {
	    assert(n1*n2 == A.size()); 
	    Vector::operator -=(A);
	  }

	inline
	void operator += (const Hom &A)
	  {
	    assert(n1 == A.n1); 
	    assert(n2 == A.n2); 
	    Vector::operator +=(A);
	  }

	inline
	void operator -= (const Hom &A)
	  {
	    assert(n1 == A.n1); 
	    assert(n2 == A.n2); 
	    Vector::operator -=(A);
	  }

	inline
	double operator () (const Hom &A)
	  {
	    assert(n1 == A.n1); 
	    assert(n2 == A.n2); 
	    double a = 0.0; 
	    const double *p;
	    const double *q;
	    for (p=Vector::head, q=A.Vector::begin(); 
		 p<Vector::tail; p++, q++) a += (*q)*(*p);
	    return a;
	  }

	inline
	const Vector & operator [] ( unsigned int i) const
	  {
	    assert(i < n2); 
	    return *L[i];
	  }

	inline
	Vector & operator [] ( unsigned int i)
	  {
	    assert(i < n2); 
	    return *L[i];
	  }
	
	inline
	double operator()(unsigned int i, unsigned int j) const { //ith column and jth row
	  assert(i<n2 && j<n1);
	  return *(head + i * n1 + j);
	}

	inline
	double & operator()(unsigned int i, unsigned int j) { //ith column and jth row
	  assert(i<n2 && j<n1);
	  return *(head + i * n1 + j);	
	}

	void print(ostream *) const;

	inline
	unsigned int size1() const {return n1;}

	inline
	unsigned int size2() const {return n2;}	

	void operator += (const HomZero &);
	void operator -= (const HomZero &);
	double operator () (const HomZero &);
	void operator += (const HomId &);
	void operator -= (const HomId &);
	double operator () (const HomId &);

	void write(ostream & os) const;
	void read(istream & is);

protected:

	unsigned int n1;
	unsigned int n2;
	Vector **L;
};

//////////////////////////////////////////////////////////////////////
// Class HomZero
//////////////////////////////////////////////////////////////////////

class HomZero : public Hom
{
public:

	HomZero();
	HomZero(const unsigned int &, const unsigned int &);
	HomZero(const unsigned int &, const unsigned int &, double * const &);
	HomZero(const unsigned int &);
	HomZero(const unsigned int &, double * const &);
	virtual ~HomZero();
	HomZero(const Hom &);
	HomZero(const HomZero &);
	HomZero & operator = (const HomZero &);
	HomZero & operator = (const VectorZero &);
	void operator += (const HomZero &);
	void operator -= (const HomZero &);
	void operator *= (const double &);
	void operator /= (const double &);
	VectorZero operator () (const Vector &) const;
	double operator () (const Vector &A, const Vector &) const;
	double operator () (const Hom &);
	double operator () (const HomId &);
};

//////////////////////////////////////////////////////////////////////
// Class HomId
//////////////////////////////////////////////////////////////////////

class HomId : public Hom
{
public:

	HomId();
	HomId(const unsigned int &);
	HomId(const unsigned int &, double * const &);
	virtual ~HomId();
	HomId(const HomId &);
	HomId & operator = (const HomId &);
	void operator += (const HomZero &);
	void operator -= (const HomZero &);
	Vector operator () (const Vector &) const;
	double operator () (const Vector &, const Vector &) const;
	double operator () (const Hom &);
	double operator () (const HomZero &);
};

}

}

//////////////////////////////////////////////////////////////////////
// Class Hom
//////////////////////////////////////////////////////////////////////

inline
Set::VectorSpace::Hom 
operator + (const Set::VectorSpace::Hom &A, 
			const Set::VectorSpace::Hom &B)
{
	Set::VectorSpace::Hom C=A; C += B; return C;
}

inline
Set::VectorSpace::Hom 
operator - (const Set::VectorSpace::Hom &A, 
			const Set::VectorSpace::Hom &B)
{
	Set::VectorSpace::Hom C=A; C -= B; return C;
}

inline
Set::VectorSpace::Hom 
operator - (const Set::VectorSpace::Hom &A)
{
	Set::VectorSpace::Hom B=A; double *p; 
	for (p=B.Vector::begin(); p!=B.Vector::end(); p++) *p = -(*p);
	return B;
}

inline
void Null(Set::VectorSpace::Hom &A)
{
	double *p; 
	for (p=A.Vector::begin(); p!=A.Vector::end(); p++) *p = 0.0;
}

inline
Set::VectorSpace::Hom
operator * (const Set::VectorSpace::Hom &A, 
			const double &p)
{
	Set::VectorSpace::Hom B=A; B *= p; return B;
}

inline
Set::VectorSpace::Hom 
operator * (const double &p, 
			const Set::VectorSpace::Hom &A)
{
	Set::VectorSpace::Hom B=A; B *= p; return B;
}

inline
Set::VectorSpace::Hom 
operator / (const Set::VectorSpace::Hom &A, 
			const double &p)
{
	Set::VectorSpace::Hom B=A; B /= p; return B;
}

inline
Set::VectorSpace::Vector
operator * (const Set::VectorSpace::Hom &A, 
			const Set::VectorSpace::Vector &B)
{
	return A(B);
}

Set::VectorSpace::Hom
operator * (const Set::VectorSpace::Hom &, 
			const Set::VectorSpace::Hom &);
void Identity(Set::VectorSpace::Hom &);
Set::VectorSpace::Hom Inverse(const Set::VectorSpace::Hom &);
double Norm(const Set::VectorSpace::Hom &);
double VectorNorm(const Set::VectorSpace::Hom &);
Set::VectorSpace::Hom Adjoint(const Set::VectorSpace::Hom &);
double Jacobian(const Set::VectorSpace::Hom &);

inline
double 
Trace(const Set::VectorSpace::Hom &A)
{
	assert (A.size1() == A.size2());
	unsigned int n1 = A.size1()+1; 
	const double *p; double a = 0.0;
	for (p=A.begin(); p<A.end(); p+=n1) a += *p;
	return a;
}

inline
Set::VectorSpace::Hom
Dyadic(const Set::VectorSpace::Vector &a, 
	   const Set::VectorSpace::Vector &b)
{
	unsigned int i;
	Set::VectorSpace::Hom A(b.size(),a.size());
	for (i=0; i<a.size(); i++) A[i] = a[i]*b;
	return A;
}

inline
Set::VectorSpace::Hom &
Dyadic(const Set::VectorSpace::Vector &a, 
       const Set::VectorSpace::Vector &b,
       Set::VectorSpace::Hom & A)
{       
        assert( A.size1() == b.size() && A.size2() == a.size() );
	
	double *head = A.begin();
	const double *pa;
	const double *pb;
	for (pa=a.begin(); pa!=a.end(); pa++) {
	  for (pb=b.begin(); pb!=b.end(); pb++) {
	    *(head++) = *pb * *pa;
	  }
	}

	return A;
}

inline
Set::VectorSpace::Hom &
DyadicSum(const Set::VectorSpace::Vector &a, 
	  const Set::VectorSpace::Vector &b,
	  Set::VectorSpace::Hom & A)
{       
        assert( A.size1() == b.size() && A.size2() == a.size() );
	
	double *head = A.begin();
	const double *pa;
	const double *pb;
	for (pa=a.begin(); pa!=a.end(); pa++) {
	  for (pb=b.begin(); pb!=b.end(); pb++) {
	    *(head++) += *pb * *pa;
	  }
	}

	return A;
}

inline
Set::VectorSpace::Hom
Adjoint(const Set::VectorSpace::Hom &A)
{
	unsigned int n1=A.size1(), n2=A.size2(); 
	Set::VectorSpace::Hom B(n2,n1);

	if (n1 == n2 && n1 == 3) {
	  double * Bh = B.begin();
	  double * Ah = A.begin();
	  *Bh = *Ah;
	  *(Bh+1) = *(Ah+3);
	  *(Bh+2) = *(Ah+6);
	  *(Bh+3) = *(Ah+1);
	  *(Bh+4) = *(Ah+4);
	  *(Bh+5) = *(Ah+7);
	  *(Bh+6) = *(Ah+2);
	  *(Bh+7) = *(Ah+5);
	  *(Bh+8) = *(Ah+8);
	}
	else {
	  unsigned int i1, i2;
	  for (i2=0; i2<n2; i2++) 
	    for (i1=0; i1<n1; i1++) 
	      B[i1][i2] = A[i2][i1];
	}
	return B;
}

inline
void
Adjoint(const Set::VectorSpace::Hom &A, Set::VectorSpace::Hom & B)
{
	unsigned int n1=A.size1(), n2=A.size2(); 
	assert(n1 == B.size2() && n2 == B.size1() );

	if (n1 == n2 && n1 == 3) {
	  double * Bh = B.begin();
	  double * Ah = A.begin();
	  *Bh     = *Ah;
	  *(Bh+1) = *(Ah+3);
	  *(Bh+2) = *(Ah+6);
	  *(Bh+3) = *(Ah+1);
	  *(Bh+4) = *(Ah+4);
	  *(Bh+5) = *(Ah+7);
	  *(Bh+6) = *(Ah+2);
	  *(Bh+7) = *(Ah+5);
	  *(Bh+8) = *(Ah+8);
	}
	else {
	  unsigned int i1, i2;
	  for (i2=0; i2<n2; i2++) 
	    for (i1=0; i1<n1; i1++) 
	      B[i1][i2] = A[i2][i1];
	}

	return;
}

//////////////////////////////////////////////////////////////////////
// Class HomZero
//////////////////////////////////////////////////////////////////////

Set::VectorSpace::HomZero 
operator + (const Set::VectorSpace::HomZero &, 
			const Set::VectorSpace::HomZero &);
Set::VectorSpace::Hom 
operator + (const Set::VectorSpace::Hom &, 
			const Set::VectorSpace::HomZero &);
Set::VectorSpace::Hom 
operator + (const Set::VectorSpace::HomZero &, 
			const Set::VectorSpace::Hom &);
Set::VectorSpace::HomZero 
operator - (const Set::VectorSpace::HomZero &, 
			const Set::VectorSpace::HomZero &);
Set::VectorSpace::Hom 
operator - (const Set::VectorSpace::Hom &, 
			const Set::VectorSpace::HomZero &);
Set::VectorSpace::Hom 
operator - (const Set::VectorSpace::HomZero &, 
			const Set::VectorSpace::Hom &);
Set::VectorSpace::HomZero 
operator - (const Set::VectorSpace::HomZero &);
Set::VectorSpace::HomZero 
operator * (const Set::VectorSpace::HomZero &,
			const double &);
Set::VectorSpace::HomZero 
operator * (const double &,
			const Set::VectorSpace::HomZero &);
Set::VectorSpace::HomZero 
operator / (const Set::VectorSpace::HomZero &, 
			const double &);
Set::VectorSpace::VectorZero
operator * (const Set::VectorSpace::HomZero &, 
			const Set::VectorSpace::Vector &);
Set::VectorSpace::HomZero
operator * (const Set::VectorSpace::HomZero &, 
			const Set::VectorSpace::HomZero &);
Set::VectorSpace::HomZero
operator * (const Set::VectorSpace::HomZero &, 
			const Set::VectorSpace::Hom &);
Set::VectorSpace::HomZero
operator * (const Set::VectorSpace::Hom &, 
			const Set::VectorSpace::HomZero &);
double Norm(const Set::VectorSpace::HomZero &);
double VectorNorm(const Set::VectorSpace::HomZero &);
Set::VectorSpace::HomZero
Adjoint(const Set::VectorSpace::HomZero &);
ostream & 
operator<<(ostream &os, 
		   const Set::VectorSpace::HomZero &);
double Trace(const Set::VectorSpace::HomZero &);

//////////////////////////////////////////////////////////////////////
// Class HomId
//////////////////////////////////////////////////////////////////////

Set::VectorSpace::Hom
operator + (const Set::VectorSpace::HomId &, 
			const Set::VectorSpace::HomId &);
Set::VectorSpace::HomZero 
operator - (const Set::VectorSpace::HomId &, 
			const Set::VectorSpace::HomId &);
Set::VectorSpace::Hom
operator * (const Set::VectorSpace::HomId &, 
			const double &);
Set::VectorSpace::Hom
operator * (const double &, 
			const Set::VectorSpace::HomId &);
Set::VectorSpace::Hom
operator / (const Set::VectorSpace::HomId &, 
			const double &);
Set::VectorSpace::HomId 
operator + (const Set::VectorSpace::HomId &, 
			const Set::VectorSpace::HomZero &);
Set::VectorSpace::HomId 
operator + (const Set::VectorSpace::HomZero &, 
			const Set::VectorSpace::HomId &);
Set::VectorSpace::HomId 
operator - (const Set::VectorSpace::HomId &, 
			const Set::VectorSpace::HomZero &);
Set::VectorSpace::Hom
operator + (const Set::VectorSpace::HomId &, 
			const Set::VectorSpace::Hom &);
Set::VectorSpace::Hom
operator + (const Set::VectorSpace::Hom &, 
			const Set::VectorSpace::HomId &);
Set::VectorSpace::Hom
operator - (const Set::VectorSpace::HomId &, 
			const Set::VectorSpace::Hom &);
Set::VectorSpace::Hom
operator - (const Set::VectorSpace::Hom &, 
			const Set::VectorSpace::HomId &);
Set::VectorSpace::Vector
operator * (const Set::VectorSpace::HomId &, 
			const Set::VectorSpace::Vector &);
Set::VectorSpace::HomId
operator * (const Set::VectorSpace::HomId &, 
			const Set::VectorSpace::HomId &);
Set::VectorSpace::Hom
operator * (const Set::VectorSpace::HomId &, 
			const Set::VectorSpace::Hom &);
Set::VectorSpace::Hom
operator * (const Set::VectorSpace::Hom &, 
			const Set::VectorSpace::HomId &);
double Norm(const Set::VectorSpace::HomId &);
double VectorNorm(const Set::VectorSpace::HomId &);
double Jacobian(const Set::VectorSpace::HomId &);
double Trace(const Set::VectorSpace::HomId &);

#endif // !defined(SET_VECTORSPACE_CATEGORY_H__INCLUDED_)
