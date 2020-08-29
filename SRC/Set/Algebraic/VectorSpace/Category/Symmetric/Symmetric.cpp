// Symmetric.cpp: Implementation of the Sym class. 
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <cassert>
#include "Symmetric.h"

namespace Set
{
namespace VectorSpace
{
//////////////////////////////////////////////////////////////////////
// Class Sym
//////////////////////////////////////////////////////////////////////

//
// Notes:
// Hom Strain=[e11, e21, e31, e12, e22, e32, e13, e23, e33] is a symmetric 3x3 matrix
// <---> Sym E=[e11, e22, e33, 2e12, 2e13, 2e23]
// note the differences to the standard voigt notation.
//

Sym::Sym(const Vector &A) : Vector(A)
{
	unsigned int n = A.size();
	float x = (float)(1 + 8*n); 
	float y = sqrt(x);
	n1 = (-1 + (unsigned int)y)/2;
	assert (n == n1*(n1+1)/2);
}

// Sym::Sym(const Hom &A) : 
// 	Vector(A.size1()*(A.size1()+1)/2), n1(A.size1())
// {
// 	assert (A.size1() == A.size2());
// 	unsigned int i, j; double *p, *q; 
// 	Set::VectorSpace::Hom B=A;
// 	for (j=0, p=this->head; j<n1; p++, j++) *p = *(B.begin()+j*n1+j);
// 	if (n1 == 1) return;
// 	for (j=0, p=this->head+n1; j<n1-1; j++)
// 		for (i=j+1, q=B.begin()+j*n1+i; i<n1; i++, p++, q++) *p = *q;
// 	Set::VectorSpace::Hom C=Adjoint(A);
// 	for (j=0, p=this->head+n1; j<n1-1; j++)
// 		for (i=j+1, q=C.begin()+j*n1+i; i<n1; i++, p++, q++) *p += *q;
// }

Sym::Sym(const Hom &A) : 
	Vector(A.size1()*(A.size1()+1)/2), n1(A.size1())
{
	assert (A.size1() == A.size2());

        switch(n1) {
            default:
            {
	        unsigned int i, j;  double *p, *q;
		for (j = 0, p = this->head; j < n1; p++, j++) *p = *(A.begin() + j * n1 + j);		
		for (j = 0, p = this->head + n1; j < n1 - 1; j++)
		  for (i = j + 1, q = A.begin() + j * n1 + i; i < n1; i++, p++, q++) *p = *q;
		Set::VectorSpace::Hom C = Adjoint(A);
		for (j = 0, p = this->head + n1; j < n1 - 1; j++)
		  for (i = j + 1, q = C.begin() + j * n1 + i; i < n1; i++, p++, q++) *p += *q;
            }
                    return;
            case 1:
                *head = *A.begin();
                return;
            case 2:
            {
                double * q = A.begin();
                *head = *q;
                *(head + 1) = *(q + 3);
                *(head + 2) = *(q + 1) + *(q + 2);
            }
                return;
            case 3:
            {
                double * q = A.begin();
                *head = *q;
                *(head + 1) = *(q + 4);
                *(head + 2) = *(q + 8);
                *(head + 3) = *(q + 1) + *(q + 3);
                *(head + 4) = *(q + 2) + *(q + 6);
                *(head + 5) = *(q + 5) + *(q + 7);
            }
                return;
        }

}


Sym::Sym(const SymDual &A) : Vector(A), n1(A.size1())
{
	double *p;
	for (p=this->head+n1; p<this->tail; p++) *p *= 2.0;
}

Sym & Sym::operator = (const Vector &A)
{
	assert(Vector::n == A.size()); 
	if (this == &A) return *this;
	Vector::operator = (A); 
	return *this;
}

Sym & Sym::operator = (const Sym &A)
{
	assert(n1 == A.n1); 
	if (this == &A) return *this;
	Vector::operator = (A); 
	return *this;
}

Sym & Sym::operator = (const Hom &A)
{
	assert (n1 = A.size1());
	assert (n1 = A.size2());

        switch(n1) {
            default:
            {
	        unsigned int i, j;  double *p, *q;
		for (j = 0, p = this->head; j < n1; p++, j++) *p = *(A.begin() + j * n1 + j);		
		for (j = 0, p = this->head + n1; j < n1 - 1; j++)
		  for (i = j + 1, q = A.begin() + j * n1 + i; i < n1; i++, p++, q++) *p = *q;
		Set::VectorSpace::Hom C = Adjoint(A);
		for (j = 0, p = this->head + n1; j < n1 - 1; j++)
		  for (i = j + 1, q = C.begin() + j * n1 + i; i < n1; i++, p++, q++) *p += *q;
            }
	    return *this;

            case 1:
                *head = *A.begin();
                return *this;

            case 2:
            {
                double * q = A.begin();
                *head = *q;
                *(head + 1) = *(q + 3);
                *(head + 2) = *(q + 1) + *(q + 2);
            }
	    return *this;

            case 3:
            {
                double * q = A.begin();
                *head = *q;
                *(head + 1) = *(q + 4);
                *(head + 2) = *(q + 8);
                *(head + 3) = *(q + 1) + *(q + 3);
                *(head + 4) = *(q + 2) + *(q + 6);
                *(head + 5) = *(q + 5) + *(q + 7);
            }
	    return *this;
        }
}

double Sym::operator () (const Vector &A) const
{
	assert (Vector::n == A.size());
	Hom B = Embed();
	Sym C(A);
	return B(C.Embed());
}

double Sym::operator () (const Sym &A) const
{
	assert (n1 == A.n1);
	Hom B = Embed();
	Hom C = A.Embed();
	return B(C);
}

double Sym::operator () (const SymDual &A) const
{
	assert (A.size() == this->size());
	double *p, *q; double a=0.0;
	for (p=A.begin(), q=this->begin(); p!=A.end(); p++, q++) 
		a += (*p)*(*q); return a;
}

double Sym::operator () (const Hom &A) const
{
	assert (n1 == A.size1());
	assert (n1 == A.size2());
	Hom B = Embed();
	return B(A);
}

double Sym::operator () (const Vector &A, 
						 const Vector &B) const
{
	return Embed()(A,B);
}

Hom Sym::Embed() const
{
  switch(n1) {
  default:
    {
	unsigned int i, j; double *p, *q; 
	Set::VectorSpace::Hom A(n1);
	for (j=0, p=this->head; j<n1; p++, j++) *(A.begin()+j*n1+j) = *p;
	if (n1 == 1) return A;
	for (j=0, p=this->head+n1; j<n1-1; j++)
		for (i=j+1, q=A.begin()+j*n1+i; i<n1; i++, p++, q++) *q = 0.5*(*p);
	Set::VectorSpace::Hom B=Adjoint(A);
	for (j=0, p=this->head+n1; j<n1-1; j++)
		for (i=j+1, q=B.begin()+j*n1+i; i<n1; i++, p++, q++) *q = 0.5*(*p);
	return B;
    }
  case 1:
    { 
      Set::VectorSpace::Hom A(1);
      double * pA = A.begin();
      *pA = *head;
      return A;
    }
  case 2:
    {
      Set::VectorSpace::Hom A(2);
      double *pA = A.begin();
      *pA = *head;
      *(pA+3) = *(head+1);
      *(pA+1) = *(pA+2) = *(head+2) * 0.5;     
      return A;
    }
  case 3:
    {
      Set::VectorSpace::Hom A(3);
      double *pA = A.begin();
      *pA = *head;
      *(pA+4) = *(head+1);
      *(pA+8) = *(head+2);
      *(pA+1) = *(pA+3) = *(head+3) * 0.5;
      *(pA+2) = *(pA+6) = *(head+4) * 0.5;
      *(pA+5) = *(pA+7) = *(head+5) * 0.5;
      return A;
    }
  }
}

void Sym::write(ostream & os) const 
{
  os.write((char*)&n1, sizeof(int));
  Array::write(os);
}


void Sym::read(istream & is)
{
  is.read((char*)&n1, sizeof(int));
  Array::read(is);
}

//////////////////////////////////////////////////////////////////////
// Class SymDual
//////////////////////////////////////////////////////////////////////

//
// Notes:
// Hom Stress=[s11, s21, s31, s12, s22, s32, s13, s23, s33] is a symmetric 3x3 matrix
// <---> Sym S=[s11, s22, s33, s12, s13, s23]
// note the differences to the standard voigt notation.
//

SymDual::SymDual() : Vector(), n1(0){}

SymDual::SymDual(const unsigned int &m1) : 
	Vector(m1*(m1+1)/2), n1(m1){}

SymDual::SymDual(const unsigned int &m1, double * const &u) : 
	Vector(m1*(m1+1)/2,u), n1(m1){}

SymDual::~SymDual(){}

SymDual::SymDual(const Vector &A) : Vector(A)
{
	unsigned int n = A.size();
	float x = (float)(1 + 8*n); 
	float y = sqrt(x);
	n1 = (-1 + (unsigned int)y)/2;
	assert (n == n1*(n1+1)/2);
}

SymDual::SymDual(const SymDual &A) : Vector(A), n1(A.n1){}

SymDual::SymDual(const Hom &A) : 
	Vector(A.size1()*(A.size1()+1)/2), n1(A.size1())
{
	assert (A.size1() == A.size2());

        switch(n1) {
            default:
            {
	        unsigned int i, j; double *p, *q; 
		Set::VectorSpace::Hom B=A;
		for (j=0, p=this->head; j<n1; p++, j++) *p = *(B.begin()+j*n1+j);
		if (n1 == 1) return;
		for (j=0, p=this->head+n1; j<n1-1; j++)
		  for (i=j+1, q=B.begin()+j*n1+i; i<n1; i++, p++, q++) *p = 0.5*(*q);
		Set::VectorSpace::Hom C=Adjoint(A);
		for (j=0, p=this->head+n1; j<n1-1; j++)
		  for (i=j+1, q=C.begin()+j*n1+i; i<n1; i++, p++, q++) *p += 0.5*(*q);
            }
                    return;
            case 1:
                *head = *A.begin();
                return;
            case 2:
            {
                double * q = A.begin();
                *head = *q;
                *(head + 1) = *(q + 3);
                *(head + 2) = *(q + 1);
            }
                return;
            case 3:
            {
                double * q = A.begin();
                *head = *q;
                *(head + 1) = *(q + 4);
                *(head + 2) = *(q + 8);
                *(head + 3) = *(q + 1);
                *(head + 4) = *(q + 2);
                *(head + 5) = *(q + 5);
            }
                return;
        }		
}

SymDual::SymDual(const Sym &A) : Vector(A), n1(A.size1())
{
	double *p;
	for (p=this->head+n1; p<this->tail; p++) *p /= 2.0;
}

SymDual & SymDual::operator = (const Vector &A)
{
	assert(Vector::n == A.size()); 
	if (this == &A) return *this;
	Vector::operator = (A); 
	return *this;
}

SymDual & SymDual::operator = (const SymDual &A)
{
	assert(n1 == A.n1); 
	if (this == &A) return *this;
	Vector::operator = (A); 
	return *this;
}

SymDual & SymDual::operator = (const Hom &A)
{
	assert (n1 = A.size1());
	assert (n1 = A.size2());
        switch(n1) {
            default:
            {
	        unsigned int i, j; double *p, *q; 
		Set::VectorSpace::Hom B=A;
		for (j=0, p=this->head; j<n1; p++, j++) *p = *(B.begin()+j*n1+j);
		if (n1 == 1) return *this;
		for (j=0, p=this->head+n1; j<n1-1; j++)
		  for (i=j+1, q=B.begin()+j*n1+i; i<n1; i++, p++, q++) *p = 0.5*(*q);
		Set::VectorSpace::Hom C=Adjoint(A);
		for (j=0, p=this->head+n1; j<n1-1; j++)
		  for (i=j+1, q=C.begin()+j*n1+i; i<n1; i++, p++, q++) *p += 0.5*(*q);
            }
	    return *this;        
           case 1:
                *head = *A.begin();
		return *this;
            case 2:
            {
                double * q = A.begin();
                *head = *q;
                *(head + 1) = *(q + 3);
                *(head + 2) = *(q + 1);
            }
	    return *this;
            case 3:
            {
                double * q = A.begin();
                *head = *q;
                *(head + 1) = *(q + 4);
                *(head + 2) = *(q + 8);
                *(head + 3) = *(q + 1);
                *(head + 4) = *(q + 2);
                *(head + 5) = *(q + 5);
            }
	    return *this;
        }		
	
}

double SymDual::operator () (const Vector &A) const
{
	assert (Vector::n == A.size());
	Hom B = Embed();
	SymDual C(A);
	return B(C.Embed());
}

double SymDual::operator () (const SymDual &A) const
{
	assert (n1 == A.n1);
	Hom B = Embed();
	Hom C = A.Embed();
	return B(C);
}

double SymDual::operator () (const Sym &A) const
{
	assert (A.size() == this->size());
	double *p, *q; double a=0.0;
	for (p=A.begin(), q=this->begin(); p!=A.end(); p++, q++) 
		a += (*p)*(*q); return a;
}

double SymDual::operator () (const Hom &A) const
{
	assert (n1 == A.size1());
	assert (n1 == A.size2());
	Hom B = Embed();
	return B(A);
}

double SymDual::operator () (const Vector &A, 
						 const Vector &B) const
{
	return Embed()(A,B);
}

Hom SymDual::Embed() const
{
  switch(n1) {
  default:
    {
	unsigned int i, j; double *p, *q; 
	Set::VectorSpace::Hom A(n1);
	for (j=0, p=this->head; j<n1; p++, j++) *(A.begin()+j*n1+j) = *p;
	if (n1 == 1) return A;
	for (j=0, p=this->head+n1; j<n1-1; j++)
		for (i=j+1, q=A.begin()+j*n1+i; i<n1; i++, p++, q++) *q = *p;
	Set::VectorSpace::Hom B=Adjoint(A);
	for (j=0, p=this->head+n1; j<n1-1; j++)
		for (i=j+1, q=B.begin()+j*n1+i; i<n1; i++, p++, q++) *q = *p;
	return B;
    }
  case 1:
    { 
      Set::VectorSpace::Hom A(1);
      double * pA = A.begin();
      *pA = *head;
      return A;
    }
  case 2:
    {
      Set::VectorSpace::Hom A(2);
      double *pA = A.begin();
      *pA = *head;
      *(pA+3) = *(head+1);
      *(pA+1) = *(pA+2) = *(head+2);
      return A;
    }
  case 3:
    {
      Set::VectorSpace::Hom A(3);
      double *pA = A.begin();
      *pA = *head;
      *(pA+4) = *(head+1);
      *(pA+8) = *(head+2);
      *(pA+1) = *(pA+3) = *(head+3);
      *(pA+2) = *(pA+6) = *(head+4);
      *(pA+5) = *(pA+7) = *(head+5);
      return A;
    }
  }
}

unsigned int SymDual::size1() const{return n1;}

unsigned int SymDual::size2() const{return n1;}	

void SymDual::write(ostream & os) const 
{
  os.write((char*)&n1, sizeof(int));
  Array::write(os);
}


void SymDual::read(istream & is)
{
  is.read((char*)&n1, sizeof(int));
  Array::read(is);
}

//////////////////////////////////////////////////////////////////////
// Class SymZero
//////////////////////////////////////////////////////////////////////

SymZero::SymZero() : VectorZero(), n1(0){}

SymZero::SymZero(const unsigned int &m1) : 
	VectorZero(m1*(m1+1)/2), n1(m1){}

SymZero::SymZero(const unsigned int &m1, double * const &u) : 
	VectorZero(m1*(m1+1)/2,u), n1(m1){}

SymZero::~SymZero(){}

SymZero::SymZero(const VectorZero &A) : VectorZero(A)
{
	unsigned int n = A.size();
	float x = (float)(1 + 8*n); 
	float y = sqrt(x);
	n1 = (-1 + (unsigned int)y)/2;
	assert (n == n1*(n1+1)/2);
}

SymZero::SymZero(const SymZero &A) : 
	VectorZero(A), n1(A.n1){}

SymZero::SymZero(const HomZero &A) : 
	VectorZero(A.size1()*(A.size1()+1)/2), n1(A.size1()){}

SymZero & SymZero::operator = (const VectorZero &A)
{
	assert(Vector::n == A.size()); 
	if (this == &A) return *this;
	return *this;
}

SymZero & SymZero::operator = (const SymZero &A)
{
	assert(n1 == A.n1); 
	return *this;
}

SymZero & SymZero::operator = (const HomZero &A)
{
	assert (n1 = A.size1());
	assert (n1 = A.size2());
	return *this;
}

double SymZero::operator () (const Vector &A) const
{
	assert (Vector::n == A.size());
	return 0.0;
}

double SymZero::operator () (const Sym &A) const
{
	assert (n1 == A.size1());
	return 0.0;
}

double SymZero::operator () (const SymDual &A) const
{
	assert (n1 == A.size1());
	return 0.0;
}

double SymZero::operator () (const Hom &A) const
{
	assert (n1 == A.size1());
	assert (n1 == A.size2());
	return 0.0;
}

double SymZero::operator () (const Vector &, 
							 const Vector &) const
{
	return 0.0;
}

HomZero SymZero::Embed() const
{
	HomZero A(n1);
	return A;
}

unsigned int SymZero::size1() const{return n1;}

unsigned int SymZero::size2() const{return n1;}	

//////////////////////////////////////////////////////////////////////
// Class SymId
//////////////////////////////////////////////////////////////////////

SymId::SymId() : Vector(), n1(0){}

SymId::SymId(const unsigned int &m1) : 
	Vector(m1*(m1+1)/2), n1(m1)
{
	double *p; 
	for (p=this->head; p<this->head+n1; p++) *p=1.0;
}

SymId::SymId(const unsigned int &m1, double * const &u) : 
	Vector(m1*(m1+1)/2,u), n1(m1)
{
	double *p; 
	for (p=this->head; p<this->head+n1; p++) *p=1.0;
}

SymId::~SymId(){}

SymId::SymId(const SymId &A) : 
	Vector(A), n1(A.n1){}

SymId::SymId(const HomId &A) : 
	Vector(A.size1()*(A.size1()+1)/2), n1(A.size1())
{
	double *p; 
	for (p=this->head; p<this->head+n1; p++) *p=1.0;
}

SymId & SymId::operator = (const SymId &A)
{
	assert(n1 == A.n1); 
	return *this;
}

SymId & SymId::operator = (const HomId &A)
{
	assert (n1 = A.size1());
	assert (n1 = A.size2());
	return *this;
}

double SymId::operator () (const Vector &A) const
{
	assert (Vector::n == A.size());
	Sym B(A); return Trace(B.Embed());
}

double SymId::operator () (const Sym &A) const
{
	assert (n1 == A.size1());
	double *p; double t = 0.0;
	for (p=A.begin(); p<A.begin()+n1; p++) t+=*p;
	return t;
}

double SymId::operator () (const SymDual &A) const
{
	assert (n1 == A.size1());
	double *p; double t = 0.0;
	for (p=A.begin(); p<A.begin()+n1; p++) t+=*p;
	return t;
}

double SymId::operator () (const Hom &A) const
{
	assert (n1 == A.size1());
	assert (n1 == A.size2());
	return Trace(A);
}

double SymId::operator () (const Vector &A, 
						   const Vector &B) const
{
	return A(B);
}

HomId SymId::Embed() const
{
	HomId A(n1);
	return A;
}

unsigned int SymId::size1() const{return n1;}

unsigned int SymId::size2() const{return n1;}	

//////////////////////////////////////////////////////////////////////
// Class SymSub<0>
//////////////////////////////////////////////////////////////////////

SymSub<0>::SymSub() : d(0) {}

SymSub<0>::SymSub(const unsigned int &d_) : d(d_) {}

SymSub<0>::~SymSub() {}

Set::VectorSpace::Sym 
SymSub<0>::operator () (const Set::VectorSpace::Hom &A)
{
	assert(A.size1() == A.size2());
	assert (d == A.size1());	
	return Set::VectorSpace::Sym(A);
}

//////////////////////////////////////////////////////////////////////
// Class SymSub<1>
//////////////////////////////////////////////////////////////////////

SymSub<1>::SymSub() : d(0)  {}

SymSub<1>::SymSub(const unsigned int &d_) : d(d_) {}

SymSub<1>::~SymSub() {}

// Set::VectorSpace::Hom 
// SymSub<1>::operator () (const Set::VectorSpace::Hom &A)
// {
// 	unsigned int i, j, m; 
// 	unsigned int n1 = d*(d+1)/2;
// 	unsigned int m1 = A.size2();
// 	assert(m1 == A.size1()); 
// 	assert(d*d == m1);

// 	double ***C = Indexing::New(A.begin(),m1,d,d);
// 	Set::VectorSpace::Hom B(n1); SymSub<0> g(d);

// 	for (m=0; m<d; m++) 
// 	{
// 		Set::VectorSpace::Hom Aloc(d,C[m][m]);
// 		B[m] = g(Aloc);
// 	}

// 	if (d == 1) 
// 	{
// 		Indexing::Delete(C,m1,d,d);
// 		return B;
// 	}

// 	for (j=0, m=d; j<d-1; j++) 
// 	{
// 		for (i=j+1; i<d; i++, m++) 
// 		{
// 			Set::VectorSpace::Hom Aloc(d,C[j][i]);
// 			Set::VectorSpace::Hom Bloc(d,C[i][j]);
// 			B[m] = 0.5*(g(Aloc) + g(Bloc));
// 		}
// 	}

// 	Indexing::Delete(C,m1,d,d);
// 	return B;
// }

Set::VectorSpace::Hom 
SymSub<1>::operator () (const Set::VectorSpace::Hom &A)
{
	unsigned int i, j, m; 
	unsigned int n1 = d*(d+1)/2;
	unsigned int m1 = A.size2();
	assert(m1 == A.size1()); 
	assert(d*d == m1);

	double ***C = Indexing::New(A.begin(),m1,d,d);
	Set::VectorSpace::Hom B(n1); 
	Set::VectorSpace::Sym g(d), gT(d);
	const double * pg;
	const double * pgT;
	double * head;

	for (m=0; m<d; m++) 
	{
		g = Set::VectorSpace::Hom(d,C[m][m]);
		head = B[m].begin();
		for ( pg = g.begin(); pg != g.end(); pg++ )
		  *(head++) = *pg;
	}

	if (d == 1) 
	{
		Indexing::Delete(C,m1,d,d);
		return B;
	}

	for (j=0, m=d; j<d-1; j++) 
	{
		for (i=j+1; i<d; i++, m++) 
		{
			g  = Set::VectorSpace::Hom(d,C[j][i]);
			gT = Set::VectorSpace::Hom(d,C[i][j]);
			head = B[m].begin();			
			for(pg = g.begin(), pgT = gT.begin(); pg != g.end(); pg++, pgT++)
			  *(head++) = 0.5 * (*pg + *pgT);			
		}
	}

	Indexing::Delete(C,m1,d,d);
	return B;
}

void
SymSub<1>::operator3D(double * Ah, Set::VectorSpace::Hom & B)
{
	unsigned int i, j, m; 
	assert( d == 3 );
	assert(B.size1() == 6 && B.size2() == 6);

	double ***C = Indexing::New(Ah, 9, 3, 3);

	Set::VectorSpace::Sym g(3), gT(3);
	const double * pg;
	const double * pgT;
	double * head;

	//m=0
	g = Set::VectorSpace::Hom(3,C[0][0]);
	head = B[0].begin();
	for ( pg = g.begin(); pg != g.end(); pg++ ) *(head++) = *pg;

	//m=1
	g = Set::VectorSpace::Hom(3,C[1][1]);
	head = B[1].begin();
	for ( pg = g.begin(); pg != g.end(); pg++ ) *(head++) = *pg;

	//m=2
	g = Set::VectorSpace::Hom(3,C[2][2]);
	head = B[2].begin();
	for ( pg = g.begin(); pg != g.end(); pg++ ) *(head++) = *pg;

	//j=0
	////i=1
	g  = Set::VectorSpace::Hom(3,C[0][1]);
	gT = Set::VectorSpace::Hom(3,C[1][0]);
	head = B[3].begin();			
	for(pg = g.begin(), pgT = gT.begin(); pg != g.end(); pg++, pgT++)
	  *(head++) = 0.5 * (*pg + *pgT);			
	
	////i=2
	g  = Set::VectorSpace::Hom(3,C[0][2]);
	gT = Set::VectorSpace::Hom(3,C[2][0]);
	head = B[4].begin();			
	for(pg = g.begin(), pgT = gT.begin(); pg != g.end(); pg++, pgT++)
	  *(head++) = 0.5 * (*pg + *pgT);	
	
	//j=1
	////i=2
	g  = Set::VectorSpace::Hom(3,C[1][2]);
	gT = Set::VectorSpace::Hom(3,C[2][1]);
	head = B[5].begin();			
	for(pg = g.begin(), pgT = gT.begin(); pg != g.end(); pg++, pgT++)
	  *(head++) = 0.5 * (*pg + *pgT);			

	Indexing::Delete(C,9,3,3);

	return;
}

//////////////////////////////////////////////////////////////////////
// Class SymSub<2>
//////////////////////////////////////////////////////////////////////

SymSub<2>::SymSub() : d(0)  {}

SymSub<2>::SymSub(const unsigned int &d_) : d(d_) {}

SymSub<2>::~SymSub() {}

Set::VectorSpace::Hom 
SymSub<2>::operator () (const Set::VectorSpace::Hom &A)
{
	unsigned int i, j, m; 
	unsigned int n1 = d*(d+1)/2;
	unsigned int m1 = A.size2();
	assert(m1*m1 == A.size1()); 
	assert(d*d == m1);

	double ***C = Indexing::New(A.begin(),m1*m1,d,d);
	Set::VectorSpace::Hom B(n1*n1,n1); SymSub<1> dg(d);

	for (m=0; m<d; m++) 
	{
		Set::VectorSpace::Hom Aloc(m1,C[m][m]);
		B[m] = dg(Aloc);
	}

	if (d == 1) 
	{
		Indexing::Delete(C,m1*m1,d,d);
		return B;
	}

	for (j=0, m=d; j<d-1; j++) 
	{
		for (i=j+1; i<d; i++, m++) 
		{
			Set::VectorSpace::Hom Aloc(m1,C[j][i]);
			Set::VectorSpace::Hom Bloc(m1,C[i][j]);
			B[m] = 0.5*(dg(Aloc) + dg(Bloc));
		}
	}

	Indexing::Delete(C,m1*m1,d,d);
	return B;
}

//////////////////////////////////////////////////////////////////////
// Class SymEmb<0>
//////////////////////////////////////////////////////////////////////

SymEmb<0>::SymEmb() : d(0) {}

SymEmb<0>::SymEmb(const unsigned int &d_) : d(d_) {}

SymEmb<0>::~SymEmb() {}

Set::VectorSpace::Hom
SymEmb<0>::operator () (const Set::VectorSpace::Sym &B)
{
	assert(B.size1() == B.size2());
	assert(d == B.size1());
	return B.Embed();
}

//////////////////////////////////////////////////////////////////////
// Class SymEmb<1>
//////////////////////////////////////////////////////////////////////

SymEmb<1>::SymEmb() : d(0) {}

SymEmb<1>::SymEmb(const unsigned int &d_) : d(d_) {}

SymEmb<1>::~SymEmb() {}

Set::VectorSpace::Hom 
SymEmb<1>::operator () (const Set::VectorSpace::Hom &B)
{
	unsigned int i, j, m;
	unsigned int m1 = d*d;
	unsigned int n1 = d*(d+1)/2;
	assert(n1 == B.size1());
	assert(n1 == B.size2());
	Set::VectorSpace::Hom A(m1); SymEmb<0> f(d);

	double ***C = Indexing::New(A.begin(),m1,d,d);

	for (m=0; m<d; m++) 
	{
		Set::VectorSpace::Vector Aloc(m1,C[m][m]);
		Set::VectorSpace::Sym Bloc(d,B[m].begin());
		Aloc = f(Bloc);
	}

	if (d == 1) 
	{
		Indexing::Delete(C,m1,d,d); 
		return A;
	}

	for (j=0, m=d; j<d-1; j++) 
	{
		for (i=j+1; i<d; i++, m++) 
		{
			Set::VectorSpace::Vector Aloc(m1,C[j][i]);
			Set::VectorSpace::Vector Bloc(m1,C[i][j]);
			Set::VectorSpace::Sym Cloc(d,B[m].begin());
			Aloc = Bloc = f(Cloc);
		}
	}

	Indexing::Delete(C,m1,d,d); 
	return A;
}

//////////////////////////////////////////////////////////////////////
// Class SymEmb<2>
//////////////////////////////////////////////////////////////////////

SymEmb<2>::SymEmb() : d(0) {}

SymEmb<2>::SymEmb(const unsigned int &d_) : d(d_) {}

SymEmb<2>::~SymEmb() {}

Set::VectorSpace::Hom 
SymEmb<2>::operator () (const Set::VectorSpace::Hom &B)
{
	unsigned int i, j, m;
	unsigned int m1 = d*d;
	unsigned int n1 = d*(d+1)/2;
	assert(n1*n1 == B.size1());
	assert(n1 == B.size2());
	Set::VectorSpace::Hom A(m1*m1,m1); SymEmb<1> df(d);

	double ***C = Indexing::New(A.begin(),m1*m1,d,d);

	for (m=0; m<d; m++) 
	{
		Set::VectorSpace::Hom Aloc(m1,C[m][m]);
		Set::VectorSpace::Hom Bloc(n1,B[m].begin());
		Aloc = df(Bloc);
	}

	if (d == 1) 
	{
		Indexing::Delete(C,m1*m1,d,d); 
		return A;
	}

	for (j=0, m=d; j<d-1; j++) 
	{
		for (i=j+1; i<d; i++, m++) 
		{
			Set::VectorSpace::Hom Aloc(m1,C[j][i]);
			Set::VectorSpace::Hom Bloc(m1,C[i][j]);
			Set::VectorSpace::Hom Cloc(n1,B[m].begin());
			Aloc = Bloc = df(Cloc);
		}
	}

	Indexing::Delete(C,m1*m1,d,d); 
	return A;
}

}

}

//////////////////////////////////////////////////////////////////////
// Sym Operators
//////////////////////////////////////////////////////////////////////

Set::VectorSpace::Sym 
operator + (const Set::VectorSpace::Sym &A, 
			const Set::VectorSpace::Sym &B)
{
	Set::VectorSpace::Sym C=A; C += B; return C;
}

Set::VectorSpace::Sym 
operator + (const Set::VectorSpace::Sym &A, 
			const Set::VectorSpace::SymId &B)
{
	Set::VectorSpace::Sym C=A; double *p; 
	for (p=C.begin(); p<C.begin()+C.size1(); p++) *p += 1.0;
	return C;
}

Set::VectorSpace::Sym 
operator + (const Set::VectorSpace::SymId &A, 
			const Set::VectorSpace::Sym &B)
{
	Set::VectorSpace::Sym C=B; double *p; 
	for (p=C.begin(); p<C.begin()+C.size1(); p++) *p += 1.0;
	return C;
}

Set::VectorSpace::Sym 
operator + (const Set::VectorSpace::Sym &A, 
			const Set::VectorSpace::SymZero &B)
{
	return A;
}

Set::VectorSpace::Sym 
operator + (const Set::VectorSpace::SymZero &A, 
			const Set::VectorSpace::Sym &B)
{
	return B;
}

Set::VectorSpace::Sym 
operator - (const Set::VectorSpace::Sym &A, 
			const Set::VectorSpace::Sym &B)
{
	Set::VectorSpace::Sym C=A; C -= B; return C;
}

Set::VectorSpace::Sym 
operator - (const Set::VectorSpace::Sym &A, 
			const Set::VectorSpace::SymId &B)
{
	Set::VectorSpace::Sym C=A; double *p; 
	for (p=C.begin(); p<C.begin()+C.size1(); p++) *p -= 1.0;
	return C;
}

Set::VectorSpace::Sym 
operator - (const Set::VectorSpace::SymId &A, 
			const Set::VectorSpace::Sym &B)
{
	Set::VectorSpace::Sym C=-B; double *p; 
	for (p=C.begin(); p<C.begin()+C.size1(); p++) *p += 1.0;
	return C;
}

Set::VectorSpace::Sym 
operator - (const Set::VectorSpace::Sym &A, 
			const Set::VectorSpace::SymZero &B)
{
	return A;
}

Set::VectorSpace::Sym 
operator - (const Set::VectorSpace::SymZero &A, 
			const Set::VectorSpace::Sym &B)
{
	return B;
}

Set::VectorSpace::Sym 
operator - (const Set::VectorSpace::Sym &A)
{
	Set::VectorSpace::Sym B=A; double *p; 
	for (p=B.begin(); p!=B.end(); p++) *p = -(*p);
	return B;
}

void Null(Set::VectorSpace::Sym &A)
{
	double *p; 
	for (p=A.begin(); p!=A.end(); p++) *p = 0.0;
}

Set::VectorSpace::Sym
operator * (const Set::VectorSpace::Sym &A, 
			const double &p)
{
	Set::VectorSpace::Sym B=A; B *= p; return B;
}

Set::VectorSpace::Sym 
operator * (const double &p, 
			const Set::VectorSpace::Sym &A)
{
	Set::VectorSpace::Sym B=A; B *= p; return B;
}

Set::VectorSpace::Sym 
operator / (const Set::VectorSpace::Sym &A, 
			const double &p)
{
	Set::VectorSpace::Sym B=A; B /= p; return B;
}

Set::VectorSpace::Vector
operator * (const Set::VectorSpace::Sym &A, 
			const Set::VectorSpace::Vector &B)
{
	return A.Embed()(B);
}

Set::VectorSpace::Sym
operator * (const Set::VectorSpace::Sym &A, 
			const Set::VectorSpace::Sym &B)
{
	return A.Embed()*B.Embed();
}

void Identity(Set::VectorSpace::Sym &A)
{
	double *p; 
	for (p=A.begin(); p<A.begin()+A.size1(); p++) *p = 1.0;
}

Set::VectorSpace::SymDual
Inverse(const Set::VectorSpace::Sym &A)
{
	return Inverse(A.Embed());
}

double Norm(const Set::VectorSpace::Sym &A)
{
	return Norm(A.Embed());
}

Set::VectorSpace::Sym
Adjoint(const Set::VectorSpace::Sym &A)
{
	return A;
}

double
Jacobian(const Set::VectorSpace::Sym &A)
{
	return Jacobian(A.Embed());
}

double 
Trace(const Set::VectorSpace::Sym &A)
{
	double *p; double a = 0.0;
	for (p=A.begin(); p<A.begin()+A.size1(); p++) a += *p;
	return a;
}

Set::VectorSpace::Sym
Deviator(const Set::VectorSpace::Sym &A)
{
	double *p; double a = 0.0;
	for (p=A.begin(); p<A.begin()+A.size1(); p++) a += *p;
	a /= 3.0; Set::VectorSpace::Sym B=A;
	for (p=B.begin(); p<B.begin()+B.size1(); p++) *p -= a;
	return B;
}

//////////////////////////////////////////////////////////////////////
// SymDual Operators
//////////////////////////////////////////////////////////////////////

Set::VectorSpace::SymDual 
operator + (const Set::VectorSpace::SymDual &A, 
			const Set::VectorSpace::SymDual &B)
{
	Set::VectorSpace::SymDual C=A; C += B; return C;
}

Set::VectorSpace::SymDual 
operator + (const Set::VectorSpace::SymDual &A, 
			const Set::VectorSpace::SymId &B)
{
	Set::VectorSpace::SymDual C=A; double *p; 
	for (p=C.begin(); p<C.begin()+C.size1(); p++) *p += 1.0;
	return C;
}

Set::VectorSpace::SymDual 
operator + (const Set::VectorSpace::SymId &A, 
			const Set::VectorSpace::SymDual &B)
{
	Set::VectorSpace::SymDual C=B; double *p; 
	for (p=C.begin(); p<C.begin()+C.size1(); p++) *p += 1.0;
	return C;
}

Set::VectorSpace::SymDual 
operator + (const Set::VectorSpace::SymDual &A, 
			const Set::VectorSpace::SymZero &B)
{
	return A;
}

Set::VectorSpace::SymDual 
operator + (const Set::VectorSpace::SymZero &A, 
			const Set::VectorSpace::SymDual &B)
{
	return B;
}

Set::VectorSpace::SymDual 
operator - (const Set::VectorSpace::SymDual &A, 
			const Set::VectorSpace::SymDual &B)
{
	Set::VectorSpace::SymDual C=A; C -= B; return C;
}

Set::VectorSpace::SymDual 
operator - (const Set::VectorSpace::SymDual &A, 
			const Set::VectorSpace::SymId &B)
{
	Set::VectorSpace::SymDual C=A; double *p; 
	for (p=C.begin(); p<C.begin()+C.size1(); p++) *p -= 1.0;
	return C;
}

Set::VectorSpace::SymDual
operator - (const Set::VectorSpace::SymId &A, 
			const Set::VectorSpace::SymDual &B)
{
	Set::VectorSpace::SymDual C=-B; double *p; 
	for (p=C.begin(); p<C.begin()+C.size1(); p++) *p += 1.0;
	return C;
}

Set::VectorSpace::SymDual 
operator - (const Set::VectorSpace::SymDual &A, 
			const Set::VectorSpace::SymZero &B)
{
	return A;
}

Set::VectorSpace::SymDual
operator - (const Set::VectorSpace::SymZero &A, 
			const Set::VectorSpace::SymDual &B)
{
	return B;
}

Set::VectorSpace::SymDual 
operator - (const Set::VectorSpace::SymDual &A)
{
	Set::VectorSpace::SymDual B=A; double *p; 
	for (p=B.begin(); p!=B.end(); p++) *p = -(*p);
	return B;
}

void Null(Set::VectorSpace::SymDual &A)
{
	double *p; 
	for (p=A.begin(); p!=A.end(); p++) *p = 0.0;
}

Set::VectorSpace::SymDual
operator * (const Set::VectorSpace::SymDual &A, 
			const double &p)
{
	Set::VectorSpace::SymDual B=A; B *= p; return B;
}

Set::VectorSpace::SymDual 
operator * (const double &p, 
			const Set::VectorSpace::SymDual &A)
{
	Set::VectorSpace::SymDual B=A; B *= p; return B;
}

Set::VectorSpace::SymDual 
operator / (const Set::VectorSpace::SymDual &A, 
			const double &p)
{
	Set::VectorSpace::SymDual B=A; B /= p; return B;
}

Set::VectorSpace::Vector
operator * (const Set::VectorSpace::SymDual &A, 
			const Set::VectorSpace::Vector &B)
{
	return A.Embed()(B);
}

Set::VectorSpace::SymDual
operator * (const Set::VectorSpace::SymDual &A, 
			const Set::VectorSpace::SymDual &B)
{
	return A.Embed()*B.Embed();
}

void Identity(Set::VectorSpace::SymDual &A)
{
	double *p; 
	for (p=A.begin(); p<A.begin()+A.size1(); p++) *p = 1.0;
}

Set::VectorSpace::Sym
Inverse(const Set::VectorSpace::SymDual &A)
{
	return Inverse(A.Embed());
}

double Norm(const Set::VectorSpace::SymDual &A)
{
	return Norm(A.Embed());
}

Set::VectorSpace::SymDual
Adjoint(const Set::VectorSpace::SymDual &A)
{
	return A;
}

double
Jacobian(const Set::VectorSpace::SymDual &A)
{
	return Jacobian(A.Embed());
}

double 
Trace(const Set::VectorSpace::SymDual &A)
{
	double *p; double a = 0.0;
	for (p=A.begin(); p<A.begin()+A.size1(); p++) a += *p;
	return a;
}

Set::VectorSpace::SymDual
Deviator(const Set::VectorSpace::SymDual &A)
{
	double *p; double a = 0.0;
	for (p=A.begin(); p<A.begin()+A.size1(); p++) a += *p;
	a /= 3.0; Set::VectorSpace::SymDual B=A;
	for (p=B.begin(); p<B.begin()+B.size1(); p++) *p -= a;
	return B;
}

//////////////////////////////////////////////////////////////////////
// SymZero Operators
//////////////////////////////////////////////////////////////////////

Set::VectorSpace::SymZero 
operator + (const Set::VectorSpace::SymZero &A, 
			const Set::VectorSpace::SymZero &B)
{
	return A;
}

Set::VectorSpace::SymZero 
operator - (const Set::VectorSpace::SymZero &A, 
			const Set::VectorSpace::SymZero &B)
{
	return A;
}
