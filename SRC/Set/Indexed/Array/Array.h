// Array.h: interface for the Array class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(SET_INDEXED_ARRAY_H__INCLUDED_)
#define SET_INDEXED_ARRAY_H__INCLUDED_

#pragma once

#include <cstdlib>
#include <cassert>
#include <math.h>
#include <iostream>

using namespace std;

namespace Set
{
//////////////////////////////////////////////////////////////////////
// Class Array
//////////////////////////////////////////////////////////////////////

class Array
{
public:

        inline
        Array() : sub(false), n(0), head(0), tail(0){}

	inline
	Array(unsigned int n0): 
	  sub(false), n(n0), head(new double [n0]), tail(head+n0) {
	  double *p; for (p=head; p<tail; p++) *p=0.0;
	}

	inline
	Array(unsigned int n0, double u):
	  sub(false), n(n0), head(new double [n0]), tail(head+n0) {
	  double *p; for (p=head; p<tail; p++) *p=u;
	}

	inline
	Array(unsigned int n0, double * pu) : sub(true), n(n0), head(pu), tail(head+n0){}

	inline
        virtual ~Array() { if (!sub) delete [] head; }

	inline
	Array(const Array &A): sub(false), n(A.n), head(new double [n]), tail(head+n) {
	    double *p, *q; 
	    for (p=head, q=A.head; p<tail; p++, q++) *p = *q;
	}

	inline
	Array & operator = (const Array & A) {
	  if (this == &A) return *this;
	  assert(n == A.n); double *p, *q; 
	  for (p=head, q=A.head; p<tail; p++, q++) *p = *q;
	  return *this;
	}

        inline
	double * const & begin() const { 
	  return head;
	}

        inline
	double * const & end() const { 
	  return tail; 
	}

        inline
	const double & operator [] (unsigned int i) const {
	  assert(i<n);
	  return *(head + i);
	}

        inline
	double & operator [] (unsigned int i) {
	  assert(i<n);
	  return *(head + i);
	}

        inline
	void Randomize() {
	  for (double *p=head; p<tail; p++) 
	    *p = 2.0*((double)rand()/(double)RAND_MAX) - 1.0;
	}

        inline
	void print(ostream * os) const {
	  for (double *p=head; p<tail; p++) *os << *p << ", ";
	}

	inline unsigned int size() const { 
	  return n; 
	}

	inline
	virtual void write(ostream & os) const {
	  os.write((char*)&n, sizeof(int));
	  os.write((char*)head, n * sizeof(double));
	}

	inline
	virtual void read(istream & is){
	  sub = false;
	  is.read((char*)&n, sizeof(int));	  
	  head = new double[n]; tail = head + n;
	  is.read((char*)head, n * sizeof(double));
	}

protected:

	bool sub;
	unsigned int n;
	double *head;
	double *tail;
};

}

//////////////////////////////////////////////////////////////////////
// Set operators
//////////////////////////////////////////////////////////////////////

bool operator != (const Set::Array &A, const Set::Array &B);
bool operator == (const Set::Array &A, const Set::Array &B);
bool operator < (const Set::Array &A, const Set::Array &B);
void Random(Set::Array &A);

//////////////////////////////////////////////////////////////////////
// Printing
//////////////////////////////////////////////////////////////////////

ostream & operator<<(ostream &os, const Set::Array &A);

#endif // !defined(SET_INDEXED_ARRAY_H__INCLUDED_)
