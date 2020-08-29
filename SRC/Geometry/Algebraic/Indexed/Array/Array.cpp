// Array.cpp: implementation the Array class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <cassert>
#include <iostream>
#include "./Array.h"

using namespace std;

namespace Geometry
{
namespace Algebraic
{
//////////////////////////////////////////////////////////////////////
// Member methods
//////////////////////////////////////////////////////////////////////

Array::Array() : sub(false), n(0), head(0), tail(0){}

Array::Array(const unsigned int &n0) : 
	sub(false), n(n0), head(new int [n]), tail(head+n)
{
	int *p; for (p=head; p<tail; p++) *p=0;
}

Array::Array(const unsigned int &n0, const int &u) : 
	sub(false), n(n0), head(new int [n]), tail(head+n)
{
	int *p; for (p=head; p<tail; p++) *p=u;
}

Array::Array(const unsigned int &n0, int * const &u) : 
	sub(true), n(n0), head(u), tail(head+n){}

Array::~Array(){if(!sub){delete [] head;};}

Array::Array(const Array &A) : 
	sub(false), n(A.n), head(new int [n]), tail(head+n)
{
	int *p, *q; 
	for (p=head, q=A.head; p<tail; p++, q++) *p = *q;
}

Array & 
Array::operator = (const Array &A)
{
	if (this == &A) return *this;
    assert(n == A.n); int *p, *q; 
	for (p=head, q=A.head; p<tail; p++, q++) *p = *q;
	return *this;
}

int * const & 
Array::begin() const{return head;}

int * const & 
Array::end() const{return tail;}

const int & 
Array::operator [] (const unsigned int &i) const
{
	assert(i < n); return *(head+i);
}

int & 
Array::operator [] (const unsigned int &i)
{
	assert(i < n); return *(head+i);
}

void 
Array::Randomize()
{
	for (int *p=head; p<tail; p++) 
		*p =(int)rand();
}

void 
Array::print(ostream *os)
{
	for (int *p=head; p<tail; p++) *os << *p << ", ";
}

unsigned int 
Array::size() const{return n;}

}

}

//////////////////////////////////////////////////////////////////////
// Set operators
//////////////////////////////////////////////////////////////////////

bool 
operator != (
	const Geometry::Algebraic::Array &A, 
	const Geometry::Algebraic::Array &B)
{
	int *p, *q;
	for (p=A.begin(),q=B.begin(); p<A.end(); p++,q++) 
		if ((*p) != (*q)) return true;
	return false;
}

bool 
operator == (
	const Geometry::Algebraic::Array &A, 
	const Geometry::Algebraic::Array &B)
{
	int *p, *q;
	for (p=A.begin(),q=B.begin(); p<A.end(); p++,q++) 
		if ((*p) != (*q)) return false;
	return true;
}

bool 
operator < (
	const Geometry::Algebraic::Array &A, 
	const Geometry::Algebraic::Array &B)
{
	const int *p, *q, *pm, *qm;
	for (p=A.end(), pm=p, pm--, q=B.end(), qm=q, qm--; 
		p!=A.begin(); p--, pm--, q--, qm--)
	{
		if (*pm < *qm) return true;
		if (*pm > *qm) return false;
	}
	return false;
}

void 
Random(Geometry::Algebraic::Array &A)
{
	A.Randomize();
}

//////////////////////////////////////////////////////////////////////
// Printing
//////////////////////////////////////////////////////////////////////

ostream & 
operator<<(ostream &os, Geometry::Algebraic::Array &A)
{
	A.print(&os); return os;
}
