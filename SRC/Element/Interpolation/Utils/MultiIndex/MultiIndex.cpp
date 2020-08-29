// MultiIndex.cpp: implementation of the MultiIndex class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

//#include <cassert>
#include "./MultiIndex.h"

namespace Element
{
namespace Interpolation
{
namespace Utils
{
//////////////////////////////////////////////////////////////////////
// Class MultiIndex
//////////////////////////////////////////////////////////////////////

MultiIndex::MultiIndex() : n(0), head(0), tail(0){}

MultiIndex::MultiIndex(const unsigned int &n_) : 
	n(n_), head(new unsigned int [n]), tail(head+n)
{
	unsigned int *p; for (p=head; p<tail; p++) *p=0;
}

MultiIndex::~MultiIndex(){delete [] head;}

MultiIndex::MultiIndex(const MultiIndex &A) : 
	n(A.n), head(new unsigned int [n]), tail(head+n)
{
	unsigned int *p, *q; 
	for (p=head, q=A.head; p<tail; p++, q++) *p = *q;
}

MultiIndex & 
MultiIndex::operator = (const MultiIndex &A)
{
	if (this == &A) return *this;
    assert(n == A.n); unsigned int *p, *q; 
	for (p=head, q=A.head; p<tail; p++, q++) *p = *q;
	return *this;
}

unsigned int * const & 
MultiIndex::begin() const{return head;}

unsigned int * const & 
MultiIndex::end() const{return tail;}

const unsigned int & 
MultiIndex::operator [] (const unsigned int &i) const
{
	assert(i < n); return *(head+i);
}

unsigned int & 
MultiIndex::operator [] (const unsigned int &i)
{
	assert(i < n); return *(head+i);
}

bool 
MultiIndex::operator < (const MultiIndex &A) const 
{
	assert (A.size() == n);
	unsigned int *p, *q;
	for (p=head, q=A.begin(); p!=tail; p++, q++)
	{
		if (*p < *q) return true;
		if (*p > *q) return false;
	}
	return false;
}

bool 
MultiIndex::operator == (const MultiIndex &A) const 
{
	assert (A.size() == n);
	unsigned int *p, *q;
	for (p=head, q=A.begin(); p!=tail; p++, q++)
		if (*p != *q) return false;
	return true;
}

bool 
MultiIndex::operator != (const MultiIndex &A) const 
{
	assert (A.size() == n);
	unsigned int *p, *q;
	for (p=head, q=A.begin(); p!=tail; p++, q++)
		if (*p != *q) return true;
	return false;
}

void 
MultiIndex::print(ostream *os) const
{
	unsigned int *p;
	for (p=head; p<tail; p++) *os << *p << ", ";
}

unsigned int 
MultiIndex::size() const
{
	return n;
}

unsigned int 
MultiIndex::degree() const
{
	unsigned int k=0; unsigned int *p;
	for (p=head; p<tail; p++) k += *p;
	return k;
}

//////////////////////////////////////////////////////////////////////
// Utilities
//////////////////////////////////////////////////////////////////////

MultiIndex
Reduce(const MultiIndex &b)
{
	MultiIndex a = b;
	unsigned int i, m=a.degree();
	if (m == 0) return a;
	for (i=0; i<a.size(); i++)
		if ((a[i] != 0) && (a[i] < m)) m = a[i];
	if (m == 1) return a;
	for (i=0; i<a.size(); i++)
		if (a[i]%m != 0) return a;
	for (i=0; i<a.size(); i++) a[i] /= m;
	return a;
}

}

}

}

//////////////////////////////////////////////////////////////////////
// Class MultiIndex: operators
//////////////////////////////////////////////////////////////////////

ostream & operator<<(ostream &os, 
	const Element::Interpolation::Utils::MultiIndex &A)
{
	A.print(&os); return os;
}
