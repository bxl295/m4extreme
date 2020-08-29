// Table.cpp: implementation of the Table class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <cassert>
#include "./Table.h"

namespace Geometry
{
namespace Algebraic
{
//////////////////////////////////////////////////////////////////////
// Member methods
//////////////////////////////////////////////////////////////////////

Table::Table() : n1(0), n2(0), L(0){}

Table::Table(
	const unsigned int &m1, 
	const unsigned int &m2) : 
	Geometry::Algebraic::Array(m1*m2), 
	n1(m1), n2(m2), 
	L(new Geometry::Algebraic::Array * [n2])
{
	unsigned int i; int *p;
	for (i=0, p=this->head; i<n2; i++, p+=n1) 
		L[i] = new Geometry::Algebraic::Array(n1,p);
}

Table::Table(
	const unsigned int &m1, 
	const unsigned int &m2, 
	int * const &u) : 
	Geometry::Algebraic::Array(m1*m2,u), 
	n1(m1), n2(m2), 
	L(new Geometry::Algebraic::Array * [n2])
{
	unsigned int i; int *p;
	for (i=0, p=this->head; i<n2; i++, p+=n1) 
		L[i] = new Geometry::Algebraic::Array(n1,p);
}

Table::~Table()
{
	for (unsigned int i=0; i<n2; i++) delete L[i];
	delete [] L; 
}

Table::Table(const Table &A) : 
	Geometry::Algebraic::Array(A), 
	n1(A.n1), n2(A.n2), 
	L(new Geometry::Algebraic::Array * [n2])
{
	unsigned int i; int *p;
	for (i=0, p=this->head; i<n2; i++, p+=n1) 
		L[i] = new Geometry::Algebraic::Array(n1,p);
}

Table & 
Table::operator = (const Table &A)
{
	if (this == &A) return *this;
	assert(n1 == A.n1); assert(n2 == A.n2); 
	Geometry::Algebraic::Array::operator = (A);
	return *this;
}

Table & 
Table::operator = (const Geometry::Algebraic::Array &A)
{
	assert(n == A.size()); 
	Geometry::Algebraic::Array::operator = (A);
	return *this;
}

const Geometry::Algebraic::Array & 
Table::operator [] (const unsigned int &i) const
{
	assert(i < n2); return *L[i];
}

Geometry::Algebraic::Array & 
Table::operator [] (const unsigned int &i)
{
	assert(i < n2); return *L[i];
}

void 
Table::print(ostream *os)
{
	for (unsigned int i=0; i<n2; i++) *os << *L[i];
}

unsigned int 
Table::size1() const
{
	return n1;
}

unsigned int 
Table::size2() const
{
	return n2;
}

}

}

//////////////////////////////////////////////////////////////////////
// Printing
//////////////////////////////////////////////////////////////////////

ostream & 
operator<<(ostream &os, Geometry::Algebraic::Table &A)
{
	A.print(&os); return os;
}
