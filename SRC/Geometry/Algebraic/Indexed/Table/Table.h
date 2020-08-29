// Table.h: interface for the Table class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(GEOMETRY_ALGEBRAIC_TABLE_H__INCLUDED_)
#define GEOMETRY_ALGEBRAIC_TABLE_H__INCLUDED_

#pragma once

#include "../Array/Array.h"

namespace Geometry
{
namespace Algebraic
{
//////////////////////////////////////////////////////////////////////
// Class Table
//////////////////////////////////////////////////////////////////////

class Table : public Geometry::Algebraic::Array
{
public:

	Table();
	Table(const unsigned int &, const unsigned int &);
	Table(const unsigned int &, const unsigned int &, int * const &);
	virtual ~Table();
	Table(const Table &);
	Table & operator = (const Table &);
	Table & operator = (const Geometry::Algebraic::Array &);
	const Geometry::Algebraic::Array & operator [] (const unsigned int &) const;
	Geometry::Algebraic::Array & operator [] (const unsigned int &);
	void print(ostream *);
	unsigned int size1() const;
	unsigned int size2() const;

protected:

	unsigned int n1;
	unsigned int n2;
	Geometry::Algebraic::Array **L;
};

}

}

//////////////////////////////////////////////////////////////////////
// Printing
//////////////////////////////////////////////////////////////////////

template<unsigned int n, class T>
ostream & operator<<(ostream &os, Geometry::Algebraic::Table &A)
{
	A.print(&os); return os;
}

#endif // !defined(GEOMETRY_ALGEBRAIC_TABLE_H__INCLUDED_)
