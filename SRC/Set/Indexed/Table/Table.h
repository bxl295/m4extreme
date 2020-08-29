// Table.h: interface for the Table class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(SET_INDEXED_TABLE_H__INCLUDED_)
#define SET_INDEXED_TABLE_H__INCLUDED_

#pragma once

#include "../Array/Array.h"

namespace Set
{
//////////////////////////////////////////////////////////////////////
// Class Table
//////////////////////////////////////////////////////////////////////

class Table : public Set::Array
{
public:

	Table();
	Table(const unsigned int &, const unsigned int &);
	Table(const unsigned int &, const unsigned int &, double * const &);
	virtual ~Table();
	Table(const Table &);
	Table & operator = (const Table &);
	Table & operator = (const Set::Array &);
	const Set::Array & operator [] (const unsigned int &) const;
	Set::Array & operator [] (const unsigned int &);
	void print(ostream *);
	unsigned int size1() const;
	unsigned int size2() const;

protected:

	unsigned int n1;
	unsigned int n2;
	Set::Array **L;
};

}

//////////////////////////////////////////////////////////////////////
// Printing
//////////////////////////////////////////////////////////////////////

template<unsigned int n, class T>
ostream & operator<<(ostream &os, Set::Table &A)
{
	A.print(&os); return os;
}

#endif // !defined(SET_INDEXED_TABLE_H__INCLUDED_)
