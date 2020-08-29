// SparseTable.h: interface for the SparseTable class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(GEOMETRY_ALGEBRAIC_SPARSETABLE_H__INCLUDED_)
#define GEOMETRY_ALGEBRAIC_SPARSETABLE_H__INCLUDED_

#pragma once

#include <map>
#include <iostream>
#include "../SparseArray/SparseArray.h"

namespace Geometry
{
namespace Algebraic
{
//////////////////////////////////////////////////////////////////////
// Class SparseTable
//////////////////////////////////////////////////////////////////////

class SparseTable : public map<unsigned int, Geometry::Algebraic::SparseArray>
{
public:

	SparseTable();
	virtual ~SparseTable();
	SparseTable(const SparseTable &);
	SparseTable & operator = (const SparseTable &);
	void print(ostream *);
};

}

}

//////////////////////////////////////////////////////////////////////
// Set operators
//////////////////////////////////////////////////////////////////////

void Min(
	const Geometry::Algebraic::SparseTable::iterator &, 
	const Geometry::Algebraic::SparseTable::iterator &, 
	unsigned int &, unsigned int &, int &, int &);

Geometry::Algebraic::SparseTable
Transpose(const Geometry::Algebraic::SparseTable &); 

void
Null(Geometry::Algebraic::SparseTable &); 

void
Identity(Geometry::Algebraic::SparseTable &); 

//////////////////////////////////////////////////////////////////////
// Elementary operations
//////////////////////////////////////////////////////////////////////

void Row(
	Geometry::Algebraic::SparseTable &, 
	Geometry::Algebraic::SparseTable &, 
	unsigned int, unsigned int, int);
void Row(
	Geometry::Algebraic::SparseTable &, 
	Geometry::Algebraic::SparseTable &, 
	unsigned int, unsigned int);
void Row(
	Geometry::Algebraic::SparseTable &, 
	Geometry::Algebraic::SparseTable &, 
	unsigned int);
void Col(
	Geometry::Algebraic::SparseTable &, 
	Geometry::Algebraic::SparseTable &, 
	unsigned int, unsigned int, int);
void Col(
	Geometry::Algebraic::SparseTable &, 
	Geometry::Algebraic::SparseTable &, 
	unsigned int, unsigned int);
void Col(
	Geometry::Algebraic::SparseTable &, 
	Geometry::Algebraic::SparseTable &, 
	unsigned int);

//////////////////////////////////////////////////////////////////////
// Printing
//////////////////////////////////////////////////////////////////////

ostream & operator<<(ostream &os, Geometry::Algebraic::SparseTable &);

#endif // !defined(GEOMETRY_ALGEBRAIC_SPARSETABLE_H__INCLUDED_)
