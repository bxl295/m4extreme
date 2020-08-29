// CellComplex.h: interface for the CellComplex class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(GEOMETRY_CELLCOMPLEX_H__INCLUDED_)
#define GEOMETRY_CELLCOMPLEX_H__INCLUDED_

#pragma once

#include <set>
#include <vector>
#include "../ChainComplex/ChainComplex.h"

using namespace std;

namespace Geometry
{
//////////////////////////////////////////////////////////////////////
// Class SubComplex
//////////////////////////////////////////////////////////////////////

class SubComplex : public vector<set<Cell *> >
{
public:

	SubComplex();
	SubComplex(const unsigned int &);
	virtual ~SubComplex();
	SubComplex(const SubComplex &);
	SubComplex(const vector<set<Cell *> > &);
	SubComplex & operator = (const SubComplex &);
	void SetBoundaries();
	void SetCoboundaries();
	const vector<set<Cell *> > Boundary(Cell * const) const;
	const vector<set<Cell *> > Incidence(Cell * const) const;
	const set<Cell *> Incidence(Cell * const, const int &) const;
	map<Cell *, Cell *> Bind(const SubComplex &) const;
	void Bind(const SubComplex &, map<Cell *, Cell *> &) const;
	const int dim() const;
	Chain RandomChain(const unsigned int &) const;
	Cochain RandomCochain(const unsigned int &) const;
	Chain NullChain(const unsigned int &) const;
	Cochain NullCochain(const unsigned int &) const;

	private:

	bool Contained(const set<Cell *> &, const vector<set<Cell *> > &) const;
};

//////////////////////////////////////////////////////////////////////
// Class CellComplex
//////////////////////////////////////////////////////////////////////

class CellComplex : public SubComplex
{
public:

	CellComplex();
	CellComplex(const unsigned int &);
	virtual ~CellComplex();
	CellComplex(const SubComplex &);

private:

	CellComplex(const CellComplex &);
	CellComplex & operator = (const CellComplex &);	
};

}

#endif // !defined(GEOMETRY_CELLCOMPLEX_H__INCLUDED_)
