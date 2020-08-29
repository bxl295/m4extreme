// SimplicialComplex.h: interface for the SimplicialComplex class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(GEOMETRY_SIMPLICIALCOMPLEX_H__INCLUDED_)
#define GEOMETRY_SIMPLICIALCOMPLEX_H__INCLUDED_

#pragma once

#include <set>
#include <vector>
#include "../CellComplex/CellComplex.h"
#include "../ChainComplex/ChainComplex.h"

using namespace std;

namespace Geometry
{
//////////////////////////////////////////////////////////////////////
// Class SimplicialComplex
//////////////////////////////////////////////////////////////////////

class SimplicialComplex: public CellComplex
{
public:

	SimplicialComplex();
	SimplicialComplex(const unsigned int &);
	virtual ~SimplicialComplex();
	void ReadTable(const set<set<Cell *> > &);
	vector<Cell *> ReadTable(const vector<vector<unsigned int> > &);
	void Reconnect(const vector<vector<unsigned int> > &, const vector<Cell *> &);

private:

	SimplicialComplex(const SimplicialComplex &);
	SimplicialComplex & operator = (const SimplicialComplex &);
};

}

#endif // !defined(GEOMETRY_SIMPLICIALCOMPLEX_H__INCLUDED_)
