// Mesh.h: interface for the Mesh class.
// Copyright (c) 2006 by Bo Li - All rights reserved
//////////////////////////////////////////////////////////////////////

#if !defined(GEOMETRY_MESH_H__INCLUDED_)
#define GEOMETRY_MESH_H__INCLUDED_

#pragma once

#include <map>
#include <vector>
#include "../../Algebraic/CellComplex/CellComplex.h"
#include "../../../Set/Manifold/Euclidean/Orthonormal/Orthonormal.h"

using namespace std;

namespace Geometry
{
class Mesh
{
public:

	Mesh(){}
	virtual ~Mesh(){}
	virtual const CellComplex & GetCellComplex() const=0;
	virtual CellComplex & GetCellComplex()=0;
	virtual const map<Cell *, Cell *> & GetBinding() const=0;
	virtual map<Cell *, Cell *> & GetBinding()=0;
	virtual const map<Cell *, Set::Euclidean::Orthonormal::Point> & GetPointSet() const=0;
	virtual map<Cell *, Set::Euclidean::Orthonormal::Point> & GetPointSet()=0;

private:

	Mesh(const Mesh &);
	Mesh & operator = (const Mesh &);	
};

}

#endif // !defined(GEOMETRY_MESH_H__INCLUDED_)
