// SparseNormalForm.h: interface for the SparseNormalForm class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(GEOMETRY_ALGEBRAIC_UTILS_NORMALFORM_H__INCLUDED_)
#define GEOMETRY_ALGEBRAIC_UTILS_NORMALFORM_H__INCLUDED_

#pragma once

#include <map>
#include "../../Indexed/SparseTable/SparseTable.h"

using namespace std;

namespace Geometry
{
namespace Algebraic
{
namespace Utils
{
//////////////////////////////////////////////////////////////////////
// Class SparseNormalForm
//////////////////////////////////////////////////////////////////////

class SparseNormalForm 
{
public:

	SparseNormalForm();
	virtual ~SparseNormalForm();
	void operator () (
		Geometry::Algebraic::SparseTable &);
	void operator () (
		Geometry::Algebraic::SparseTable &,
		Geometry::Algebraic::SparseTable &,
		Geometry::Algebraic::SparseTable &);

private:

	SparseNormalForm(const SparseNormalForm &);
	SparseNormalForm & operator = (const SparseNormalForm &);	
};

}

}

}

#endif // !defined(GEOMETRY_ALGEBRAIC_UTILS_NORMALFORM_H__INCLUDED_)
