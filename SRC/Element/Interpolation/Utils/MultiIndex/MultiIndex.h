// MultiIndex.h: interface for the MultiIndex class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(ELEMENT_INTERPOLATION_UTILS_MULTIINDEX_H__INCLUDED_)
#define ELEMENT_INTERPOLATION_UTILS_MULTIINDEX_H__INCLUDED_

#pragma once

#include <cassert>
#include <iostream>

using namespace std;

namespace Element
{
namespace Interpolation
{
namespace Utils
{
//////////////////////////////////////////////////////////////////////
// Class MultiIndex
//////////////////////////////////////////////////////////////////////

class MultiIndex
{
public:

	MultiIndex();
	MultiIndex(const unsigned int &);
	virtual ~MultiIndex();
	MultiIndex(const MultiIndex &A);
	MultiIndex & operator = (const MultiIndex &);
	unsigned int * const & begin() const;
	unsigned int * const & end() const;
	const unsigned int & operator [] (const unsigned int &) const;
	unsigned int & operator [] (const unsigned int &);
	bool operator <  (const MultiIndex &) const;
	bool operator == (const MultiIndex &) const;
	bool operator != (const MultiIndex &) const;
	void print(ostream *) const;
	unsigned int size() const;
	unsigned int degree() const;

protected:

	unsigned int n;
	unsigned int *head;
	unsigned int *tail;
};

//////////////////////////////////////////////////////////////////////
// Utilities
//////////////////////////////////////////////////////////////////////

MultiIndex Reduce(const MultiIndex &);

}

}

}

//////////////////////////////////////////////////////////////////////
// Class MultiIndex: Operators
//////////////////////////////////////////////////////////////////////

ostream & operator<<(ostream &, 
	const Element::Interpolation::Utils::MultiIndex &);

#endif // !defined(ELEMENT_INTERPOLATION_UTILS_MULTIINDEX_H__INCLUDED_

