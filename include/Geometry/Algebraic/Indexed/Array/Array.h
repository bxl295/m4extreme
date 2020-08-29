// Array.h: interface for the Array class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(GEOMETRY_ALGEBRAIC_ARRAY_H__INCLUDED_)
#define GEOMETRY_ALGEBRAIC_ARRAY_H__INCLUDED_

#pragma once

#include <stdlib.h>
#include <iostream>
using namespace std;

namespace Geometry
{
namespace Algebraic
{
//////////////////////////////////////////////////////////////////////
// Class Array
//////////////////////////////////////////////////////////////////////

class Array
{
public:

	Array();
	Array(const unsigned int &);
	Array(const unsigned int &, const int &);
	Array(const unsigned int &, int * const &u);
	virtual ~Array();
	Array(const Array &A);
	Array & operator = (const Array &);
	int * const & begin() const;
	int * const & end() const;
	const int & operator [] (const unsigned int &) const;
	int & operator [] (const unsigned int &);
	void Randomize();
	void print(ostream *);
	unsigned int size() const;

protected:

	bool sub;
	unsigned int n;
	int *head;
	int *tail;
};

}

}

bool operator != (const Geometry::Algebraic::Array &, const Geometry::Algebraic::Array &);
bool operator == (const Geometry::Algebraic::Array &, const Geometry::Algebraic::Array &);
bool operator <  (const Geometry::Algebraic::Array &, const Geometry::Algebraic::Array &);
void Random(Geometry::Algebraic::Array &);
ostream & operator<<(ostream &, Geometry::Algebraic::Array &);

#endif // !defined(GEOMETRY_ALGEBRAIC_ARRAY_H__INCLUDED_)
