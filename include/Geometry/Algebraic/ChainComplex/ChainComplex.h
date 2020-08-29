// ChainComplex.h: interface for the ChainComplex class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(GEOMETRY_CHAINCOMPLEX_H__INCLUDED_)
#define GEOMETRY_CHAINCOMPLEX_H__INCLUDED_

#pragma once

#include <stdlib.h>
#include <map>
#include <iostream>

using namespace std;

namespace Geometry
{
class Cell;
class Chain;
class Cochain;

//////////////////////////////////////////////////////////////////////
// Class Chain
//////////////////////////////////////////////////////////////////////

class Chain : public map<Cell *,int>
{
public:

	Chain();
	Chain(const int &);
	virtual ~Chain();
	Chain(const Chain &);
	Chain & operator = (const Chain &);
	void operator += (const Chain &);
	void operator -= (const Chain &);
	void operator *= (const int &);
	int operator () (const Cochain &);
	Chain Boundary() const;
	void print(ostream *);
	const int & dim() const;

private:

	int n;
};

//////////////////////////////////////////////////////////////////////
// Class Cochain
//////////////////////////////////////////////////////////////////////

class Cochain : public map<Cell *,int>
{
public:

	Cochain();
	Cochain(const int &);
	virtual ~Cochain();
	Cochain(const Cochain &);
	Cochain & operator = (const Cochain &);
	void operator += (const Cochain &);
	void operator -= (const Cochain &);
	void operator *= (const int &);
	int operator () (const Chain &);
	Cochain Coboundary() const;
	void print(ostream *);
	const int & dim() const;

private:

	int n;
};

//////////////////////////////////////////////////////////////////////
// Class Cell
//////////////////////////////////////////////////////////////////////

class Cell 
{
public:

	Cell();
	Cell(const int &);
	virtual ~Cell();
	Cell(const Cell &);
	Cell & operator = (const Cell &);
	const Chain & Boundary() const;
	Chain & Boundary();
	const Cochain & Coboundary() const;
	Cochain & Coboundary();
	void print(ostream *);
	const int & dim() const;
	void Bo2Co();
	void Co2Bo();

private:

	int n;
	Chain bo;
	Cochain co;
};

}

Geometry::Chain Boundary(const Geometry::Cell &);
Geometry::Cochain Coboundary(const Geometry::Cell &);
ostream & operator<<(ostream &, Geometry::Cell &);
Geometry::Chain operator - (const Geometry::Chain &);
Geometry::Chain operator + (const Geometry::Chain &, const Geometry::Chain &);
Geometry::Chain operator - (const Geometry::Chain &, const Geometry::Chain &);
Geometry::Chain operator * (const Geometry::Chain &, const int &);
Geometry::Chain operator * (const int &, const Geometry::Chain &);
Geometry::Chain Boundary(Geometry::Chain &);
void Random(Geometry::Chain &);
ostream & operator<<(ostream &, Geometry::Chain &);
Geometry::Cochain operator - (const Geometry::Cochain &);
Geometry::Cochain operator + (const Geometry::Cochain &, const Geometry::Cochain &);
Geometry::Cochain operator - (const Geometry::Cochain &, const Geometry::Cochain &);
Geometry::Cochain operator * (const Geometry::Cochain &, const int &);
Geometry::Cochain operator * (const int &, const Geometry::Cochain &);
Geometry::Cochain Coboundary(Geometry::Cochain &);
void Random(Geometry::Cochain &);
ostream & operator<<(ostream &, Geometry::Cochain &);

#endif // !defined(GEOMETRY_CHAINCOMPLEX_H__INCLUDED_)
