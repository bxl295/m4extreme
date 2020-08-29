// NormalForm.h: interface for the NormalForm class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(GEOMETRY_ALGEBRAIC_NORMALFORM_H__INCLUDED_)
#define GEOMETRY_ALGEBRAIC_NORMALFORM_H__INCLUDED_

#pragma once

#include "../../Indexed/Array/Array.h"
#include "../../Indexed/Table/Table.h"

using namespace std;

namespace Geometry
{
namespace Algebraic
{
namespace Utils
{
//////////////////////////////////////////////////////////////////////
// Class SubTable
//////////////////////////////////////////////////////////////////////

class SubTable
{
public:

	SubTable();
	SubTable(
		const Geometry::Algebraic::Table &, 
		const unsigned int &);
	virtual ~SubTable();
	const Geometry::Algebraic::Array & 
		operator [] (const unsigned int &) const;
	Geometry::Algebraic::Array & 
		operator [] (const unsigned int &);
	void print(ostream *);
	unsigned int size1() const;
	unsigned int size2() const;	

protected:

	unsigned int n1;
	unsigned int n2;
	Geometry::Algebraic::Array **L;	

private:

	SubTable(const SubTable &);
	SubTable & operator = (const SubTable &);
};

//////////////////////////////////////////////////////////////////////
// Class NormalForm
//////////////////////////////////////////////////////////////////////

class NormalForm 
{
public:

	NormalForm();
	virtual ~NormalForm();
	void operator () (Geometry::Algebraic::Table &);
	void operator () (
		Geometry::Algebraic::Table &,
		Geometry::Algebraic::Table &,
		Geometry::Algebraic::Table &);

private:

	void Row(SubTable &, unsigned int, unsigned int);
	void Row(SubTable &, unsigned int);
	void Row(SubTable &, unsigned int, unsigned int, unsigned int);
	void Col(SubTable &, unsigned int, unsigned int);
	void Col(SubTable &, unsigned int);
	void Col(SubTable &, unsigned int, unsigned int, unsigned int);
	void Min(SubTable &, unsigned int &, unsigned int &, unsigned int &);

private:

	NormalForm(const NormalForm &);
	NormalForm & operator = (const NormalForm &);	
};

}

}

}

ostream & operator<<(ostream &os, Geometry::Algebraic::Utils::SubTable &);
void Null(Geometry::Algebraic::Table &);
void Identity(Geometry::Algebraic::Table &);

#endif // !defined(GEOMETRY_ALGEBRAIC_NORMALFORM_H__INCLUDED_)
