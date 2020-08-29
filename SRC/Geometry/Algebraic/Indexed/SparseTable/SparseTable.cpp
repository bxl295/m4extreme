// SparseTable.cpp: implementation of the SparseTable class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <cassert>
#include "./SparseTable.h"

namespace Geometry
{
namespace Algebraic
{
//////////////////////////////////////////////////////////////////////
// Member methods
//////////////////////////////////////////////////////////////////////

SparseTable::SparseTable() {}

SparseTable::~SparseTable() {}

SparseTable::SparseTable(const SparseTable &T) : 
	map<unsigned int, Geometry::Algebraic::SparseArray>(T) {}

SparseTable & SparseTable::operator = (const SparseTable &T)
{
	if (this == &T) return *this;
	map<unsigned int, Geometry::Algebraic::SparseArray>::operator = (T);
	return *this;
}

void SparseTable::print(ostream *os)
{
	map<unsigned int, Geometry::Algebraic::SparseArray>::iterator pT;
	for (pT=this->begin(); pT!=this->end(); pT++) *os << pT->second;
}

}

}

//////////////////////////////////////////////////////////////////////
// Set operators
//////////////////////////////////////////////////////////////////////

void Min(
	const Geometry::Algebraic::SparseTable::iterator &pA_beg, 
	const Geometry::Algebraic::SparseTable::iterator &pA_end, 
	unsigned int &imin, unsigned int &jmin, int &alpha, int &beta)
{
	unsigned int i, j;
	int aij, gamma;

	Geometry::Algebraic::SparseTable::iterator pA;
	Geometry::Algebraic::SparseArray::iterator ppA;

	alpha = 0;

	for (pA=pA_beg; pA!=pA_end; pA++)
	{
		j = pA->first;
		for (ppA=pA->second.begin(); ppA!=pA->second.end(); ppA++)
		{
			i = ppA->first; 
			aij = ppA->second; 
			gamma = abs(aij);
			if (alpha == 0)
			{
				if (gamma > 0) 
				{
					alpha = gamma; 
					beta = aij; 
					imin = i; 
					jmin = j;
				}
			}
			else 
			{
				if ((gamma > 0) && (gamma < alpha)) 
				{
					alpha = gamma; 
					beta = aij; 
					imin = i; 
					jmin = j;
				}
			}
		}
	}
}

Geometry::Algebraic::SparseTable 
Transpose(const Geometry::Algebraic::SparseTable &S)
{
	Geometry::Algebraic::SparseTable T;
	Geometry::Algebraic::SparseTable::const_iterator pS;
	Geometry::Algebraic::SparseArray::const_iterator ppS;
	for (pS=S.begin(); pS!=S.end(); pS++)
		for (ppS=pS->second.begin(); ppS!=pS->second.end(); ppS++)
			T[ppS->first][pS->first] = ppS->second;
	return T;
}

//////////////////////////////////////////////////////////////////////
// Elementary operations
//////////////////////////////////////////////////////////////////////

void Row(
	Geometry::Algebraic::SparseTable &S,
	Geometry::Algebraic::SparseTable &T,
	unsigned int row1, unsigned int row2, int q)
{
	assert(row1 != row2); T[row1] += q*T[row2];
	Geometry::Algebraic::SparseArray::iterator ppT;
	for (ppT=T[row1].begin(); ppT!=T[row1].end(); ppT++)
		S[ppT->first][row1] = ppT->second;
}

void Row(
	Geometry::Algebraic::SparseTable &S,
	Geometry::Algebraic::SparseTable &T,
	unsigned int row1, unsigned int row2)
{
	if (row1 == row2) return;
	Geometry::Algebraic::SparseArray::const_iterator ppT;
	for (ppT=T[row1].begin(); ppT!=T[row1].end(); ppT++)
		S[ppT->first].erase(row1);
	for (ppT=T[row2].begin(); ppT!=T[row2].end(); ppT++)
		S[ppT->first].erase(row2);
	Geometry::Algebraic::SparseArray A = T[row1];
	T[row1] = T[row2]; T[row2] = A;
	for (ppT=T[row1].begin(); ppT!=T[row1].end(); ppT++)
		S[ppT->first][row1] = ppT->second;
	for (ppT=T[row2].begin(); ppT!=T[row2].end(); ppT++)
		S[ppT->first][row2] = ppT->second;
}

void Row(
	Geometry::Algebraic::SparseTable &S,
	Geometry::Algebraic::SparseTable &T,
	unsigned int row)
{
	Geometry::Algebraic::SparseArray::iterator ppT;
	for (ppT=T[row].begin(); ppT!=T[row].end(); ppT++)
	{
		ppT->second = -ppT->second;
		S[ppT->first][row] = ppT->second;
	}
}

void Col(
	Geometry::Algebraic::SparseTable &S,
	Geometry::Algebraic::SparseTable &T,
	unsigned int col1, unsigned int col2, int q)
{
	assert(col1 != col2); S[col1] += q*S[col2];
	Geometry::Algebraic::SparseArray::iterator ppS;
	for (ppS=S[col1].begin(); ppS!=S[col1].end(); ppS++)
		T[ppS->first][col1] = ppS->second;
}

void Col(
	Geometry::Algebraic::SparseTable &S,
	Geometry::Algebraic::SparseTable &T,
	unsigned int col1, unsigned int col2)
{
	if (col1 == col2) return;
	Geometry::Algebraic::SparseArray::const_iterator ppS;
	for (ppS=S[col1].begin(); ppS!=S[col1].end(); ppS++)
		T[ppS->first].erase(col1);
	for (ppS=S[col2].begin(); ppS!=S[col2].end(); ppS++)
		T[ppS->first].erase(col2);
	Geometry::Algebraic::SparseArray A = S[col1];
	S[col1] = S[col2]; S[col2] = A;
	for (ppS=S[col1].begin(); ppS!=S[col1].end(); ppS++)
		T[ppS->first][col1] = ppS->second;
	for (ppS=S[col2].begin(); ppS!=S[col2].end(); ppS++)
		T[ppS->first][col2] = ppS->second;
}

void Col(
	Geometry::Algebraic::SparseTable &S,
	Geometry::Algebraic::SparseTable &T,
	unsigned int col)
{
	Geometry::Algebraic::SparseArray::iterator ppS;
	for (ppS=S[col].begin(); ppS!=S[col].end(); ppS++)
	{
		ppS->second = -ppS->second;
		T[ppS->first][col] = ppS->second;
	}
}


//////////////////////////////////////////////////////////////////////
// Printing
//////////////////////////////////////////////////////////////////////

ostream & 
operator<<(ostream &os, Geometry::Algebraic::SparseTable &T)
{
	T.print(&os); return os;
}
