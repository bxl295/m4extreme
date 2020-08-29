// SparseNormalForm.cpp: Implementation of the SparseNormalForm class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include "SparseNormalForm.h"

namespace Geometry
{
namespace Algebraic
{
namespace Utils
{
/////////////////////////////////////////////////////////////////////
// Class SparseNormalForm
//////////////////////////////////////////////////////////////////////

//	Reduction algorithm
//	Munkres, J.R., "Elements of Algebric Topology", 
//	Perseus Publishing, 1984, p.56.

SparseNormalForm::SparseNormalForm(){}

SparseNormalForm::~SparseNormalForm() {}

void
SparseNormalForm::operator () (
	Geometry::Algebraic::SparseTable &A)
{
	unsigned int imin, jmin;
	int alpha, beta, gamma;
	
	Geometry::Algebraic::SparseTable AT = Transpose(A);

	Geometry::Algebraic::SparseTable::iterator pA;
	Geometry::Algebraic::SparseArray::iterator ppA;

	Geometry::Algebraic::SparseTable::iterator pA_beg;
	Geometry::Algebraic::SparseTable::iterator pAT_beg;

	for (pA_beg=A.begin(), pAT_beg=AT.begin(); 
		(pA_beg!=A.end() && pAT_beg!=A.end()) ; pA_beg++, pAT_beg++)
	{
		while (true)
		{
			unsigned int i, j;
			Min(pA_beg,A.end(),imin,jmin,alpha,beta);
			if (alpha == 0)
			{
				for (pA=pA_beg; pA!=A.end(); pA++)
					pA->second.clear(); 
				for (pA=pAT_beg; pA!=AT.end(); pA++)
					pA->second.clear(); 
				return;
			}
			for (pA=pA_beg; pA!=A.end(); pA++)
			{
				j = pA->first;
				for (ppA=pA->second.begin(); ppA!=pA->second.end(); ppA++)
				{
					i = ppA->first; gamma = ppA->second;
					if (abs(gamma)%alpha != 0) goto divide;
				}
			}
			break; 
divide: 
			if (j == jmin)
			{
				int q = -gamma/beta;
				Row(A,AT,i,imin,q);
			}
			else if (i == imin)
			{
				int q = -gamma/beta;
				Col(A,AT,j,jmin,q);
			}
			else 
			{
				int q = -A[jmin][i]/beta;
				Row(A,AT,i,imin,q);
				Row(A,AT,imin,i,1);
			}
		}
        if (alpha == 0) break;
		unsigned int l = pA_beg->first;
		Row(A,AT,l,imin);
		Col(A,AT,l,jmin);
		int diag = A[l][l];
		if (diag < 0) Row(A,AT,l);
		diag = abs(diag);
		if (pA_beg->second.size() > 1)
			for (ppA=++pA_beg->second.begin(); ppA!=pA_beg->second.end(); ppA++)
			{
				int r = - ppA->second/diag;
				Row(A,AT,ppA->first,l,r); 
				AT[ppA->first].erase(l);
			}
		A[l].clear(); A[l][l] = diag;
		if (pAT_beg->second.size() > 1)
			for (ppA=++pAT_beg->second.begin(); ppA!=pAT_beg->second.end(); ppA++)
			{
				int r = - ppA->second/diag;
				Col(A,AT,ppA->first,l,r); 
				A[ppA->first].erase(l);
			}
		AT[l].clear(); AT[l][l] = diag;
	}
}

void
SparseNormalForm::operator () (
	Geometry::Algebraic::SparseTable &A,
	Geometry::Algebraic::SparseTable &S,
	Geometry::Algebraic::SparseTable &T)
{
	unsigned int imin, jmin;
	int alpha, beta, gamma;
	
	Geometry::Algebraic::SparseTable AT = Transpose(A);
	Geometry::Algebraic::SparseTable ST = Transpose(S);
	Geometry::Algebraic::SparseTable TT = Transpose(T);

	Geometry::Algebraic::SparseTable::iterator pA;
	Geometry::Algebraic::SparseArray::iterator ppA;

	Geometry::Algebraic::SparseTable::iterator pA_beg;
	Geometry::Algebraic::SparseTable::iterator pAT_beg;

	for (pA_beg=A.begin(), pAT_beg=AT.begin(); 
		 pA_beg!=A.end(); pA_beg++, pAT_beg++)
	{
		while (true)
		{
			unsigned int i, j;
			Min(pA_beg,A.end(),imin,jmin,alpha,beta);
			if (alpha == 0)
			{
				for (pA=pA_beg; pA!=A.end(); pA++)
					pA->second.clear(); 
				for (pA=pAT_beg; pA!=AT.end(); pA++)
					pA->second.clear(); 
				return;
			}
			for (pA=pA_beg; pA!=A.end(); pA++)
			{
				j = pA->first;
				for (ppA=pA->second.begin(); ppA!=pA->second.end(); ppA++)
				{
					i = ppA->first; gamma = ppA->second;
					if (abs(gamma)%alpha != 0) goto divide;
				}
			}
			break; divide: 
			if (j == jmin)
			{
				int q = -gamma/beta;
				Row(A,AT,i,imin,q); Row(S,ST,i,imin,q);
			}
			else if (i == imin)
			{
				int q = -gamma/beta;
				Col(A,AT,j,jmin,q); Col(T,TT,j,jmin,q);
			}
			else 
			{
				int q = -A[jmin][i]/beta;
				Row(A,AT,i,imin,q); Row(S,ST,i,imin,q);
				Row(A,AT,imin,i,1); Row(S,ST,imin,i,1);
			}
		}
        if (alpha == 0) break;
		unsigned int l = pA_beg->first;
		Row(A,AT,l,imin); Row(S,ST,l,imin);
		Col(A,AT,l,jmin); Col(T,TT,l,jmin);
		int diag = A[l][l];
		if (diag < 0) {Row(A,AT,l); Row(S,ST,l);}
		diag = abs(diag);
		if (pA_beg->second.size() > 1)
			for (ppA=++pA_beg->second.begin(); ppA!=pA_beg->second.end(); ppA++)
			{
				int r = - ppA->second/diag;
				Row(A,AT,ppA->first,l,r); 
				Row(S,ST,ppA->first,l,r);
				AT[ppA->first].erase(l);
			}
		A[l].clear(); A[l][l] = diag;
		if (pAT_beg->second.size() > 1)
			for (ppA=++pAT_beg->second.begin(); ppA!=pAT_beg->second.end(); ppA++)
			{
				int r = - ppA->second/diag;
				Col(A,AT,ppA->first,l,r); 
				Col(T,TT,ppA->first,l,r);
				A[ppA->first].erase(l);
			}
		AT[l].clear(); AT[l][l] = diag;
	}
}

}

}

}
