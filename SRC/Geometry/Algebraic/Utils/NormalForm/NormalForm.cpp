// NormalForm.cpp: Implementation of the NormalForm class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <cassert>
#include "NormalForm.h"

#define MIN(a,b) ((a < b) ? a : b)

namespace Geometry
{
namespace Algebraic
{
namespace Utils
{
//////////////////////////////////////////////////////////////////////
// Class SubTable
//////////////////////////////////////////////////////////////////////

SubTable::SubTable() {}

SubTable::SubTable(
	const Geometry::Algebraic::Table &A, 
	const unsigned int &n) : 
	n1(A.size1()-n), n2(A.size2()-n),
	L(new Geometry::Algebraic::Array * [n2])
{
	unsigned int m1 = A.size1();
	unsigned int m2 = A.size2();
	assert (n <= MIN(m1,m2));
	unsigned int i; int *p;
	for (i=0, p=A.begin()+n*m1+n; i<n2; i++, p+=m1) 
		L[i] = new Geometry::Algebraic::Array(n1,p);
}

SubTable::~SubTable()
{
	for (unsigned int i=0; i<n2; i++) delete L[i];
	delete [] L; 
}

const Geometry::Algebraic::Array & 
SubTable::operator [] (const unsigned int &i) const
{
	assert(i < n2); return *L[i];
}

Geometry::Algebraic::Array & 
SubTable::operator [] (const unsigned int &i)
{
	assert(i < n2); return *L[i];
}

void 
SubTable::print(ostream *os)
{
	for (unsigned int i=0; i<n2; i++) *os << *L[i];
}

unsigned int 
SubTable::size1() const
{
	return n1;
}

unsigned int 
SubTable::size2() const
{
	return n2;
}

/////////////////////////////////////////////////////////////////////
// Class NormalForm
//////////////////////////////////////////////////////////////////////

//	Reduction algorithm
//	Munkres, J.R., "Elements of Algebric Topology", 
//	Perseus Publishing, 1984, p.56.

NormalForm::NormalForm(){}

NormalForm::~NormalForm() {}

void
NormalForm::operator () (Geometry::Algebraic::Table &A0)
{
	unsigned int i, j, l;
	unsigned int nrow = A0.size1();
	unsigned int ncol = A0.size2();
	unsigned int min = MIN(nrow,ncol);

	for (l=0; l<min; l++)
	{
		SubTable A(A0,l);
		unsigned int imin, jmin, alpha;
		while (true)
		{
			unsigned int i, j;
			Min(A,imin,jmin,alpha);
			if (alpha == 0) return;
			for (j=0; j<A.size2(); j++)
				for (i=0; i<A.size1(); i++)
					if (abs(A[j][i])%alpha != 0) goto divide;
			break; divide: 
			if (j == jmin)
			{
				int q = -A[j][i]/A[jmin][imin];
				Row(A,i,imin,q);
			}
			else if (i == imin)
			{
				int q = -A[j][i]/A[jmin][imin];
				Col(A,j,jmin,q);
			}
			else 
			{
				int q = -A[jmin][i]/A[jmin][imin];
				Row(A,i,imin,q);
				Row(A,imin,i,1);
			}
		}
        if (alpha == 0) break;
		Row(A,0,imin);
		Col(A,0,jmin);
		if (A[0][0] < 0) Row(A,0);

		for (i=1; i<A.size1(); i++)
		{
			int r = - A[0][i]/A[0][0];
			Row(A,i,0,r);
		}
		for (j=1; j<A.size2(); j++)
		{
			int r = - A[j][0]/A[0][0];
			Col(A,j,0,r);
		}
	}
}

void
NormalForm::operator () (
	Geometry::Algebraic::Table &A0,
	Geometry::Algebraic::Table &S0,
	Geometry::Algebraic::Table &T0)
{
	unsigned int i, j, l;
	unsigned int nrow = A0.size1();
	unsigned int ncol = A0.size2();
	unsigned int min = MIN(nrow,ncol);
	SubTable S(S0,0); SubTable T(T0,0);

	for (l=0; l<min; l++)
	{
		SubTable A(A0,l);
		unsigned int imin, jmin, alpha;
		while (true)
		{
			unsigned int i, j;
			Min(A,imin,jmin,alpha);
			if (alpha == 0) return;
			for (j=0; j<A.size2(); j++)
				for (i=0; i<A.size1(); i++)
					if (abs(A[j][i])%alpha != 0) goto divide;
			break; divide: 
			if (j == jmin)
			{
				int q = -A[j][i]/A[jmin][imin];
				Row(A,i,imin,q); Row(S,l+i,l+imin,q);
			}
			else if (i == imin)
			{
				int q = -A[j][i]/A[jmin][imin];
				Col(A,j,jmin,q); Col(T,l+j,l+jmin,q);
			}
			else 
			{
				int q = -A[jmin][i]/A[jmin][imin];
				Row(A,i,imin,q); Row(S,l+i,l+imin,q);
				Row(A,imin,i,1); Row(S,l+imin,l+i,1);
			}
		}
        if (alpha == 0) break;
		Row(A,0,imin); Row(S,l,l+imin);
		Col(A,0,jmin); Col(T,l,l+jmin);
		if (A[0][0] < 0) {Row(A,0); Row(S,l);}

		for (i=1; i<A.size1(); i++)
		{
			int r = - A[0][i]/A[0][0];
			Row(A,i,0,r); Row(S,l+i,l,r);
		}
		for (j=1; j<A.size2(); j++)
		{
			int r = - A[j][0]/A[0][0];
			Col(A,j,0,r); Col(T,l+j,l,r);
		}
	}
}

void 
NormalForm::Row(
	SubTable &A, 
	unsigned int row1,
	unsigned int row2)
{
	int save; 
	for (unsigned int j=0; j<A.size2(); j++)
	{
		save = A[j][row1]; 
		A[j][row1] = A[j][row2]; 
		A[j][row2] = save;
	}
}

void 
NormalForm::Row(
	SubTable &A, 
	unsigned int row)
{
	for (unsigned int j=0; j<A.size2(); j++)
		A[j][row] = - A[j][row]; 
}

void 
NormalForm::Row(
	SubTable &A, 
	unsigned int row1,
	unsigned int row2,
	unsigned int q)
{
	for (unsigned int j=0; j<A.size2(); j++)
		A[j][row1] = A[j][row1] + q*A[j][row2]; 
}

void 
NormalForm::Col(
	SubTable &A, 
	unsigned int col1,
	unsigned int col2)
{
	int save;
	for (unsigned int i=0; i<A.size1(); i++)
	{
		save = A[col1][i]; 
		A[col1][i] = A[col2][i]; 
		A[col2][i] = save;
	}
}

void 
NormalForm::Col(
	SubTable &A, 
	unsigned int col)
{
	for (unsigned int i=0; i<A.size1(); i++)
		A[col][i] = - A[col][i]; 
}

void 
NormalForm::Col(
	SubTable &A, 
	unsigned int col1,
	unsigned int col2,
	unsigned int q)
{
	for (unsigned int i=0; i<A.size1(); i++)
		A[col1][i] = A[col1][i] + q*A[col2][i]; 
}

void 
NormalForm::Min(
	SubTable &A, 
	unsigned int &imin,
	unsigned int &jmin,
	unsigned int &alpha)
{
	alpha = 0;
	for (unsigned int j=0; j<A.size2(); j++)
	{
		for (unsigned int i=0; i<A.size1(); i++)
		{
			unsigned int a = abs(A[j][i]);
			if (alpha == 0){if (a > 0) {alpha = a; imin = i; jmin = j;}}
			else {if ((a != 0) && (a < alpha)) {alpha = a; imin = i; jmin = j;}}
		}
	}
}

}

}

}

//////////////////////////////////////////////////////////////////////
// Printing
//////////////////////////////////////////////////////////////////////

ostream & 
operator<<(ostream &os, Geometry::Algebraic::Utils::SubTable &A)
{
	A.print(&os); return os;
}

void 
Null(Geometry::Algebraic::Table &A)
{
	for (unsigned int j=0; j<A.size2(); j++)
		for (unsigned int i=0; i<A.size1(); i++) A[j][i] = 0;
}
void 
Identity(Geometry::Algebraic::Table &A)
{
	assert(A.size1() == A.size2()); Null(A);
	for (unsigned int i=0; i<A.size1(); i++) A[i][i] = 1;
}
