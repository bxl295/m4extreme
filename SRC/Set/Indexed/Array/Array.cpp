// Array.cpp: implementation the Array class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include "Array.h"

//////////////////////////////////////////////////////////////////////
// Set operators
//////////////////////////////////////////////////////////////////////

bool operator != (const Set::Array &A, const Set::Array &B)
{
	double *p, *q; double tol=1.e-15;
	for (p=A.begin(),q=B.begin(); p<A.end(); p++,q++) 
		if (fabs((*p)-(*q)) > tol) return true;
	return false;
}

bool operator == (const Set::Array &A, const Set::Array &B)
{
	double *p, *q; double tol=1.e-15;
	for (p=A.begin(),q=B.begin(); p<A.end(); p++,q++) 
		if (fabs((*p)-(*q)) > tol) return false;
	return true;
}

bool operator < (const Set::Array &A, const Set::Array &B)
{
	const double *p, *q, *pm, *qm;
	for (p=A.end(), pm=p, pm--, q=B.end(), qm=q, qm--; 
		p!=A.begin(); p--, pm--, q--, qm--)
	{
		if (*pm < *qm) return true;
		if (*pm > *qm) return false;
	}
	return false;
}

void Random(Set::Array &A){A.Randomize();}

//////////////////////////////////////////////////////////////////////
// Printing
//////////////////////////////////////////////////////////////////////

ostream & operator<<(ostream &os, const Set::Array &A)
{
	A.print(&os); return os;
}
