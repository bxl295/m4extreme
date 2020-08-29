// MultiVector.h: implementation of the MultiVector class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <cassert>
#include "./MultiVector.h"

using namespace std;

bool 
operator != (const vector<Set::VectorSpace::Vector> &A, 
			 const vector<Set::VectorSpace::Vector> &B)
{
	vector<Set::VectorSpace::Vector>::const_iterator pA;
	vector<Set::VectorSpace::Vector>::const_iterator pB;
	for (pA=A.begin(), pB=B.begin(); pA!=A.end(); pA++, pB++)
		if (*pA == *pB) return false;
	return true;
}

bool 
operator == (const vector<Set::VectorSpace::Vector> &A, 
			 const vector<Set::VectorSpace::Vector> &B)
{
	vector<Set::VectorSpace::Vector>::const_iterator pA;
	vector<Set::VectorSpace::Vector>::const_iterator pB;
	for (pA=A.begin(), pB=B.begin(); pA!=A.end(); pA++, pB++)
		if (*pA != *pB) return false;
	return true;
}

vector<Set::VectorSpace::Vector>
operator + (const vector<Set::VectorSpace::Vector> &A, 
			const vector<Set::VectorSpace::Vector> &B)
{
	vector<Set::VectorSpace::Vector> C = A;
	vector<Set::VectorSpace::Vector>::const_iterator pB;
	vector<Set::VectorSpace::Vector>::iterator pC;
	for (pB=B.begin(), pC=C.begin(); pB!=B.end(); pB++, pC++) *pC += *pB; 
	return C;
}

vector<Set::VectorSpace::Vector>
operator - (const vector<Set::VectorSpace::Vector> &A, 
			const vector<Set::VectorSpace::Vector> &B)
{
	vector<Set::VectorSpace::Vector> C = A;
	vector<Set::VectorSpace::Vector>::const_iterator pB;
	vector<Set::VectorSpace::Vector>::iterator pC;
	for (pB=B.begin(), pC=C.begin(); pB!=B.end(); pB++, pC++) *pC -= *pB; 
	return C;
}

void
operator *= (vector<Set::VectorSpace::Vector> &A, 
			 const double &p)
{
	vector<Set::VectorSpace::Vector>::iterator pA;
	for (pA=A.begin(); pA!=A.end(); pA++) *pA *= p;
}

void
operator /= (vector<Set::VectorSpace::Vector> &A, 
			 const double &p)
{
	vector<Set::VectorSpace::Vector>::iterator pA;
	for (pA=A.begin(); pA!=A.end(); pA++) *pA /= p;
}

vector<Set::VectorSpace::Vector>
operator * (const vector<Set::VectorSpace::Vector> &A, 
			const double &p)
{
	vector<Set::VectorSpace::Vector> B = A;
	vector<Set::VectorSpace::Vector>::iterator pB;
	for (pB=B.begin(); pB!=B.end(); pB++) *pB *= p;
	return B;
}

vector<Set::VectorSpace::Vector>
operator * (const double &p,
			const vector<Set::VectorSpace::Vector> &A)
{
	vector<Set::VectorSpace::Vector> B = A;
	vector<Set::VectorSpace::Vector>::iterator pB;
	for (pB=B.begin(); pB!=B.end(); pB++) *pB *= p;
	return B;
}

vector<Set::VectorSpace::Vector>
operator / (const vector<Set::VectorSpace::Vector> &A, 
			const double &p)
{
	vector<Set::VectorSpace::Vector> B = A;
	vector<Set::VectorSpace::Vector>::iterator pB;
	for (pB=B.begin(); pB!=B.end(); pB++) *pB /= p;
	return B;
}

double 
operator * (const vector<Set::VectorSpace::Vector> &A, 
			const vector<Set::VectorSpace::Vector> &B)
{
	double p = 0.0;
	vector<Set::VectorSpace::Vector>::const_iterator pA;
	vector<Set::VectorSpace::Vector>::const_iterator pB;
	for (pA=A.begin(), pB=B.begin(); pA!=A.end(); pA++, pB++) p += (*pA)(*pB); 
	return p;
}

void 
Random(vector<Set::VectorSpace::Vector> &A)
{
	vector<Set::VectorSpace::Vector>::iterator pA;
	for (pA=A.begin(); pA!=A.end(); pA++) pA->Randomize();
}
