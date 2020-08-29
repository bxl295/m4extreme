// MultiHom.cpp: implementation of the MultiHom class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include "./MultiHom.h" 

bool 
operator != (const vector<vector<Set::VectorSpace::Hom> > &A, 
			 const vector<vector<Set::VectorSpace::Hom> > &B)
{
	vector<vector<Set::VectorSpace::Hom> >::const_iterator pA;
	vector<vector<Set::VectorSpace::Hom> >::const_iterator pB;
	for (pA=A.begin(), pB=B.begin(); pA!=A.end(); pA++, pB++)
	{
		vector<Set::VectorSpace::Hom>::const_iterator pAcol;
		vector<Set::VectorSpace::Hom>::const_iterator pBcol;
		for (pAcol=pA->begin(), pBcol=pB->begin(); pAcol!=pA->end(); pAcol++, pBcol++)
			if (*pAcol == *pBcol) return false;
	}
	return true;
}

bool 
operator == (const vector<vector<Set::VectorSpace::Hom> > &A, 
			 const vector<vector<Set::VectorSpace::Hom> > &B)
{
	vector<vector<Set::VectorSpace::Hom> >::const_iterator pA;
	vector<vector<Set::VectorSpace::Hom> >::const_iterator pB;
	for (pA=A.begin(), pB=B.begin(); pA!=A.end(); pA++, pB++)
	{
		vector<Set::VectorSpace::Hom>::const_iterator pAcol;
		vector<Set::VectorSpace::Hom>::const_iterator pBcol;
		for (pAcol=pA->begin(), pBcol=pB->begin(); pAcol!=pA->end(); pAcol++, pBcol++)
			if (*pAcol != *pBcol) return false;
	}
	return true;
}

vector<Set::VectorSpace::Vector>
operator * (const vector<vector<Set::VectorSpace::Hom> > &A,
		    const vector<Set::VectorSpace::Vector> &B)
{
	vector<Set::VectorSpace::Vector> C;
	vector<Set::VectorSpace::Hom>::const_iterator pAcol;
	for (pAcol=A.begin()->begin(); pAcol!=A.begin()->end(); pAcol++)
		C.push_back(pAcol->size1());
	vector<Set::VectorSpace::Vector>::const_iterator pB;
	vector<vector<Set::VectorSpace::Hom> >::const_iterator pA;
	for (pA=A.begin(), pB=B.begin(); pA!=A.end(); pA++, pB++)
	{
		vector<Set::VectorSpace::Vector>::iterator pC;
		vector<Set::VectorSpace::Hom>::const_iterator pAcol;
		for (pAcol=pA->begin(), pC=C.begin(); pAcol!=pA->end(); pAcol++, pC++)
			*pC += (*pAcol)(*pB);
	}
	return C;
}
