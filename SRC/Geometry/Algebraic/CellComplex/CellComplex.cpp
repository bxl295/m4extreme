// CellComplex.cpp: Implementation of the CellComplex class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include "CellComplex.h"

namespace Geometry
{
/////////////////////////////////////////////////////////////////////
// Class SubComplex
//////////////////////////////////////////////////////////////////////

SubComplex::SubComplex(){}

SubComplex::SubComplex(const unsigned int &n) 
: vector<set<Cell *> >(n+1){}

SubComplex::~SubComplex()
{
	vector<set<Cell *> >::iterator pS;
	for (pS=this->begin(); pS!=this->end(); pS++)
		pS->clear();
}

SubComplex::SubComplex(const SubComplex &S)
: vector<set<Cell *> >(S) {}

SubComplex::SubComplex(const vector<set<Cell *> > &S)
: vector<set<Cell *> >(S) {}

SubComplex & 
SubComplex::operator = (const SubComplex &S)
{
	if (this == &S) return *this;
	vector<set<Cell *> >::operator = (S);
	return *this;
}

void 
SubComplex::SetBoundaries()
{
	vector<set<Cell *> >::iterator pS;
	for (pS=this->begin(); pS!=this->end(); pS++)
	{
		set<Cell *>::iterator pE;
		for (pE=pS->begin(); pE!=pS->end(); pE++)
			(*pE)->Co2Bo();
	}
}

void 
SubComplex::SetCoboundaries()
{
	vector<set<Cell *> >::iterator pS;
	for (pS=this->begin(); pS!=this->end(); pS++)
	{
		set<Cell *>::iterator pE;
		for (pE=pS->begin(); pE!=pS->end(); pE++)
			(*pE)->Bo2Co();
	}
}

const vector<set<Cell *> > 
SubComplex::Boundary(Cell * const e) const
{
	int i, j; 
	int p = e->dim();
	int n = (int)this->dim();

	vector<set<Cell *> > I(n+1);
	I[p].insert(e); 
	if (p > 0)
	{
		for (i=p; i>0; i--)
		{
			j = i-1;
			set<Cell *>::iterator pE;
			for (pE=I[i].begin(); pE!=I[i].end(); pE++)
			{
				Chain bo = (*pE)->Boundary();
				Chain::iterator pbo; 
				for (pbo=bo.begin(); pbo!=bo.end(); pbo++)
					I[j].insert(pbo->first);
			}
		}
	}
	return I;
}
//
const vector<set<Cell *> > 
SubComplex::Incidence(Cell * const e) const
{
	int i, j; 
	int p = e->dim();
	int n = (int)this->dim();

	vector<set<Cell *> > I(n+1);
	I[p].insert(e); 
	if (p < n)
	{
		for (i=p; i<n; i++)
		{
			j = i+1;
			set<Cell *>::iterator pE;
			for (pE=I[i].begin(); pE!=I[i].end(); pE++)
			{
				Cochain co = (*pE)->Coboundary();
				Cochain::iterator pco; 
				for (pco=co.begin(); pco!=co.end(); pco++)
					I[j].insert(pco->first);
			}
		}
	}
	if (p > 0)
	{
		for (i=p; i>0; i--)
		{
			j = i-1;
			set<Cell *>::iterator pE;
			for (pE=I[i].begin(); pE!=I[i].end(); pE++)
			{
				Chain bo = (*pE)->Boundary();
				Chain::iterator pbo; 
				for (pbo=bo.begin(); pbo!=bo.end(); pbo++)
					I[j].insert(pbo->first);
			}
		}
	}
	return I;
}

const set<Cell *> 
SubComplex::Incidence(Cell * const e, const int &q) const
{
	int i, j; 
	int p = e->dim();
	int n = (int)this->dim();
	vector<set<Cell *> > I(n+1);

	I[p].insert(e); 

	if (p < q)
	{
		for (i=p; i<q; i++)
		{
			j = i+1;
			set<Cell *>::iterator pE;
			for (pE=I[i].begin(); pE!=I[i].end(); pE++)
			{
				Cochain co = (*pE)->Coboundary();
				Cochain::iterator pco; 
				for (pco=co.begin(); pco!=co.end(); pco++)
					I[j].insert(pco->first);
			}
		}
	}

	if (p > q)
	{
		for (i=p; i>q; i--)
		{
			j = i-1;
			set<Cell *>::iterator pE;
			for (pE=I[i].begin(); pE!=I[i].end(); pE++)
			{
				Chain bo = (*pE)->Boundary();
				Chain::iterator pbo; 
				for (pbo=bo.begin(); pbo!=bo.end(); pbo++)
					I[j].insert(pbo->first);
			}
		}
	}

	return I[q];
}

void 
SubComplex::Bind(
	const SubComplex &C, 
	map<Cell *, Cell *> &f) const
{	
	if (this->dim() == 0) return;

	{
		set<Cell *> e0_set;
		map<Cell *, Cell *>::const_iterator pf;
		for (pf=f.begin(); pf!=f.end(); pf++)
			e0_set.insert(pf->first);
		if (e0_set != C[0]) throw(0);
	}

	const SubComplex &B = *this;

	unsigned int p, q;
	for (p=1; p<C.size(); p++)
	{
		set<Cell *>::const_iterator pe;
		for (pe=C[p].begin(); pe!=C[p].end(); pe++)
		{
			Cell *e = *pe; f[e] = 0;
			set<Cell *> e0_set = C.Incidence(e,0);
			set<Cell *>::iterator pe0_set;
			for (q=p; q<C.size(); q++)
			{
				set<Cell *> e_set;
				for (pe0_set=e0_set.begin(); pe0_set!=e0_set.end(); pe0_set++)
					e_set.insert(f[*pe0_set]);
				set<Cell *>::const_iterator pB;
				for (pB=B[q].begin(); pB!=B[q].end(); pB++)
				{
					Cell *b = *pB;
					vector<set<Cell *> > Bb = B.Boundary(b);
					if (Contained(e_set,Bb)){f[e] = b; goto next;}
				}
			}
			throw(1);
			next: continue;
		}
	}
}

map<Cell *, Cell *> 
SubComplex::Bind(const SubComplex &C) const
{
	map<Cell *, Cell *> f;
	this->Bind(C,f); return f;
}

bool
SubComplex::Contained(
	const set<Cell *> &e_set, 
	const vector<set<Cell *> > &B) const
{
	unsigned int p; set<Cell *> B_set;
	set<Cell *>::const_iterator pE;
	for (p=0; p<B.size(); p++)
		for (pE=B[p].begin(); pE!=B[p].end(); pE++)
			B_set.insert(*pE);
	for (pE=e_set.begin(); pE!=e_set.end(); pE++)
		if (B_set.find(*pE) == B_set.end()) return false;
	return true;
}

const int
SubComplex::dim() const
{
	return (int)this->size()-1;
}
	
Chain 
SubComplex::RandomChain(const unsigned int &p) const
{
	const SubComplex &S = *this;
	Chain c(p);
	set<Cell *>::const_iterator pE;
	for (pE=S[p].begin(); pE!=S[p].end(); pE++)
		c[*pE] = 2*rand() - RAND_MAX;
	return c;
}

Cochain 
SubComplex::RandomCochain(const unsigned int &p) const
{
	const SubComplex &S = *this;
	Cochain c(p);
	set<Cell *>::const_iterator pE;
	for (pE=S[p].begin(); pE!=S[p].end(); pE++)
		c[*pE] = 2*rand() - RAND_MAX;
	return c;
}
	
Chain 
SubComplex::NullChain(const unsigned int &p) const
{
	const SubComplex &S = *this;
	Chain c(p);
	set<Cell *>::const_iterator pE;
	for (pE=S[p].begin(); pE!=S[p].end(); pE++)
		c[*pE] = 0;
	return c;
}

Cochain 
SubComplex::NullCochain(const unsigned int &p) const
{
	const SubComplex &S = *this;
	Cochain c(p);
	set<Cell *>::const_iterator pE;
	for (pE=S[p].begin(); pE!=S[p].end(); pE++)
		c[*pE] = 0;
	return c;
}

//////////////////////////////////////////////////////////////////////
// Class CellComplex
//////////////////////////////////////////////////////////////////////

CellComplex::CellComplex(){}

CellComplex::CellComplex(const unsigned int &n) : SubComplex(n){}

CellComplex::~CellComplex()
{
	vector<set<Cell *> >::iterator pS;
	for (pS=this->begin(); pS!=this->end(); pS++)
	{
		set<Cell *>::iterator pE;
		for (pE=pS->begin(); pE!=pS->end(); pE++) {
		  if ( *pE != NULL ) {
		    delete *pE;
		  }
		}

		pS->clear();
	}
}

CellComplex::CellComplex(const SubComplex &S) : SubComplex(S) {}

}
