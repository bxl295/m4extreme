// SimplicialComplex.cpp: Implementation of the SimplicialComplex class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <cassert>
#include "SimplicialComplex.h"

namespace Geometry
{
//////////////////////////////////////////////////////////////////////
// Class SimplicialComplex
//////////////////////////////////////////////////////////////////////

SimplicialComplex::SimplicialComplex(){}
	
SimplicialComplex::SimplicialComplex(const unsigned int &n) 
: CellComplex(n){}

SimplicialComplex::~SimplicialComplex(){}

void 
SimplicialComplex::ReadTable(
	const set<set<Cell *> > &CT)
{
	if (CT.size() == 0) return;

	unsigned int n = (unsigned int)CT.begin()->size();
	set<set<Cell *> >::const_iterator pCT;
	for (pCT=CT.begin(); pCT!=CT.end(); pCT++)
		if (pCT->size() > n) n = pCT->size();
	unsigned int dim = n-1;
	if (dim == 0) return;

	unsigned int p, q; 
	if (this->size() == 0)
	{
		for (p=0; p<=dim; p++)
			this->push_back(set<Cell *>());
	}
	else
	{
		assert(n == this->size());
		set<Cell *>::iterator pE;
		for (pE=(*this)[0].begin(); pE!=(*this)[0].end(); pE++) 
			(*pE)->Coboundary().clear();
		for (p=1; p<=dim; p++)
		{
			set<Cell *>::iterator pE;
			for (pE=(*this)[p].begin(); pE!=(*this)[p].end(); pE++) delete *pE;
			(*this)[p].clear();
		}
	}

	vector<set<set<Cell *> > > K; 
	for (p=0; p<=dim; p++)
		K.push_back(set<set<Cell *> >());
	for (pCT=CT.begin(); pCT!=CT.end(); pCT++)
		K[pCT->size()-1].insert(*pCT);

	vector<map<set<Cell *>, Cell *> > e_map;
	for (p=0; p<=dim; p++)
		e_map.push_back(map<set<Cell *>, Cell *>());

	set<set<Cell *> >::iterator pK;
	set<Cell *>::iterator pk;

	if (dim > 1)
	{
		for (p=dim; p>1; p--)
		{
			q = p - 1;
			for (pK=K[p].begin(); pK!=K[p].end(); pK++)
			{
				set<Cell *> k = *pK;
				for (pk=k.begin(); pk!=k.end(); pk++)
				{
					set<Cell *> m = k;
					m.erase(*pk);
					K[q].insert(m); 
				}
			}
		}
	}

	if ((*this)[0].size() == 0)
	{
		for (pK=K[dim].begin(); pK!=K[dim].end(); pK++) 
		{
			set<Cell *> k = *pK;
			for (pk=k.begin(); pk!=k.end(); pk++) 
				(*this)[0].insert(*pk);
		}
	}

	for (p=1; p<=dim; p++)
	{
		for (pK=K[p].begin(); pK!=K[p].end(); pK++) 
		{
			Cell *e = new Cell(p);
			e_map[p][*pK] = e;
			(*this)[p].insert(e);
		}
	}

	int sign;

	for (pK=K[1].begin(); pK!=K[1].end(); pK++)
	{
		set<Cell *> k = *pK;
		set<Cell *>::iterator pk;
		Cell *e = e_map[1][k];
		Chain &e_bo = e->Boundary();
		for (pk=k.begin(), sign=1; pk!=k.end(); pk++, sign=-sign)
		{
			set<Cell *> m = k;
			m.erase(*pk); 
			Cell *f = *m.begin();
			e_bo[f] = sign;
			Cochain &f_co = f->Coboundary();
			f_co[e] = sign;
		}
	}

	if (dim > 1)
	{
		for (p=2; p<=dim; p++)
		{
			q = p - 1;
			for (pK=K[p].begin(); pK!=K[p].end(); pK++)
			{
				set<Cell *> k = *pK;
				set<Cell *>::iterator pk;
				Cell *e = e_map[p][k];
				Chain &e_bo = e->Boundary();
				for (pk=k.begin(), sign=1; pk!=k.end(); pk++, sign=-sign)
				{
					set<Cell *> m = k;
					m.erase(*pk); 
					Cell *f = e_map[q][m];
					e_bo[f] = sign;
					Cochain &f_co = f->Coboundary();
					f_co[e] = sign;
				}
			}
		}
	}
}

vector<Cell *>
SimplicialComplex::ReadTable(
	const vector<vector<unsigned int> > &CA)
{
	assert (CA.size() > 0);
	unsigned int i, j;

	set<unsigned int> N;
	for (i=0; i<CA.size(); i++)
		for (j=0; j<CA[i].size(); j++)
			N.insert(CA[i][j]);

	vector<Cell *> CV;
	for (i=0; i<N.size(); i++) 
		CV.push_back(new Cell(0));

	set<set<Cell *> > CT; 
	for (i=0; i<CA.size(); i++)
	{
		set<Cell *> k;
		for (j=0; j<CA[i].size(); j++) 
			k.insert(CV[CA[i][j]]);
		CT.insert(k);
	}
	ReadTable(CT);

	return CV;
}

void
SimplicialComplex::Reconnect(
	const vector<vector<unsigned int> > &CA,
	const vector<Cell *> &CV)
{
	assert (CA.size() > 0);
	unsigned int i, j;

	set<set<Cell *> > CT; 
	for (i=0; i<CA.size(); i++)
	{
		set<Cell *> k;
		for (j=0; j<CA[i].size(); j++) 
			k.insert(CV[CA[i][j]]);
		CT.insert(k);
	}
	ReadTable(CT);
}

}
