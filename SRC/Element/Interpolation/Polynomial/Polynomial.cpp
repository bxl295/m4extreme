// Polynomial.cpp: implementation of the Polynomial class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <set>
#include <cassert>
#include "./Polynomial.h"

namespace Element
{
namespace Interpolation
{
namespace Polynomial
{
//////////////////////////////////////////////////////////////////////
// Class Shape<0>
//////////////////////////////////////////////////////////////////////

Shape<0>::Shape() {}

Shape<0>::Shape(const set<point_type *> &x) : n((unsigned int)(x.size()-1))
{
	unsigned int i, npp = n + 1;
	set<point_type *>::const_iterator px;
	for (i=0, px=x.begin(); i<npp; i++, px++)
	{
		multiindex_type a(npp); a[i] = 1;
		polynomial_type p(npp); 
		p.GetCoefficients().insert(make_pair(a,1.0));
		w.insert(make_pair(*px,p));
	}
}

Shape<0>::Shape(const map<point_type *, multiindex_type> &I) : 
	n(I.begin()->second.size()-1)
{
//	Initialization

	unsigned int i, npp = n + 1;
	map<point_type *, multiindex_type>::const_iterator pI;
	for (pI=I.begin(); pI!=I.end(); pI++)
	{
		assert(pI->second.size() == npp);
		assert(pI->second.degree() > 0);
	}

//	Shape functions

	map<point_type *, vector_type> d;

	for (pI=I.begin(); pI!=I.end(); pI++)
	{
		multiindex_type a = pI->second;
		polynomial_type m(npp); 
		m.GetCoefficients().insert(make_pair(a,1.0));
		vector_type l(npp); 
		unsigned int k = a.degree();
		for (i=0; i<npp; i++) 
			l[i] = (double)a[i]/(double)k;
		polynomial_type p = m;
		map<point_type *, polynomial_type>::iterator pw;
		map<point_type *, vector_type>::iterator pd;
		for (pw=w.begin(), pd=d.begin(); pw!=w.end(); pw++, pd++)
		{
			polynomial_type &q = pw->second;
			const vector_type &y = pd->second;
			p -= m(y)*q;
		}
		p /= p(l); 
		for (pw=w.begin(), pd=d.begin(); pw!=w.end(); pw++, pd++)
		{
			polynomial_type &q = pw->second;
			const vector_type &y = pd->second;
			q -= q(l)*p; q /= q(y);
		}
		w.insert(make_pair(pI->first,p));
		d.insert(make_pair(pI->first,l)); 
	}
}

Shape<0>::~Shape() {}

map<point_type *, double> 
Shape<0>::operator () (const vector_type &x) const
{
	assert(n == x.size());
	unsigned int i, npp = n + 1;
	vector_type l(npp); l[n] = 1.0;
	for (i=0; i<n; i++){l[i] = x[i]; l[n] -= x[i];}
	map<point_type *, double> N;
	map<point_type *, polynomial_type>::const_iterator pw;
	for (pw=w.begin(); pw!=w.end(); pw++)
		N[pw->first] = pw->second(l);

	return N;
}

//////////////////////////////////////////////////////////////////////
// Class Shape<1>
//////////////////////////////////////////////////////////////////////

Shape<1>::Shape() {}

Shape<1>::Shape(const Shape<0> &S) : n(S.n)
{
	unsigned int i, npp = n + 1;

	map<point_type *, polynomial_type>::const_iterator pw;
	for (pw=S.w.begin(); pw!=S.w.end(); pw++)
	{
		vector<polynomial_type> Dwa;
		for (i=0; i<npp; i++) 
			Dwa.push_back(pw->second.diff(i));
		Dw.insert(make_pair(pw->first,Dwa));
	}
}

Shape<1>::~Shape() {}

map<point_type *, vector_type>
Shape<1>::operator () (const vector_type &x) const
{
	assert(n == x.size());
	unsigned int i, npp = n + 1;
	vector_type l(npp); l[n] = 1.0;
	for (i=0; i<n; i++) {l[i] = x[i]; l[n] -= x[i];}
	map<point_type *, vector_type> DN;
	map<point_type *, vector<polynomial_type> >::const_iterator pDw;
	for (pDw=Dw.begin(); pDw!=Dw.end(); pDw++)
	{
		const vector<polynomial_type> &Dwa = pDw->second;
		vector_type A(npp); vector_type B(n);
		for (i=0; i<npp; i++) A[i] = Dwa[i](l); 
		for (i=0; i<n; i++) B[i] = A[i] - A[n];
		DN.insert(make_pair(pDw->first,B));
	}

	return DN;
}

//////////////////////////////////////////////////////////////////////
// Utilities
//////////////////////////////////////////////////////////////////////

vector_type 
RandomPoint(const unsigned int &n)
{
	assert(n >= 1); 
	unsigned int i, npp = n+1;
	vector_type l(npp);
	double sum = 0.0;
	for (i=0; i<npp; i++)
	{
		l[i] = (double)rand()/(double)RAND_MAX;
		sum += l[i];
	}
	for (i=0; i<npp; i++)
		l[i] /= sum;
	vector_type t(n);
	for (i=0; i<n; i++) t[i] = l[i];
	return t;
}

map<point_type *, vector_type>
NodalCoordinates(const set<point_type *> &P)
{
	unsigned int npp = (unsigned int)P.size();
	unsigned int i, n = npp - 1;
	map<point_type *, vector_type> x;
	set<point_type *>::const_iterator pP;
	for (i=0, pP=P.begin(); i<n; i++, pP++)
	{
		vector_type xloc(n); xloc[i] = 1.0;
		x.insert(make_pair(*pP,xloc));
	}
	vector_type xloc(n);
	x.insert(make_pair(*pP,xloc));
	return x;
}

map<point_type *, vector_type>
NodalCoordinates(const map<point_type *, multiindex_type> &I)
{
	unsigned int i, k, n;
	map<point_type *, vector_type> x;
	map<point_type *, multiindex_type>::const_iterator pI;
	for (pI=I.begin(); pI!=I.end(); pI++)
	{
		multiindex_type a = pI->second;
		n = a.size()-1; k = a.degree();
		vector_type xa(n);
		for (i=0; i<n; i++) 
			xa[i] = (double)a[i]/(double)k;
		x.insert(make_pair(pI->first,xa));
	}
	return x;
}

set<multiindex_type>
IndexSetReduced(const unsigned int &n, const unsigned int &k)
{
	set<multiindex_type> I0 = 
		Element::Interpolation::Utils::IndexSet(n,k);
	set<multiindex_type> I;
	set<multiindex_type>::const_iterator pI;
	for (pI=I0.begin(); pI!=I0.end(); pI++)
		I.insert(Reduce(*pI)); return I;
}

set<multiindex_type>
IndexSetInteriorReduced(const unsigned int &n, const unsigned int &k)
{
	set<multiindex_type> I1 = 
		Element::Interpolation::Utils::IndexSetInterior(n,k);
	set<multiindex_type> I;
	set<multiindex_type>::const_iterator pI;
	for (pI=I1.begin(); pI!=I1.end(); pI++)
		I.insert(Reduce(*pI)); return I;
}

}

}

}
