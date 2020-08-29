// Polynomial.cpp: implementation of the Polynomial class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include "./Polynomial.h"

namespace Element
{
namespace Interpolation
{
namespace Utils
{
//////////////////////////////////////////////////////////////////////
// Class Polynomial
//////////////////////////////////////////////////////////////////////

Polynomial::Polynomial() : n(0) {}

Polynomial::Polynomial(const unsigned int &n_) : n(n_) {}

Polynomial::~Polynomial() {}

Polynomial::Polynomial(const Polynomial &P) : n(P.n), c(P.c) {}

Polynomial & 
Polynomial::operator = (const Polynomial &P)
{
	if (&P == this) return *this;
	n = P.n; c = P.c; return *this;
}

const unsigned int &
Polynomial::dim() const
{
	return n;
}

map<MultiIndex, double> & 
Polynomial::GetCoefficients()
{
	return c;
}

const map<MultiIndex, double> & 
Polynomial::GetCoefficients() const
{
	return c;
}

void 
Polynomial::operator *= (const double &a)
{
	map<MultiIndex, double>::iterator pc;
	for (pc=c.begin(); pc!=c.end(); pc++) pc->second *= a;
}

void 
Polynomial::operator /= (const double &a)
{
	map<MultiIndex, double>::iterator pc;
	for (pc=c.begin(); pc!=c.end(); pc++) pc->second /= a;
}

void 
Polynomial::operator += (const Polynomial &P)
{
	map<MultiIndex, double>::const_iterator pc;
	for (pc=P.GetCoefficients().begin(); pc!=P.GetCoefficients().end(); pc++) 
		c[pc->first] += pc->second;
}

void 
Polynomial::operator -= (const Polynomial &P)
{
	map<MultiIndex, double>::const_iterator pc;
	for (pc=P.GetCoefficients().begin(); pc!=P.GetCoefficients().end(); pc++) 
		c[pc->first] -= pc->second;
}

double 
Polynomial::operator () (const Set::VectorSpace::Vector &x) const
{
	assert(x.size() == n); 
	unsigned int i; double p = 0.0;
	map<MultiIndex, double>::const_iterator pc;
	for (pc=c.begin(); pc!=c.end(); pc++)
	{
		MultiIndex a = pc->first;
		double monomial = 1.0;
		for (i=0; i<a.size(); i++)
			monomial *= pow(x[i],(double)a[i]);
		p += pc->second*monomial;
	}
	return p;
}

Polynomial 
Polynomial::diff(const unsigned int &i) const
{
	Polynomial q(n);
	map<MultiIndex, double> &d = q.GetCoefficients();
	map<MultiIndex, double>::const_iterator pc;
	for (pc=c.begin(); pc!=c.end(); pc++)
	{
		const MultiIndex &a = pc->first;
		if (a[i] != 0)
		{
			MultiIndex b = a; b[i] -= 1;
			d[b] = ((double)a[i])*pc->second;
		}
	}

	return q;
}

//////////////////////////////////////////////////////////////////////
// Utilities
//////////////////////////////////////////////////////////////////////

set<MultiIndex>
IndexSet(const unsigned int &n, const unsigned int &k)
{
	unsigned int i, j, npp=n+1, nmm=n-1;

	set<MultiIndex> I;

	if (n == 0)
	{
		MultiIndex a(npp); 
		a[0] = k; I.insert(a);
	}
	else
	{
		for (i=0; i<=k; i++)
		{
			set<MultiIndex> J = IndexSet(nmm,i);
			set<MultiIndex>::iterator pJ;
			for (pJ=J.begin(); pJ!=J.end(); pJ++)
			{
				MultiIndex b = *pJ;
				MultiIndex a(npp);
				for(j=0; j<n; j++) a[j] = b[j]; 
				a[n] = k - i; I.insert(a);
			}
		}
	}

	return I;
}

set<MultiIndex>
IndexSetInterior(const unsigned int &n, const unsigned int &k)
{
	unsigned int i, j, npp=n+1, nmm=n-1;

	set<MultiIndex> I;

	if (n == 0)
	{
		MultiIndex a(npp); 
		a[0] = k; I.insert(a);
	}
	else
	{
		for (i=1; i<k; i++)
		{
			set<MultiIndex> J = IndexSetInterior(nmm,i);
			set<MultiIndex>::iterator pJ;
			for (pJ=J.begin(); pJ!=J.end(); pJ++)
			{
				MultiIndex b = *pJ;
				MultiIndex a(npp);
				for(j=0; j<n; j++) a[j] = b[j]; 
				a[n] = k - i; I.insert(a);
			}
		}
	}

	return I;
}

unsigned int
IndexSetInteriorSize(const unsigned int &n, const unsigned int &k)
{
	if (n == 0) return 1;
	unsigned int i, nmm=n-1, m=0;
	for (i=1; i<k; i++) m += IndexSetInteriorSize(nmm,i);
	return m;
}

}

}

}

//////////////////////////////////////////////////////////////////////
// Class Polynomial: operators
//////////////////////////////////////////////////////////////////////

Element::Interpolation::Utils::Polynomial operator + (
	const Element::Interpolation::Utils::Polynomial &p, 
	const Element::Interpolation::Utils::Polynomial &q)
{
	Element::Interpolation::Utils::Polynomial r = p; r += q; return r;
}

Element::Interpolation::Utils::Polynomial operator - (
	const Element::Interpolation::Utils::Polynomial &p, 
	const Element::Interpolation::Utils::Polynomial &q)
{
	Element::Interpolation::Utils::Polynomial r = p; r -= q; return r;
}

Element::Interpolation::Utils::Polynomial operator * (
	const double &a, 
	const Element::Interpolation::Utils::Polynomial &p)
{
	Element::Interpolation::Utils::Polynomial q = p; q *= a; return q;
}

Element::Interpolation::Utils::Polynomial operator * (
	const Element::Interpolation::Utils::Polynomial &p,
	const double &a)
{
	Element::Interpolation::Utils::Polynomial q = p; q *= a; return q;
}

Element::Interpolation::Utils::Polynomial operator / (
	const Element::Interpolation::Utils::Polynomial &p,
	const double &a)
{
	Element::Interpolation::Utils::Polynomial q = p; q /= a; return q;
}
