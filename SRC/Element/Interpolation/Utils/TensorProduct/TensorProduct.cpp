// TensorProduct.cpp: implementation of the TensorProduct class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include "./TensorProduct.h"

namespace Element
{
namespace Interpolation
{
namespace Utils
{
//////////////////////////////////////////////////////////////////////
// Class TensorProduct
//////////////////////////////////////////////////////////////////////

TensorProduct::TensorProduct() {}

TensorProduct::~TensorProduct() {}

TensorProduct::TensorProduct(const Polynomial &p) : Polynomial(p) {}

TensorProduct::TensorProduct(const Polynomial &p, const Polynomial &q) : 
	Polynomial(p.dim()+q.dim())
{
	unsigned int i;
	const map<MultiIndex, double> &cp = p.GetCoefficients();
	const map<MultiIndex, double> &cq = q.GetCoefficients();
	map<MultiIndex, double>::const_iterator pcp;
	map<MultiIndex, double>::const_iterator pcq;
	for (pcp=cp.begin(); pcp!=cp.end(); pcp++)
	{
		const MultiIndex &aploc = pcp->first; 
		const double &cploc = pcp->second; 
		for (pcq=cq.begin(); pcq!=cq.end(); pcq++)
		{
			const MultiIndex &aqloc = pcq->first; 
			const double &cqloc = pcq->second; 
			MultiIndex a(n);
			for (i=0; i<p.dim(); i++) a[i] = aploc[i];
			for (i=0; i<q.dim(); i++) a[i+p.dim()] = aqloc[i];
			c.insert(make_pair(a,cploc*cqloc));
		}
	}
}

TensorProduct::TensorProduct(const vector<Polynomial> &p)
{
	assert(p.size() > 0); unsigned int i;
	for (i=0; i<p.size(); i++) n += p[i].dim();
	*this = p[0]; if (p.size() == 1) return;
	for (i=1; i<p.size(); i++) *this = TensorProduct(*this,p[i]);
}

TensorProduct & 
TensorProduct::operator = (const Polynomial &P)
{
	Polynomial::operator = (P);
	return *this;
}

double 
TensorProduct::operator () (const Set::VectorSpace::Vector &x) const
{
	return Polynomial::operator () (x);
}

double 
TensorProduct::operator () (const vector<Set::VectorSpace::Vector> &x) const
{
	Set::VectorSpace::Vector y(n); 
	unsigned int j, k;
	vector<Set::VectorSpace::Vector>::const_iterator px;
	for (px=x.begin(), j=0; px!=x.end(); px++)
	{
		const Set::VectorSpace::Vector &xloc = *px;
		for (k=0; k<xloc.size(); k++, j++) y[j] = xloc[k];
	}
	return this->operator () (y);
}

}

}

}
