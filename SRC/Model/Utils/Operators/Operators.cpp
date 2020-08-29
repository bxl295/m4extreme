// Operators.cpp: implementation of the Operators class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include "./Operators.h"

void 
operator += (
	set<Set::Manifold::Point *> &x,
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &u)
{
	set<Set::Manifold::Point *>::iterator px;
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator pu;
	for (px = x.begin(), pu = u.begin(); px!= x.end(); px++, pu++) 
		*(*px) += pu->second;
}

void 
operator -= (
	set<Set::Manifold::Point *> &x,
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &u)
{
	set<Set::Manifold::Point *>::iterator px;
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator pu;
	for (px = x.begin(), pu = u.begin(); px!= x.end(); px++, pu++) 
		*(*px) -= pu->second;
}

void 
operator *= (
	map<Set::Manifold::Point *, Set::VectorSpace::Vector> &u,
	const double &a)
{
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pu;
	for (pu = u.begin(); pu!= u.end(); pu++) pu->second *= a;
}

void 
operator /= (
	map<Set::Manifold::Point *, Set::VectorSpace::Vector> &u,
	const double &a)
{
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pu;
	for (pu = u.begin(); pu!= u.end(); pu++) pu->second /= a;
}

map<Set::Manifold::Point *, Set::VectorSpace::Vector>
operator + (
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &u,
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &v)
{
	map<Set::Manifold::Point *, Set::VectorSpace::Vector> w = u;
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pw;
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator pv;
	for (pw = w.begin(), pv = v.begin(); pw!= w.end(); pw++, pv++) 
		pw->second += pv->second; return w;
}

map<Set::Manifold::Point *, Set::VectorSpace::Vector>
operator - (
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &u,
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &v)
{
	map<Set::Manifold::Point *, Set::VectorSpace::Vector> w = u;
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pw;
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator pv;
	for (pw = w.begin(), pv = v.begin(); pw!= w.end(); pw++, pv++) 
		pw->second -= pv->second; return w;
}

map<Set::Manifold::Point *, Set::VectorSpace::Vector>
operator - (
	const map<Set::Manifold::Point *, Set::Euclidean::Orthonormal::Point> &y,
	const map<Set::Manifold::Point *, Set::Euclidean::Orthonormal::Point> &x)
{
	typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> vector_type;
	map<Set::Manifold::Point *, Set::VectorSpace::Vector> u;
	map<Set::Manifold::Point *, Set::Euclidean::Orthonormal::Point>::const_iterator px;
	map<Set::Manifold::Point *, Set::Euclidean::Orthonormal::Point>::const_iterator py;
	for (px = x.begin(), py = y.begin(); px != x.end(); px++, py++) 
		u.insert(vector_type::value_type(px->first,py->second-px->second));
	return u;
}

map<Set::Manifold::Point *, Set::VectorSpace::Vector>
operator - (const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &u)
{
	map<Set::Manifold::Point *, Set::VectorSpace::Vector> v = u;
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pv;
	for (pv = v.begin(); pv!= v.end(); pv++) 
		pv->second = -(pv->second); return v;
}

map<Set::Manifold::Point *, Set::VectorSpace::Vector>
operator * (
	const double &a,
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &u)
{
	map<Set::Manifold::Point *, Set::VectorSpace::Vector> w = u;
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pw;
	for (pw = w.begin(); pw!= w.end(); pw++) 
		pw->second *= a; return w;
}

map<Set::Manifold::Point *, Set::VectorSpace::Vector>
operator * (
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &u,
	const double &a)
{
	map<Set::Manifold::Point *, Set::VectorSpace::Vector> w = u;
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pw;
	for (pw = w.begin(); pw!= w.end(); pw++) 
		pw->second *= a; return w;
}

map<Set::Manifold::Point *, Set::VectorSpace::Vector>
operator / (
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &u,
	const double &a)
{
	map<Set::Manifold::Point *, Set::VectorSpace::Vector> w = u;
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pw;
	for (pw = w.begin(); pw!= w.end(); pw++) 
		pw->second /= a; return w;
}

double
operator * (
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &u,
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &v)
{
	double a = 0.0;
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator pu;
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator pv;
	for (pu = u.begin(), pv = v.begin(); pu!= u.end(); pu++, pv++) 
		a += pu->second(pv->second); return a;
}
void Null(map<Set::Manifold::Point *, Set::VectorSpace::Vector> &u)
{
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pu;
	for (pu = u.begin(); pu!= u.end(); pu++) Null(pu->second);
}

void Random(map<Set::Manifold::Point *, Set::VectorSpace::Vector> &u)
{
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pu;
	for (pu = u.begin(); pu!= u.end(); pu++) Random(pu->second);
}

double Norm(const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &u)
{
	double a = 0.0;
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator pu;
	for (pu = u.begin(); pu!= u.end(); pu++) 
		a += pu->second(pu->second); return sqrt(a);
}
