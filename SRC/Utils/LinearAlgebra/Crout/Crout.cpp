// Crout.cpp: Implementation of the Crout class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <cassert>
#include "Crout.h"

namespace LinearAlgebra  
{
Crout::Crout() : n(0), index(0){}

Crout::Crout(const unsigned int &n0) 
: n(n0), index(new unsigned int [n])
{
	for (unsigned int i=0; i<n; i++) index[i]=0;
}

Crout::~Crout()
{
	delete [] index;
}

void Crout::Decomposition(Set::VectorSpace::Hom &a)
{
	assert (n == a.size1());
	assert (n == a.size2());

	int i, j, k, imax; 
	int nn = (int) n;
	double aamax, dum, aijabs;
	double diag, save, sum;
	double *vv;

	vv = new double [n];

	sign = 1;

	for (i=0; i<nn; i++)
	{
		aamax = 0.0;
		for (j=0; j<nn; j++)
		{
			aijabs = fabs(a[j][i]);
			if (aijabs > aamax) aamax = aijabs;
		}
		if (aamax == 0.0) throw(0);
		vv[i] = 1.0/aamax;
	}

	for (j=0; j<nn; j++)
	{
		for (i=0; i<j; i++)
		{
			sum = a[j][i];
			for (k=0; k<i; k++) sum -= a[k][i]*a[j][k];
			a[j][i] = sum;
		}

		aamax = 0.0;

		for (i=j; i<nn; i++)
		{
			sum = a[j][i];
			for (k=0; k<j; k++) sum -= a[k][i]*a[j][k];
			a[j][i] = sum;
			dum = vv[i]*fabs(sum);
			if (dum >= aamax)
			{
				imax = i;
				aamax = dum;
			}
		}

		if (j != imax)
		{
			for (k=0; k<nn; k++)
			{
				save = a[k][imax];
				a[k][imax] = a[k][j];
				a[k][j] = save;
			}
			sign = -sign;
			vv[imax] = vv[j];
		}

		if (a[j][j] == 0.0) throw(0);

		index[j] = imax;

		if (j != nn-1)
		{
			diag = 1.0/a[j][j];
			for (i=j+1; i<nn; i++) a[j][i] = a[j][i]*diag;
		}
	}

	delete [] vv;
}

void Crout::Substitution(const Set::VectorSpace::Hom &a, 
						 Set::VectorSpace::Vector &b)
{
	assert (n == a.size1());
	assert (n == a.size2());
	assert (n == b.size());

	int i, ii, j, ll; 
	int nn = (int) n;
	double sum;

	ii = 0;

	for (i=0; i<nn; i++)
	{
		ll = index[i];
		sum = b[ll];
		b[ll] = b[i];

		if (ii >= 0)
		{
			for (j=ii; j<i; j++) sum -= a[j][i]*b[j];
		}
		else if (sum != 0.0) ii = i;

		b[i] = sum;
	}

	for (i=nn-1; i>=0; i--)
	{
		sum = b[i];
		for (j=i+1; j<nn; j++) sum -= a[j][i]*b[j];
		b[i] = sum/a[i][i];
	}
}

void Crout::Substitution(const Set::VectorSpace::Hom &a, 
						 const Set::VectorSpace::Vector &b, 
						 Set::VectorSpace::Vector &x)
{
	assert (n == a.size1());
	assert (n == a.size2());
	assert (n == b.size());
	assert (n == x.size());

	int i; int nn = (int) n;

	for (i=0; i<nn; i++) x[i] = b[i];

	Substitution(a,x);
}

void Crout::Solve(const Set::VectorSpace::Hom &a0, 
				  Set::VectorSpace::Vector &b)
{
	assert (n == a0.size1());
	assert (n == a0.size2());
	assert (n == b.size());

	unsigned int i; 
	Set::VectorSpace::Hom a(a0);

	for (i=0; i<n; i++) index[i] = 0;

	try{Decomposition(a);}
    catch(int flag){throw(flag);}

	Substitution(a,b);
}

void Crout::Solve(const Set::VectorSpace::Hom &a0, 
				  const Set::VectorSpace::Vector &b0, 
				  Set::VectorSpace::Vector &x)
{
	assert (n == a0.size1());
	assert (n == a0.size2());
	assert (n == b0.size());
	assert (n == x.size());

	unsigned int i; 
	Set::VectorSpace::Hom a(a0);
	Set::VectorSpace::Vector b(b0);

	for (i=0; i<n; i++) index[i] = 0;

	try{Decomposition(a);}
    catch(int flag){throw(flag);}

	Substitution(a,b,x);
}

void Crout::Invert(const Set::VectorSpace::Hom &a0, 
				   Set::VectorSpace::Hom &ainv)
{
	assert (n == a0.size1());
	assert (n == a0.size2());
	assert (n == ainv.size1());
	assert (n == ainv.size2());

	unsigned int i, j;
	Set::VectorSpace::Vector b(n);
	Set::VectorSpace::Hom a(a0);

	for (i=0; i<n; i++) index[i] = 0;

	try{Decomposition(a);}
    catch(int flag){throw(flag);}

	for (j=0; j<n; j++)
	{
		for (i=0; i<n; i++) b[i] = 0.0;
		b[j] = 1.0;
		Substitution(a,b,ainv[j]);
	}
}

double Crout::Determinant(const Set::VectorSpace::Hom &a0)
{
	assert (n == a0.size1());
	assert (n == a0.size2());

	unsigned int i; 
	double det;
	Set::VectorSpace::Hom a(a0);

	for (i=0; i<n; i++) index[i] = 0;

	try{Decomposition(a);}
    catch(int flag){if (flag==0) return 0.0;}

	det = 1.0;
	for (i=0; i<n; i++) det *= a[i][i];

	return sign*det;
}

}
