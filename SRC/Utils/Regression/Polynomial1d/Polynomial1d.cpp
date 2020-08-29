// Polynomial1d.cpp: implementation for the Polynomial1d class
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <cassert>
#include <math.h>

#include "./Polynomial1d.h"

namespace Regression1d
{
//////////////////////////////////////////////////////////////////////
// Class Polynomial
//////////////////////////////////////////////////////////////////////

Polynomial::Polynomial(): Beta(0), p(0), R2(0.0) {}

Polynomial::Polynomial(const Polynomial &rhs): 
	xi(rhs.xi), yi(rhs.yi), 
	p(rhs.p), Beta(rhs.Beta), R2(rhs.R2) {}

Polynomial::Polynomial(
	const Set::Array *xi_, 
	const Set::Array *yi_, 
	const unsigned int &p_): 
	xi(xi_), yi(yi_), p(p_), Beta(p_+1)
{
    unsigned int m = yi->size();
    unsigned int n = p+1;
    double xpow, factor;
    Set::VectorSpace::Hom X(m,n);
    Set::VectorSpace::Vector Y(m, yi->begin()), Res(m);
    Set::VectorSpace::Vector B(n, Beta.begin());

    for(unsigned int i=0; i<m; i++)
    {
        X[0][i] = 1.0;
        xpow = (*xi)[i];    factor = xpow;
        
        for(unsigned int j=1; j<n; j++)
        {
            X[j][i] = xpow;
            xpow *= factor;                
        }
    }

    Set::VectorSpace::Hom XTX = Adjoint(X)*X;
    Set::VectorSpace::Vector XTY = Adjoint(X)*Y;
    
    LinearAlgebra::Cholesky C(n);
    C.Solve(XTX, XTY, B);  
    
    double ym(0.0), RSS(0.0), SYY(0.0);
    for(unsigned int i=0; i<m; i++)  ym += Y[i];
    ym /= (double)m;
    for(unsigned int i=0; i<m; i++)  SYY += (Y[i]-ym)*(Y[i]-ym);
    RSS = (X(B)-Y)(X(B)-Y);
    R2 = 1.0 - RSS/SYY;    
}

Polynomial::~Polynomial() {}

Regression1d::Polynomial * 
Polynomial::Clone() const
{
    return new Polynomial(*this);
}

const Set::Array
Polynomial::GetBeta() const
{
    return Beta;
}

const double
Polynomial::GetR2() const
{
    return R2;
}


double
Polynomial::operator () (const double &x, const unsigned int &d)
{
    if(d > p) return 0.0;
    
    if(d == 0)
    {
        double xpow = x;
        double y = Beta[0];
        for(unsigned int i = 1; i<p+1; i++)
        {
            y += Beta[i]*xpow;
            xpow *= x;
        }
        return y;
    }
    else
    {
        double xpow = x;
        double y = Beta[d]*((double)Permutation(d,0));
        for(unsigned int i = d+1; i<p+1; i++)
        {
            y += Beta[i]*(double)Permutation(i,i-d)*xpow;
            xpow *= x;
        }            
    return y;    
    }
}
   
Polynomial &
Polynomial::operator = (const Polynomial &A)
{
    if(this == &A) return *this;
    xi = A.xi; yi = A.yi; 
    p = A.p; Beta = A.Beta;
	return *this;
}

unsigned long int
Polynomial::Permutation(const unsigned int &n, const unsigned int &r)
{
    unsigned long int prod = 1;        
    for(unsigned int i = n; i>r; i--) prod *= i;        
    return prod;
    
}
    
}
