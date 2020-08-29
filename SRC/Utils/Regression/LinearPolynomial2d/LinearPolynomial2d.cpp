// LinearPolynomial2d.cpp: implementation for LinearPolynomial2d class
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include "LinearPolynomial2d.h"
#include <cassert>
#include <math.h>

namespace Regression2d
{
//////////////////////////////////////////////////////////////////////
// Class LinearPolynomial
//////////////////////////////////////////////////////////////////////
   
LinearPolynomial::LinearPolynomial(): 
	xmin(0.0), xmax(0.0), dx(0.0), 
	PolyDegree(0), completebase(false) {}
    
LinearPolynomial::LinearPolynomial(const LinearPolynomial &rhs): 
	x(rhs.x), y(rhs.y), f(rhs.f), 
	xmin(rhs.xmin), xmax(rhs.xmax), 
	dx(rhs.dx), PolyDegree(rhs.PolyDegree), 
	RegressionMap(rhs.RegressionMap), 
	completebase(rhs.completebase) {}
    
LinearPolynomial::LinearPolynomial(
	const Set::Array *x_, 
	const Set::Array *y_, 
	const Set::Table *f_, 
	unsigned int &PolyDegree_): 
    x(x_), y(y_), f(f_), 
	PolyDegree(PolyDegree_), 
	completebase(false)
{
    int px = x->size();
    int py = y->size();
    assert(px == f->size2());
    assert(py == f->size1());
    xmin = (*x)[0];     xmax = (*x)[px-1];      
    dx = (xmax - xmin)/(double(px-1));
}
  
LinearPolynomial::~LinearPolynomial() 
{
    map<unsigned int, Regression1d::Polynomial *>::iterator it;
    for(it = RegressionMap.begin(); it!=RegressionMap.end(); it++)
        delete it->second;
}
    
Regression2d::LinearPolynomial * LinearPolynomial::Clone() const
{
    return new LinearPolynomial(*this);
}
    
double
LinearPolynomial::operator () (const double &xq, const double &yq, 
        const unsigned int &pdx, const unsigned int &pdy)
{            
    unsigned int jm1 = locate(xq);
    Regression1d::Polynomial &fjm1 = *ProvideRegression(jm1);
    Regression1d::Polynomial &fj = *ProvideRegression(jm1+1);

    if(pdx == 0)
    {
        return ( fjm1(yq,pdy)*( (*x)[jm1+1] - xq)/dx +
                        fj(yq,pdy)*( xq - (*x)[jm1] )/dx );
    }
    
    if(pdx == 1) return ( fj(yq,pdy) - fjm1(yq,pdy) )/dx;
 
    return 0.0;
}

unsigned int
LinearPolynomial::locate(const double &xq)
{
    assert(xmin<=xq && xq<=xmax);
    if(xq==xmax) return (x->size()-2);
    return (int)floor((xq - xmin)/dx);
}
    
Regression1d::Polynomial *
LinearPolynomial::ProvideRegression(const unsigned int &col)
{
    assert( col <= x->size());
    if(!completebase)
    {
        if(RegressionMap.find(col) == RegressionMap.end())
        {
            RegressionMap[col] = new Regression1d::Polynomial(y, &((*f)[col]), PolyDegree);     
            if(RegressionMap.size() == x->size()) completebase = true;           
            return RegressionMap[col];
        }
        else
            return RegressionMap[col];        
    }
    else
        return RegressionMap[col];  
}    

LinearPolynomial &
LinearPolynomial::operator = (const LinearPolynomial &A)
{
    if(this == &A) return *this;
    x = A.x; y = A.y; f = A.f;
    xmin = A.xmin; xmax = A.xmax; dx = A.dx; 
	PolyDegree = A.PolyDegree;
    RegressionMap = A.RegressionMap; 
	completebase = A.completebase;  
	return *this;
}
    
}

