// LinearPolynomial2d.h: interface for LinearPolynomial2d class
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(UTILS_REGRESSION_LINEARPOLYNOMIAL2D_INCLUDED_)
#define UTILS_REGRESSION_LINEARPOLYNOMIAL2D_INCLUDED_

#pragma once

#include <map>
#include "../Polynomial1d/Polynomial1d.h"
#include "../../../Set/Indexed/Table/Table.h"

namespace Regression2d
{
class LinearPolynomial
{
public:

    LinearPolynomial();                      
    LinearPolynomial(const LinearPolynomial &);       
    LinearPolynomial(const Set::Array *, const Set::Array *,
		const Set::Table *, unsigned int &);  
    virtual ~LinearPolynomial();
    Regression2d::LinearPolynomial *Clone() const;
    double operator () (const double &, const double &, 
		const unsigned int &, const unsigned int &);    
    
private:

    const Set::Array *x;
    const Set::Array *y;
    const Set::Table *f;

    double xmin, xmax, dx;
    
    map<unsigned int, Regression1d::Polynomial *> RegressionMap;
    unsigned int PolyDegree;
    bool completebase;
          
private:

    unsigned int locate(const double &);
    Regression1d::Polynomial * ProvideRegression(const unsigned int &);
    LinearPolynomial & operator = (const LinearPolynomial &);
};
  
}

#endif // !defined(UTILS_REGRESSION_LINEARPOLYNOMIAL2D_INCLUDED_)
