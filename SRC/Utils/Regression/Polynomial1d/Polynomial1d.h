// Polynomial1d.h: interface for the Polynomial1d class
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(UTILS_REGRESSION_POLYNOMIAL1D_INCLUDED_)
#define UTILS_REGRESSION_POLYNOMIAL1D_INCLUDED_

#pragma once

#include "../../../Set/Indexed/Table/Table.h"
#include "../../../Set/Algebraic/VectorSpace/VectorSpace.h"
#include "../../../Set/Algebraic/VectorSpace/Category/Category.h"
#include "../../../Utils/LinearAlgebra/Cholesky/Cholesky.h"

namespace Regression1d
{
class Polynomial
{
public:
    Polynomial();                      
    Polynomial(const Polynomial &);       
    Polynomial(const Set::Array *, const Set::Array *, 
		const unsigned int &);
    virtual ~Polynomial();
    Regression1d::Polynomial *Clone() const;
    const Set::Array GetBeta() const;
    const double GetR2() const;
    double operator () (const double &, const unsigned int &);
      
private:

    const Set::Array *xi;
    const Set::Array *yi;
    unsigned int p;
    Set::Array Beta;
    double R2;
         
private:

    unsigned long int Permutation(const unsigned int &, const unsigned int &);
    Polynomial & operator = (const Polynomial &);
     
};
  
}

#endif // !defined(UTILS_REGRESSION_POLYNOMIAL1D_INCLUDED_)
