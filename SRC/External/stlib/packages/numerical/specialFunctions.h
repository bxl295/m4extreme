// -*- C++ -*-

/*!
  \file numerical/specialFunctions.h
  \brief Includes the special function classes.
*/

#if !defined(__numerical_specialFunctions_h__)
#define __numerical_specialFunctions_h__

#include "specialFunctions/BinomialCoefficient.h"
#include "specialFunctions/ErrorFunction.h"
#include "specialFunctions/ExponentialForSmallArgument.h"
#include "specialFunctions/Factorial.h"
#include "specialFunctions/Gamma.h"
#include "specialFunctions/HarmonicNumber.h"
#include "specialFunctions/LogarithmOfFactorial.h"
#include "specialFunctions/LogarithmOfFactorialCached.h"

#endif

namespace numerical {

/*!
  \page numerical_specialFunctions Special Functions.

  - There are functions and functors for computing
    \ref numerical_specialFunctions_BinomialCoefficient "binomial coefficients".
  - There are functions for computing the
    \ref numerical_specialFunctions_ErrorFunction "error functions" \c erf()
    and \c erfc().
  - The ExponentialForSmallArgument functor uses truncated Taylor series to
    efficiently computes the exponential function for very small arguments.
  - There are functions and functors for computing the
    \ref numerical_specialFunctions_Factorial "factorial function"
    and the
    \ref numerical_specialFunctions_LogarithmOfFactorial "logarithm of the factorial function".
    The LogarithmOfFactorialCached functor stores function values in a table.
  - One can compute the Gamma function with LogarithmOfGamma.
  - There are functions and functors for computing the
  \ref numerical_specialFunctions_HarmonicNumber "harmonic number function".
*/

} // namespace numerical
