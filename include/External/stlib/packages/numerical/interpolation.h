// -*- C++ -*-

/*!
  \file interpolation.h
  \brief Includes the interpolation functions.
*/

/*!
  \page interpolation Interpolation Package

  - There are several functors for performing 
  \ref interpolation_PolynomialInterpolationUsingCoefficients "polynomial interpolation"
  on a 1-D regular grid. Because these store the polynomial coefficients
  they are fast, but use more memory than methods that work directly on the
  grid of function values.
  - The numerical::QuinticInterpolation2D class performs quintic
  interpolation on 2-D grids (period or plain). It stores the gradient
  and Hessian at the grid points in order to accelerate the
  interpolation.
  - There are several functors for performing linear and cubic 
  \ref interpolation_InterpolatingFunction1DRegularGrid "interpolation on a 1-D regular grid".
  - There are also functors for performing
  \ref interpolation_InterpolatingFunctionRegularGrid "interpolation on an multi-dimensional regular grid".
  Constant value,
  multi-linear, and multi-cubic interpolation is supported for plain
  and periodic data.
  - \ref interpolation_simplex
  - \ref interpolation_hermite
  - The numerical::LinInterpGrid class is a functor for performing
    linear interpolation on a regular grid.
  - The numerical::PolynomialInterpolationNonNegative class may be useful
  for interpolating non-negative, nearly singular data.

  Use the interpolation package by including the file interpolation.h.
*/

#if !defined(__numerical_interpolation_h__)
#define __numerical_interpolation_h__

#include "interpolation/hermite.h"
#include "interpolation/InterpolatingFunction1DRegularGrid.h"
#include "interpolation/InterpolatingFunction1DRegularGridReference.h"
#include "interpolation/InterpolatingFunctionMultiArrayOf1DRegularGrids.h"
#include "interpolation/InterpolatingFunctionRegularGrid.h"
#include "interpolation/InterpolatingFunctionRegularGridReference.h"
#include "interpolation/InterpolatingFunctionMultiArrayOfRegularGrids.h"
#include "interpolation/LinInterpGrid.h"
#include "interpolation/PolynomialInterpolationUsingCoefficients.h"
#include "interpolation/PolynomialInterpolationUsingCoefficientsReference.h"
#include "interpolation/PolynomialInterpolationUsingCoefficientsVector.h"
#include "interpolation/PolynomialInterpolationUsingCoefficientsMultiArray.h"
#include "interpolation/simplex.h"

#endif
