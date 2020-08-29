// -*- C++ -*-

/*!
  \file Point.h
  \brief Define functions to treat a std::tr1::array as a point or vector.
*/

#if !defined(__geom_Point_h__)
#define __geom_Point_h__

#include "../../ext/array.h"

namespace geom {

/*!
  \page point Point Functions

  We define functions for treating an std::tr1::array as a point.
  We group the functions according to:
  \ref point_math "mathematical functions",
  \ref point_angle "angle functions" and
  \ref point_rotation "rotation functions".

  \todo Check the efficiency of these functions.
*/

//-----------------------------------------------------------------------------
/*! \defgroup point_math Point: Mathematical Functions */
// @{

//! Compute an orthogonal vector.
template<typename T>
void
computeAnOrthogonalVector(std::tr1::array<T, 3> vector,
                          std::tr1::array<T, 3>* orthogonal);

// @}
//-----------------------------------------------------------------------------
/*! \defgroup point_angle Point: Angles */
// @{


//! Positive turn: return 1.  No turn: return 0.  Negative turn: return -1.
template<typename T>
int
computeSignOfTurn(const std::tr1::array<T, 2>& p,
                  const std::tr1::array<T, 2>& q,
                  const std::tr1::array<T, 2>& r);


// CONTINUE: perhaps replace this with something better.
//! Positive turn: return 1.  No turn: return 0.  Negative turn: return -1.
template<typename T>
int
computeApproximateSignOfTurn(const std::tr1::array<T, 2>& p,
                             const std::tr1::array<T, 2>& q,
                             const std::tr1::array<T, 2>& r);


//! Return the pseudo-angle between vec and the x axis.
template<typename T>
T
computePseudoAngle(const std::tr1::array<T, 2>& vec);


// CONTINUE: Check the implementation.
//! Return the angle between the two vectors.
/*!
  The angle is in the range [0..pi].
*/
template<std::size_t N, typename T>
T
computeAngle(const std::tr1::array<T, N>& a, const std::tr1::array<T, N>& b);


// @}
//-----------------------------------------------------------------------------
/*! \defgroup point_rotation Point: Rotation */
// @{


//! Rotate the vector + pi / 2.
template<typename T>
void
rotatePiOver2(std::tr1::array<T, 2>* p);


//! Rotate the vector - pi / 2.
template<typename T>
void
rotateMinusPiOver2(std::tr1::array<T, 2>& p);


// @}

} // namespace geom

#define __geom_Point_ipp__
#include "Point.ipp"
#undef __geom_Point_ipp__

#endif
