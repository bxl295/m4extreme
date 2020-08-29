// -*- C++ -*-

/*!
  \file Ball.h
  \brief A ball in N-dimensional space.
*/

#if !defined(__geom_kernel_Ball_h__)
#define __geom_kernel_Ball_h__

#include "BBox.h"

namespace geom {

//! A ball in N-dimensional space.
/*!
  \param _T is the number type.
  \param _N is the dimension.

  A ball is defined by a center and a radius.
  This class is an aggregate type. Thus it has no user-defined constructors.
*/
template<typename _T, std::size_t _N>
struct Ball {
   //
   // Constants.
   //

   //! The space dimension.
   static const std::size_t Dimension = _N;

   //
   // Types
   //

   //! The number type.
   typedef _T Number;
   //! The representation of a point.
   typedef std::tr1::array<Number, Dimension> Point;

   //
   // Data
   //

   //! The center of the ball.
   Point center;
   //! The radius of the ball.
   Number radius;
};


//! Make a ball from the center and radius.
/*! \relates Ball */
template<typename _T, std::size_t _N>
inline
Ball<_T, _N>
makeBall(const std::tr1::array<_T, _N>& center, const _T radius) {
   Ball<_T, _N> x = {center, radius};
   return x;
}


//
// Equality Operators.
//

//! Return true if the balls are equal.
/*! \relates Ball */
template<typename _T, std::size_t _N>
inline
bool
operator==(const Ball<_T, _N>& x, const Ball<_T, _N>& y) {
   return (x.center == y.center && x.radius == y.radius);
}


//! Return true if the balls are not equal.
/*! \relates Ball */
template<typename _T, std::size_t _N>
inline
bool
operator!=(const Ball<_T, _N>& x, const Ball<_T, _N>& y) {
   return !(x == y);
}


//
// File I/O Operators.
//

//! Read a ball.
/*! \relates Ball */
template<typename _T, std::size_t _N>
inline
std::istream&
operator>>(std::istream& in, Ball<_T, _N>& x) {
   return in >> x.center >> x.radius;
}

//! Write the ball.
/*! \relates Ball */
template<typename _T, std::size_t _N>
inline
std::ostream&
operator<<(std::ostream& out, const Ball<_T, _N>& x) {
   return out << x.center << " " << x.radius;
}

//
// Mathematical operations.
//


//! Return true if the point is inside this ball.
template<typename _T, std::size_t _N>
inline
bool
isInside(const Ball<_T, _N>& ball, const std::tr1::array<_T, _N>& position) {
   return squaredDistance(ball.center, position) < ball.radius * ball.radius;
}


//! Make a bounding box around the ball.
/*! \relates Ball */
template<typename _T, std::size_t _N>
inline
BBox<_T, _N>
bound(const Ball<_T, _N>& ball) {
   BBox<_T, _N> box = {ball.center - ball.radius,
                       ball.center + ball.radius};
   return box;
}


//! Return true if the two balls intersect.
/*! \relates Ball */
template<typename _T, std::size_t _N>
inline
bool
doIntersect(const Ball<_T, _N>& a, const Ball<_T, _N>& b) {
   return squaredDistance(a.center, b.center) <= 
      (a.radius + b.radius) * (a.radius + b.radius);
}

//! Calculate the signed distance to the surface.
/*! \relates Ball */
template<typename _T, std::size_t _N>
inline
_T
distance(const Ball<_T, _N>& ball, const std::tr1::array<_T, _N>& x) {
   return euclideanDistance(ball.center, x) - ball.radius;
}

//! Calculate the closest point on the surface and return the signed distance.
/*! \relates Ball */
template<typename _T, std::size_t _N>
inline
_T
closestPoint(const Ball<_T, _N>& ball, const std::tr1::array<_T, _N>& x,
             std::tr1::array<_T, _N>* closest) {
   // Start at the point.
   *closest = x;
   // Translate to the origin.
   *closest -= ball.center;
   // Move to the surface.
   const _T r = magnitude(*closest);
   // Special case that the point is at the center.
   if (r < std::numeric_limits<_T>::epsilon()) {
      // Pick an arbitrary closest point.
      std::fill(closest->begin(), closest->end(), 0);
      (*closest)[0] = ball.radius;
   }
   else {
      *closest *= (ball.radius / r);
   }
   // Translate back to the ball.
   *closest += ball.center;
   // Return the distance.
   return r - ball.radius;
}

} // namespace geom

#endif
