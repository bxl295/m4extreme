// -*- C++ -*-

/*!
  \file geom/kernel/distance.h
  \brief Define functions to compute distance.
*/
/*!
  \page geom_distance Distance

  There are functions for computing upper and lower bounds for the distance to
  object(s) contained in a bounding box.

  signed_distance_upper_bound(const BBox<T,N>& box,const typename BBox<T,N>::point_type& x)

  signed_distance_lower_bound(const BBox<T,N>& box,const typename BBox<T,N>::point_type& x)

*/

#if !defined(__geom_distance_h__)
#define __geom_distance_h__

#include "BBox.h"
#include "Point.h"
#include "content.h"

#include "../../loki/TypeManip.h"

#include <limits>

namespace geom {

// CONTINUE: continue fixing.  I think I should remove the signed distance
// functions.

//! Compute an upper bound on the squared distance between the objects bounded by the specified boxes.
/*! Bounding boxes have the property that the bounded object must touch every
  face of box. Thus one may obtain on upper bound on the distance between
  the objects by computing maximum distances between the faces. */
template<typename _T>
_T
computeUpperBoundSquaredDistance(const BBox<_T, 3>& a, const BBox<_T, 3>& b);

//! Compute a lower bound on the squared distance between the object bounded by the box and the point.
template<std::size_t _N, typename _T>
_T
computeLowerBoundSquaredDistance(const BBox<_T, _N>& box,
                                 const std::tr1::array<_T, _N>& x);

//! Compute an upper bound on the squared distance between the object bounded by the box and the point.
template<std::size_t _N, typename _T>
_T
computeUpperBoundSquaredDistance(const BBox<_T, _N>& box,
                                 const std::tr1::array<_T, _N>& x);

//---------------------------------------------------------------------------
// Bounds on the distance to objects in a bounding box.
//---------------------------------------------------------------------------


//! Return an upper bound on the signed distance from the point to the objects in the box.
/*!
  Consider some objects contained in the bounding box: \c box.  (Specifically,
  the objects form an N-D manifold.  The distance from the point \c x to
  the manifold is signed; positive outside and negative inside.)  This
  function returns an upper bound on the distance to the objects.
*/
template<std::size_t N, typename T>
inline
T
computeUpperBoundOnSignedDistance(const BBox<T, N>& box,
                                  const typename BBox<T, N>::Point& x) {
   return computeUpperBoundOnSignedDistance(box, x, Loki::Int2Type < N != 1 > ());
}


//! Return an lower bound on the signed distance from the point to the objects in the box.
/*!
  Consider some objects contained in the bounding box: \c box.  (Specifically,
  the objects form an N-D manifold.  The distance from the point \c x to
  the manifold is signed; positive outside and negative inside.)  This
  function returns an lower bound on the distance to the objects.
*/
template<std::size_t N, typename T>
inline
T
computeLowerBoundOnSignedDistance(const BBox<T, N>& box,
                                  const typename BBox<T, N>::Point& x) {
   return computeLowerBoundOnSignedDistance(box, x, Loki::Int2Type < N != 1 > ());
}


//! Return an upper bound on the signed distance.
template<std::size_t N, typename T>
T
computeUpperBoundOnUnsignedDistance(const BBox<T, N>& box,
                                    const typename BBox<T, N>::Point& x);


//! Return a lower bound on the signed distance.
template<std::size_t N, typename T>
T
computeLowerBoundOnUnsignedDistance(const BBox<T, N>& box,
                                    const typename BBox<T, N>::Point& x);

} // namespace geom

#define __geom_distance_ipp__
#include "distance.ipp"
#undef __geom_distance_ipp__

#endif
