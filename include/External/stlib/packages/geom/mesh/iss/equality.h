// -*- C++ -*-

/*!
  \file geom/mesh/iss/equality.h
  \brief Equality operators for IndSimpSet.
*/

#if !defined(__geom_mesh_iss_equality_h__)
#define __geom_mesh_iss_equality_h__

#include "IndSimpSetIncAdj.h"

namespace geom {

//-----------------------------------------------------------------------------
/*! \defgroup iss_equality Equality
  Test the equality of two indexed simplex sets.  These functions are for
  debugging purposes only.  They don't do anything fancy like check if the
  vertices or simplices are given in different order.
*/
// @{

//! Return true if the vertices and indexed simplices are equal.
/*! \relates IndSimpSet */
template < std::size_t N, std::size_t M,
         typename T >
inline
bool
operator==(const IndSimpSet<N, M, T>& x,
           const IndSimpSet<N, M, T>& y) {
   return (x.getVertices() == y.getVertices() &&
           x.getIndexedSimplices() == y.getIndexedSimplices());
}

//! Return true if the vertices and indexed simplices are not equal.
/*! \relates IndSimpSet */
template < std::size_t N, std::size_t M,
         typename T >
inline
bool
operator!=(const IndSimpSet<N, M, T>& x,
           const IndSimpSet<N, M, T>& y) {
   return !(x == y);
}

//
// Equality.
//

//! Return true if the meshes are equal.
/*! \relates IndSimpSetIncAdj */
template < std::size_t N, std::size_t M,
         typename T >
inline
bool
operator==(const IndSimpSetIncAdj<N, M, T>& x,
           const IndSimpSetIncAdj<N, M, T>& y) {
   return (x.getVertices() == y.getVertices() &&
           x.getIndexedSimplices() == y.getIndexedSimplices() &&
           x.getVertexSimplexIncidence() == y.getVertexSimplexIncidence() &&
           x.getSimplexAdjacencies() == y.getSimplexAdjacencies());
}

//! Return true if the meshes are not equal.
/*! \relates IndSimpSetIncAdj */
template < std::size_t N, std::size_t M,
         typename T >
inline
bool
operator!=(const IndSimpSetIncAdj<N, M, T>& x,
           const IndSimpSetIncAdj<N, M, T>& y) {
   return !(x == y);
}


// @}

} // namespace geom

#endif
