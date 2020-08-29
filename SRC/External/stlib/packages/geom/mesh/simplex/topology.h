// -*- C++ -*-

/*!
  \file geom/mesh/simplex/topology.h
  \brief Topology functions for simplices.
*/

#if !defined(__geom_mesh_simplex_topology_h__)
#define __geom_mesh_simplex_topology_h__

#include "../../../ext/array.h"

namespace geom {

//-----------------------------------------------------------------------------
/*! \defgroup simplex_topology Topology Functions
*/
//@{

//---------------------------------------------------------------------------
// Simplex indices.
//---------------------------------------------------------------------------


//! Compute the other indices of the simplex.
void
computeOtherIndices(std::size_t i, std::size_t j, std::size_t* a,
                    std::size_t* b);


//! Compute the other index of the simplex.
std::size_t
computeOtherIndex(std::size_t i, std::size_t j, std::size_t k);


//! Reverse the orientation of the simplex.
template<typename _T, std::size_t _N>
inline
void
reverseOrientation(std::tr1::array<_T, _N>* simplex) {
   if (simplex->size() != 1) {
      std::swap((*simplex)[0], (*simplex)[1]);
   }
}

//! Get the face obtained by removing the n_th vertex.
/*!
  For the simplex (v[0], ... v[N]) the face is
  (-1)^n (v[0], ..., v[n-1], v[n+1], ..., v[N]).
*/
template<typename _T, std::size_t _N>
inline
void
getFace(const std::tr1::array<_T, _N>& simplex, const std::size_t n,
        std::tr1::array < _T, _N - 1 > * face) {
#ifdef DEBUG_stlib
   assert(n < simplex.size());
#endif
   std::size_t j = 0;
   for (std::size_t i = 0; i != simplex.size(); ++i) {
      if (i != n) {
         (*face)[j++] = simplex[i];
      }
   }
   if (n % 2 == 1) {
      reverseOrientation(face);
   }
}

//! Return the face obtained by removing the n_th vertex.
/*!
  For the simplex (v[0], ... v[N]) return
  (-1)^n (v[0], ..., v[n-1], v[n+1], ..., v[N]).
*/
template<typename _T, std::size_t _N>
inline
std::tr1::array < _T, _N - 1 >
getFace(const std::tr1::array<_T, _N>& simplex, const std::size_t n) {
   std::tr1::array < _T, _N - 1 > f;
   getFace(simplex, n, &f);
   return f;
}

//! Return true if the two simplices have the same orientation.
/*!
  \pre \c x and \c y must have the same vertices.
*/
template<typename _T>
inline
bool
#ifdef DEBUG_stlib
haveSameOrientation(const std::tr1::array < _T, 0 + 1 >& x,
                    const std::tr1::array < _T, 0 + 1 >& y) {
   assert(x[0] == y[0]);
   return true;
}
#else
haveSameOrientation(const std::tr1::array < _T, 0 + 1 >& /*x*/,
                    const std::tr1::array < _T, 0 + 1 >& /*y*/) {
   return true;
}
#endif

//! Return true if the two simplices have the same orientation.
/*!
  \pre \c x and \c y must have the same vertices.
*/
template<typename _T>
inline
bool
haveSameOrientation(const std::tr1::array < _T, 1 + 1 > & x,
                    const std::tr1::array < _T, 1 + 1 > & y) {
#ifdef DEBUG_stlib
   assert(hasElement(x, y[0]) && hasElement(x, y[1]));
#endif
   return x[0] == y[0];
}

//! Return true if the two simplices have the same orientation.
/*!
  \pre \c x and \c y must have the same vertices.
*/
template<typename _T>
inline
bool
haveSameOrientation(const std::tr1::array < _T, 2 + 1 > & x,
                    const std::tr1::array < _T, 2 + 1 > & y) {
#ifdef DEBUG_stlib
   assert(hasElement(x, y[0]) && hasElement(x, y[1]) && hasElement(x, y[2]));
#endif
   if (x[0] == y[0]) {
      return x[1] == y[1];
   }
   else if (x[0] == y[1]) {
      return x[1] == y[2];
   }
   // else x[0] == y[2]
   return x[1] == y[0];
}

//! Return true if the N-simplex has the specified (N-1)-face.
/*!
  If true, set the face index.
*/
template<typename _T, std::size_t _N, std::size_t _M>
inline
bool
hasFace(const std::tr1::array<_T, _N>& simplex,
        const std::tr1::array<_T, _M>& face, std::size_t* faceIndex) {
   // Loop over the vertices of the face.
   for (std::size_t i = 0; i != face.size(); ++i) {
      // If the vertex of the face is not in the simplex.
      if (! hasElement(simplex, face[i])) {
         *faceIndex = i;
         return false;
      }
   }
   // If we get here then all the vertices in the face are in the simplex.
   return true;
}

//! Return true if the simplex has the given face as a sub-simplex.
/*!
  This function does not check the orientation of the face.  It returns
  true if the simplex has each of the vertices in the face.
*/
template<typename _T, std::size_t _N, std::size_t _M>
inline
bool
hasFace(const std::tr1::array<_T, _N>& simplex,
        const std::tr1::array<_T, _M>& face) {
   std::size_t faceIndex;
   return hasFace(simplex, face, &faceIndex);
}

//! Return true if the 3-simplex has the face specified by the three vertices.
template<typename _T>
inline
bool
hasFace(const std::tr1::array < _T, 3 + 1 > & simplex,
        const _T& x, const _T& y, const _T& z) {
   return hasElement(simplex, x) && hasElement(simplex, y) &&
          hasElement(simplex, z);
}

//! Return true if the 3-simplex has the face specified by the three vertices.
/*!
  Set the face index.
*/
template<typename _T>
inline
bool
hasFace(const std::tr1::array < _T, 3 + 1 > & simplex,
        const _T& x, const _T& y, const _T& z,
        std::size_t* faceIndex) {
   if (hasFace(simplex, x, y, z)) {
      *faceIndex = computeOtherIndex(index(simplex, x),
                                     index(simplex, y),
                                     index(simplex, z));
      return true;
   }
   return false;
}

//@}

} // namespace geom

#define __geom_mesh_simplex_topology_ipp__
#include "topology.ipp"
#undef __geom_mesh_simplex_topology_ipp__

#endif
