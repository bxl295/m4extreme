// -*- C++ -*-

/*!
  \file geom/mesh/simplex/geometry.h
  \brief Geometric functions for simplices.
*/

#if !defined(__geom_mesh_simplex_geometry_h__)
#define __geom_mesh_simplex_geometry_h__

#include "topology.h"

#include "../../kernel/BBox.h"
#include "../../kernel/Plane.h"

#include "../../../numerical/constants.h"

namespace geom {

//-----------------------------------------------------------------------------
/*! \defgroup simplex_geometry Geometric Functions
*/
//@{

//! Calculate a bounding box around the simplex.
template<typename _T, std::size_t _N, typename _Number>
inline
void
computeBBox(const std::tr1::array<_T, _N>& simplex,
            BBox<_Number, _N - 1> * bb) {
   bb->bound(simplex.begin(), simplex.end());
}

//! Calculate the centroid of the simplex.
template<typename _T, std::size_t _N>
inline
void
computeCentroid(const std::tr1::array<_T, _N>& simplex, _T* centroid) {
   typedef typename _T::value_type Number;
   std::fill(centroid->begin(), centroid->end(), 0);
   for (std::size_t i = 0; i != simplex.size(); ++i) {
      *centroid += simplex[i];
   }
   *centroid /= Number(simplex.size());
}

//---------------------------------------------------------------------------
// Angles
//---------------------------------------------------------------------------

//! The dihedral angle between two faces.
template<typename _T>
_T
computeAngle(const std::tr1::array < std::tr1::array<_T, 3>, 3 + 1 > & s, std::size_t a,
             std::size_t b);


//! The solid angle at a vertex.
template<typename _T>
_T
computeAngle(const std::tr1::array < std::tr1::array<_T, 3>, 3 + 1 > & s,
             std::size_t n);

//! The interior angle at a vertex.
template<typename _T>
_T
computeAngle(const std::tr1::array < std::tr1::array<_T, 2>, 2 + 1 > & s,
             std::size_t n);

//! The interior angle at a vertex is 1.
template<typename _T>
_T
computeAngle(const std::tr1::array < std::tr1::array<_T, 1>, 1 + 1 > & s,
             std::size_t n);


//---------------------------------------------------------------------------
// Project
//---------------------------------------------------------------------------


//! Project the simplex to a lower dimension.
template<typename _T>
void
projectToLowerDimension(const std::tr1::array < std::tr1::array<_T, 2>, 1 + 1 > & s,
                        std::tr1::array < std::tr1::array<_T, 1>, 1 + 1 > * t);


//! Project the simplex to a lower dimension.
template<typename _T>
void
projectToLowerDimension(const std::tr1::array < std::tr1::array<_T, 3>, 1 + 1 > & s,
                        std::tr1::array < std::tr1::array<_T, 1>, 1 + 1 > * t);


//! Project the simplex to a lower dimension.
template<typename _T>
void
projectToLowerDimension(const std::tr1::array < std::tr1::array<_T, 3>, 2 + 1 > & s,
                        std::tr1::array < std::tr1::array<_T, 2>, 2 + 1 > * t);

//@}

} // namespace geom

#define __geom_mesh_simplex_geometry_ipp__
#include "geometry.ipp"
#undef __geom_mesh_simplex_geometry_ipp__

#endif
