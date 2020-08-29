// -*- C++ -*-

/*!
  \file surfaceArea/sphereMeshing.h
  \brief Calculate solvent-accessible surface areas considering pairwise interactions.
*/

#if !defined(__surfaceArea_sphereMeshing_h__)
#define __surfaceArea_sphereMeshing_h__

#include "pointsOnSphere.h"

#include "../geom/kernel/BBox.h"
#include "../ext/vector.h"

#include <vector>

namespace surfaceArea {

//! Clip the points on a sphere with another sphere.
/*! The input is a set of points along with the number of active points.
 The active points are in the range [0..numActive). The inactive points
 are in the range [numActive..points->size()). The active points are
 clipped by a sphere (defined by the specified center and radius).
 Clipped points are moved into the inactive range. The number of active
 points after clipping is returned. */
template<typename _Number>
std::size_t
clip(std::vector<std::tr1::array<_Number, 3> >* points,
     std::size_t numActive, const std::tr1::array<_Number, 3>& center,
     const _Number radius);

} // namespace surfaceArea

#define __surfaceArea_sphereMeshing_ipp__
#include "sphereMeshing.ipp"
#undef __surfaceArea_sphereMeshing_ipp__

#endif
