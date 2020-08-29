// -*- C++ -*-

/*!
  \file PointsOnManifold.h
  \brief Represent the features of a mesh.
*/

#if !defined(__geom_PointsOnManifold_h__)
#define __geom_PointsOnManifold_h__

#include "build.h"
#include "geometry.h"
#include "quality.h"

#include "../simplicial/SimpMeshRed.h"
#include "../simplicial/geometry.h"
#include "../../kernel/SegmentMath.h"

#include "../../../ads/algorithm/OrderedPair.h"
#include "../../../array/MultiArray.h"
#include "../../../numerical/constants.h"

#include <set>
#include <map>

namespace geom {

//! Local closest point to a simplicial complex.
/*!
  \param N is the space dimension.
  \param M is the simplex dimension.
  \param SD is the spline degree.
  \param T is the number type.  By default it is double.
*/
template < std::size_t N, std::size_t M, std::size_t SD, typename T = double >
class PointsOnManifold;

} // namespace geom

// Include the implementations for 3-2 and N-1 meshes.
#include "PointsOnManifoldN11.h"
#include "PointsOnManifold321.h"

#endif
