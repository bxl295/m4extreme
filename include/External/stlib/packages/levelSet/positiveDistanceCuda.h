/* -*- C++ -*- */

#if !defined(__levelSet_positiveDistanceCuda_h__)
#define __levelSet_positiveDistanceCuda_h__

#include "cuda.h"

#ifndef __CUDA_ARCH__

#include "Grid.h"
#include "../cuda/check.h"
#include "../geom/kernel/Ball.h"

#include <vector>

namespace levelSet {

//! Construct a level set for a union of balls.
/*!
  The level set function is the positive distance from the surface for points
  outside of the union of balls. For points inside the balls the function
  has negative values, but due to intersecting characteristics it is not
  necessarily the signed distance.
*/
void
positiveDistanceCuda(Grid<float, 3, PatchExtent>* grid,
                     const std::vector<geom::Ball<float, 3> >& balls,
                     float offset = 0, float maxDistance = 0);


//! Use this version when the device memory for the grid is reused.
void
positiveDistanceCuda
(std::size_t numRefined,
 float* patchesDev,
 uint3* indicesDev,
 float3 lowerCorner,
 float spacing,
 const std::vector<geom::Ball<float, 3> >& balls,
 float offset,
 const array::StaticArrayOfArrays<std::size_t>& dependencies);


inline
void
positiveDistanceDependencies
(const Grid<float, 3, PatchExtent>& grid,
 const std::vector<geom::Ball<float, 3> >& balls,
 const float offset, const float maxDistance,
 array::StaticArrayOfArrays<std::size_t>* dependencies) {
   // Offset the balls' radii to include the volume of calculated distance.
   std::vector<geom::Ball<float, 3> > offsetBalls(balls);
   for (std::size_t i = 0; i != offsetBalls.size(); ++i) {
      offsetBalls[i].radius += offset + maxDistance;
   }
   // Calculate the dependencies.
   patchDependencies(grid, offsetBalls.begin(), offsetBalls.end(),
                     dependencies);
}
                             

} // namespace levelSet

#endif

#endif
