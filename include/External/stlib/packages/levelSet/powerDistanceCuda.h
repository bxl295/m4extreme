/* -*- C++ -*- */

#if !defined(__levelSet_powerDistanceCuda_h__)
#define __levelSet_powerDistanceCuda_h__

#ifndef __CUDA_ARCH__

#include "cuda.h"
#include "Grid.h"
#include "../geom/kernel/BallSquared.h"
#include "../geom/kernel/Ball.h"

#include <vector>

namespace levelSet {

//! Construct a level set for the power distance to a set of balls.
/*! For each grid point, compute the distance to each of the balls. 
  Translate the data to CUDA format and call powerDistanceKernel(). */
void
powerDistanceCuda(array::EquilateralArray<float, 3, PatchExtent>* patch,
                  const std::tr1::array<float, 3>& lowerCorner,
                  float spacing,
                  const std::vector<geom::BallSquared<float, 3> >& balls);

//! Construct a level set for the power distance to a set of balls.
void
negativePowerDistanceCuda(Grid<float, 3, PatchExtent>* grid,
                          const std::vector<geom::Ball<float, 3> >& balls);

} // namespace levelSet

#endif

#endif
