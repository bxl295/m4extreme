/* -*- C++ -*- */

#ifndef __levelSet_gridCuda_h__
#define __levelSet_gridCuda_h__

#include "cuda.h"
#include "../cuda/check.h"

namespace levelSet {

//! Allocate device memory for the grid.
void
allocateGridCuda(const Grid<float, 3, PatchExtent>& grid, float** patchesDev,
                 uint3** indicesDev);

} // namespace levelSet

#define __levelSet_gridCuda_ipp__
#include "gridCuda.ipp"
#undef __levelSet_gridCuda_ipp__

#endif
