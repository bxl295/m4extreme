/* -*- C++ -*- */

#if !defined(__levelSet_outsideCuda_h__)
#define __levelSet_outsideCuda_h__

#ifndef __CUDA_ARCH__

#include "../ext/array.h"

#include <vector>

namespace levelSet {

void
markOutsideAsNegativeInf
(const std::tr1::array<std::size_t, 3>& gridExtents,
 std::size_t numRefined,
 float* patchesDev,
 const uint3* indicesDev,
 const std::vector<std::tr1::array<std::size_t, 3> >& negativePatches,
 std::vector<bool>* outsideAtLowerCorners);

} // namespace levelSet

#endif

#endif
