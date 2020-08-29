/* -*- C++ -*- */

#ifndef __levelSet_ballsCuda_h__
#define __levelSet_ballsCuda_h__

#include "cuda.h"
#include "../cuda/check.h"

namespace levelSet {

//! Allocate device memory for the ball indices.
void
allocateBallIndicesCuda
(const array::StaticArrayOfArrays<std::size_t>& dependencies,
 std::size_t** ballIndexOffsetsDev,
 std::size_t** packedBallIndicesDev);

                        

} // namespace levelSet

#define __levelSet_ballsCuda_ipp__
#include "ballsCuda.ipp"
#undef __levelSet_ballsCuda_ipp__

#endif
