/* -*- C++ -*- */

#ifndef __levelSet_solventAccessibleCavitiesCuda_h__
#define __levelSet_solventAccessibleCavitiesCuda_h__

#include "positiveDistanceCuda.h"
#include "powerDistanceCuda.h"
#include "gridCuda.h"
#include "outsideCuda.h"
#include "solventExcluded.h"

#include <vector>

#include <cassert>

namespace levelSet {

//! Construct a level set for the solvent-accessible surface for solvents inside the protein.
void
solventAccessibleCavitiesCuda(Grid<float, 3, PatchExtent>* grid,
                              const std::vector<geom::Ball<float, 3> >& balls,
                              float probeRadius);

} // namespace levelSet

#define __levelSet_solventAccessibleCavitiesCuda_ipp__
#include "solventAccessibleCavitiesCuda.ipp"
#undef __levelSet_solventAccessibleCavitiesCuda_ipp__

#endif
