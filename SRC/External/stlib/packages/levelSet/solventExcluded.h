// -*- C++ -*-

#if !defined(__levelSet_solventExcluded_h__)
#define __levelSet_solventExcluded_h__

#include "boolean.h"
#include "countGrid.h"
#include "negativeDistance.h"
#include "positiveDistance.h"

#include "../hj/hj.h"

namespace levelSet {

/*! \defgroup levelSetSolventExcluded Solvent-Excluded Surface
These functions calculate level sets related to the solvent-excluded surface.
*/
//@{

//! Subtract the balls from the level set.
template<typename _T, std::size_t _D, std::size_t _N, typename _Base>
void
subtract(array::EquilateralArrayImp<_T, _D, _N, _Base>* patch,
         const std::tr1::array<_T, _D>& lowerCorner, _T spacing,
         const std::vector<geom::Ball<_T, _D> >& balls);



//! Construct a level set for the solvent-excluded surface.
template<typename _T, std::size_t _D>
void
solventExcluded(GridUniform<_T, _D>* grid,
                const std::vector<geom::Ball<_T, _D> >& balls,
                _T probeRadius);

//! Construct a level set for the solvent-excluded surface.
template<typename _T, std::size_t _D, std::size_t _PatchExtent>
void
solventExcluded(Grid<_T, _D, _PatchExtent>* grid,
                const std::vector<geom::Ball<_T, _D> >& balls,
                _T probeRadius);

//! Construct a level set for the solvent-accessible surface for solvents inside the protein.
template<typename _T, std::size_t _D>
void
solventAccessibleCavities(array::SimpleMultiArrayRef<_T, _D>* grid,
                          const geom::BBox<_T, _D>& domain,
                          const std::vector<geom::Ball<_T, _D> >& balls,
                          _T probeRadius);


//! Construct a level set for the solvent-accessible surface for solvents inside the protein.
template<typename _T, std::size_t _D>
inline
void
solventAccessibleCavities(GridUniform<_T, _D>* grid,
                          const std::vector<geom::Ball<_T, _D> >& balls,
                          const _T probeRadius) {
   solventAccessibleCavities(grid, grid->domain(), balls, probeRadius);
}
   

//! Construct a level set for the solvent-excluded cavities.
template<typename _T, std::size_t _D>
void
solventExcludedCavities(GridUniform<_T, _D>* grid,
                        const std::vector<geom::Ball<_T, _D> >& balls,
                        _T probeRadius);


//@}

} // namespace levelSet

#define __levelSet_solventExcluded_ipp__
#include "solventExcluded.ipp"
#undef __levelSet_solventExcluded_ipp__

#endif
