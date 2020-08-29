// -*- C++ -*-

#if !defined(__levelSet_powerDistance_h__)
#define __levelSet_powerDistance_h__

#include "Grid.h"

#include "../geom/kernel/Ball.h"
#include "../geom/kernel/BallSquared.h"

namespace levelSet {

/*! \defgroup levelSetPowerDistance Power Distance
  These functions calculate the power distance to a set of balls.
*/
//@{

//! Construct a level set for the power distance to a set of balls.
/*! For each grid point, compute the distance to each of the balls. */
template<typename _T, std::size_t _D, std::size_t _N, typename _Base>
void
powerDistance(array::EquilateralArrayImp<_T, _D, _N, _Base>* patch,
              const std::tr1::array<_T, _D>& lowerCorner,
              _T spacing, const std::vector<geom::BallSquared<_T, _D> >& balls);


//! Construct a level set for the power distance to a set of balls.
template<typename _T, std::size_t _D, std::size_t _N>
void
negativePowerDistance(Grid<_T, _D, _N>* grid,
                      const std::vector<geom::Ball<_T, _D> >& balls);

//@}

} // namespace levelSet

#define __levelSet_powerDistance_ipp__
#include "powerDistance.ipp"
#undef __levelSet_powerDistance_ipp__

#endif
