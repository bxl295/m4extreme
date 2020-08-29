// -*- C++ -*-

#if !defined(__levelSet_negativeDistance_h__)
#define __levelSet_negativeDistance_h__

#include "GridUniform.h"
#include "Grid.h"
#include "geometry.h"
#include "IntersectionPoint.h"
#include "IntersectionCircle.h"

#include "../geom/grid/SimpleRegularGrid.h"

namespace levelSet {


/*! \defgroup levelSetNegativeDistance Negative Distance
These functions calculate the Negative distance for the union of a set of balls.
*/
//@{

//! Compute the negative distances for a union of balls.
/*!
  The distance is correctly computed up to the minimum radius of the set 
  of balls. Set the positive distances to NaN. Set the negative far-away
  distances to \f$-\infty\f$.
*/
template<typename _T, std::size_t _D>
void
negativeDistance(GridUniform<_T, _D>* grid,
                 const std::vector<geom::Ball<_T, _D> >& balls);


//! Compute the negative distances for a union of balls.
/*!
  The distance is correctly computed up to the minimum radius of the set 
  of balls. Set the positive distances to NaN. Set the negative far-away
  distances to \f$-\infty\f$.
*/
template<typename _T, std::size_t _D, std::size_t _N>
void
negativeDistance(Grid<_T, _D, _N>* grid,
                 const std::vector<geom::Ball<_T, _D> >& balls);


//@}

} // namespace levelSet

#define __levelSet_negativeDistance_ipp__
#include "negativeDistance.ipp"
#undef __levelSet_negativeDistance_ipp__

#endif
