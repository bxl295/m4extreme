// -*- C++ -*-

#if !defined(__levelSet_positiveDistance_h__)
#define __levelSet_positiveDistance_h__

#include "GridUniform.h"
#include "Grid.h"

#include "../geom/grid/SimpleRegularGrid.h"
#include "../geom/kernel/Ball.h"

namespace levelSet {


/*! \defgroup levelSetPositiveDistance Positive Distance
These functions calculate the positive distance for the union of a set of balls.
*/
//@{

//! Construct a level set for a union of balls.
/*!
  The level set function is the positive distance from the surface for points
  outside of the union of balls. For points inside the balls the function
  has negative values, but due to intersecting characteristics it is not
  necessarily the signed distance.
*/
template<typename _T, std::size_t _D>
void
positiveDistance(array::SimpleMultiArrayRef<_T, _D>* grid,
                 const geom::BBox<_T, _D>& domain,
                 const std::vector<geom::Ball<_T, _D> >& balls,
                 _T offset = 0, _T maxDistance = 0);


//! Construct a level set for a union of balls.
/*!
  The level set function is the positive distance from the surface for points
  outside of the union of balls. For points inside the balls the function
  has negative values, but due to intersecting characteristics it is not
  necessarily the signed distance.
*/
template<typename _T, std::size_t _D>
inline
void
positiveDistance(GridUniform<_T, _D>* grid,
                 const std::vector<geom::Ball<_T, _D> >& balls,
                 _T offset = 0, _T maxDistance = 0) {
   positiveDistance(grid, grid->domain(), balls, offset, maxDistance);
}


//! Construct a level set for a union of balls.
/*!
  The level set function is the positive distance from the surface for points
  outside of the union of balls. For points inside the balls the function
  has negative values, but due to intersecting characteristics it is not
  necessarily the signed distance.
*/
template<typename _T, std::size_t _D, std::size_t _N>
void
positiveDistance(Grid<_T, _D, _N>* grid,
                 const std::vector<geom::Ball<_T, _D> >& balls,
                 _T offset = 0, _T maxDistance = 0);


//@}

} // namespace levelSet

#define __levelSet_positiveDistance_ipp__
#include "positiveDistance.ipp"
#undef __levelSet_positiveDistance_ipp__

#endif
