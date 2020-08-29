// -*- C++ -*-

#if !defined(__levelSet_Grid_h__)
#define __levelSet_Grid_h__

#include "Patch.h"
#include "GridUniform.h"

#include "../array/EquilateralArray.h"
#include "../array/SimpleMultiArray.h"
#include "../array/SimpleMultiIndexRangeIterator.h"
#include "../array/StaticArrayOfArrays.h"
#include "../geom/kernel/BBox.h"
#include "../numerical/constants/Logarithm.h"

#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace levelSet {

//! A grid is a multi-array of patches.
/*!
  \param _T The value type.
  \param _D The dimension
  \param _N The extent in each dimension for a patch.
  \param _R The real number type. By default it is _T.
*/
template<typename _T, std::size_t _D, std::size_t _N, typename _R = _T>
class Grid :
      public array::SimpleMultiArray<Patch<_T, _D, _N>, _D> {
   //
   // Constants.
   //
public:

   //! The space dimension.
   static const std::size_t Dimension = _D;
   //! The number of vertices (and also the number of voxels) per patch.
   static const std::size_t NumVerticesPerPatch =
      numerical::Exponentiation<std::size_t, _N, _D>::Result;
   //! The logarithm (base 2) of the patch extent.
   static const std::size_t LogarithmOfExtent =
      numerical::Logarithm<std::size_t, 2, _N>::Result;

   //
   // Types.
   //
private:

   //! A multi-array of patches.
   typedef array::SimpleMultiArray<Patch<_T, _D, _N>, _D> Base;

public:
   
   //! A patch without ghost elements.
   /*! This is used for working with vertices. */
   typedef Patch<_T, _D, _N> VertexPatch;
   //! A patch that has been padded on the upper sides with one layer of ghost elements.
   /*! This is used for working with voxels. */
   typedef array::EquilateralArray<_T, _D, _N+1> VoxelPatch;
   //! A Cartesian point.
   typedef std::tr1::array<_R, Dimension> Point;
   //! A bounding box.
   typedef geom::BBox<_R, Dimension> BBox;
   //! A pair of grid/patch indices.
   typedef std::pair<typename Base::IndexList, typename Base::IndexList>
   DualIndices;

   //
   // Member data.
   //
public:

   //! The Cartesian coordinates of the lower corner of the grid.
   const Point lowerCorner;
   //! The grid spacing for the patches.
   const _R spacing;

protected:

   //! The data for all of the patches is contiguous in memory.
   std::vector<_T> _data;

   //--------------------------------------------------------------------------
   /*! \name Constructors etc.
     We use the default destructor. The copy constructor and assignment
     operator are disabled.
   */
   // @{
public:

   //! Construct from the Cartesian domain and the suggested grid patch spacing.
   /*!
     The grid spacing will be no greater than the suggested grid spacing and
     is the same in all dimensions. The domain will be expanded in the upper
     limits to exactly accomodate the grid.
   */
   Grid(const BBox& domain, const _R targetSpacing);

   //! Construct from the Cartesian domain and the grid extents (number of patches).
   Grid(const BBox& domain, const typename Base::IndexList& extents);

private:

   //! Calculate the extents.
   static
   typename Base::IndexList
   calculateExtents(const BBox& domain, const _R targetSpacing);

   //! The copy constructor is not implemented.
   Grid(const Grid&);

   //! The assignment operator is not implemented.
   Grid&
   operator=(const Grid&);

   // @}
   //--------------------------------------------------------------------------
   //! \name Accessors.
   // @{
public:

   //! The data for all of the patches is contiguous in memory.
   const _T*
   data() const {
      return &_data[0];
   }

   using Base::operator();

   //! Return the value at the specified grid point.
   _T
   operator()(const typename Base::IndexList& patch,
              const typename Base::IndexList& index) const;

   //! Return the value at the specified grid point.
   _T
   operator()(const DualIndices& pair) const {
      return operator()(pair.first, pair.second);
   }

   //! Return the Cartesian domain spanned by the grid.
   BBox
   domain() const;

   //! Return the length of a side of a voxel patch.
   _R
   getVoxelPatchLength() const {
      return spacing * _N;
   }

   //! Return the domain for the specified vertex patch.
   BBox
   getVertexPatchDomain(const typename Base::IndexList& i) const;

   //! Return the lower corner of the Cartesian domain for the specified patch.
   /*! Note that it doesn't matter whether one is referring to a vertex patch
    or a voxel patch. */
   Point
   getPatchLowerCorner(const typename Base::IndexList& i) const;

   //! Get the specified voxel patch.
   /*! The values on the upper boundaries are copied from the adjacent, refined 
    vertex patches. If there is no refined, adjacent vertex patch, the
    nearest value in the vertex patch is used.

    \pre The requested patch must be refined.
   */
   void
   getVoxelPatch(const typename Base::IndexList& i, VoxelPatch* patch) const;
      
   //! Return the Cartesian position of the specified vertex.
   Point
   indexToLocation(const typename Base::IndexList& patch,
                   const typename Base::IndexList& index) const;

   //! Get the number of refined patches.
   std::size_t
   numRefined() const;

   //! Get the total number of vertices in the refined patches.
   std::size_t
   numVertices() const {
      return _data.size();
   }

   //! Report the specified set of grid points as patch/grid multi-index pairs.
   template<typename _OutputIterator>
   void
   report(const typename Base::IndexList& patch,
          const typename Base::Range& range,
          _OutputIterator neighbors) const;

   //! Report the set of adjacent neighbors in the specified direction.
   template<typename _OutputIterator>
   void
   adjacentNeighbors(const DualIndices& pair,
                     const std::size_t direction,
                     _OutputIterator neighbors) const {
      adjacentNeighbors(pair.first, pair.second, direction, neighbors);
   }

   //! Report the set of adjacent neighbors in the specified direction.
   template<typename _OutputIterator>
   void
   adjacentNeighbors(const typename Base::IndexList& patch,
                     const typename Base::IndexList& index,
                     const std::size_t direction,
                     _OutputIterator neighbors) const;

   //! Report the set of all adjacent neighbors.
   template<typename _OutputIterator>
   void
   adjacentNeighbors(const DualIndices& pair, _OutputIterator neighbors) const;

   //! Report the set of all neighbors.
   /*! For the refined case, there are 3^_D - 1 neighbors. */
   template<typename _OutputIterator>
   void
   allNeighbors(const DualIndices& pair, _OutputIterator neighbors) const;

   //! Return true if the grid is valid.
   bool
   isValid() const;

   // @}
   //--------------------------------------------------------------------------
   //! \name Manipulators.
   // @{
public:

   //! The data for all of the patches is contiguous in memory.
   _T*
   data() {
      return &_data[0];
   }

   //! Return a reference to the value at the specified grid point.
   _T&
   operator()(const typename Base::IndexList& patch,
              const typename Base::IndexList& index);

   //! Return the value at the specified grid point.
   _T&
   operator()(const DualIndices& pair) {
      return operator()(pair.first, pair.second);
   }

   //! Make all of the grids unrefined.
   /*! The grid extents and domain are not altered. */
   void
   clear();

   //! Refine the specified patches.
   void
   refine(const std::vector<std::size_t>& indices);

   //! Refine the patches that have one or more dependencies.
   void
   refine(const array::StaticArrayOfArrays<std::size_t>& dependencies);

   //! Remove unecessary refinement.
   void
   coarsen();

   // @}
};


//! Add a constant to all vertices and fill values.
template<typename _T, std::size_t _D, std::size_t _N, typename _R>
Grid<_T, _D, _N, _R>&
operator+=(Grid<_T, _D, _N, _R>& grid, _T x);


//! Subtract a constant from all vertices and fill values.
template<typename _T, std::size_t _D, std::size_t _N, typename _R>
Grid<_T, _D, _N, _R>&
operator-=(Grid<_T, _D, _N, _R>& grid, _T x);


//! Determine the objects whose bounding boxes overlap each patch.
/*! \relates Grid */
template<typename _T, std::size_t _D, std::size_t _N, typename _R,
         typename _InputIterator>
void
patchDependencies(const Grid<_T, _D, _N, _R>& grid, _InputIterator begin,
                  _InputIterator end, 
                  array::StaticArrayOfArrays<std::size_t>* dependencies);


//! Determine the patches that intersect the bounding box.
/*! \relates Grid */
template<typename _T, std::size_t _D, std::size_t _PatchExtent, typename _R,
         typename _OutputIterator>
void
getIntersectingPatches(const Grid<_T, _D, _PatchExtent, _R>& grid,
                       geom::BBox<_R, _D> box, _OutputIterator indices);


//! Write the grid in VTK XML format.
/*! \relates Grid */
template<typename _T, std::size_t _N, typename _R>
void
writeVtkXml(const Grid<_T, 3, _N, _R>& grid, std::ostream& out);


//! Write the grid in VTK XML format.
/*! \relates Grid */
template<typename _T, std::size_t _N, typename _R>
void
writeVtkXml(const Grid<_T, 2, _N, _R>& grid, std::ostream& out);


//! Print information about the grid.
/*! \relates Grid */
template<typename _T, std::size_t _D, std::size_t _N, typename _R>
void
printInfo(const Grid<_T, _D, _N, _R>& grid, std::ostream& out);


} // namespace levelSet

#define __levelSet_Grid_ipp__
#include "Grid.ipp"
#undef __levelSet_Grid_ipp__

#endif
