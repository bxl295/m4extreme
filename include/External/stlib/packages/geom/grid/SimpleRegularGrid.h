// -*- C++ -*-

/*!
  \file RegularGrid.h
  \brief A class for a regular grid.
*/

#if !defined(__geom_SimpleRegularGrid_h__)
#define __geom_SimpleRegularGrid_h__

#include "../kernel/BBox.h"
#include "../../array/SimpleMultiIndexRange.h"

#include <iostream>
#include <limits>

namespace geom {

//! A regular grid in N-D.
/*!
  \param _T is the number type.
  \param _Dimension is the dimension.
*/
template<typename _T, std::size_t _Dimension>
class SimpleRegularGrid {
   //
   // Friends.
   //

   template<typename T_, std::size_t Dimension_>
   friend
   std::istream&
   operator>>(std::istream&, SimpleRegularGrid<T_, Dimension_>&);

   //
   // Types.
   //
public:

   //! The number type.
   typedef _T Number;
   //! An array index type is \c std::size_t.
   typedef std::size_t Index;
   //! A list of indices.
   typedef std::tr1::array<std::size_t, _Dimension> IndexList;
   //! A multi-index range.
   typedef array::SimpleMultiIndexRange<_Dimension> Range;
   //! The point type.
   typedef std::tr1::array<Number, _Dimension> Point;
   //! A bounding box.
   typedef geom::BBox<Number, _Dimension> BBox;

   //
   // Data
   //
protected:

   //! Extents of the grid.
   IndexList _extents;
   //! The lower corner of the domain spanned by the grid.
   Point _lower;
   //! Grid spacing.
   Point _delta;
   //! The inverse of the grid spacing.
   /*! Store this redundant data to avoid costly divisions. */
   Point _inverseDelta;


public:

   //--------------------------------------------------------------------------
   /*! \name Constructors etc.
     We use the synthesized copy constructor, assignment operator, 
     and destructor.
   */
   // @{

   //! Default constructor.  Uninitialized memory.
   SimpleRegularGrid() :
      _extents(),
      _lower(),
      _delta(),
      _inverseDelta() {
   }

   //! Construct from grid dimensions and a Cartesian domain.
   /*!
     Construct a regular grid given the grid extents and the Cartesian domain
     that the grid spans.

     \param extents the number of grid points in each direction.
     \param domain the Cartesian domain spanned by the grid.
   */
   SimpleRegularGrid(const IndexList& extents, const BBox& domain);

   // @}
   //--------------------------------------------------------------------------
   //! \name Accesors.
   // @{

   //! Return the grid dimensions.
   const IndexList&
   getExtents() const {
      return _extents;
   }

   //! Return the lower corner of the domain spanned by the grid.
   const Point&
   getLower() const {
      return _lower;
   }

   //! Return the grid spacings.
   const Point&
   getDelta() const {
      return _delta;
   }

   // @}
   //--------------------------------------------------------------------------
   //! \name Mathematical functions.
   // @{

   //! Convert a Cartesian coordinate to a grid index coordinate.
   void
   locationToIndex(Point* p) const {
      for (std::size_t i = 0; i != p->size(); ++i) {
         (*p)[i] = ((*p)[i] - _lower[i]) * _inverseDelta[i];
      }
   }

   //! Convert a vector in Cartesian coordinates to index coordinates.
   void
   vectorToIndex(Point* p) const {
      *p *= _inverseDelta;
   }

   //! Convert a grid index coordinate to a Cartesian coordinate.
   void
   indexToLocation(Point* p) const {
      for (std::size_t i = 0; i != p->size(); ++i) {
         (*p)[i] = _lower[i] + (*p)[i] * _delta[i];
      }
   }

   //! Convert a grid index coordinate to a Cartesian coordinate.
   Point
   indexToLocation(const IndexList& index) const {
      Point p = _lower;
      for (std::size_t i = 0; i != p.size(); ++i) {
         p[i] += index[i] * _delta[i];
      }
      return p;
   }

   //! Calculate the index range that contains all of the indices in the box.
   /*! The window is first converted to index coordinates and then rounded. */
   Range
   computeRange(const BBox& box) const {
      IndexList lower, upper;

      Point p = box.lower;
      locationToIndex(&p);
      for (std::size_t i = 0; i != lower.size(); ++i) {
         // Round up to the closed lower bound.
         // Use only points in the grid range.
         lower[i] = Index(std::max(std::ceil(p[i]), Number(0.)));
      }

      p = box.upper;
      locationToIndex(&p);
      for (std::size_t i = 0; i != upper.size(); ++i) {
         // Round up to get an open upper bound.
         // Use only points in the grid range.
         // Indicate an empty range by setting the upper bound to the lower
         // bound.
         const Number u = std::ceil(p[i]);
         // Check the special case that the coordinate is negative. Otherwise
         // we would cast it to an unsigned integer.
         if (u < 0) {
            upper[i] = lower[i];
         }
         else {
            upper[i] = std::max(std::min(Index(u), _extents[i]), lower[i]);
         }
      }

      Range range = {upper - lower, lower};
      return range;
   }

   //! Convert a set of Cartesian coordinates to grid index coordinates.
   template<typename OutputIterator>
   void
   locationsToIndices(OutputIterator begin, OutputIterator end) const {
      for (; begin != end; ++begin) {
         locationToIndex(&*begin);
      }
   }

   //! Convert a set of grid index coordinates to Cartesian coordinates.
   template<typename OutputIterator>
   void
   indicesToLocations(OutputIterator begin, OutputIterator end) const {
      for (; begin != end; ++begin) {
         indexToLocation(&*begin);
      }
   }

   //! Convert a set of vectors in Cartesian coordinates to index coordinates.
   template<typename OutputIterator>
   void
   vectorsToIndices(OutputIterator begin, OutputIterator end) const {
      for (; begin != end; ++begin) {
         vectorToIndex(&*begin);
      }
   }

   //! Convert the Cartesian coordinates of the bounding box to a grid index coordinates.
   void
   locationsToIndices(BBox* box) const {
      convertLocationToIndex(&box->lower);
      convertLocationToIndex(&box->upper);
   }

   //! Convert the grid index coordinates of a bounding box to Cartesian coordinates.
   void
   indicesToLocations(BBox* box) const {
      convertIndexToLocation(&box->lower);
      convertIndexToLocation(&box->upper);
   }

   //! Convert the grid index coordinates of the first bounding box to Cartesian coordinates in the second.
   void
   indicesToLocations(const geom::BBox<Index, _Dimension>& indexBox,
                      BBox* cartesianBox) const {
      // Lower corner.
      Point p;
      // First convert to the number type.
      for (std::size_t i = 0; i != p.size(); ++i) {
         p[i] = indexBox.lower[i];
      }
      convertIndexToLocation(&p);
      cartesianBox->lower = p;
      // Upper corner.
      for (std::size_t i = 0; i != p.size(); ++i) {
         p[i] = indexBox.upper[i];
      }
      convertIndexToLocation(&p);
      cartesianBox->upper = p;
   }

   // @}
};


//
// Equality operators.
//


//! Return true if the true SimpleRegularGrid's are equal.
/*! \relates SimpleRegularGrid */
template<typename _T, std::size_t _Dimension>
bool
operator==(const SimpleRegularGrid<_T, _Dimension>& a,
           const SimpleRegularGrid<_T, _Dimension>& b);


//! Return true if the true SimpleRegularGrid's are not equal.
/*! \relates SimpleRegularGrid */
template<typename _T, std::size_t _Dimension>
inline
bool
operator!=(const SimpleRegularGrid<_T, _Dimension>& a,
           const SimpleRegularGrid<_T, _Dimension>& b) {
   return !(a == b);
}


//
// File I/O
//


//! Write to a file stream.
/*! \relates SimpleRegularGrid */
template<typename _T, std::size_t _Dimension>
std::ostream&
operator<<(std::ostream& out, const SimpleRegularGrid<_T, _Dimension>& grid);


//! Read from a file stream.
/*! \relates SimpleRegularGrid */
template<typename _T, std::size_t _Dimension>
std::istream&
operator>>(std::istream& in, SimpleRegularGrid<_T, _Dimension>& grid);


} // namespace geom

#define __geom_SimpleRegularGrid_ipp__
#include "SimpleRegularGrid.ipp"
#undef __geom_SimpleRegularGrid_ipp__

#endif
