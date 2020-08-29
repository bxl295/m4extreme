// -*- C++ -*-

/*!
  \file BBox.h
  \brief Implements a class for an axis-oriented bounding box in N dimensions.
*/

#if !defined(__geom_BBox_h__)
#define __geom_BBox_h__

#include "bound.h"

#include "../../ext/array.h"

#include <vector>

#include <cmath>

namespace geom {

//! An axis-oriented bounding box in the specified dimension.
/*!
  \param _T is the number type.
  \param _Dim is the dimension.
  
  This class is an aggregate type. Thus it has no user-defined constructors.
*/
template<typename _T, std::size_t _Dim>
struct BBox {
   //
   // Types.
   //

   //! The number type.
   typedef _T Number;
   //! The point type.
   typedef std::tr1::array<_T, _Dim> Point;

   //
   // Constants.
   //

   //! The space dimension.
   static const std::size_t Dimension = _Dim;

   //
   // Member data
   //

   //! The lower corner.
   Point lower;
   //! The upper corner.
   Point upper;

   //
   // Member functions.
   //

   //! Return true if the BBox is empty.
   /*! An invalid bounding box has NaN for its coordinates. */
   bool
   isValid() const {
      return lower[0] == lower[0];
   }

   //! Invalidate the box by setting the coordinates to NaN.
   void
   invalidate() {
      for (std::size_t i = 0; i != Dimension; ++i) {
         lower[i] = std::numeric_limits<Number>::quiet_NaN();
         upper[i] = std::numeric_limits<Number>::quiet_NaN();
      }
   }

   // CONTINUE: See if I can remove this and use isValid instead.
   //! Return true if the BBox is empty.
   bool
   isEmpty() const;

   //! Offset (expand or contract) by the specified amount.
   void
   offset(const Number value) {
      lower -= value;
      upper += value;
   }

   //! Return true if the object is contained in this bounding box.
   template<typename _Object>
   bool
   isIn(const _Object& x) const {
      // Check if the bounding box around the object is contained in this 
      // bounding box.
      return isIn(geom::bound(x));
   }

   //! Return true if the point is in this bounding box.
   bool
   isIn(const Point& p) const;

   //! Return true if the bounding box \c x is in this bounding box.
   bool
   isIn(const BBox& x) const {
      return isIn(x.lower) && isIn(x.upper);
   }

   //! Make the bounding box expand to contain the object.
   template<typename _Object>
   void
   add(const _Object& x) {
      add(geom::bound(x));
   }

   //! Make the bounding box expand to contain the point.
   void
   add(const Point& p);

   //! Make the bounding box expand to contain the box.
   void
   add(const BBox& x) {
      add(x.lower);
      add(x.upper);
   }

   //! Bound the sequence of objects.
   template<typename _InputIterator>
   void
   bound(_InputIterator begin, _InputIterator end);
};


//
// Equality Operators.
//


//! Equality.
/*! \relates BBox */
template<typename _T, std::size_t _Dim>
inline
bool
operator==(const BBox<_T, _Dim>& a, const BBox<_T, _Dim>& b) {
   return (a.lower == b.lower && a.upper == b.upper);
}


//! Inequality.
/*! \relates BBox */
template<typename _T, std::size_t _Dim>
inline
bool
operator!=(const BBox<_T, _Dim>& a, const BBox<_T, _Dim>& b) {
   return !(a == b);
}


//
// File I/O.
//


//! Read the bounding box.
/*! \relates BBox */
template<typename _T, std::size_t _Dim>
inline
std::istream&
operator>>(std::istream& in, BBox<_T, _Dim>& x) {
   return in >> x.lower >> x.upper;
}


//! Write the bounding box.
/*! \relates BBox */
template<typename _T, std::size_t _Dim>
inline
std::ostream&
operator<<(std::ostream& out, const BBox<_T, _Dim>& x) {
   return out << x.lower << ' ' << x.upper;
}


//! Print in a nice format.
/*! \relates BBox */
template<typename _T, std::size_t _Dim>
void
printFormatted(std::ostream& out, const BBox<_T, _Dim>& x);


//
// Mathematical Functions.
//


//! Make a bounding box from the lower and upper corners.
template<typename _T, std::size_t _Dim>
inline
BBox<_T, _Dim>
makeBBox(const std::tr1::array<_T, _Dim>& lower,
         const std::tr1::array<_T, _Dim>& upper) {
   BBox<_T, _Dim> x = {lower, upper};
   return x;
}


//! Return the content (length, area, volume, etc.) of the box.
template<typename _T, std::size_t _Dim>
inline
_T
content(const BBox<_T, _Dim>& box) {
   _T x = 1;
   for (std::size_t n = 0; n != _Dim; ++n) {
      x *= box.upper[n] - box.lower[n];
   }
   return x;
}


//! Put a bounding box around a point.
template<typename _T, std::size_t _Dim>
inline
BBox<_T, _Dim>
bound(const std::tr1::array<_T, _Dim>& p) {
   BBox<_T, _Dim> box = {p, p};
   return box;
}


//! Put a bounding box around two points.
template<typename _T, std::size_t _Dim>
inline
BBox<_T, _Dim>
bound(const std::tr1::array<_T, _Dim>& p, const std::tr1::array<_T, _Dim>& q) {
   BBox<_T, _Dim> box;
   for (std::size_t i = 0; i != _Dim; ++i) {
      box.lower[i] = std::min(p[i], q[i]);
      box.upper[i] = std::max(p[i], q[i]);
   }
   return box;
}


//! Put a bounding box around three points.
template<typename _T, std::size_t _Dim>
inline
BBox<_T, _Dim>
bound(const std::tr1::array<_T, _Dim>& p, const std::tr1::array<_T, _Dim>& q,
      const std::tr1::array<_T, _Dim>& r) {
   BBox<_T, _Dim> box;
   for (std::size_t i = 0; i != _Dim; ++i) {
      box.lower[i] = std::min(p[i], std::min(q[i], r[i]));
      box.upper[i] = std::max(p[i], std::max(q[i], r[i]));
   }
   return box;
}


//! Put a bounding box around a box.
template<typename _T, std::size_t _Dim>
inline
const BBox<_T, _Dim>&
bound(const BBox<_T, _Dim>& x) {
   return x;
}

// CONTINUE: REMOVE
inline
BBox<int, 3>
bound(const BBox<int, 3>& x) {
   assert(false);
   return x;
}

// CONTINUE: REMOVE
#if 0
//! Compute a bounding box around the sequence of objects.
/*!
  The element type must be able to be bounded by the function
  <code>bound(const Element& x, BBox<_Dim, _T>* box)</code>.
  \pre The sequence must not be empty.
*/
template<typename _InputIterator, std::size_t _Dim, typename _T>
void
bound(_InputIterator begin, _InputIterator end, BBox<_Dim, _T>* box);
#endif


//! Return the squared distance between two 1-D intervals.
template<typename _T>
_T
squaredDistanceBetweenIntervals(const _T lower1, const _T upper1,
                                const _T lower2, const _T upper2);


//! Return the squared distance between two bounding boxes.
template<typename _T, std::size_t _Dim>
_T
squaredDistance(const BBox<_T, _Dim>& x, const BBox<_T, _Dim>& y);


//! Return true if the bounding boxes overlap.
/*! \relates BBox */
template<typename _T, std::size_t _Dim>
bool
doOverlap(const BBox<_T, _Dim>& a, const BBox<_T, _Dim>& b);


// CONTINUE
#if 0
//! Return true if the domains overlap.
/*! \relates BBox */
template<typename _T, std::size_t _Dim>
bool
doOverlap(const BBox<_T, _Dim>& a, const SemiOpenInterval<_T, _Dim>& b);
#endif


//! Return the intersection of the bounding boxes.
/*! \relates BBox */
template<typename _T, std::size_t _Dim>
BBox<_T, _Dim>
intersection(const BBox<_T, _Dim>& a, const BBox<_T, _Dim>& b);


//! Scan convert the index bounding box.
/*!
  \relates BBox

  \param indices is an output iterator for multi-indices.  The value type
  must be std::tr1::array<Index,3> assignable where Index is the integer
  index type.
  \param box describes the range of indices.  It is a bounding box of
  some floating point number type.  This box will be converted to an
  integer bounding box.  Then the below scan conversion function is used.
*/
template<typename Index, typename MultiIndexOutputIterator>
void
scanConvert(MultiIndexOutputIterator indices, const BBox<double, 3>& box);


//! Scan convert the index bounding box.
/*!
  \relates BBox

  \param indices is an output iterator for multi-indices.  The value type
  must be std::tr1::array<Index,3> assignable where Index is the integer
  index type.
  \param box describes the range of indices.
*/
template<typename MultiIndexOutputIterator, typename Index>
void
scanConvertIndex(MultiIndexOutputIterator indices, const BBox<Index, 3>& box);


//! Scan convert the index bounding box on the specified index domain.
/*!
  \relates BBox

  \param indices is an output iterator for multi-indices.  The value type
  must be std::tr1::array<Index,3> assignable where Index is the integer
  index type.
  \param box describes the range of indices.  It is a bounding box of
  some floating point number type.  This box will be converted to an
  integer bounding box.  Then the below scan conversion function is used.
  \param domain is the closed range of indices on which to perform the
  scan conversion.
*/
template<typename MultiIndexOutputIterator, typename Index>
void
scanConvert(MultiIndexOutputIterator indices, const BBox<double, 3>& box,
            const BBox<Index, 3>& domain);


//! Scan convert the index bounding box on the specified index domain.
/*!
  \relates BBox

  \param indices is an output iterator for multi-indices.  The value type
  must be std::tr1::array<Index,3> assignable where Index is the integer
  index type.
  \param box is the closed range of indices.
  \param domain is the closed range of indices on which to perform the
  scan conversion.
*/
template<typename MultiIndexOutputIterator, typename Index>
inline
void
scanConvert(MultiIndexOutputIterator indices, const BBox<Index, 3>& box,
            const BBox<Index, 3>& domain);

// CONTINUE: get rid of this
#if 0
//! Scan convert the BBox.  Return the number of points scan converted.
/*! \relates BBox */
template<typename T, typename _IndexType>
int
scanConvert(std::vector<_IndexType>& indexSet, const BBox<_T, 3>& box);
#endif


// CONTINUE: get rid of this
#if 0
//! Scan convert the BBox for a grid on the domain.
/*! \relates BBox */
template<typename T, typename _IndexType>
int
scanConvert(std::vector<_IndexType>& indexSet,
            const BBox<_T, 3>& box, const BBox<_T, 3>& domain);
#endif


} // namespace geom

#define __geom_BBox_ipp__
#include "BBox.ipp"
#undef __geom_BBox_ipp__

#endif
