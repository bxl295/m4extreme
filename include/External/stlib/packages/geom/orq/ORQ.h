// -*- C++ -*-

/*!
  \file geom/orq/ORQ.h
  \brief A base class for a data structure for doing orthogonal range queries
  in N-D.
*/

#if !defined(__geom_ORQ_h__)
#define __geom_ORQ_h__

#include "../kernel/BBox.h"

#include <limits>

namespace geom {

//! Base class for an orthogonal range query data structure in N-D.
/*!
  \param N The space dimension.
  \param _Record The record type.  All the derived data structures hold
  records.  The record type is most likely a pointer to a class or an iterator
  into a container.
  \param _MultiKey An N-tuple of the key type.  The multi-key type must be
  subscriptable.  For example, the multi-key type could be
  \c std::tr1::array<double,3> or \c double*.
  \param _Key The number type.
  \param _MultiKeyAccessor A functor that takes the record type as an argument
  and returns the multi-key for the record.  If possible, it should return
  a constant reference to the multi-key.
*/
template < std::size_t N,
         typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor >
class ORQ {
public:

   //
   // Public types.
   //

   //! The record type.
   typedef _Record Record;
   //! The multi-key type.
   typedef _MultiKey MultiKey;
   //! The key type.
   typedef _Key Key;
   //! The multy-key accessor.
   typedef _MultiKeyAccessor MultiKeyAccessor;

   //! A Cartesian point.
   typedef std::tr1::array<Key, N> Point;
   //! Bounding box.
   typedef geom::BBox<Key, N> BBox;

private:

   //! Number of records in the data structure.
   std::size_t _size;

protected:

   //! The multi-key accessor.
   MultiKeyAccessor _multiKeyAccessor;

public:

   //--------------------------------------------------------------------------
   //! \name Constructors.
   // @{

   //! Default constructor.
   ORQ() :
      _size(0),
      _multiKeyAccessor() {
   }

   //! Copy constructor.
   ORQ(const ORQ& other) :
      _size(other._size),
      _multiKeyAccessor(other._multiKeyAccessor) {}

   // @}
   //--------------------------------------------------------------------------
   //! \name Assignment operator.
   // @{

   //! Assignment operator.
   ORQ&
   operator=(const ORQ& other) {
      if (this != &other) {
         _size = other._size;
         _multiKeyAccessor = other._multiKeyAccessor;
      }
      return *this;
   }

   // @}
   //--------------------------------------------------------------------------
   //! \name Accessors.
   // @{

   //! Return the number of records.
   std::size_t
   getSize() const {
      return _size;
   }

   //! Return true if the grid is empty.
   bool
   isEmpty() const {
      return _size == 0;
   }

   // @}
   //--------------------------------------------------------------------------
   //! \name File I/O.
   // @{

   //! Print the number of records.
   void
   put(std::ostream& out) const {
      out << getSize();
   }

   // @}

protected:

   //! Increment the number of records.
   void
   incrementSize() {
      ++_size;
   }

   //! Decrement the number of records.
   void
   decrementSize() {
      --_size;
   }

   //! Set the number of records.
   void
   setSize(const std::size_t size) {
      _size = size;
   }

   //! Determine an appropriate domain to contain the records.
   /*!
     \param first The first record.
     \param last The last record.
     \param domain The domain that contains the records.

     \note This function assumes that the records are iterators.
   */
   template<class InputIterator>
   void
   computeDomain(InputIterator first, InputIterator last,
                 BBox* domain) const;

};


//! Write to a file stream.
/*! \relates ORQ */
template < std::size_t N, typename _Record, typename _MultiKey, typename _Key,
         typename _MultiKeyAccessor >
inline
std::ostream&
operator<<(std::ostream& out,
           const ORQ<N, _Record, _MultiKey, _Key, _MultiKeyAccessor>& x) {
   x.put(out);
   return out;
}


// Determine an appropriate domain to contain the records.
template < std::size_t N,
         typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor >
template<class InputIterator>
void
ORQ<N, _Record, _MultiKey, _Key, _MultiKeyAccessor>::
computeDomain(InputIterator first, InputIterator last, BBox* domain) const {
   // Check the special case that there are no records.
   if (first == last) {
      // Unit box.
      domain->lower = ext::filled_array<MultiKey>(0);
      domain->upper = ext::filled_array<MultiKey>(1);
      return;
   }
   // The number of records must be non-zero.
   assert(first != last);
   // Compute the domain.
   domain->lower = _multiKeyAccessor(first);
   domain->upper = _multiKeyAccessor(first);
   for (++first; first != last; ++first) {
      domain->add(_multiKeyAccessor(first));
   }
   // Because the upper sides are open, expand the domain.
   const Key Epsilon = std::sqrt(std::numeric_limits<Key>::epsilon());
   Point diagonal = domain->upper;
   diagonal -= domain->lower;
   Point upperCorner = domain->upper;
   for (std::size_t n = 0; n != N; ++n) {
      if (diagonal[n] == 0 && upperCorner[n] == 0) {
         upperCorner[n] = Epsilon;
      }
      else {
         upperCorner[n] += Epsilon *
                           std::max(diagonal[n], std::abs(upperCorner[n]));
      }
   }
   domain->upper = upperCorner;
}


} // namespace geom

#endif
