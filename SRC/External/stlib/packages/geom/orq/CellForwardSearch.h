// -*- C++ -*-

/*!
  \file CellForwardSearch.h
  \brief A class for a cell array coupled with a forward search.
*/

#if !defined(__geom_CellForwardSearch_h__)
#define __geom_CellForwardSearch_h__

#include "CellSearch.h"

#include "../../ads/functor/Dereference.h"

namespace geom {


//! Data structure for performing forward searching on records.
template < std::size_t N,
         typename _Record,
         typename _Key,
         typename _MultiKeyAccessor >
class ForwardSearch :
   public Search<N, _Record, _MultiKeyAccessor> {
private:

   typedef Search<N, _Record, _MultiKeyAccessor> Base;
   typedef _MultiKeyAccessor MultiKeyAccessor;

public:

   //
   // Public types.
   //

   //! The record type.
   typedef _Record Record;
   //! The key type.
   typedef _Key Key;

   //! The record type.
   typedef typename Base::value_type value_type;
   //! A pointer to the value type.
   typedef typename Base::pointer pointer;
   //! A const pointer to the value type.
   typedef typename Base::const_pointer const_pointer;
   //! An iterator to the value type.
   typedef typename Base::iterator iterator;
   //! A const iterator to the value type.
   typedef typename Base::const_iterator const_iterator;
   //! A reference to the value type.
   typedef typename Base::reference reference;
   //! A const reference to the value type.
   typedef typename Base::const_reference const_reference;
   //! The size type.
   typedef typename Base::size_type size_type;

private:

   //
   // Member data
   //

   MultiKeyAccessor _multiKeyAccessor;
   mutable iterator _current;

   //
   // Not implemented.
   //

   //! Assignment operator not implemented.
   ForwardSearch&
   operator=(const ForwardSearch&);

public:

   //--------------------------------------------------------------------------
   //! \name Constructors.
   // @{

   //! Default constructor.
   ForwardSearch() :
      Base(),
      _multiKeyAccessor(),
      _current(0) {}

   //! Construct and reserve memory for n elements.
   explicit
   ForwardSearch(const size_type size) :
      Base(size),
      _multiKeyAccessor(),
      _current(0) {}

   //! Copy constructor.
   ForwardSearch(const ForwardSearch& other) :
      Base(other),
      _multiKeyAccessor(other._multiKeyAccessor),
      _current(other._current) {}

   //! Construct from a range.
   template<typename InputIterator>
   ForwardSearch(InputIterator first, InputIterator last) :
      Base(first, last),
      _multiKeyAccessor(),
      _current(0) {}

   // @}
   //--------------------------------------------------------------------------
   /*! \name Accessors.
     Functionality inherited from std::vector.
   */
   // @{

   //! Return the number of records.
   using Base::size;

   //! Return the beginning of the range of records.
   using Base::begin;

   //! Return the end of the range of records.
   using Base::end;

   // @}
   //--------------------------------------------------------------------------
   /*! \name Manipulators.
     Functionality inherited from std::vector.
   */
   // @{

   // Return the beginning of the range of records.
   //using Base::begin;

   // Return the end of the range of records.
   //using Base::end;

   // @}
   //--------------------------------------------------------------------------
   //! \name Memory usage.
   // @{

   //! Return the memory usage.
   size_type
   getMemoryUsage() const {
      return (sizeof(ForwardSearch) + size() * sizeof(value_type));
   }

   // @}
   //--------------------------------------------------------------------------
   //! \name Searching.
   // @{

   //! Initialize for a set of queries.
   void
   initialize() {
      _current = begin();
   }

   //! Do a forward search to find the record.
   /*!
     \param x is a key in the final coordinate.
   */
   iterator
   search(const Key x) {
      // Move the current pointer to the specified value.
      for (; _current != end() && _multiKeyAccessor(*_current)[N-1] < x;
            ++_current)
         ;
      return _current;
   }

   //! Do a forward search to find the record.
   /*!
     \param x is a key in the final coordinate.
   */
   const_iterator
   search(const Key x) const {
      // Move the current pointer to the specified value.
      for (; _current != end() && _multiKeyAccessor(*_current)[N-1] < x;
            ++_current)
         ;
      return _current;
   }
};









//! A cell array in combined with a forward search in the final coordinate.
/*!
  The cell array spans N-1 coordinates.  Record access
  is accomplished with array indexing in these coordinates and a
  forward search of a sorted vector in the final coordinate.
*/
template < std::size_t N,
         typename _Record,
         typename _MultiKey =
         typename std::iterator_traits<_Record>::value_type,
         typename _Key = typename _MultiKey::value_type,
         typename _MultiKeyAccessor = ads::Dereference<_Record> >
class CellForwardSearch :
   public CellSearch < N, _Record, _MultiKey, _Key, _MultiKeyAccessor,
      ForwardSearch<N, _Record, _Key, _MultiKeyAccessor> > {
private:

   //
   // Private types.
   //

   typedef CellSearch < N, _Record, _MultiKey, _Key, _MultiKeyAccessor,
           ForwardSearch<N, _Record, _Key, _MultiKeyAccessor> > Base;

protected:

   //
   // Protected types.
   //

   //! The cell type.
   typedef typename Base::Cell Cell;
   //! A single index.
   typedef typename Base::Index Index;
   //! A multi-index into the cell array.
   typedef typename Base::IndexList IndexList;

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

   //! A Cartesian point.
   typedef typename Base::Point Point;
   //! Bounding box.
   typedef typename Base::BBox BBox;

protected:

   using Base::_multiKeyAccessor;

private:

   //
   // Not implemented
   //

   //! Default constructor not implemented.
   CellForwardSearch();

   //! Copy constructor not implemented
   CellForwardSearch(const CellForwardSearch&);

   //! Assignment operator not implemented
   CellForwardSearch&
   operator=(const CellForwardSearch&);

public:

   //--------------------------------------------------------------------------
   //! \name Constructors.
   //@{

   //! Construct from the size of a cell and a Cartesian domain.
   /*!
     Construct given the cell size and the Cartesian domain that
     contains the records.

     \param delta the suggested size of a cell.
     \param domain the Cartesian domain that contains the records.
   */
   CellForwardSearch(const Point& delta, const BBox& domain) :
      Base(delta, domain) {}

   //! Construct from the Cartesian domain and a range of records.
   /*!
     Construct given the cell size, the Cartesian domain that
     contains the records and a range of records.

     \param delta the suggested size of a cell.
     \param domain the Cartesian domain that contains the records.
     \param first points to the begining of the range of records.
     \param last points to the end of the semi-open range.
   */
   template<typename InputIterator>
   CellForwardSearch(const Point& delta, const BBox& domain,
                     InputIterator first, InputIterator last) :
      Base(delta, domain, first, last) {}

   //! Construct from a range of records.
   /*!
     Construct given the cell size a range of records. Compute an appropriate
     domain.

     \param delta the suggested size of a cell.
     \param first points to the begining of the range of records.
     \param last points to the end of the semi-open range.
   */
   template<typename ForwardIterator>
   CellForwardSearch(const Point& delta,
                     ForwardIterator first, ForwardIterator last) :
      Base(delta, first, last) {}

   //@}
   //--------------------------------------------------------------------------
   /*! \name Accesors.
     Functionality inherited from ORQ.
    */
   // @{

   //! Return the number of records.
   using Base::getSize;

   //! Return true if the grid is empty.
   using Base::isEmpty;

   // @}
   //--------------------------------------------------------------------------
   /*! \name Insert records.
     Functionality inherited from CellSearch.
   */
   //@{

   //! Insert a single record or a range of records.
   using Base::insert;

   //@}
   //--------------------------------------------------------------------------
   /*! \name Memory usage.
     Functionality inherited from CellSearch.
   */
   //@{

   //! Return the memory usage.
   using Base::getMemoryUsage;

   //@}
   //--------------------------------------------------------------------------
   //! \name Window queries.
   //@{

   //! Get the records in the window.  Return the # of records inside.
   template<typename OutputIterator>
   std::size_t
   computeWindowQuery(OutputIterator iter, const BBox& window) const;

   //@}
   //--------------------------------------------------------------------------
   /*! \name Sorting.
     Functionality inherited from CellSearch.
   */
   //@{

   //! Sort the records.
   using Base::sort;

   //! Initialize for a set of queries.
   using Base::initialize;

   //@}
   //--------------------------------------------------------------------------
   /*! \name File I/O
     Functionality inherited from CellSearch.
   */
   //@{

   //! Print the data structure.
   using Base::put;

   //@}

private:

   //! Convert the multikey to a cell array index.
   using Base::convertMultiKeyToIndices;

   //! Return the extents of the cell array.
   using Base::getCellArrayExtents;

   //! Return the specified cell.
   using Base::getCell;
};


//
// File I/O
//


//! Write to a file stream.
/*! \relates CellForwardSearch */
template < std::size_t N,
         typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor >
inline
std::ostream&
operator<<
(std::ostream& out,
 const CellForwardSearch<N, _Record, _MultiKey, _Key, _MultiKeyAccessor>& x) {
   x.put(out);
   return out;
}


} // namespace geom

#define __geom_CellForwardSearch_ipp__
#include "CellForwardSearch.ipp"
#undef __geom_CellForwardSearch_ipp__

#endif
