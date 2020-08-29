// -*- C++ -*-

/*!
  \file CellSearch.h
  \brief A base class for a cell array with searching in N-D.
*/

#if !defined(__geom_CellSearch_h__)
#define __geom_CellSearch_h__

#include "ORQ.h"

#include "../../array/MultiArray.h"

#include <vector>
#include <algorithm>

namespace geom {

//! Base class for a search structure in the final coordinate.
template < std::size_t N,
         typename _Record,
         typename _MultiKeyAccessor >
class Search :
   public std::vector<_Record> {
private:

   typedef std::vector<_Record> Base;

   // CONTINUE: With the Intel compiler, private members are not accessible
   // in nested classes.
#ifdef __INTEL_COMPILER
public:
#else
private:
#endif

   typedef _MultiKeyAccessor MultiKeyAccessor;

public:

   //
   // Public types.
   //

   //! The record type.
   typedef _Record Record;

   //! The record type.
   typedef typename Base::value_type value_type;
   //! A pointer to value_type.
   typedef typename Base::pointer pointer;
   //! A pointer to const value_type.
   typedef typename Base::const_pointer const_pointer;
   //! An iterator on value_type.
   typedef typename Base::iterator iterator;
   //! An iterator on const value_type.
   typedef typename Base::const_iterator const_iterator;
   //! A reference to value_type.
   typedef typename Base::reference reference;
   //! A reference to const value_type.
   typedef typename Base::const_reference const_reference;
   //! The size type.
   typedef typename Base::size_type size_type;

private:

   //
   // Functors.
   //

   //! Less than comparison in the final coordinate.
   class LessThanCompare :
      public std::binary_function<Record, Record, bool> {
   private:

      MultiKeyAccessor _f;

   public:

      //! Less than comparison in the final coordinate.
      bool
      operator()(const Record x, const Record y) const {
         return _f(x)[N-1] < _f(y)[N-1];
      }
   };

   //
   // Member data.
   //


   //
   // Not implemented.
   //

   // Assignment operator not implemented.
   Search&
   operator=(const Search&);

public:

   //--------------------------------------------------------------------------
   //! \name Constructors.
   // @{

   //! Default constructor.
   Search() :
      Base() {}

   //! Construct and reserve memory for n elements.
   explicit
   Search(const size_type size) :
      Base() {
      reserve(size);
   }

   //! Copy constructor.
   Search(const Search& other) :
      Base(other) {}

   //! Construct from a range.
   template<typename InputIterator>
   Search(InputIterator first, InputIterator last) :
      Base(first, last) {}

   // @}
   //--------------------------------------------------------------------------
   /*! \name Accesors.
     Functionality inherited from std::vector.
   */
   // @{

   //! Return the number of records.
   using Base::size;

   //! Return the beginning of the range of record pointers.
   using Base::begin;

   //! Return the end of the range of record pointers.
   using Base::end;

   // @}
   //--------------------------------------------------------------------------
   /*! \name Manipulators.
     Functionality inherited from std::vector.
   */
   // @{

   // Return the beginning of the range of record pointers.
   //using Base::begin;

   // Return the end of the range of record pointers.
   //using Base::end;

   // @}
   //--------------------------------------------------------------------------
   //! \name Memory usage.
   // @{

   //! Return the memory usage.
   size_type
   getMemoryUsage() const {
      return (sizeof(Search) + size() * sizeof(value_type));
   }

   // @}
   //--------------------------------------------------------------------------
   //! \name Sorting and searching.
   // @{

   //! Sort the records.
   void
   sort() {
      // CONTINUE: This causes warnings with the PGI compiler.
      LessThanCompare compare;
      std::sort(begin(), end(), compare);
   }

   //! Initialize for a set of queries.
   void
   initialize() {}

   //@}
};









//! Base class for a cell array combined with a search structure.
/*!
  A base class for a cell array in N-D.
  This class implements the common functionality of data structures
  which have cell arrays in the first N-1 coordinates and have a search
  data structure in the final coordinate.  It does not store pointers to
  the records.  Instead it has info on the number and size of the cells.
*/
template < std::size_t N,
         typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor,
         typename _Cell >
class CellSearch :
   public ORQ<N, _Record, _MultiKey, _Key, _MultiKeyAccessor> {
private:

   //
   // Private types.
   //

   typedef ORQ<N, _Record, _MultiKey, _Key, _MultiKeyAccessor> Base;

public:

   //
   // Public types.
   //

   //! A pointer to the record type.
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

   //
   // Protected types.
   //

   //! The cell type.
   typedef _Cell Cell;
   //! The cell array type.
   typedef array::MultiArray < Cell, N - 1 > DenseArray;
   //! A single index.
   typedef typename DenseArray::Index Index;
   //! A multi-index into the cell array.
   typedef typename DenseArray::IndexList IndexList;
   //! A type to describe the extents of the cell array.
   typedef typename DenseArray::SizeList SizeList;
   //! The (N-1)-dimensional domain spanned by the cell array.
   typedef geom::BBox<Key, N - 1> Domain;

private:

   //
   // Member data.
   //

   //! The semi-open domain spanned by the grid.
   Domain _domain;

   //! The (N-1)-D array of cells that span the first N-1 coordinates.
   DenseArray _cellArray;

   //! The inverse cell sizes.
   std::tr1::array < Key, N - 1 > _inverseCellSizes;

protected:

   using Base::_multiKeyAccessor;

private:

   //
   // Not implemented
   //

   //! Default constructor not implemented.
   CellSearch();

   //! Copy constructor not implemented
   CellSearch(const CellSearch&);

   //! Assignment operator not implemented
   CellSearch&
   operator=(const CellSearch&);

public:

   //--------------------------------------------------------------------------
   //! \name Constructors.
   //@{

   //! Construct from the size of a cell and a Cartesian domain.
   /*!
     Construct given the cell size and the Cartesian domain
     that contains the records.

     \param delta is the suggested size of a cell.  The final coordinate
     is ignored.
     \param domain is the Cartesian domain spanned by the records.
     The final coordinate is ignored.
   */
   CellSearch(const Point& delta, const BBox& domain);

   //! Construct from the cell size, the Cartision domain and a range of records.
   /*!
     Construct given the cell size, the Cartesian domain and a range of records.

     \param delta is the suggested size of a cell.
     \param domain is the Cartesian domain that contains the records.
     \param first points to the begining of the range of records.
     \param last points to the end of the range.
   */
   template<typename InputIterator>
   CellSearch(const Point& delta, const BBox& domain,
              InputIterator first, InputIterator last);

   //! Construct from the cell size and a range of records.
   /*!
     Construct given the cell size and a range of records. An appropriate
     domain will be computed.

     \param delta is the suggested size of a cell.
     \param first points to the begining of the range of records.
     \param last points to the end of the range.
   */
   template<typename ForwardIterator>
   CellSearch(const Point& delta, ForwardIterator first, ForwardIterator last);

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
   //! \name Insert records.
   //@{

   //! Insert a single record.
   void
   insert(const Record record) {
      getCell(_multiKeyAccessor(record)).push_back(record);
      incrementSize();
   }

   //! Insert a range of records.
   /*!
     The input iterators are to a container of records.
   */
   template<typename InputIterator>
   void
   insert(InputIterator first, InputIterator last) {
      while (first != last) {
         insert(first);
         ++first;
      }
   }

   //@}
   //--------------------------------------------------------------------------
   //! \name Memory usage.
   //@{

   //! Return the memory usage.
   std::size_t
   getMemoryUsage() const;

   //@}
   //--------------------------------------------------------------------------
   //! \name Sorting and searching.
   //@{

   //! Sort the records.
   void
   sort();

   //! Initialize for a set of queries.
   void
   initialize();

   //@}
   //--------------------------------------------------------------------------
   //! \name File I/O
   //@{

   //! Print the data structure.
   void
   put(std::ostream& out) const;

   //@}

protected:

   //! Increment the number of records.
   using Base::incrementSize;

   //! Determine an appropriate domain to contain the records.
   template<class InputIterator>
   void
   computeDomain(InputIterator first, InputIterator last) {
      BBox cartesianDomain;
      Base::computeDomain(first, last, &cartesianDomain);
      setDomain(cartesianDomain);
   }

   //! Convert the multikey to a cell array index.
   template<typename AnyMultiKeyType>
   void
   convertMultiKeyToIndices(const AnyMultiKeyType& multiKey, IndexList* mi)
   const {
      for (std::size_t n = 0; n != N - 1; ++n) {
         (*mi)[n] = Index((multiKey[n] - _domain.lower[n]) *
                          _inverseCellSizes[n]);
      }
   }

   //! Return the extents of the cell array.
   const SizeList&
   getCellArrayExtents() const {
      return _cellArray.extents();
   }

   //! Return a reference to the specified cell.
   Cell&
   getCell(const IndexList& index) {
      return _cellArray(index);
   }

   //! Return a const reference to the specified cell.
   const Cell&
   getCell(const IndexList& index) const {
      return _cellArray(index);
   }

   //! Return a const reference to the search structure that would hold the point.
   template<typename AnyMultiKeyType>
   const Cell&
   getCell(const AnyMultiKeyType& multiKey) const {
      IndexList mi;
      convertMultiKeyToIndices(multiKey, &mi);
      return _cellArray(mi);
   }

   //! Return a reference to the search structure that would hold the point.
   template<typename AnyMultiKeyType>
   Cell&
   getCell(const AnyMultiKeyType& multiKey) {
      IndexList mi;
      convertMultiKeyToIndices(multiKey, &mi);
      return _cellArray(mi);
   }

private:

   //! Compute the array extents and the sizes for the cells.
   void
   computeExtentsAndSizes(const Point& suggestedCellSize);

   //! Set the domain. Ignore the last coordinate.
   void
   setDomain(const BBox& domain);
};


//
// File I/O
//


//! Write to a file stream.
/*! \relates CellSearch */
template < std::size_t N,
         typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor,
         typename _Cell >
inline
std::ostream&
operator<<
(std::ostream& out,
 const CellSearch<N, _Record, _MultiKey, _Key, _MultiKeyAccessor, _Cell>& x) {
   x.put(out);
   return out;
}


} // namespace geom

#define __geom_CellSearch_ipp__
#include "CellSearch.ipp"
#undef __geom_CellSearch_ipp__

#endif
