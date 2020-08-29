// -*- C++ -*-

/*!
  \file CellArrayStatic.h
  \brief A class for a static cell array in N-D.
*/

#if !defined(__geom_CellArrayStatic_h__)
#define __geom_CellArrayStatic_h__

#include "CellArrayBase.h"

#include "../../array/StaticArrayOfArrays.h"

namespace geom {

//! A static cell array in N-D.
/*!
  A static, dense cell array in N-D.
*/
template<std::size_t N,
         typename _Record,
         typename _MultiKey =
         typename std::iterator_traits<_Record>::value_type,
         typename _Key = typename _MultiKey::value_type,
         typename _MultiKeyAccessor = ads::Dereference<_Record> >
class CellArrayStatic :
   public CellArrayBase<N, _Record, _MultiKey, _Key, _MultiKeyAccessor> {
private:

   //
   // Private types.
   //

   typedef CellArrayBase<N, _Record, _MultiKey, _Key, _MultiKeyAccessor> Base;

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

   //! A single index.
   typedef typename Base::Index Index;
   //! A multi-index.
   typedef typename Base::IndexList IndexList;

private:

   //
   // Data
   //

   // The strides are used for array indexing.
   IndexList _strides;
   // The array of cells.
   array::StaticArrayOfArrays<Record> _cellArray;

protected:

   using Base::_multiKeyAccessor;

private:

   //
   // Not implemented
   //

   // Default constructor not implemented.
   CellArrayStatic();

   // Copy constructor not implemented
   CellArrayStatic(const CellArrayStatic&);

   // Assignment operator not implemented
   CellArrayStatic&
   operator=(const CellArrayStatic&);

   //--------------------------------------------------------------------------
   //! \name Constructors.
   // @{
public:

   //! Construct from the cell size and the domain.
   /*!
     Construct a cell grid given the cell size and the Cartesian domain.

     \param delta The suggested size of a cell.
     \param domain The Cartesian domain that contains the records.
   */
   CellArrayStatic(const Point& delta, const BBox& domain) :
      Base(delta, domain),
      _strides(),
      _cellArray() {
   }

   //! Construct from a range of records.
   /*!
     Construct a cell grid given the cell size, the Cartesian domain
     and a range of records.

     \pre The records must lie in the specified domain.
     \pre There must be a non-zero number of records.

     \param delta The suggested size of a cell.
     \param domain The Cartesian domain that contains the records.
     \param first The first record.
     \param last The last record.

     \note This function assumes that the records are iterators.
   */
   template<class ForwardIterator>
   CellArrayStatic(const Point& delta, const BBox& domain,
                   ForwardIterator first, ForwardIterator last) :
      Base(delta, domain),
      _strides(),
      _cellArray() {
      // Insert the grid elements in the range.
      insert(first, last);
   }

   //! Construct from the cell size and a range of records.
   /*!
     \pre There must be a non-zero number of records.

     \param delta The suggested size of a cell.
     \param first The first record.
     \param last The last record.

     \c first and \c last must be forward iterators because there is one pass
     to compute an appropriate domain and one pass to insert the records.

     \note This function assumes that the records are iterators.
   */
   template<class ForwardIterator>
   CellArrayStatic(const Point& delta, ForwardIterator first,
                   ForwardIterator last) :
      Base(delta, first, last),
      _strides(),
      _cellArray() {
      // Insert the grid elements in the range.
      insert(first, last);
   }

   //! Construct from the number of cells and a range of records.
   /*!
     \pre There must be a non-zero number of records.

     \param numberOfCells the suggested number of cells.
     \param first The first record.
     \param last The last record.

     \c first and \c last must be forward iterators because there is one pass
     to compute an appropriate domain and one pass to insert the records.

     \note This function assumes that the records are iterators.
   */
   template<class ForwardIterator>
   CellArrayStatic(const std::size_t numberOfCells, ForwardIterator first,
                   ForwardIterator last) :
      Base(numberOfCells, first, last),
      _strides(),
      _cellArray() {
      // Insert the grid elements in the range.
      insert(first, last);
   }

   //! Construct from a range of records.
   /*!
     \pre There must be a non-zero number of records.

     \param first The first record.
     \param last The last record.

     \c first and \c last must be forward iterators because there is one pass
     to compute an appropriate domain and one pass to insert the records.
     The number of cells will be approximately the number of records.

     \note This function assumes that the records are iterators.
   */
   template<class ForwardIterator>
   CellArrayStatic(ForwardIterator first, ForwardIterator last) :
      Base(std::distance(first, last), first, last),
      _strides(),
      _cellArray() {
      // Insert the grid elements in the range.
      insert(first, last);
   }

   // @}
   //--------------------------------------------------------------------------
   /*! \name Accesors.
     Functionality inherited from CellArrayBase.
   */
   // @{
public:

   //! Return the number of records.
   using Base::getSize;

   //! Return true if the grid is empty.
   using Base::isEmpty;

   //! Return the domain spanned by the grid.
   using Base::getDomain;

   //! Return the number of cells in each dimension.
   using Base::getExtents;

protected:

   //! Get the inverse cell sizes.
   using Base::getInverseCellSizes;

   // @}
   //--------------------------------------------------------------------------
   //! \name Insert/Erase records.
   // @{
public:

   //! Insert a number of records.
   /*!
     \param first The first record.
     \param last The last record.
     \pre The set of records must be empty.
     \note This function assumes that the records are iterators.
   */
   template<typename ForwardIterator>
   void
   insert(ForwardIterator first, ForwardIterator last) {
      assert(isEmpty());
      // Set the size.
      setSize(std::distance(first, last));
      // Compute the strides.
      computeStrides();
      // Determine where each record should be placed and the sizes for the
      // cells.
      std::vector<std::size_t> indices(getSize());
      std::vector<std::size_t> sizes(product(getExtents()), std::size_t(0));
      Record record = first;
      for (std::size_t i = 0; i != indices.size(); ++i) {
         indices[i] = computeIndex(_multiKeyAccessor(record++));
         ++sizes[indices[i]];
      }
      // Allocate the static array of arrays that represent the cell array.
      _cellArray.rebuild(sizes.begin(), sizes.end());
      // Copy the records into the array.
      std::vector<typename array::StaticArrayOfArrays<Record>::iterator >
      positions(sizes.size());
      for (std::size_t i = 0; i != positions.size(); ++i) {
         positions[i] = _cellArray.begin(i);
      }
      record = first;
      for (std::size_t i = 0; i != indices.size(); ++i) {
         *positions[indices[i]]++ = record++;
      }
   }

   //! Clear all records.
   void
   clear() {
      _cellArray.clear();
      // There are now no records.
      setSize(0);
   }

   //! Rebuild the data structure for a new set of records or following modifications to the record multi-keys.
   template<class ForwardIterator>
   void
   rebuild(ForwardIterator first, ForwardIterator last) {
      clear();
      Base::rebuild(std::distance(first, last), first, last);
      // Insert the grid elements in the range.
      insert(first, last);
   }

   // @}
   //--------------------------------------------------------------------------
   //! \name Memory usage.
   // @{
public:

   //! Return the memory usage.
   std::size_t
   getMemoryUsage() const {
      return _cellArray.getMemoryUsage();
   }

   // @}
   //--------------------------------------------------------------------------
   //! \name Window Queries.
   // @{
public:

   //! Get the records in the window.  Return the # of records inside.
   template<typename OutputIterator>
   std::size_t
   computeWindowQuery(OutputIterator iter, const BBox& window) const;

   // @}
   //--------------------------------------------------------------------------
   //! \name File I/O.
   // @{
public:

   //! Print the records.
   void
   put(std::ostream& out) const;

   // @}

protected:

   //! Increment the number of records.
   using Base::incrementSize;

   //! Set the number of records.
   using Base::setSize;

   //! Convert the multikey to a cell index.
   template<typename AnyMultiKeyType>
   std::size_t
   computeIndex(const AnyMultiKeyType& multiKey) const {
      std::size_t index = 0;
      for (std::size_t n = 0; n != N; ++n) {
         index += _strides[n] * Index((multiKey[n] -
                                       getDomain().lower[n]) *
                                      getInverseCellSizes()[n]);
      }
      return index;
   }

private:

   // Compute the strides from the extents.
   void
   computeStrides() {
      _strides[0] = 1;
      for (std::size_t n = 1; n != N; ++n) {
         _strides[n] = _strides[n-1] * getExtents()[n-1];
      }
   }

};

//
// File I/O
//

//! Write to a file stream.
/*! \relates CellArrayStatic */
template<std::size_t N,
         typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor>
inline
std::ostream&
operator<<(std::ostream& out,
           const CellArrayStatic<N, _Record, _MultiKey, _Key,
           _MultiKeyAccessor>& x) {
   x.put(out);
   return out;
}

} // namespace geom

#define __geom_CellArrayStatic_ipp__
#include "CellArrayStatic.ipp"
#undef __geom_CellArrayStatic_ipp__

#endif
