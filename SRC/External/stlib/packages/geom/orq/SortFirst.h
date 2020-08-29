// -*- C++ -*-

/*!
  \file SortFirst.h
  \brief A class for a vector of records in N-D sorted in each dimension.
*/

#if !defined(__geom_SortFirst_h__)
#define __geom_SortFirst_h__

#include "ORQ.h"

#include "../../ads/algorithm/is_sorted.h"
#include "../../ads/functor/Dereference.h"

#include <iostream>
#include <vector>
#include <algorithm>

namespace geom {


//! A sorted vector of records in N-D.
/*!
  Dimension sorted vectors of records.
*/
template < std::size_t N,
         typename _Record,
         typename _MultiKey =
         typename std::iterator_traits<_Record>::value_type,
         typename _Key = typename _MultiKey::value_type,
         typename _MultiKeyAccessor = ads::Dereference<_Record> >
class SortFirst :
   public ORQ<N, _Record, _MultiKey, _Key, _MultiKeyAccessor> {
private:

   //
   // Private types.
   //

   typedef ORQ<N, _Record, _MultiKey, _Key, _MultiKeyAccessor> Base;
   //! The container of record iterators.
   typedef std::vector<_Record> Container;
   //! A const iterator.
   typedef typename Container::const_iterator ConstIterator;

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

   using Base::_multiKeyAccessor;

protected:

   //
   // Comparison functors.
   //

   //! Less than comparison in a specified coordinate.
   class LessThanCompare :
      public std::binary_function<Record, Record, bool> {
   private:

      MultiKeyAccessor _f;

   public:

      //! Less than comparison in the first coordinate.
      bool
      operator()(const Record x, const Record y) const {
         return _f(x)[0] < _f(y)[0];
      }
   };


   //! Less than comparison in the first coordinate for a record iterator and a multi-key.
   class LessThanCompareValueAndMultiKey :
      public std::binary_function<Record, MultiKey, bool> {
   private:

      MultiKeyAccessor _f;

   public:

      //! Less than comparison in the first coordinate.
      bool
      operator()(const Record value, const MultiKey& multiKey) const {
         return _f(value)[0] < multiKey[0];
      }
   };


   //! Less than comparison in the first coordinate for a multi-key and a record iterator.
   class LessThanCompareMultiKeyAndValue :
      public std::binary_function<MultiKey, Record, bool> {
   private:

      MultiKeyAccessor _f;

   public:

      //! Less than comparison in the first coordinate.
      bool
      operator()(const MultiKey& multiKey, const Record value) const {
         return  multiKey[0] < _f(value)[0];
      }
   };

   //
   // Member data.
   //

private:

   //! Records sorted by the first coordinate.
   Container _sorted;
   mutable LessThanCompare _lessThanCompare;
   mutable LessThanCompareValueAndMultiKey _lessThanCompareValueAndMultiKey;
   mutable LessThanCompareMultiKeyAndValue _lessThanCompareMultiKeyAndValue;

private:

   //
   // Not implemented
   //

   //! Copy constructor not implemented
   SortFirst(const SortFirst&);

   //! Assignment operator not implemented
   SortFirst&
   operator=(const SortFirst&);

public:

   //-------------------------------------------------------------------------
   //! \name Constructors.
   // @{

   //! Default constructor.
   SortFirst() :
      Base(),
      _sorted() {}

   //! Reserve storage for \c size records.
   explicit
   SortFirst(const std::size_t size) :
      Base(),
      _sorted() {
      _sorted.reserve(size);
   }

   //! Construct from a range of records.
   /*!
     \param first the beginning of a range of records.
     \param last the end of a range of records.
   */
   template<class InputIterator>
   SortFirst(InputIterator first, InputIterator last) :
      Base(),
      _sorted() {
      insert(first, last);
      sort();
   }

   // @}
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
   //-------------------------------------------------------------------------
   //! \name Insert records.
   // @{

   //! Insert a single record.
   void
   insert(const Record record) {
      _sorted.push_back(record);
      incrementSize();
   }


   //! Insert a range of records.
   template<class InputIterator>
   void
   insert(InputIterator first, InputIterator last) {
      while (first != last) {
         insert(first);
         ++first;
      }
   }

   // @}
   //-------------------------------------------------------------------------
   //! \name Memory usage.
   // @{

   //! Return the total memory usage.
   std::size_t
   getMemoryUsage() const {
      return sizeof(SortFirst) + getSize() * sizeof(Record);
   }

   // @}
   //-------------------------------------------------------------------------
   //! \name Window Queries.
   // @{

   //! Sort the records by x coordinate.
   void
   sort() {
      std::sort(_sorted.begin(), _sorted.end(), _lessThanCompare);
   }

   //! Get the records in the window.  Return the # of records inside.
   template<typename OutputIterator>
   std::size_t
   computeWindowQuery(OutputIterator iter, const BBox& window) const;

   // @}
   //-------------------------------------------------------------------------
   //! \name File I/O.
   // @{

   //! Print the records.
   void
   put(std::ostream& out) const;

   // @}
   //-------------------------------------------------------------------------
   //! \name Validity check.
   // @{

   //! Check the validity of the data structure.
   bool
   isValid() const;

   //@}

protected:

   //! Increment the number of records.
   using Base::incrementSize;
};


//! Write to a file stream.
/*! \relates SortFirst */
template < std::size_t N,
         typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor >
inline
std::ostream&
operator<<(std::ostream& out,
           const SortFirst<N, _Record, _MultiKey, _Key, _MultiKeyAccessor>& x) {
   x.put(out);
   return out;
}


} // namespace geom

#define __geom_SortFirst_ipp__
#include "SortFirst.ipp"
#undef __geom_SortFirst_ipp__

#endif
