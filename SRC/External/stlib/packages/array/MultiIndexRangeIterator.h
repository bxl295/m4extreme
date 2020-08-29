// -*- C++ -*-

/*!
  \file array/MultiIndexRangeIterator.h
  \brief An index range iterator.
*/

#if !defined(__array_MultiIndexRangeIterator_h__)
#define __array_MultiIndexRangeIterator_h__

#include "MultiIndexRange.h"

#include <iterator>

namespace array {

//! An index range iterator.
template<std::size_t _Dimension>
class
      MultiIndexRangeIterator {
   //
   // Constants.
   //
public:

   //! The number of dimensions.
   static const std::size_t Dimension = _Dimension;

   //
   // Types.
   //
private:

   typedef MultiIndexTypes<_Dimension> Types;

public:

   //! An index range.
   typedef MultiIndexRange<_Dimension> Range;
   //! The size type.
   typedef typename Types::size_type size_type;
   //! An %array index is a signed integer.
   typedef typename Types::Index Index;
   //! A list of indices.
   typedef typename Types::IndexList IndexList;
   //! A list of sizes.
   typedef typename Types::SizeList SizeList;

   // Iterator types.

   //! Random access iterator category.
   typedef std::random_access_iterator_tag iterator_category;
   //! Value type.
   typedef IndexList value_type;
   //! Pointer difference type.
   typedef typename Types::difference_type difference_type;
   //! Const reference to the value type.
   typedef const value_type& reference;
   //! Const pointer to the value type.
   typedef const value_type* pointer;

   //
   // Member data.
   //
private:

   //! An index list.
   IndexList _indexList;
   //! The rank of the index list.
   Index _rank;
   //! The index range.
   Range _range;

   //--------------------------------------------------------------------------
   //! \name Constructors etc.
   //@{
public:

   //! Return an iterator to the beginning of the index range.
   static
   MultiIndexRangeIterator
   begin(const Range& range);

   //! Return an iterator to the beginning of the index range.
   /*! Since the argument is the range extents, the bases are zero. */
   static
   MultiIndexRangeIterator
   begin(const SizeList& extents);

   //! Return an iterator to the end of the index range.
   static
   MultiIndexRangeIterator
   end(const Range& range);

   //! Return an iterator to the end of the index range.
   /*! Since the argument is the range extents, the bases are zero. */
   static
   MultiIndexRangeIterator
   end(const SizeList& extents);

   //! Copy constructor.
   MultiIndexRangeIterator(const MultiIndexRangeIterator& other);

   //! Assignment operator.
   MultiIndexRangeIterator&
   operator=(const MultiIndexRangeIterator& other);

private:

   //! Default constructor. Uninitialized memory.
   MultiIndexRangeIterator() {
   }

   //@}
   //--------------------------------------------------------------------------
   //! \name Accessors.
   //@{
public:

   //! The index range.
   const Range&
   range() const {
      return _range;
   }

   //@}
   //--------------------------------------------------------------------------
   //! \name Validity.
   //@{
private:

   //! Return true if the iterator is valid.
   /*!
     It's valid if it is in the range [begin(), end()).
   */
   bool
   isValid() const;

   //! Return true if the iterator is at the beginning.
   bool
   isBegin() const;

   //! Return true if the iterator is at the end.
   bool
   isEnd() const;

   //@}
   //--------------------------------------------------------------------------
   //! \name Forward iterator requirements.
   //@{
public:

   reference
   operator*() const {
      return _indexList;
   }

   pointer
   operator->() const {
      return &_indexList;
   }

   //! Pre-increment.
   MultiIndexRangeIterator&
   operator++();

   //! Post-increment.
   /*!
     \warning This function is inefficient. Use pre-increment instead.
   */
   MultiIndexRangeIterator
   operator++(int);

   //@}
   //--------------------------------------------------------------------------
   //! \name Bidirectional iterator requirements.
   //@{
public:

   //! Pre-decrement.
   MultiIndexRangeIterator&
   operator--();

   //! Post-decrement.
   /*!
     \warning This function is inefficient. Use pre-increment instead.
   */
   MultiIndexRangeIterator
   operator--(int);

   //@}
   //--------------------------------------------------------------------------
   //! \name Random access iterator requirements.
   //@{
public:

   //! Iterator indexing.
   /*!
     \warning This function is inefficient.
   */
   value_type
   operator[](const difference_type n) const {
      return *(*this + n);
   }

   MultiIndexRangeIterator&
   operator+=(const difference_type n) {
      _rank += n;
      calculateIndexList();
      return *this;
   }

   MultiIndexRangeIterator
   operator+(const difference_type n) const {
      MultiIndexRangeIterator tmp(*this);
      tmp += n;
      return tmp;
   }

   MultiIndexRangeIterator&
   operator-=(const difference_type n) {
      _rank -= n;
      calculateIndexList();
      return *this;
   }

   MultiIndexRangeIterator
   operator-(const difference_type n) const {
      MultiIndexRangeIterator tmp(*this);
      tmp -= n;
      return tmp;
   }

   Index
   base() const {
      return _rank;
   }

private:

   //! Calculate the index list from the rank.
   void
   calculateIndexList();

   //@}
};

//---------------------------------------------------------------------------
// Equality.

//! Return true if the iterators are equal.
/*! \relates MultiIndexRangeIterator */
template<std::size_t _Dimension>
inline
bool
operator==(const MultiIndexRangeIterator<_Dimension>& x,
           const MultiIndexRangeIterator<_Dimension>& y) {
#ifdef DEBUG_stlib
   // The must be iterators over the same index range.
   assert(x.range() == y.range());
#endif
   return x.base() == y.base();
}

//! Return true if they are not equal.
/*! \relates MultiIndexRangeIterator */
template<std::size_t _Dimension>
inline
bool
operator!=(const MultiIndexRangeIterator<_Dimension>& x,
           const MultiIndexRangeIterator<_Dimension>& y) {
   return !(x == y);
}


//! Return true if the first precedes the second.
/*! \relates MultiIndexRangeIterator */
template<std::size_t _Dimension>
inline
bool
operator<(const MultiIndexRangeIterator<_Dimension>& x,
          const MultiIndexRangeIterator<_Dimension>& y) {
#ifdef DEBUG_stlib
   // The must be iterators over the same index range.
   assert(x.range() == y.range());
#endif
   return x.base() < y.base();
}

//! Return y < x.
/*! \relates MultiIndexRangeIterator */
template<std::size_t _Dimension>
inline
bool
operator>(const MultiIndexRangeIterator<_Dimension>& x,
          const MultiIndexRangeIterator<_Dimension>& y) {
   return y < x;
}

//! Return !(y < x).
/*! \relates MultiIndexRangeIterator */
template<std::size_t _Dimension>
inline
bool
operator<=(const MultiIndexRangeIterator<_Dimension>& x,
           const MultiIndexRangeIterator<_Dimension>& y) {
   return !(y < x);
}

//! Return !(x < y).
/*! \relates MultiIndexRangeIterator */
template<std::size_t _Dimension>
inline
bool
operator>=(const MultiIndexRangeIterator<_Dimension>& x,
           const MultiIndexRangeIterator<_Dimension>& y) {
   return !(x < y);
}

//! Return the difference between the two iterators.
/*! \relates MultiIndexRangeIterator */
template<std::size_t _Dimension>
inline
typename MultiIndexRangeIterator<_Dimension>::difference_type
operator-(const MultiIndexRangeIterator<_Dimension>& x,
          const MultiIndexRangeIterator<_Dimension>& y) {
   return x.base() - y.base();
}

//! Advance the iterator.
/*! \relates MultiIndexRangeIterator */
template<std::size_t _Dimension>
inline
MultiIndexRangeIterator<_Dimension>
operator+(const typename MultiIndexRangeIterator<_Dimension>::difference_type& n,
          const MultiIndexRangeIterator<_Dimension>& x) {
   return x + n;
}

} // namespace array

#define __array_MultiIndexRangeIterator_ipp__
#include "MultiIndexRangeIterator.ipp"
#undef __array_MultiIndexRangeIterator_ipp__

#endif
