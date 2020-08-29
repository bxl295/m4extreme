// -*- C++ -*-

/*!
  \file array/SimpleMultiIndexExtentsIterator.h
  \brief An index range iterator.
*/

#if !defined(__array_SimpleMultiIndexExtentsIterator_h__)
#define __array_SimpleMultiIndexExtentsIterator_h__

#include "../ext/array.h"

#include <iterator>

namespace array {

//! An index range iterator.
template<std::size_t _Dimension>
class SimpleMultiIndexExtentsIterator {
   //
   // Constants.
   //
public:

   //! The number of dimensions.
   static const std::size_t Dimension = _Dimension;

   //
   // Types.
   //
public:

   //! The size type.
   typedef std::size_t size_type;
   //! An array index is the same as the size type.
   typedef size_type Index;
   //! A list of indices.
   typedef std::tr1::array<size_type, Dimension> IndexList;

   // Iterator types.

   //! Random access iterator category.
   typedef std::random_access_iterator_tag iterator_category;
   //! Value type.
   typedef IndexList value_type;
   //! Pointer difference type.
   typedef std::ptrdiff_t difference_type;
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
   //! The extents of the index range.
   IndexList _extents;

   //--------------------------------------------------------------------------
   /*! \name Constructors etc.
     We use the default copy constructor, assignment operator, and destructor.
     The default constructor is private.
   */
   //@{
public:

   //! Return an iterator to the beginning of the index range.
   static
   SimpleMultiIndexExtentsIterator
   begin(const IndexList& extents);

   //! Return an iterator to the end of the index range.
   static
   SimpleMultiIndexExtentsIterator
   end(const IndexList& extents);

private:

   //! Default constructor. Uninitialized memory.
   SimpleMultiIndexExtentsIterator() {
   }

   //@}
   //--------------------------------------------------------------------------
   //! \name Accessors.
   //@{
public:

   //! The extents of the index range.
   const IndexList&
   extents() const {
      return _extents;
   }

   //@}
   //--------------------------------------------------------------------------
   //! \name Validity.
   //@{
private:

   //! Return true if the iterator is valid.
   /*! It's valid if it is in the range [begin(), end()). */
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
   SimpleMultiIndexExtentsIterator&
   operator++();

   //! Post-increment.
   /*!
     \warning This function is inefficient. Use pre-increment instead.
   */
   SimpleMultiIndexExtentsIterator
   operator++(int);

   //@}
   //--------------------------------------------------------------------------
   //! \name Bidirectional iterator requirements.
   //@{
public:

   //! Pre-decrement.
   SimpleMultiIndexExtentsIterator&
   operator--();

   //! Post-decrement.
   /*!
     \warning This function is inefficient. Use pre-increment instead.
   */
   SimpleMultiIndexExtentsIterator
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

   SimpleMultiIndexExtentsIterator&
   operator+=(const difference_type n) {
      _rank += n;
      calculateIndexList();
      return *this;
   }

   SimpleMultiIndexExtentsIterator
   operator+(const difference_type n) const {
      SimpleMultiIndexExtentsIterator tmp(*this);
      tmp += n;
      return tmp;
   }

   SimpleMultiIndexExtentsIterator&
   operator-=(const difference_type n) {
      _rank -= n;
      calculateIndexList();
      return *this;
   }

   SimpleMultiIndexExtentsIterator
   operator-(const difference_type n) const {
      SimpleMultiIndexExtentsIterator tmp(*this);
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
/*! \relates SimpleMultiIndexExtentsIterator */
template<std::size_t _Dimension>
inline
bool
operator==(const SimpleMultiIndexExtentsIterator<_Dimension>& x,
           const SimpleMultiIndexExtentsIterator<_Dimension>& y) {
#ifdef DEBUG_stlib
   // The must be iterators over the same index range.
   assert(x.extents() == y.extents());
#endif
   return x.base() == y.base();
}

//! Return true if they are not equal.
/*! \relates SimpleMultiIndexExtentsIterator */
template<std::size_t _Dimension>
inline
bool
operator!=(const SimpleMultiIndexExtentsIterator<_Dimension>& x,
           const SimpleMultiIndexExtentsIterator<_Dimension>& y) {
   return !(x == y);
}


//! Return true if the first precedes the second.
/*! \relates SimpleMultiIndexExtentsIterator */
template<std::size_t _Dimension>
inline
bool
operator<(const SimpleMultiIndexExtentsIterator<_Dimension>& x,
          const SimpleMultiIndexExtentsIterator<_Dimension>& y) {
#ifdef DEBUG_stlib
   // The must be iterators over the same index range.
   assert(x.extents() == y.extents());
#endif
   return x.base() < y.base();
}

//! Return y < x.
/*! \relates SimpleMultiIndexExtentsIterator */
template<std::size_t _Dimension>
inline
bool
operator>(const SimpleMultiIndexExtentsIterator<_Dimension>& x,
          const SimpleMultiIndexExtentsIterator<_Dimension>& y) {
   return y < x;
}

//! Return !(y < x).
/*! \relates SimpleMultiIndexExtentsIterator */
template<std::size_t _Dimension>
inline
bool
operator<=(const SimpleMultiIndexExtentsIterator<_Dimension>& x,
           const SimpleMultiIndexExtentsIterator<_Dimension>& y) {
   return !(y < x);
}

//! Return !(x < y).
/*! \relates SimpleMultiIndexExtentsIterator */
template<std::size_t _Dimension>
inline
bool
operator>=(const SimpleMultiIndexExtentsIterator<_Dimension>& x,
           const SimpleMultiIndexExtentsIterator<_Dimension>& y) {
   return !(x < y);
}

//! Return the difference between the two iterators.
/*! \relates SimpleMultiIndexExtentsIterator */
template<std::size_t _Dimension>
inline
typename SimpleMultiIndexExtentsIterator<_Dimension>::difference_type
operator-(const SimpleMultiIndexExtentsIterator<_Dimension>& x,
          const SimpleMultiIndexExtentsIterator<_Dimension>& y) {
   typedef typename SimpleMultiIndexExtentsIterator<_Dimension>::difference_type
      difference_type;
   return difference_type(x.base()) - difference_type(y.base());
}

//! Advance the iterator.
/*! \relates SimpleMultiIndexExtentsIterator */
template<std::size_t _Dimension>
inline
SimpleMultiIndexExtentsIterator<_Dimension>
operator+(const typename
          SimpleMultiIndexExtentsIterator<_Dimension>::difference_type& n,
          const SimpleMultiIndexExtentsIterator<_Dimension>& x) {
   return x + n;
}

} // namespace array

#define __array_SimpleMultiIndexExtentsIterator_ipp__
#include "SimpleMultiIndexExtentsIterator.ipp"
#undef __array_SimpleMultiIndexExtentsIterator_ipp__

#endif
