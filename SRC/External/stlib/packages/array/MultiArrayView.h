// -*- C++ -*-

/*!
  \file array/MultiArrayView.h
  \brief Multi-dimensional view of an %array.
*/

#if !defined(__array_MultiArrayView_h__)
#define __array_MultiArrayView_h__

#include "MultiArrayConstView.h"
#include "MultiIndexRangeIterator.h"

namespace array {

// Forward declaration for assignment operator.
template<typename _T, std::size_t _Dimension>
class MultiArrayConstRef;

//! Multi-dimensional view of an %array.
/*!
  <b>Constructors, etc.</b>

  Since this %array references externally allocated memory, there is
  no default constructor. This class uses the automatically-generated
  copy constructor; the %array data is referenced. Ordinarilly one would create
  an instance of this class with the view() member function. However,
  one can also use the constructors.

  You can construct an %array from a pointer to the data and its index extents.
  Below we make a 2x4x8 %array with index range [0..1]x[0..3]x[0..7]
  \code
  double data[2 * 4 * 8];
  array::MultiArrayConstView<double, 3>::SizeList extents(2, 4, 8)
  array::MultiArrayConstView<double, 3> a(data, extents);
  \endcode

  You can also specify the index bases. Below we make a
  2x4x8 %array with index range [-1..0]x[2..5]x[-3..4]
  \code
  double data[2 * 4 * 8];
  array::MultiArrayConstView<double, 3>::SizeList extents(2, 4, 8)
  array::MultiArrayConstView<double, 3>::IndexList bases(-1, 2, -3)
  array::MultiArrayConstView<double, 3> a(data, extents, bases);
  \endcode

  The copy constructors create shallow copies of the argument, i.e. the
  array data is referenced.
  \code
  array::MultiArray<int, 3> a(extents);
  array::MultiArrayView<int, 3> b(a);
  \endcode
  The argument may be a MultiArray, a MultiArrayRef, or a MultiArrayView.
  The dimension and value type must be the same.

  The assignment operators copy the element values. The argument must have
  the same index ranges as the %array, though they can differ in the value
  type.
  \code
  array::MultiArray<int, 3> a(extents);
  {
    int* data = new int[product(extents)];
    array::MultiArrayRef<int, 3> b(data, extents);
    b = a;
  }
  {
    double* data = new double[product(extents)];
    array::MultiArray<double, 3> c(data, extents);
    c = a;
  }
  \endcode
  The argument may be any of the multidimensional %array types.

  <b>Container Member Functions</b>

  MultiArrayView inherits the following functionality for treating the %array as
  a random access container.

  - MultiArrayBase::empty()
  - MultiArrayBase::size()
  - MultiArrayBase::max_size()

  It defines the following functions.

  - begin()
  - end()
  - rbegin()
  - rend()
  - fill()

  <b>%Array Indexing Member Functions</b>

  MultiArrayView inherits the following %array indexing functionality.

  - MultiArrayBase::extents()
  - MultiArrayBase::bases()
  - MultiArrayBase::setBases()
  - MultiArrayBase::range()
  - MultiArrayBase::storage()
  - MultiArrayBase::strides()
  - MultiArrayBase::offset()

  It defines the following functions.

  - operator()()
  - data()
  - view()

  <b>Free Functions</b>

  - \ref MultiArrayViewAssignmentOperatorsScalar
*/
template<typename _T, std::size_t _Dimension>
class MultiArrayView :
      virtual public MultiArrayConstView<_T, _Dimension> {
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

   typedef MultiArrayTypes<_T, _Dimension> Types;
   typedef MultiArrayConstView<_T, _Dimension> Base;

public:

   // Types for STL compliance.

   //! The element type of the %array.
   typedef typename Types::value_type value_type;
   //! A pointer to an %array element.
   typedef typename Types::pointer pointer;
   //! A pointer to a constant %array element.
   typedef typename Types::const_pointer const_pointer;
   //! A iterator on elements in the %array.
   typedef MultiViewIterator<MultiArrayView, false> iterator;
   //! A reverse iterator on elements in the %array.
   typedef std::reverse_iterator<iterator> reverse_iterator;
   //! A iterator on constant elements in the %array.
   typedef MultiViewIterator<MultiArrayView, true> const_iterator;
   //! A reverse iterator on constant elements in the %array.
   typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
   //! A reference to an %array element.
   typedef typename Types::reference reference;
   //! A reference to a constant %array element.
   typedef typename Types::const_reference const_reference;
   //! The size type.
   typedef typename Types::size_type size_type;
   //! Pointer difference type.
   typedef typename Types::difference_type difference_type;

   // Other types.

   //! The parameter type.
   /*! This is used for passing the value type as an argument. */
   typedef typename Types::Parameter Parameter;
   //! An %array index is a signed integer.
   typedef typename Types::Index Index;
   //! A list of indices.
   typedef typename Types::IndexList IndexList;
   //! A list of sizes.
   typedef typename Types::SizeList SizeList;
   //! The storage order.
   typedef typename Types::Storage Storage;
   //! An index range.
   typedef typename Base::Range Range;
   //! A constant view of this %array.
   typedef typename Base::ConstView ConstView;
   //! A view of this %array.
   typedef MultiArrayView View;

   //
   // Member data.
   //
protected:

   //! Pointer to the data.
   pointer _data;

   //--------------------------------------------------------------------------
   //! \name Constructors etc.
   //@{
public:

   // The default copy constructor is fine.

   //! Construct from a pointer to the memory, the %array extents, the index bases, the storage order, and the strides.
   MultiArrayView(pointer data, const SizeList& extents,
                  const IndexList& bases, const Storage& storage,
                  const IndexList& strides) :
      Base(data, extents, bases, storage, strides),
      _data(data) {
   }

   //! Assignment operator for other %array views.
   /*! \pre The arrays must have the same index range. */
   template<typename _T2>
   MultiArrayView&
   operator=(const MultiArrayConstView<_T2, _Dimension>& other);

   //! Assignment operator for arrays with contiguous memory.
   /*!
     \pre The arrays must have the same index range.
     \note This version is faster than the assignment operator that takes a
     MultiArrayConstView as an argument because arrays with contiguous memory
     have faster iterators.
   */
   template<typename _T2>
   MultiArrayView&
   operator=(const MultiArrayConstRef<_T2, _Dimension>& other);

   //! Assignment operator.
   /*! \pre The arrays must have the same index range. */
   MultiArrayView&
   operator=(const MultiArrayView& other);

   //! Destructor does not deallocate memory.
   virtual
   ~MultiArrayView() {
   }

protected:

   //! Rebuild the data structure.
   void
   rebuild(pointer data, const SizeList& extents, const IndexList& bases,
           const Storage& storage, const IndexList& strides) {
      Base::rebuild(data, extents, bases, storage, strides);
      _data = data;
   }

private:

   //! Default constructor not implemented.
   MultiArrayView();

   //@}
   //--------------------------------------------------------------------------
   //! \name Random Access Container.
   //@{
public:

   using Base::empty;
   using Base::size;
   using Base::max_size;

   //! Return a const iterator to the first value.
   const_iterator
   begin() const {
      return const_iterator::begin(*this);
   }

   //! Return an iterator to the first value.
   iterator
   begin() {
      return iterator::begin(*this);
   }

   //! Return a const iterator to one past the last value.
   const_iterator
   end() const {
      return const_iterator::end(*this);
   }

   //! Return an iterator to one past the last value.
   iterator
   end() {
      return iterator::end(*this);
   }

   //! Return a const reverse iterator to the end of the sequence.
   const_reverse_iterator
   rbegin() const {
      return const_reverse_iterator(end());
   }

   //! Return a reverse iterator to the end of the sequence.
   reverse_iterator
   rbegin() {
      return reverse_iterator(end());
   }

   //! Return a const reverse iterator to the beginning of the sequence.
   const_reverse_iterator
   rend() const {
      return const_reverse_iterator(begin());
   }

   //! Return a reverse iterator to the beginning of the sequence.
   reverse_iterator
   rend() {
      return reverse_iterator(begin());
   }

   //! Fill the %array with the specified value.
   template<typename _T2>
   void
   fill(const _T2& value) {
      std::fill(begin(), end(), value);
   }

   //@}
   //--------------------------------------------------------------------------
   //! \name Array indexing.
   //@{
public:

   using Base::extents;
   using Base::bases;
   using Base::setBases;
   using Base::range;
   using Base::storage;
   using Base::strides;
   using Base::offset;
   using Base::operator();
   using Base::data;
   using Base::view;

   //! Array indexing.
   reference
   operator()(const IndexList& indices) {
#ifdef DEBUG_stlib
      for (size_type n = 0; n != Dimension; ++n) {
         assert(bases()[n] <= indices[n] &&
                indices[n] < bases()[n] + Index(extents()[n]));
      }
#endif
      return _data[arrayIndex(indices)];
   }

   //! Array indexing.
   /*! \note The array dimension must be one in order to use this function. */
   reference
   operator()(const Index i0) {
#ifdef DEBUG_stlib
      assert(bases()[0] <= i0 && i0 < bases()[0] + Index(extents()[0]));
#endif
      return _data[arrayIndex(i0)];
   }

   //! Array indexing.
   /*! \note The array dimension must be two in order to use this function. */
   reference
   operator()(const Index i0, const Index i1) {
#ifdef DEBUG_stlib
      assert(bases()[0] <= i0 && i0 < bases()[0] + Index(extents()[0]) &&
             bases()[1] <= i1 && i1 < bases()[1] + Index(extents()[1]));
#endif
      return _data[arrayIndex(i0, i1)];
   }

   //! Array indexing.
   /*! \note The array dimension must be three in order to use this function. */
   reference
   operator()(const Index i0, const Index i1, const Index i2) {
#ifdef DEBUG_stlib
      assert(bases()[0] <= i0 && i0 < bases()[0] + Index(extents()[0]) &&
             bases()[1] <= i1 && i1 < bases()[1] + Index(extents()[1]) &&
             bases()[2] <= i2 && i2 < bases()[2] + Index(extents()[2]));
#endif
      return _data[arrayIndex(i0, i1, i2)];
   }

   //! Return a pointer to the beginning of the data.
   pointer
   data() {
      return _data;
   }

   //! Make a sub-array view with the index range and optionally the specified bases.
   /*! The bases for the view are the same as that for the index range. */
   View
   view(const Range& range) {
      // Note: Don't use the function call operator because the range 
      // may be empty. In this case array indexing is not valid.
      return View(this->data() + this->arrayIndex(range.bases()),
                  range.extents(), range.bases(),
                  storage(), strides() * range.steps());
   }

   using Base::arrayIndex;

protected:

   //! Set the data pointer.
   void
   setData(pointer data) {
      Base::setData(data);
      _data = data;
   }

   //@}
};

//----------------------------------------------------------------------------
// Free functions.
//----------------------------------------------------------------------------

#if 0
//! For the overlapping elements, set this %array to the unary function of the second.
/*!
  \relates MultiArrayView
  \pre The data for the arrays should not overlap.

  For each overlapping index list <code>i</code>, perform the assignment
  <code>(*a)(i) = f(b(i))</code>.
*/
template < typename _T1, std::size_t _Dimension, typename _T2,
         typename _UnaryFunction >
void
applyUnaryToOverlap
(MultiArrayConstView<_T1, _Dimension>* a,
 const MultiArrayConstView<_T2, _Dimension>& b,
 _UnaryFunction f = std::identity<_T2>());
#endif

//----------------------------------------------------------------------------
//! \defgroup MultiArrayViewAssignmentOperatorsScalar Assignment Operators with Scalar Operand
//@{

//! Array-scalar addition.
/*! \relates MultiArrayView */
template<typename _T, std::size_t _Dimension>
MultiArrayView<_T, _Dimension>&
operator+=(MultiArrayView<_T, _Dimension>& x,
           typename Loki::TypeTraits<_T>::ParameterType value);

//! Array-scalar subtraction.
/*! \relates MultiArrayView */
template<typename _T, std::size_t _Dimension>
MultiArrayView<_T, _Dimension>&
operator-=(MultiArrayView<_T, _Dimension>& x,
           typename Loki::TypeTraits<_T>::ParameterType value);

//! Array-scalar multiplication.
/*! \relates MultiArrayView */
template<typename _T, std::size_t _Dimension>
MultiArrayView<_T, _Dimension>&
operator*=(MultiArrayView<_T, _Dimension>& x,
           typename Loki::TypeTraits<_T>::ParameterType value);

//! Array-scalar division.
/*!
  \relates MultiArrayView
  \note This does not check for division by zero as the value type may not be
  as number type.
*/
template<typename _T, std::size_t _Dimension>
MultiArrayView<_T, _Dimension>&
operator/=(MultiArrayView<_T, _Dimension>& x,
           typename Loki::TypeTraits<_T>::ParameterType value);

//! Array-scalar modulus.
/*!
  \relates MultiArrayView
  \note This does not check for division by zero as the value type may not be
  as number type.
*/
template<typename _T, std::size_t _Dimension>
MultiArrayView<_T, _Dimension>&
operator%=(MultiArrayView<_T, _Dimension>& x,
           typename Loki::TypeTraits<_T>::ParameterType value);

//! Left shift.
/*! \relates MultiArrayView */
template<typename _T, std::size_t _Dimension>
MultiArrayView<_T, _Dimension>&
operator<<=(MultiArrayView<_T, _Dimension>& x, const int offset);

//! Right shift.
/*! \relates MultiArrayView */
template<typename _T, std::size_t _Dimension>
MultiArrayView<_T, _Dimension>&
operator>>=(MultiArrayView<_T, _Dimension>& x, const int offset);

//@}

//---------------------------------------------------------------------------
// File I/O.

// CONTINUE: Add input.

} // namespace array

#define __array_MultiArrayView_ipp__
#include "MultiArrayView.ipp"
#undef __array_MultiArrayView_ipp__

#endif
