// -*- C++ -*-

/*!
  \file array/EquilateralArrayImp.h
  \brief A multi-array that has equal extents in each dimension.
*/

#if !defined(__array_EquilateralArrayImp_h__)
#define __array_EquilateralArrayImp_h__

#include "../ext/array.h"
#include "../numerical/constants/Exponentiation.h"
#include "../loki/static_check.h"
#include "../loki/TypeTraits.h"

namespace array {

//! A multi-array that has equal extents in each dimension.
/*!
  \param _T The value type.
  \param _D The dimension
  \param _N The extent in each dimension.
  \param _Base The base class that implements the std::tr1::array interface.
*/
template<typename _T, std::size_t _D, std::size_t _N, typename _Base>
class EquilateralArrayImp :
   public _Base {
   //
   // Constants.
   //
public:

   //! The total number of elements.
   static const std::size_t Size =
      numerical::Exponentiation<std::size_t, _N, _D>::Result;

   //
   // Types.
   //
private:

   typedef _Base Base;

public:

   //! An array index is the same as the size type.
   typedef typename Base::size_type Index;
   //! A list of indices.
   typedef std::tr1::array<Index, _D> IndexList;

   //--------------------------------------------------------------------------
   /*! \name Constructors etc.
     We use the default constructor, copy constructor, assignment operator, and
     destructor. */
   // @{
public:

   // @}
   //--------------------------------------------------------------------------
   //! \name Accessors.
   // @{
public:

   // GCC 4.0 does not properly implement the data() member for array.
#if __GNUC__ == 4 && __GNUC_MINOR__ == 0
   //! Return a const iterator to the data.
   typename Base::const_iterator
   data() const {
      return Base::begin();
   }
   //! Return an iterator to the data.
   typename Base::iterator
   data() {
      return Base::begin();
   }
#endif

   //! Return the extents.
   std::tr1::array<typename Base::size_type, _D>
   extents() const {
      return ext::filled_array<std::tr1::array<typename Base::size_type, _D> >
         (_N);
   }

   //! Return a const reference to the specified element in the 0-D array.
   typename Base::const_reference
   operator()() const {
      LOKI_STATIC_CHECK(_D == 0, BadDimension);
      return Base::operator[](0);
   }

   //! Return a const reference to the specified element in the 0-D array.
   template<typename _Integer>
   typename Base::const_reference
   operator()(const std::tr1::array<_Integer, 0>& /*i*/) const {
      return operator()();
   }

   //! Return a const reference to the specified element in the 1-D array.
   typename Base::const_reference
   operator()(const typename Base::size_type i0) const {
      LOKI_STATIC_CHECK(_D == 1, BadDimension);
      return operator[](i0);
   }

   //! Return a const reference to the specified element in the 1-D array.
   template<typename _Integer>
   typename Base::const_reference
   operator()(const std::tr1::array<_Integer, 1>& i) const {
      return operator()(i[0]);
   }

   //! Return a const reference to the specified element in the 2-D array.
   typename Base::const_reference
   operator()(const typename Base::size_type i0,
              const typename Base::size_type i1) const {
      LOKI_STATIC_CHECK(_D == 2, BadDimension);
      return operator[](i0 + i1 * _N);
   }

   //! Return a const reference to the specified element in the 2-D array.
   template<typename _Integer>
   typename Base::const_reference
   operator()(const std::tr1::array<_Integer, 2>& i) const {
      return operator()(i[0], i[1]);
   }

   //! Return a const reference to the specified element in the 3-D array.
   typename Base::const_reference
   operator()(const typename Base::size_type i0,
              const typename Base::size_type i1,
              const typename Base::size_type i2)
   const {
      LOKI_STATIC_CHECK(_D == 3, BadDimension);
      return operator[](i0 + i1 * _N + i2 * _N * _N);
   }

   //! Return a const reference to the specified element in the 3-D array.
   template<typename _Integer>
   typename Base::const_reference
   operator()(const std::tr1::array<_Integer, 3>& i) const {
      return operator()(i[0], i[1], i[2]);
   }

   //! Return a const reference to the specified element in the 4-D array.
   typename Base::const_reference
   operator()(const typename Base::size_type i0,
              const typename Base::size_type i1,
              const typename Base::size_type i2,
              const typename Base::size_type i3) const {
      LOKI_STATIC_CHECK(_D == 4, BadDimension);
      return operator[](i0 + i1 * _N + i2 * _N * _N + i3 * _N * _N * _N);
   }

   //! Return a const reference to the specified element in the 4-D array.
   template<typename _Integer>
   typename Base::const_reference
   operator()(const std::tr1::array<_Integer, 4>& i) const {
      return operator()(i[0], i[1], i[2], i[3]);
   }

   //! Return a const reference to the specified element in the 5-D array.
   typename Base::const_reference
   operator()(const typename Base::size_type i0,
              const typename Base::size_type i1,
              const typename Base::size_type i2,
              const typename Base::size_type i3,
              const typename Base::size_type i4) const {
      LOKI_STATIC_CHECK(_D == 5, BadDimension);
      return operator[](i0 + i1 * _N + i2 * _N * _N + i3 * _N * _N * _N +
                        i4 * _N * _N * _N * _N);
   }

   //! Return a const reference to the specified element in the 5-D array.
   template<typename _Integer>
   typename Base::const_reference
   operator()(const std::tr1::array<_Integer, 5>& i) const {
      return operator()(i[0], i[1], i[2], i[3], i[4]);
   }

   //! Return a const reference to the specified element in the 6-D array.
   typename Base::const_reference
   operator()(const typename Base::size_type i0,
              const typename Base::size_type i1,
              const typename Base::size_type i2,
              const typename Base::size_type i3,
              const typename Base::size_type i4,
              const typename Base::size_type i5) const {
      LOKI_STATIC_CHECK(_D == 6, BadDimension);
      return operator[](i0 + i1 * _N + i2 * _N * _N + i3 * _N * _N * _N +
                        i4 * _N * _N * _N * _N + i5 * _N * _N * _N * _N * _N);
   }

   //! Return a const reference to the specified element in the 6-D array.
   template<typename _Integer>
   typename Base::const_reference
   operator()(const std::tr1::array<_Integer, 6>& i) const {
      return operator()(i[0], i[1], i[2], i[3], i[4], i[5]);
   }

   //! Return a const reference to the specified element in the 7-D array.
   typename Base::const_reference
   operator()(const typename Base::size_type i0,
              const typename Base::size_type i1,
              const typename Base::size_type i2,
              const typename Base::size_type i3,
              const typename Base::size_type i4,
              const typename Base::size_type i5,
              const typename Base::size_type i6) const {
      LOKI_STATIC_CHECK(_D == 7, BadDimension);
      return operator[](i0 + i1 * _N + i2 * _N * _N + i3 * _N * _N * _N +
                        i4 * _N * _N * _N * _N + i5 * _N * _N * _N * _N * _N +
                        i6 * _N * _N * _N * _N * _N * _N);
   }

   //! Return a const reference to the specified element in the 7-D array.
   template<typename _Integer>
   typename Base::const_reference
   operator()(const std::tr1::array<_Integer, 7>& i) const {
      return operator()(i[0], i[1], i[2], i[3], i[4], i[5], i[6]);
   }

   //! Return a const reference to the specified element in the 8-D array.
   typename Base::const_reference
   operator()(const typename Base::size_type i0,
              const typename Base::size_type i1,
              const typename Base::size_type i2,
              const typename Base::size_type i3,
              const typename Base::size_type i4,
              const typename Base::size_type i5,
              const typename Base::size_type i6,
              const typename Base::size_type i7) const {
      LOKI_STATIC_CHECK(_D == 8, BadDimension);
      return operator[](i0 + i1 * _N + i2 * _N * _N + i3 * _N * _N * _N +
                        i4 * _N * _N * _N * _N + i5 * _N * _N * _N * _N * _N +
                        i6 * _N * _N * _N * _N * _N * _N +
                        i7 * _N * _N * _N * _N * _N * _N * _N);
   }

   //! Return a const reference to the specified element in the 8-D array.
   template<typename _Integer>
   typename Base::const_reference
   operator()(const std::tr1::array<_Integer, 8>& i) const {
      return operator()(i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7]);
   }

   //! Return a const reference to the specified element in the 9-D array.
   typename Base::const_reference
   operator()(const typename Base::size_type i0,
              const typename Base::size_type i1,
              const typename Base::size_type i2,
              const typename Base::size_type i3,
              const typename Base::size_type i4,
              const typename Base::size_type i5,
              const typename Base::size_type i6,
              const typename Base::size_type i7,
              const typename Base::size_type i8) const {
      LOKI_STATIC_CHECK(_D == 9, BadDimension);
      return operator[](i0 + i1 * _N + i2 * _N * _N + i3 * _N * _N * _N +
                        i4 * _N * _N * _N * _N + i5 * _N * _N * _N * _N * _N +
                        i6 * _N * _N * _N * _N * _N * _N +
                        i7 * _N * _N * _N * _N * _N * _N * _N +
                        i8 * _N * _N * _N * _N * _N * _N * _N * _N);
   }

   //! Return a const reference to the specified element in the 9-D array.
   template<typename _Integer>
   typename Base::const_reference
   operator()(const std::tr1::array<_Integer, 9>& i) const {
      return operator()(i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8]);
   }

   //! Return a const reference to the specified element in the 10-D array.
   typename Base::const_reference
   operator()(const typename Base::size_type i0,
              const typename Base::size_type i1,
              const typename Base::size_type i2,
              const typename Base::size_type i3,
              const typename Base::size_type i4,
              const typename Base::size_type i5,
              const typename Base::size_type i6,
              const typename Base::size_type i7,
              const typename Base::size_type i8,
              const typename Base::size_type i9) const {
      LOKI_STATIC_CHECK(_D == 10, BadDimension);
      return operator[](i0 + i1 * _N + i2 * _N * _N + i3 * _N * _N * _N +
                        i4 * _N * _N * _N * _N + i5 * _N * _N * _N * _N * _N +
                        i6 * _N * _N * _N * _N * _N * _N +
                        i7 * _N * _N * _N * _N * _N * _N * _N +
                        i8 * _N * _N * _N * _N * _N * _N * _N * _N +
                        i9 * _N * _N * _N * _N * _N * _N * _N * _N * _N);
   }

   //! Return a const reference to the specified element in the 10-D array.
   template<typename _Integer>
   typename Base::const_reference
   operator()(const std::tr1::array<_Integer, 10>& i) const {
      return operator()(i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8],
                        i[9]);
   }

   // @}
   //--------------------------------------------------------------------------
   //! \name Manipulators.
   // @{
public:

   //! Return a reference to the specified element in the 0-D array.
   typename Base::reference
   operator()() {
      LOKI_STATIC_CHECK(_D == 0, BadDimension);
      return Base::operator[](0);
   }

   //! Return a reference to the specified element in the 0-D array.
   template<typename _Integer>
   typename Base::reference
   operator()(const std::tr1::array<_Integer, 0>& /*i*/) {
      return operator()();
   }

   //! Return a reference to the specified element in the 1-D array.
   typename Base::reference
   operator()(const typename Base::size_type i0) {
      LOKI_STATIC_CHECK(_D == 1, BadDimension);
      return operator[](i0);
   }

   //! Return a reference to the specified element in the 1-D array.
   template<typename _Integer>
   typename Base::reference
   operator()(const std::tr1::array<_Integer, 1>& i) {
      return operator()(i[0]);
   }

   //! Return a reference to the specified element in the 2-D array.
   typename Base::reference
   operator()(const typename Base::size_type i0,
              const typename Base::size_type i1) {
      LOKI_STATIC_CHECK(_D == 2, BadDimension);
      return operator[](i0 + i1 * _N);
   }

   //! Return a reference to the specified element in the 2-D array.
   template<typename _Integer>
   typename Base::reference
   operator()(const std::tr1::array<_Integer, 2>& i) {
      return operator()(i[0], i[1]);
   }

   //! Return a reference to the specified element in the 3-D array.
   typename Base::reference
   operator()(const typename Base::size_type i0,
              const typename Base::size_type i1,
              const typename Base::size_type i2) {
      LOKI_STATIC_CHECK(_D == 3, BadDimension);
      return operator[](i0 + i1 * _N + i2 * _N * _N);
   }

   //! Return a reference to the specified element in the 3-D array.
   template<typename _Integer>
   typename Base::reference
   operator()(const std::tr1::array<_Integer, 3>& i) {
      return operator()(i[0], i[1], i[2]);
   }

   //! Return a reference to the specified element in the 4-D array.
   typename Base::reference
   operator()(const typename Base::size_type i0,
              const typename Base::size_type i1,
              const typename Base::size_type i2,
              const typename Base::size_type i3) {
      LOKI_STATIC_CHECK(_D == 4, BadDimension);
      return operator[](i0 + i1 * _N + i2 * _N * _N + i3 * _N * _N * _N);
   }

   //! Return a reference to the specified element in the 4-D array.
   template<typename _Integer>
   typename Base::reference
   operator()(const std::tr1::array<_Integer, 4>& i) {
      return operator()(i[0], i[1], i[2], i[3]);
   }

   //! Return a reference to the specified element in the 5-D array.
   typename Base::reference
   operator()(const typename Base::size_type i0,
              const typename Base::size_type i1,
              const typename Base::size_type i2,
              const typename Base::size_type i3,
              const typename Base::size_type i4) {
      LOKI_STATIC_CHECK(_D == 5, BadDimension);
      return operator[](i0 + i1 * _N + i2 * _N * _N + i3 * _N * _N * _N +
                        i4 * _N * _N * _N * _N);
   }

   //! Return a reference to the specified element in the 5-D array.
   template<typename _Integer>
   typename Base::reference
   operator()(const std::tr1::array<_Integer, 5>& i) {
      return operator()(i[0], i[1], i[2], i[3], i[4]);
   }

   //! Return a reference to the specified element in the 6-D array.
   typename Base::reference
   operator()(const typename Base::size_type i0,
              const typename Base::size_type i1,
              const typename Base::size_type i2,
              const typename Base::size_type i3,
              const typename Base::size_type i4,
              const typename Base::size_type i5) {
      LOKI_STATIC_CHECK(_D == 6, BadDimension);
      return operator[](i0 + i1 * _N + i2 * _N * _N + i3 * _N * _N * _N +
                        i4 * _N * _N * _N * _N + i5 * _N * _N * _N * _N * _N);
   }

   //! Return a reference to the specified element in the 6-D array.
   template<typename _Integer>
   typename Base::reference
   operator()(const std::tr1::array<_Integer, 6>& i) {
      return operator()(i[0], i[1], i[2], i[3], i[4], i[5]);
   }

   //! Return a reference to the specified element in the 7-D array.
   typename Base::reference
   operator()(const typename Base::size_type i0,
              const typename Base::size_type i1,
              const typename Base::size_type i2,
              const typename Base::size_type i3,
              const typename Base::size_type i4,
              const typename Base::size_type i5,
              const typename Base::size_type i6) {
      LOKI_STATIC_CHECK(_D == 7, BadDimension);
      return operator[](i0 + i1 * _N + i2 * _N * _N + i3 * _N * _N * _N +
                        i4 * _N * _N * _N * _N + i5 * _N * _N * _N * _N * _N +
                        i6 * _N * _N * _N * _N * _N * _N);
   }

   //! Return a reference to the specified element in the 7-D array.
   template<typename _Integer>
   typename Base::reference
   operator()(const std::tr1::array<_Integer, 7>& i) {
      return operator()(i[0], i[1], i[2], i[3], i[4], i[5], i[6]);
   }

   //! Return a reference to the specified element in the 8-D array.
   typename Base::reference
   operator()(const typename Base::size_type i0,
              const typename Base::size_type i1,
              const typename Base::size_type i2,
              const typename Base::size_type i3,
              const typename Base::size_type i4,
              const typename Base::size_type i5,
              const typename Base::size_type i6,
              const typename Base::size_type i7) {
      LOKI_STATIC_CHECK(_D == 8, BadDimension);
      return operator[](i0 + i1 * _N + i2 * _N * _N + i3 * _N * _N * _N +
                        i4 * _N * _N * _N * _N + i5 * _N * _N * _N * _N * _N +
                        i6 * _N * _N * _N * _N * _N * _N +
                        i7 * _N * _N * _N * _N * _N * _N * _N);
   }

   //! Return a reference to the specified element in the 8-D array.
   template<typename _Integer>
   typename Base::reference
   operator()(const std::tr1::array<_Integer, 8>& i) {
      return operator()(i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7]);
   }

   //! Return a reference to the specified element in the 9-D array.
   typename Base::reference
   operator()(const typename Base::size_type i0,
              const typename Base::size_type i1,
              const typename Base::size_type i2,
              const typename Base::size_type i3,
              const typename Base::size_type i4,
              const typename Base::size_type i5,
              const typename Base::size_type i6,
              const typename Base::size_type i7,
              const typename Base::size_type i8) {
      LOKI_STATIC_CHECK(_D == 9, BadDimension);
      return operator[](i0 + i1 * _N + i2 * _N * _N + i3 * _N * _N * _N +
                        i4 * _N * _N * _N * _N + i5 * _N * _N * _N * _N * _N +
                        i6 * _N * _N * _N * _N * _N * _N +
                        i7 * _N * _N * _N * _N * _N * _N * _N +
                        i8 * _N * _N * _N * _N * _N * _N * _N * _N);
   }

   //! Return a reference to the specified element in the 9-D array.
   template<typename _Integer>
   typename Base::reference
   operator()(const std::tr1::array<_Integer, 9>& i) {
      return operator()(i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8]);
   }

   //! Return a reference to the specified element in the 10-D array.
   typename Base::reference
   operator()(const typename Base::size_type i0,
              const typename Base::size_type i1,
              const typename Base::size_type i2,
              const typename Base::size_type i3,
              const typename Base::size_type i4,
              const typename Base::size_type i5,
              const typename Base::size_type i6,
              const typename Base::size_type i7,
              const typename Base::size_type i8,
              const typename Base::size_type i9) {
      LOKI_STATIC_CHECK(_D == 10, BadDimension);
      return operator[](i0 + i1 * _N + i2 * _N * _N + i3 * _N * _N * _N +
                        i4 * _N * _N * _N * _N + i5 * _N * _N * _N * _N * _N +
                        i6 * _N * _N * _N * _N * _N * _N +
                        i7 * _N * _N * _N * _N * _N * _N * _N +
                        i8 * _N * _N * _N * _N * _N * _N * _N * _N +
                        i9 * _N * _N * _N * _N * _N * _N * _N * _N * _N);
   }

   //! Return a reference to the specified element in the 10-D array.
   template<typename _Integer>
   typename Base::reference
   operator()(const std::tr1::array<_Integer, 10>& i) {
      return operator()(i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8],
                        i[9]);
   }

   // @}
};


//----------------------------------------------------------------------------
//! \defgroup EquilateralArrayAssignmentOperatorsScalar Assignment Operators with Scalar Operand
//@{


//! Array-scalar addition.
/*! \relates EquilateralArrayImp */
template<typename _T, std::size_t _D, std::size_t _N, typename _B>
inline
EquilateralArrayImp<_T, _D, _N, _B>&
operator+=(EquilateralArrayImp<_T, _D, _N, _B>& x,
           typename Loki::TypeTraits<_T>::ParameterType value) {
   typedef typename EquilateralArrayImp<_T, _D, _N, _B>::iterator iterator;
   for (iterator i = x.begin(); i != x.end(); ++i) {
      *i += value;
   }
   return x;
}


//! Array-scalar subtraction.
/*! \relates EquilateralArrayImp */
template<typename _T, std::size_t _D, std::size_t _N, typename _B>
inline
EquilateralArrayImp<_T, _D, _N, _B>&
operator-=(EquilateralArrayImp<_T, _D, _N, _B>& x,
           typename Loki::TypeTraits<_T>::ParameterType value) {
   typedef typename EquilateralArrayImp<_T, _D, _N, _B>::iterator iterator;
   for (iterator i = x.begin(); i != x.end(); ++i) {
      *i -= value;
   }
   return x;
}


//! Array-scalar multiplication.
/*! \relates EquilateralArrayImp */
template<typename _T, std::size_t _D, std::size_t _N, typename _B>
inline
EquilateralArrayImp<_T, _D, _N, _B>&
operator*=(EquilateralArrayImp<_T, _D, _N, _B>& x,
           typename Loki::TypeTraits<_T>::ParameterType value) {
   typedef typename EquilateralArrayImp<_T, _D, _N, _B>::iterator iterator;
   for (iterator i = x.begin(); i != x.end(); ++i) {
      *i *= value;
   }
   return x;
}


//! Array-scalar division.
/*!
  \relates EquilateralArrayImp
  \note This does not check for division by zero as the value type may not be
  as number type.
*/
template<typename _T, std::size_t _D, std::size_t _N, typename _B>
inline
EquilateralArrayImp<_T, _D, _N, _B>&
operator/=(EquilateralArrayImp<_T, _D, _N, _B>& x,
           typename Loki::TypeTraits<_T>::ParameterType value) {
   typedef typename EquilateralArrayImp<_T, _D, _N, _B>::iterator iterator;
   for (iterator i = x.begin(); i != x.end(); ++i) {
      *i /= value;
   }
   return x;
}


//! Array-scalar modulus.
/*!
  \relates EquilateralArrayImp
  \note This does not check for division by zero as the value type may not be
  as number type.
*/
template<typename _T, std::size_t _D, std::size_t _N, typename _B>
inline
EquilateralArrayImp<_T, _D, _N, _B>&
operator%=(EquilateralArrayImp<_T, _D, _N, _B>& x,
           typename Loki::TypeTraits<_T>::ParameterType value) {
   typedef typename EquilateralArrayImp<_T, _D, _N, _B>::iterator iterator;
   for (iterator i = x.begin(); i != x.end(); ++i) {
      *i %= value;
   }
   return x;
}


//! Left shift.
/*! \relates EquilateralArrayImp */
template<typename _T, std::size_t _D, std::size_t _N, typename _B>
inline
EquilateralArrayImp<_T, _D, _N, _B>&
operator<<=(EquilateralArrayImp<_T, _D, _N, _B>& x, const int offset) {
   typedef typename EquilateralArrayImp<_T, _D, _N, _B>::iterator iterator;
   for (iterator i = x.begin(); i != x.end(); ++i) {
      *i <<= offset;
   }
   return x;
}


//! Right shift.
/*! \relates EquilateralArrayImp */
template<typename _T, std::size_t _D, std::size_t _N, typename _B>
inline
EquilateralArrayImp<_T, _D, _N, _B>&
operator>>=(EquilateralArrayImp<_T, _D, _N, _B>& x, const int offset) {
   typedef typename EquilateralArrayImp<_T, _D, _N, _B>::iterator iterator;
   for (iterator i = x.begin(); i != x.end(); ++i) {
      *i >>= offset;
   }
   return x;
}


//@}


} // namespace array

#endif
