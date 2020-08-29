// -*- C++ -*-

/*!
  \file cuda/equality.h
  \brief Floating-point approximate equality.
*/

#if !defined(__cuda_equality_h__)
#define __cuda_equality_h__

#include "../ext/array.h"
#include "limits.h"
#include "iterator.h"

//#include <iterator>
//#include <vector>

//#include <cmath>

namespace numerical {

/*! \defgroup numericalEquality Floating-point Approximate Equality */
//@{

//! Return true if the two numbers are approximately equal.
/*!
  \return True if |x-y| <= t eps |x| and |x-y| <= t eps |y| where
  t is the \c toleranceFactor and eps is the machine precision.

  The default tolerance factor is two.
*/
bool
areEqual(double x, double y, double toleranceFactor = 2);

//! Return true if the two numbers are approximately equal.
/*!
  \return True if |x-y| <= t eps |x| and |x-y| <= t eps |y| where
  t is the \c toleranceFactor and eps is the machine precision.

  The default tolerance factor is two.
*/
bool
areEqual(float x, float y, float toleranceFactor = 2);

//! Return true if the two arrays are approximately equal.
/*! The default tolerance factor is two. */
template<typename _T, std::size_t _N>
bool
areEqual(const std::tr1::array<_T, _N>& x, const std::tr1::array<_T, _N>& y,
         _T toleranceFactor = 2);

//! Return true if the two sequences are approximately equal.
/*!
  \return True if each pair is approximately equal.

  The default tolerance factor is two.
*/
template<typename _InputIterator1, typename _InputIterator2>
bool
areSequencesEqual(_InputIterator1 begin1, _InputIterator1 end1,
                  _InputIterator2 begin2,
                  typename std::iterator_traits<_InputIterator1>::value_type
                  toleranceFactor = 2);


//! Return true if the two numbers are approximately equal. Use the absolute error.
/*!
  \return True if |x-y| <= s t eps, where \e s is the scale, \e t is
  the \c toleranceFactor, and eps is the machine precision.

  The default scale is one and the default tolerance factor is two.
*/
bool
areEqualAbs(double x, double y, double scale = 1, double toleranceFactor = 2);

//! Return true if the two numbers are approximately equal. Use the absolute error.
/*!
  \return True if |x-y| <= s t eps, where \e s is the scale, \e t is
  the \c toleranceFactor, and eps is the machine precision.

  The default scale is one and the default tolerance factor is two.
*/
bool
areEqualAbs(float x, float y, float scale = 1, float toleranceFactor = 2);

//! Return true if the two arrays are approximately equal. Use the absolute error.
/*! The default scale is one and the default tolerance factor is two. */
template<typename _T, std::size_t _N>
bool
areEqualAbs(const std::tr1::array<_T, _N>& x, const std::tr1::array<_T, _N>& y,
            _T scale = 1, _T toleranceFactor = 2);

//! Return true if the two sequences are approximately equal. Use the absolute error.
/*!
  \return True if each pair is approximately equal.

  The default scale is one and the default tolerance factor is two.
*/
template<typename _InputIterator1, typename _InputIterator2>
bool
areSequencesEqualAbs(_InputIterator1 begin1, _InputIterator1 end1,
                     _InputIterator2 begin2,
                     typename std::iterator_traits<_InputIterator1>::value_type
                     scale = 1,
                     typename std::iterator_traits<_InputIterator1>::value_type
                     toleranceFactor = 2);


//! Return true if the number is small.
/*!
  \return True if |x| <= t * eps where
  t is the \c toleranceFactor and eps is the machine precision.

  The default tolerance factor is unity.
*/
bool
isSmall(double x, double toleranceFactor = 1);

//! Return true if the number is small.
/*!
  \return True if |x| <= t * eps where
  t is the \c toleranceFactor and eps is the machine precision.

  The default tolerance factor is unity.
*/
bool
isSmall(float x, float toleranceFactor = 1);

//@}

} // namespace numerical

#define __cuda_equality_ipp__
#include "equality.ipp"
#undef __cuda_equality_ipp__

#endif
