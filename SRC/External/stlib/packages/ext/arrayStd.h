// -*- C++ -*-

/*!
  \file ext/array.h
  \brief Functions for array.
*/

#if !defined(__ext_arrayStd_h__)
#define __ext_arrayStd_h__

#include "../loki/Typelist.h"

#if _MSC_VER >= 1600
#include <array>
#else
#include <tr1/array>
#endif

#include <algorithm>
#include <numeric>
#include <iostream>
#include <iterator>
#include <limits>

#include <cmath>
#include <cassert>

#ifdef BOOST
namespace boost {
#else
namespace std {
namespace tr1 {
#endif

/*!
\page extArray Extensions to std::tr1::array

Here we provide functions to extend the functionality of the std::tr1::array
class [\ref extBecker2007 "Becker, 2007"]. The functions are grouped into
the following categories:
- \ref extArrayAssignmentScalar
- \ref extArrayAssignmentArray
- \ref extArrayUnary
- \ref extArrayBinary
- \ref extArrayFile
- \ref extArrayMathematical
- \ref extArrayMake
*/

//-------------------------------------------------------------------------
/*! \defgroup extArrayAssignmentScalar Assignment Operators with a Scalar Operand.

  These functions apply an assignment operation to each element of the array.
  Thus you can add, subtract, etc. a scalar from each element.
  \verbatim
  std::tr1::array<double, 3> x = {{0, 0, 0}};
  x += 2;
  x -= 3;
  x *= 5;
  x /= 7; \endverbatim

  \verbatim
  std::tr1::array<unsigned, 3> x = {{0, 0, 0}};
  x %= 2;
  x <<= 3;
  x >>= 5; \endverbatim
*/
//@{

//! Array-scalar addition.
template<typename _T, size_t _Dim, typename _T2>
inline
array<_T, _Dim>&
operator+=(array<_T, _Dim>& x, const _T2& value) {
   for (typename array<_T, _Dim>::iterator i = x.begin(); i != x.end(); ++i) {
      *i += value;
   }
   return x;
}

//! Array-scalar subtraction.
template<typename _T, size_t _Dim, typename _T2>
inline
array<_T, _Dim>&
operator-=(array<_T, _Dim>& x, const _T2& value) {
   for (typename array<_T, _Dim>::iterator i = x.begin(); i != x.end(); ++i) {
      *i -= value;
   }
   return x;
}

//! Array-scalar multiplication.
template<typename _T, size_t _Dim, typename _T2>
inline
array<_T, _Dim>&
operator*=(array<_T, _Dim>& x, const _T2& value) {
   for (typename array<_T, _Dim>::iterator i = x.begin(); i != x.end(); ++i) {
      *i *= value;
   }
   return x;
}

//! Array-scalar division.
template<typename _T, size_t _Dim, typename _T2>
inline
array<_T, _Dim>&
operator/=(array<_T, _Dim>& x, const _T2& value) {
#ifdef DEBUG_stlib
   assert(value != 0);
#endif
   for (typename array<_T, _Dim>::iterator i = x.begin(); i != x.end(); ++i) {
      *i /= value;
   }
   return x;
}

//! Array-scalar modulus.
template<typename _T, size_t _Dim, typename _T2>
inline
array<_T, _Dim>&
operator%=(array<_T, _Dim>& x, const _T2& value) {
#ifdef DEBUG_stlib
   assert(value != 0);
#endif
   for (typename array<_T, _Dim>::iterator i = x.begin(); i != x.end(); ++i) {
      *i %= value;
   }
   return x;
}

//! Left shift.
template<typename _T, size_t _Dim>
inline
array<_T, _Dim>&
operator<<=(array<_T, _Dim>& x, const int offset) {
   for (typename array<_T, _Dim>::iterator i = x.begin(); i != x.end(); ++i) {
      *i <<= offset;
   }
   return x;
}

//! Right shift.
template<typename _T, size_t _Dim>
inline
array<_T, _Dim>&
operator>>=(array<_T, _Dim>& x, const int offset) {
   for (typename array<_T, _Dim>::iterator i = x.begin(); i != x.end(); ++i) {
      *i >>= offset;
   }
   return x;
}

//@}
//-------------------------------------------------------------------------
/*! \defgroup extArrayAssignmentArray Assignment Operators with an Array Operand.

  These functions define assignment operations between arrays. The assignment
  operation is applied element-wise.
  \verbatim
  std::tr1::array<double, 3> x, y;
  ...
  x += y;
  x -= y;
  x *= y;
  x /= y; \endverbatim
  \verbatim
  std::tr1::array<unsigned, 3> x, y;
  ...
  x %= y;
  x <<= y;
  x >>= y; \endverbatim
*/
//@{

//! Array-array addition.
template<typename _T, size_t _Dim, typename _T2>
inline
array<_T, _Dim>&
operator+=(array<_T, _Dim>& x, const array<_T2, _Dim>& y) {
   for (size_t n = 0; n != x.size(); ++n) {
      x[n] += y[n];
   }
   return x;
}

//! Array-array subtraction.
template<typename _T, size_t _Dim, typename _T2>
inline
array<_T, _Dim>&
operator-=(array<_T, _Dim>& x, const array<_T2, _Dim>& y) {
   for (size_t n = 0; n != x.size(); ++n) {
      x[n] -= y[n];
   }
   return x;
}

//! Array-array multiplication.
template<typename _T, size_t _Dim, typename _T2>
inline
array<_T, _Dim>&
operator*=(array<_T, _Dim>& x, const array<_T2, _Dim>& y) {
   for (size_t n = 0; n != x.size(); ++n) {
      x[n] *= y[n];
   }
   return x;
}

//! Array-array division.
template<typename _T, size_t _Dim, typename _T2>
inline
array<_T, _Dim>&
operator/=(array<_T, _Dim>& x, const array<_T2, _Dim>& y) {
   for (size_t n = 0; n != x.size(); ++n) {
#ifdef DEBUG_stlib
      assert(y[n] != 0);
#endif
      x[n] /= y[n];
   }
   return x;
}

//! Array-array modulus.
template<typename _T, size_t _Dim, typename _T2>
inline
array<_T, _Dim>&
operator%=(array<_T, _Dim>& x, const array<_T2, _Dim>& y) {
   for (size_t n = 0; n != x.size(); ++n) {
#ifdef DEBUG_stlib
      assert(y[n] != 0);
#endif
      x[n] %= y[n];
   }
   return x;
}

//! Array-array left shift.
template<typename _T1, size_t _Dim, typename _T2>
inline
array<_T1, _Dim>&
operator<<=(array<_T1, _Dim>& x, const array<_T2, _Dim>& y) {
   for (size_t n = 0; n != x.size(); ++n) {
      x[n] <<= y[n];
   }
   return x;
}

//! Array-array right shift.
template<typename _T1, size_t _Dim, typename _T2>
inline
array<_T1, _Dim>&
operator>>=(array<_T1, _Dim>& x, const array<_T2, _Dim>& y) {
   for (size_t n = 0; n != x.size(); ++n) {
      x[n] >>= y[n];
   }
   return x;
}

//@}
//-------------------------------------------------------------------------
/*! \defgroup extArrayUnary Unary Operators

  These functions define unary operations for arrays.
  \verbatim
  std::tr1::array<double, 3> x, y;
  ...
  x = +y;
  x = -y; \endverbatim
*/
//@{

//! Unary positive operator.
template<typename _T, size_t _Dim>
inline
const array<_T, _Dim>&
operator+(const array<_T, _Dim>& x) {
   return x;
}

//! Unary negative operator.
template<typename _T, size_t _Dim>
inline
array<_T, _Dim>
operator-(const array<_T, _Dim>& x) {
   array<_T, _Dim> y(x);
   for (size_t n = 0; n != _Dim; ++n) {
      y[n] = -y[n];
   }
   return y;
}

//@}
//-------------------------------------------------------------------------
/*! \defgroup extArrayBinary Binary Operators
  These functions define binary operators for arrays. The operands may be
  arrays or scalars.

  \verbatim
  std::tr1::array<double, 3> x, y;
  ...
  // Addition
  x = x + y;
  x = x + 2.;
  x = 2. + x;
  // Subtraction.
  x = x - y;
  x = x - 2.;
  x = 2. - x;
  // Multiplication.
  x = x * y;
  x = x * 2.;
  x = 2. * x;
  // Division.
  x = x / y;
  x = x / 2.;
  x = 2. / x; \endverbatim
  \verbatim
  std::tr1::array<unsigned, 3> x, y;
  ...
  // Modulus.
  x = x % y;
  x = x % 2;
  x = 2 % x; \endverbatim

  \note Because these function instantiate std::tr1::array objects,
  they are not as efficient
  as their corresponding \ref extArrayAssignmentArray "assignment operators."
  For example, the following
  \verbatim
  std::tr1::array<double, 3> center;
  double radius;
  ...
  center += radius; \endverbatim
  is more efficient than
  \verbatim
  center = center + radius; \endverbatim
*/
//@{

//-------------------------------------------------------------------------
// These do not use the typeof operator.


//! The result type for an arithmetic expression.
/*!
  \todo This is a dirty hack and does not work for arbitrary types.
  Use the typeof versions when it makes its way into the C++ standard.
*/
template<typename _XCoord, typename _YCoord>
struct ArithmeticResult {
   //! An ordered list of numeric types.
   typedef LOKI_TYPELIST_10(double, float, unsigned long, long, unsigned, int,
                            unsigned short, short, unsigned char, signed char)
   OrderedTypes;
   //! The result type.
   typedef typename Loki::Select < int(Loki::TL::IndexOf<OrderedTypes, _XCoord>::value)
   < int(Loki::TL::IndexOf<OrderedTypes, _YCoord>::value),
   _XCoord, _YCoord >::Result Type;
};


//! Array-scalar addition.
template<typename _T, typename _Value, size_t _Dim>
inline
array<typename ArithmeticResult<_T, _Value>::Type, _Dim>
operator+(const array<_T, _Dim>& x, const _Value& y) {
   array<typename ArithmeticResult<_T, _Value>::Type, _Dim> z;
   for (std::size_t i = 0; i != _Dim; ++i) {
      z[i] = x[i] + y;
   }
   return z;
}

//! Scalar-Array addition.
template<typename _T, typename _Value, size_t _Dim>
inline
array<typename ArithmeticResult<_T, _Value>::Type, _Dim>
operator+(const _T& x, const array<_Value, _Dim>& y) {
   array<typename ArithmeticResult<_T, _Value>::Type, _Dim> z;
   for (std::size_t i = 0; i != _Dim; ++i) {
      z[i] = x + y[i];
   }
   return z;
}

//! Array-array addition.
template<typename _T, typename _Value, size_t _Dim>
inline
array<typename ArithmeticResult<_T, _Value>::Type, _Dim>
operator+(const array<_T, _Dim>& x, const array<_Value, _Dim>& y) {
   array<typename ArithmeticResult<_T, _Value>::Type, _Dim> z;
   for (std::size_t i = 0; i != _Dim; ++i) {
      z[i] = x[i] + y[i];
   }
   return z;
}

//! Array-scalar subtraction.
template<typename _T, typename _Value, size_t _Dim>
inline
array<typename ArithmeticResult<_T, _Value>::Type, _Dim>
operator-(const array<_T, _Dim>& x, const _Value& y) {
   array<typename ArithmeticResult<_T, _Value>::Type, _Dim> z;
   for (std::size_t i = 0; i != _Dim; ++i) {
      z[i] = x[i] - y;
   }
   return z;
}

//! Scalar-Array subtraction.
template<typename _T, typename _Value, size_t _Dim>
inline
array<typename ArithmeticResult<_T, _Value>::Type, _Dim>
operator-(const _T& x, const array<_Value, _Dim>& y) {
   array<typename ArithmeticResult<_T, _Value>::Type, _Dim> z;
   for (std::size_t i = 0; i != _Dim; ++i) {
      z[i] = x - y[i];
   }
   return z;
}

//! Array-array subtraction.
template<typename _T, typename _Value, size_t _Dim>
inline
array<typename ArithmeticResult<_T, _Value>::Type, _Dim>
operator-(const array<_T, _Dim>& x, const array<_Value, _Dim>& y) {
   array<typename ArithmeticResult<_T, _Value>::Type, _Dim> z;
   for (std::size_t i = 0; i != _Dim; ++i) {
      z[i] = x[i] - y[i];
   }
   return z;
}

//! Array-scalar multiplication.
template<typename _T, typename _Value, size_t _Dim>
inline
array<typename ArithmeticResult<_T, _Value>::Type, _Dim>
operator*(const array<_T, _Dim>& x, const _Value& y) {
   array<typename ArithmeticResult<_T, _Value>::Type, _Dim> z;
   for (std::size_t i = 0; i != _Dim; ++i) {
      z[i] = x[i] * y;
   }
   return z;
}

//! Scalar-Array multiplication.
template<typename _T, typename _Value, size_t _Dim>
inline
array<typename ArithmeticResult<_T, _Value>::Type, _Dim>
operator*(const _T& x, const array<_Value, _Dim>& y) {
   array<typename ArithmeticResult<_T, _Value>::Type, _Dim> z;
   for (std::size_t i = 0; i != _Dim; ++i) {
      z[i] = x * y[i];
   }
   return z;
}

//! Array-array multiplication.
template<typename _T, typename _Value, size_t _Dim>
inline
array<typename ArithmeticResult<_T, _Value>::Type, _Dim>
operator*(const array<_T, _Dim>& x, const array<_Value, _Dim>& y) {
   array<typename ArithmeticResult<_T, _Value>::Type, _Dim> z;
   for (std::size_t i = 0; i != _Dim; ++i) {
      z[i] = x[i] * y[i];
   }
   return z;
}

//! Array-scalar division.
template<typename _T, typename _Value, size_t _Dim>
inline
array<typename ArithmeticResult<_T, _Value>::Type, _Dim>
operator/(const array<_T, _Dim>& x, const _Value& y) {
#ifdef DEBUG_stlib
   assert(y != 0);
#endif
   array<typename ArithmeticResult<_T, _Value>::Type, _Dim> z;
   for (std::size_t i = 0; i != _Dim; ++i) {
      z[i] = x[i] / y;
   }
   return z;
}

//! Scalar-Array division.
template<typename _T, typename _Value, size_t _Dim>
inline
array<typename ArithmeticResult<_T, _Value>::Type, _Dim>
operator/(const _T& x, const array<_Value, _Dim>& y) {
   array<typename ArithmeticResult<_T, _Value>::Type, _Dim> z;
   for (std::size_t i = 0; i != _Dim; ++i) {
#ifdef DEBUG_stlib
      assert(y[i] != 0);
#endif
      z[i] = x / y[i];
   }
   return z;
}

//! Array-array division.
template<typename _T, typename _Value, size_t _Dim>
inline
array<typename ArithmeticResult<_T, _Value>::Type, _Dim>
operator/(const array<_T, _Dim>& x, const array<_Value, _Dim>& y) {
   array<typename ArithmeticResult<_T, _Value>::Type, _Dim> z;
   for (std::size_t i = 0; i != _Dim; ++i) {
#ifdef DEBUG_stlib
      assert(y[i] != 0);
#endif
      z[i] = x[i] / y[i];
   }
   return z;
}

//! Array-scalar modulus.
template<typename _T, typename _Value, size_t _Dim>
inline
array<typename ArithmeticResult<_T, _Value>::Type, _Dim>
operator%(const array<_T, _Dim>& x, const _Value& y) {
#ifdef DEBUG_stlib
   assert(y != 0);
#endif
   array<typename ArithmeticResult<_T, _Value>::Type, _Dim> z;
   for (std::size_t i = 0; i != _Dim; ++i) {
      z[i] = x[i] % y;
   }
   return z;
}

//! Scalar-Array modulus.
template<typename _T, typename _Value, size_t _Dim>
inline
array<typename ArithmeticResult<_T, _Value>::Type, _Dim>
operator%(const _T& x, const array<_Value, _Dim>& y) {
   array<typename ArithmeticResult<_T, _Value>::Type, _Dim> z;
   for (std::size_t i = 0; i != _Dim; ++i) {
#ifdef DEBUG_stlib
      assert(y[i] != 0);
#endif
      z[i] = x % y[i];
   }
   return z;
}

//! Array-array modulus.
template<typename _T, typename _Value, size_t _Dim>
inline
array<typename ArithmeticResult<_T, _Value>::Type, _Dim>
operator%(const array<_T, _Dim>& x, const array<_Value, _Dim>& y) {
   array<typename ArithmeticResult<_T, _Value>::Type, _Dim> z;
   for (std::size_t i = 0; i != _Dim; ++i) {
#ifdef DEBUG_stlib
      assert(y[i] != 0);
#endif
      z[i] = x[i] % y[i];
   }
   return z;
}






// CONTINUE: Use this when typeof becomes standard.
#if 0
//-------------------------------------------------------------------------
// These use the typeof operator.

//! Array-scalar addition.
template<typename _T, typename _Value, size_t _Dim>
inline
array < typeof(_T() + _Value()), _Dim >
operator+(const array<_T, _Dim>& x, const _Value& y) {
   array < typeof(_T() + _Value()), _Dim > z;
   for (std::size_t i = 0; i != _Dim; ++i) {
      z[i] = x[i] + y;
   }
   return z;
}

//! Scalar-Array addition.
template<typename _T, typename _Value, size_t _Dim>
inline
array < typeof(_T() + _Value()), _Dim >
operator+(const _T& x, const array<_Value, _Dim>& y) {
   array < typeof(_T() + _Value()), _Dim > z;
   for (std::size_t i = 0; i != _Dim; ++i) {
      z[i] = x + y[i];
   }
   return z;
}

//! Array-array addition.
template<typename _T, typename _Value, size_t _Dim>
inline
array < typeof(_T() + _Value()), _Dim >
operator+(const array<_T, _Dim>& x, const array<_Value, _Dim>& y) {
   array < typeof(_T() + _Value()), _Dim > z;
   for (std::size_t i = 0; i != _Dim; ++i) {
      z[i] = x[i] + y[i];
   }
   return z;
}

//! Array-scalar subtraction.
template<typename _T, typename _Value, size_t _Dim>
inline
array < typeof(_T() - _Value()), _Dim >
operator-(const array<_T, _Dim>& x, const _Value& y) {
   array < typeof(_T() - _Value()), _Dim > z;
   for (std::size_t i = 0; i != _Dim; ++i) {
      z[i] = x[i] - y;
   }
   return z;
}

//! Scalar-Array subtraction.
template<typename _T, typename _Value, size_t _Dim>
inline
array < typeof(_T() - _Value()), _Dim >
operator-(const _T& x, const array<_Value, _Dim>& y) {
   array < typeof(_T() - _Value()), _Dim > z;
   for (std::size_t i = 0; i != _Dim; ++i) {
      z[i] = x - y[i];
   }
   return z;
}

//! Array-array subtraction.
template<typename _T, typename _Value, size_t _Dim>
inline
array < typeof(_T() - _Value()), _Dim >
operator-(const array<_T, _Dim>& x, const array<_Value, _Dim>& y) {
   array < typeof(_T() - _Value()), _Dim > z;
   for (std::size_t i = 0; i != _Dim; ++i) {
      z[i] = x[i] - y[i];
   }
   return z;
}

//! Array-scalar multiplication.
template<typename _T, typename _Value, size_t _Dim>
inline
array < typeof(_T()*_Value()), _Dim >
operator*(const array<_T, _Dim>& x, const _Value& y) {
   array < typeof(_T() * _Value()), _Dim > z;
   for (std::size_t i = 0; i != _Dim; ++i) {
      z[i] = x[i] * y;
   }
   return z;
}

//! Scalar-Array multiplication.
template<typename _T, typename _Value, size_t _Dim>
inline
array < typeof(_T()*_Value()), _Dim >
operator*(const _T& x, const array<_Value, _Dim>& y) {
   array < typeof(_T() * _Value()), _Dim > z;
   for (std::size_t i = 0; i != _Dim; ++i) {
      z[i] = x * y[i];
   }
   return z;
}

//! Array-array multiplication.
template<typename _T, typename _Value, size_t _Dim>
inline
array < typeof(_T()*_Value()), _Dim >
operator*(const array<_T, _Dim>& x, const array<_Value, _Dim>& y) {
   array < typeof(_T() * _Value()), _Dim > z;
   for (std::size_t i = 0; i != _Dim; ++i) {
      z[i] = x[i] * y[i];
   }
   return z;
}

//! Array-scalar division.
template<typename _T, typename _Value, size_t _Dim>
inline
array < typeof(_T() / _Value()), _Dim >
operator/(const array<_T, _Dim>& x, const _Value& y) {
#ifdef DEBUG_stlib
   assert(y != 0);
#endif
   array < typeof(_T() / _Value()), _Dim > z;
   for (std::size_t i = 0; i != _Dim; ++i) {
      z[i] = x[i] / y;
   }
   return z;
}

//! Scalar-Array division.
template<typename _T, typename _Value, size_t _Dim>
inline
array < typeof(_T() / _Value()), _Dim >
operator/(const _T& x, const array<_Value, _Dim>& y) {
   array < typeof(_T() / _Value()), _Dim > z;
   for (std::size_t i = 0; i != _Dim; ++i) {
#ifdef DEBUG_stlib
      assert(y[i] != 0);
#endif
      z[i] = x / y[i];
   }
   return z;
}

//! Array-array division.
template<typename _T, typename _Value, size_t _Dim>
inline
array < typeof(_T() / _Value()), _Dim >
operator/(const array<_T, _Dim>& x, const array<_Value, _Dim>& y) {
   array < typeof(_T() / _Value()), _Dim > z;
   for (std::size_t i = 0; i != _Dim; ++i) {
#ifdef DEBUG_stlib
      assert(y[i] != 0);
#endif
      z[i] = x[i] / y[i];
   }
   return z;
}

//! Array-scalar modulus.
template<typename _T, typename _Value, size_t _Dim>
inline
array < typeof(_T() % _Value()), _Dim >
operator%(const array<_T, _Dim>& x, const _Value& y) {
#ifdef DEBUG_stlib
   assert(y != 0);
#endif
   array < typeof(_T() % _Value()), _Dim > z;
   for (std::size_t i = 0; i != _Dim; ++i) {
      z[i] = x[i] % y;
   }
   return z;
}

//! Scalar-Array modulus.
template<typename _T, typename _Value, size_t _Dim>
inline
array < typeof(_T() % _Value()), _Dim >
operator%(const _T& x, const array<_Value, _Dim>& y) {
   array < typeof(_T() % _Value()), _Dim > z;
   for (std::size_t i = 0; i != _Dim; ++i) {
#ifdef DEBUG_stlib
      assert(y[i] != 0);
#endif
      z[i] = x % y[i];
   }
   return z;
}

//! Array-array modulus.
template<typename _T, typename _Value, size_t _Dim>
inline
array < typeof(_T() % _Value()), _Dim >
operator%(const array<_T, _Dim>& x, const array<_Value, _Dim>& y) {
   array < typeof(_T() % _Value()), _Dim > z;
   for (std::size_t i = 0; i != _Dim; ++i) {
#ifdef DEBUG_stlib
      assert(y[i] != 0);
#endif
      z[i] = x[i] % y[i];
   }
   return z;
}

#endif

#if 0
//-------------------------------------------------------------------------
// Old implementation.

//! Array-array addition.
template<typename _T, size_t _Dim, typename _T2>
inline
array<_T, _Dim>
operator+(const array<_T, _Dim>& x, const array<_T2, _Dim>& y) {
   array<_T, _Dim> z = x;
   return z += y;
}

//! Array-scalar addition.
template<typename _T, size_t _Dim>
inline
array<_T, _Dim>
operator+(const array<_T, _Dim>& x, const _T& y) {
   array<_T, _Dim> z = x;
   return z += y;
}

//! Scalar-Array addition.
template<typename _T, size_t _Dim>
inline
array<_T, _Dim>
operator+(const _T& y, const array<_T, _Dim>& x) {
   return x + y;
}

//! Array-array subtraction.
template<typename _T, size_t _Dim, typename _T2>
inline
array<_T, _Dim>
operator-(const array<_T, _Dim>& x, const array<_T2, _Dim>& y) {
   array<_T, _Dim> z = x;
   return z -= y;
}

//! Array-scalar subtraction.
template<typename _T, size_t _Dim>
inline
array<_T, _Dim>
operator-(const array<_T, _Dim>& x, const _T& y) {
   array<_T, _Dim> z = x;
   return z -= y;
}

//! Scalar-Array subtraction.
template<typename _T, size_t _Dim>
inline
array<_T, _Dim>
operator-(const _T& y, const array<_T, _Dim>& x) {
   return - x + y;
}

//! Array-array multiplication.
template<typename _T, size_t _Dim, typename _T2>
inline
array<_T, _Dim>
operator*(const array<_T, _Dim>& x, const array<_T2, _Dim>& y) {
   array<_T, _Dim> z = x;
   return z *= y;
}

//! Array-scalar multiplication.
template<typename _T, size_t _Dim>
inline
array<_T, _Dim>
operator*(const array<_T, _Dim>& x, const _T& y) {
   array<_T, _Dim> z = x;
   return z *= y;
}

//! Scalar-Array multiplication.
template<typename _T, size_t _Dim>
inline
array<_T, _Dim>
operator*(const _T& y, const array<_T, _Dim>& x) {
   return x * y;
}

//! Array-array division.
template<typename _T, size_t _Dim, typename _T2>
inline
array<_T, _Dim>
operator/(const array<_T, _Dim>& x, const array<_T2, _Dim>& y) {
   array<_T, _Dim> z = x;
   return z /= y;
}

//! Array-scalar division.
template<typename _T, size_t _Dim>
inline
array<_T, _Dim>
operator/(const array<_T, _Dim>& x, const _T& y) {
   array<_T, _Dim> z = x;
   return z /= y;
}

//! Scalar-Array division.
template<typename _T, size_t _Dim>
inline
array<_T, _Dim>
operator/(const _T& y, const array<_T, _Dim>& x) {
   array<_T, _Dim> z;
   fill(z.begin(), z.end(), y);
   return z /= x;
}

//! Array-array modulus.
template<typename _T, size_t _Dim, typename _T2>
inline
array<_T, _Dim>
operator%(const array<_T, _Dim>& x, const array<_T2, _Dim>& y) {
   array<_T, _Dim> z = x;
   return z %= y;
}

//! Array-scalar modulus.
template<typename _T, size_t _Dim>
inline
array<_T, _Dim>
operator%(const array<_T, _Dim>& x, const _T& y) {
   array<_T, _Dim> z = x;
   return z %= y;
}

//! Scalar-Array modulus.
template<typename _T, size_t _Dim>
inline
array<_T, _Dim>
operator%(const _T& y, const array<_T, _Dim>& x) {
   array<_T, _Dim> z;
   fill(z.begin(), z.end(), y);
   return z %= x;
}

#endif

//@}
//-------------------------------------------------------------------------
/*! \defgroup extArrayFile File I/O

  These functions write and read std::tr1::array's in ascii and binary format.
  The file format is simply the sequence of elements, the number of elements
  is not read or written.
  \verbatim
  std::tr1::array<double, 3> x;
  ...
  // Ascii.
  std::cin >> x;
  std::cout << x;
  // Binary.
  std::ifstream in("input.bin");
  read(in, &x);
  std::ofstream out("output.bin");
  write(out, x); \endverbatim
*/
//@{

//! Write the space-separated elements.
/*!
  Format:
  \verbatim
  x[0] x[1] x[2] ... \endverbatim
*/
template<typename _T, size_t _Dim>
inline
std::ostream&
operator<<(std::ostream& out, const array<_T, _Dim>& x) {
   std::copy(x.begin(), x.end(), std::ostream_iterator<_T>(out, " "));
   return out;
}

//! Read the elements.
template<typename _T, size_t _Dim>
inline
std::istream&
operator>>(std::istream& in, array<_T, _Dim>& x) {
   for (size_t n = 0; n != x.size(); ++n) {
      in >> x[n];
   }
   return in;
}

//! Write the elements in binary format.
template<typename _T, size_t _Dim>
inline
void
write(const array<_T, _Dim>& x, std::ostream& out) {
   out.write(reinterpret_cast<const char*>(&x), sizeof(array<_T, _Dim>));
}

//! Read the elements in binary format.
template<typename _T, size_t _Dim>
inline
void
read(array<_T, _Dim>* x, std::istream& in) {
   in.read(reinterpret_cast<char*>(x), sizeof(array<_T, _Dim>));
}

//@}
//-------------------------------------------------------------------------
/*! \defgroup extArrayMathematical Mathematical Functions
  These functions define some common mathematical operations on
  std::tr1::array's. There are utility functions for the sum, product,
  minimum, maximum, etc.
  \verbatim
  std::tr1::array<double, 3> x, y, z;
  ...
  // Sum.
  const double total = sum(x);
  // Product.
  const double volume = product(x);
  // Minimum.
  const double minValue = min(x);
  // Maximum.
  const double maxValue = max(x);
  // Element-wise minimum.
  z = min(x, y);
  // Element-wise maximum.
  z = min(x, y);
  // Existence of a value.
  const bool hasNull = hasElement(x, 0);
  // Existence and index of a value.
  std::size_t i;
  const bool hasOne = hasElement(x, 1, &i);
  // Index of a value.
  i = index(x); \endverbatim

  There are also function for treating a std::tr1::array as a Cartesian
  point or vector.
  \verbatim
  std::tr1::array<double, 3> x, y, z;
  ...
  // Dot product.
  const double d = dot(x, y);
  // Cross product.
  z = cross(x, y);
  // Cross product that avoids constructing an array.
  cross(x, y, &z);
  // Triple product.
  const double volume = tripleProduct(x, y, z);
  // Discriminant.
  std::tr1::array<double, 2> a, b;
  const double disc = discriminant(a, b);
  // Squared magnitude.
  const double sm = squaredMagnitude(x);
  // Magnitude.
  const double m = magnitude(x);
  // Normalize a vector.
  normalize(&x);
  // Negate a vector.
  negateElements(&x);
  // Squared distance.
  const double d2 = squaredDistance(x, y);
  // Distance.
  const double d = euclideanDistance(x, y); \endverbatim
*/
//@{

//! Return the sum of the components.
template<typename _T, size_t _Dim>
inline
_T
sum(const array<_T, _Dim>& x) {
   return std::accumulate(x.begin(), x.end(), _T(0));
}

//! Return the product of the components.
template<typename _T, size_t _Dim>
inline
_T
product(const array<_T, _Dim>& x) {
   return std::accumulate(x.begin(), x.end(), _T(1), std::multiplies<_T>());
}

//! Return the minimum component.  Use < for comparison.
template<typename _T, size_t _Dim>
inline
_T
min(const array<_T, _Dim>& x) {
#ifdef DEBUG_stlib
   assert(x.size() != 0);
#endif
   return *std::min_element(x.begin(), x.end());
}

//! Return the maximum component.  Use > for comparison.
template<typename _T, size_t _Dim>
inline
_T
max(const array<_T, _Dim>& x) {
#ifdef DEBUG_stlib
   assert(x.size() != 0);
#endif
   return *std::max_element(x.begin(), x.end());
}

//! Return an array that is element-wise the minimum of the two.
template<typename _T, size_t _Dim>
inline
array<_T, _Dim>
min(const array<_T, _Dim>& x, const array<_T, _Dim>& y) {
   array<_T, _Dim> z;
   for (size_t n = 0; n != _Dim; ++n) {
      z[n] = std::min(x[n], y[n]);
   }
   return z;
}

//! Return an array that is element-wise the maximum of the two.
template<typename _T, size_t _Dim>
inline
array<_T, _Dim>
max(const array<_T, _Dim>& x, const array<_T, _Dim>& y) {
   array<_T, _Dim> z;
   for (size_t n = 0; n != _Dim; ++n) {
      z[n] = std::max(x[n], y[n]);
   }
   return z;
}

//! Return true if the array has the specified element.
template<typename _T, size_t _Dim, typename _Comparable>
inline
bool
hasElement(const array<_T, _Dim>& x, const _Comparable& a) {
   return std::count(x.begin(), x.end(), a);
}

//! Return true if the array has the specified element.
/*!
  If true, compute the index of the elementt.
*/
template<typename _T, size_t _Dim, typename _Comparable>
inline
bool
hasElement(const array<_T, _Dim>& x, const _Comparable& a, std::size_t* i) {
   for (*i = 0; *i != x.size(); ++*i) {
      if (a == x[*i]) {
         return true;
      }
   }
   return false;
}

//! Return the index of the specified element. Return std::numeric_limits<std::size_t>::max() if the element is not in the array.
template<typename _T, size_t _Dim, typename _Comparable>
inline
std::size_t
index(const array<_T, _Dim>& x, const _Comparable& a) {
   for (std::size_t i = 0; i != x.size(); ++i) {
      if (a == x[i]) {
         return i;
      }
   }
   return std::numeric_limits<std::size_t>::max();
}

//! Return the dot product of the two arrays.
template<typename _T, size_t _Dim>
inline
_T
dot(const array<_T, _Dim>& x, const array<_T, _Dim>& y) {
   // Clean version:
   // return std::inner_product(x.begin(), x.end(), y.begin(), _T(0));
   // More efficient because of loop unrolling:
   _T p = 0;
   for (std::size_t i = 0; i != _Dim; ++i) {
      p += x[i] * y[i];
   }
   return p;
}

//! Return the dot product of the two arrays.
/*! This specialization is a little faster than the dimension-general code. */
template<typename _T>
inline
_T
dot(const array<_T, 3>& x, const array<_T, 3>& y) {
   return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}

//! Return the cross product of the two arrays.
template<typename _T>
inline
array<_T, 3>
cross(const array<_T, 3>& x, const array<_T, 3>& y) {
   array<_T, 3> result = {{
         x[1]*y[2] - y[1]*x[2],
         y[0]*x[2] - x[0]*y[2],
         x[0]*y[1] - y[0]*x[1]
      }
   };
   return result;
}

//! Return the cross product and derivative of cross product of the two arrays.
template<typename _T>
inline
array<_T, 3>
cross(const array<_T, 3>& x, const array<_T, 3>& y,
      array<array<_T, 3>, 3>* dx, array<array<_T, 3>, 3>* dy) {
   (*dx)[0][0] = 0.0;
   (*dx)[0][1] = y[2];
   (*dx)[0][2] = -y[1];
   (*dx)[1][0] = -y[2];
   (*dx)[1][1] = 0.0;
   (*dx)[1][2] = y[0];
   (*dx)[2][0] = y[1];
   (*dx)[2][1] = -y[0];
   (*dx)[2][2] = 0.0;
   (*dy)[0][0] = 0.0;
   (*dy)[0][1] = -x[2];
   (*dy)[0][2] = x[1];
   (*dy)[1][0] = x[2];
   (*dy)[1][1] = 0.0;
   (*dy)[1][2] = -x[0];
   (*dy)[2][0] = -x[1];
   (*dy)[2][1] = x[0];
   (*dy)[2][2] = 0.0;
   //*dx = {{ {{ 0.0,   y[2],  -y[1] }},
   //         {{ -y[2], 0.0,   y[0]  }},
   //         {{ y[1],  -y[0], 0.0   }} }};
   //*dy = {{ {{ 0.0,   -x[2], x[1]  }},
   //         {{ x[2],  0.0,   -x[0] }},
   //         {{ -x[1], x[0],  0.0   }} }};
   array<_T, 3> result = {{
         x[1]*y[2] - y[1]*x[2],
         y[0]*x[2] - x[0]*y[2],
         x[0]*y[1] - y[0]*x[1] }};
   return result;
}

//! Compute the cross product of the two arrays.
template<typename _T>
inline
void
cross(const array<_T, 3>& x, const array<_T, 3>& y, array<_T, 3>* result) {
   (*result)[0] = x[1] * y[2] - y[1] * x[2];
   (*result)[1] = y[0] * x[2] - x[0] * y[2];
   (*result)[2] = x[0] * y[1] - y[0] * x[1];
}

//! Return the triple product of the three arrays.
template<typename _T>
inline
_T
tripleProduct(const array<_T, 3>& x, const array<_T, 3>& y,
              const array<_T, 3>& z) {
   array<_T, 3> t;
   cross(y, z, &t);
   return dot(x, t);
}

//! Return the discriminant of the two arrays.
template<typename _T>
inline
_T
discriminant(const array<_T, 2>& x, const array<_T, 2>& y) {
   return x[0] * y[1] - x[1] * y[0];
}

//! Return the squared magnitude.
template<typename _T, size_t _Dim>
inline
_T
squaredMagnitude(const array<_T, _Dim>& x) {
   return dot(x, x);
}

//! Return the magnitude.
template<typename _T, size_t _Dim>
inline
_T
magnitude(const array<_T, _Dim>& x) {
   return sqrt(squaredMagnitude(x));
}

//! Normalize the vector to have unit magnitude.
template<typename _T, size_t _Dim>
inline
void
normalize(array<_T, _Dim>* x) {
   const _T mag = magnitude(*x);
   if (mag != 0) {
      *x /= mag;
   }
   else {
      // If the vector has zero length, choose the unit vector whose first
      // coordinate is 1.
      std::fill(x->begin(), x->end(), 0);
      (*x)[0] = 1;
   }
}

//! Negate the vector.
template<typename _T, size_t _Dim>
inline
void
negateElements(array<_T, _Dim>* x) {
   for (std::size_t n = 0; n != _Dim; ++n) {
      (*x)[n] = - (*x)[n];
   }
}

//! Negate the vector.
template<size_t _Dim>
inline
void
negateElements(array<bool, _Dim>* x) {
   for (std::size_t n = 0; n != _Dim; ++n) {
      (*x)[n] = !(*x)[n];
   }
}

//! Return the squared distance between the two points.
template<typename _T, size_t _Dim>
inline
_T
squaredDistance(const array<_T, _Dim>& x, const array<_T, _Dim>& y) {
   _T d = 0;
   for (std::size_t n = 0; n != _Dim; ++n) {
      d += (x[n] - y[n]) * (x[n] - y[n]);
   }
   return d;
}

//! Return the squared distance between the two points.
/*! This specialization is a little faster than the dimension-general code. */
template<typename _T>
inline
_T
squaredDistance(const array<_T, 3>& x, const array<_T, 3>& y) {
   return (x[0] - y[0]) * (x[0] - y[0]) +
      (x[1] - y[1]) * (x[1] - y[1]) +
      (x[2] - y[2]) * (x[2] - y[2]);
}

//! Return the Euclidean distance between the two points.
/*!
  \note distance would not be a good name for this function because
  std::distance() calculates the distance between iterators.
*/
template<typename _T, size_t _Dim>
inline
_T
euclideanDistance(const array<_T, _Dim>& x, const array<_T, _Dim>& y) {
   return sqrt(squaredDistance(x, y));
}

//@}

#ifdef BOOST
}
#else
}
}
#endif

#endif
