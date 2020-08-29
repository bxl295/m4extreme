// -*- C++ -*-

#if !defined(__ext_array_h__)
#define __ext_array_h__

#ifdef __CUDA_ARCH__
#include "arrayCuda.h"
#else
#include "arrayStd.h"
#endif

#include "../cuda/compatible.h"

namespace ext {
//-------------------------------------------------------------------------
/*! \defgroup extArrayMake Convenience Constructor Functions.

  These functions are useful for constructing std::tr1::array's.
  \verbatim
  typedef std::tr1::array<double, 3> Point;
  Point x;
  ...
  // Convert to an array with a different value type.
  std::tr1::array<int, 3> a = ext::convert_array<int>(x);
  // Make an array filled with a specified value.
  x = ext::filled_array<Point>(0.);
  // Copy from a different kind of array.
  double y[3] = {2, 3, 5};
  x = ext::copy_array(y);
  // Make an array by specifying the elements.
  std::tr1::array<double, 1> x1;
  x1 = ext::make_array(2);
  std::tr1::array<double, 2> x2;
  x2 = ext::make_array(2, 3);
  std::tr1::array<double, 3> x3;
  x3 = ext::make_array(2, 3, 5);
  std::tr1::array<double, 4> x4;
  x4 = ext::make_array(2, 3, 5, 7); \endverbatim
*/
//@{


//! Convert to an array with different value type.
template<typename _Target, typename _Source, std::size_t _Dim>
inline
CUDA_COMPATIBLE
std::tr1::array<_Target, _Dim>
convert_array(const std::tr1::array<_Source, _Dim>& x) {
   std::tr1::array<_Target, _Dim> result;
   for (std::size_t i = 0; i != _Dim; ++i) {
      result[i] = x[i];
   }
   return result;
}

//! Return an array filled with the specified value.
template<typename _Array>
inline
CUDA_COMPATIBLE
_Array
filled_array(const typename _Array::value_type& value) {
   _Array x;
   for (std::size_t i = 0; i != x.size(); ++i) {
      x[i] = value;
   }
   return x;
}

//! Copy from the input iterator to make an array.
template<typename _Array, typename _InputIterator>
inline
CUDA_COMPATIBLE
_Array
copy_array(_InputIterator input) {
   _Array x;
   for (std::size_t i = 0; i != x.size(); ++i) {
      x[i] = *input++;
   }
   return x;
}

//! Return an array with the specified components.
template<typename _T>
inline
CUDA_COMPATIBLE
std::tr1::array<_T, 1>
make_array(const _T& x0) {
   std::tr1::array<_T, 1> x = {{x0}};
   return x;
}

//! Return an array with the specified components.
template<typename _T>
inline
CUDA_COMPATIBLE
std::tr1::array<_T, 2>
make_array(const _T& x0, const _T& x1) {
   std::tr1::array<_T, 2> x = {{x0, x1}};
   return x;
}

//! Return an array with the specified components.
template<typename _T>
inline
CUDA_COMPATIBLE
std::tr1::array<_T, 3>
make_array(const _T& x0, const _T& x1, const _T& x2) {
   std::tr1::array<_T, 3> x = {{x0, x1, x2}};
   return x;
}

//! Return an array with the specified components.
template<typename _T>
inline
CUDA_COMPATIBLE
std::tr1::array<_T, 4>
make_array(const _T& x0, const _T& x1, const _T& x2, const _T& x3) {
   std::tr1::array<_T, 4> x = {{x0, x1, x2, x3}};
   return x;
}

//! Return an array with the specified components.
template<typename _T>
inline
CUDA_COMPATIBLE
std::tr1::array<_T, 5>
make_array(const _T& x0, const _T& x1, const _T& x2, const _T& x3,
const _T& x4) {
   std::tr1::array<_T, 5> x = {{x0, x1, x2, x3, x4}};
   return x;
}

//! Return an array with the specified components.
template<typename _T>
inline
CUDA_COMPATIBLE
std::tr1::array<_T, 6>
make_array(const _T& x0, const _T& x1, const _T& x2, const _T& x3,
const _T& x4, const _T& x5) {
   std::tr1::array<_T, 6> x = {{x0, x1, x2, x3, x4, x5}};
   return x;
}

//! Return an array with the specified components.
template<typename _T>
inline
CUDA_COMPATIBLE
std::tr1::array<_T, 7>
make_array(const _T& x0, const _T& x1, const _T& x2, const _T& x3,
const _T& x4, const _T& x5, const _T& x6) {
   std::tr1::array<_T, 7> x = {{x0, x1, x2, x3, x4, x5, x6}};
   return x;
}

//! Return an array with the specified components.
template<typename _T>
inline
CUDA_COMPATIBLE
std::tr1::array<_T, 8>
make_array(const _T& x0, const _T& x1, const _T& x2, const _T& x3,
const _T& x4, const _T& x5, const _T& x6, const _T& x7) {
   std::tr1::array<_T, 8> x = {{x0, x1, x2, x3, x4, x5, x6, x7}};
   return x;
}

//! Return an array with the specified components.
template<typename _T>
inline
CUDA_COMPATIBLE
std::tr1::array<_T, 9>
make_array(const _T& x0, const _T& x1, const _T& x2, const _T& x3,
const _T& x4, const _T& x5, const _T& x6, const _T& x7,
const _T& x8) {
   std::tr1::array<_T, 9> x = {{x0, x1, x2, x3, x4, x5, x6, x7, x8}};
   return x;
}

//! Return an array with the specified components.
template<typename _T>
inline
CUDA_COMPATIBLE
std::tr1::array<_T, 10>
make_array(const _T& x0, const _T& x1, const _T& x2, const _T& x3,
const _T& x4, const _T& x5, const _T& x6, const _T& x7,
const _T& x8, const _T& x9) {
   std::tr1::array<_T, 10> x = {{x0, x1, x2, x3, x4, x5, x6, x7, x8, x9}};
   return x;
}

//@}

} // namespace ext


#endif
