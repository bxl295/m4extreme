// -*- C++ -*-

/*!
  \file bound.h
  \brief Forward declarations for objects that can be bounded.
*/

#if !defined(__geom_kernel_bound_h__)
#define __geom_kernel_bound_h__

#include <tr1/array>

// CONTINUE REMOVE
#include <cassert>

namespace geom {


template<typename _T, std::size_t _Dim>
struct BBox;

// std::tr1::array

template<typename _T, std::size_t _Dim>
BBox<_T, _Dim>
bound(const std::tr1::array<_T, _Dim>& p);

template<typename _T, std::size_t _Dim>
BBox<_T, _Dim>
bound(const std::tr1::array<_T, _Dim>& p, const std::tr1::array<_T, _Dim>& q);

template<typename _T, std::size_t _Dim>
inline
BBox<_T, _Dim>
bound(const std::tr1::array<_T, _Dim>& p, const std::tr1::array<_T, _Dim>& q,
      const std::tr1::array<_T, _Dim>& r);

// BBox

template<typename _T, std::size_t _Dim>
const BBox<_T, _Dim>&
bound(const BBox<_T, _Dim>& x);


// CONTINUE: REMOVE
#if 0
BBox<int, 3>
bound(const BBox<int, 3>& x);
#endif


// Ball

template<typename _T, std::size_t _Dim>
class Ball;

template<typename _T, std::size_t _Dim>
BBox<_T, _Dim>
bound(const Ball<_T, _Dim>& ball);


} // namespace geom

#endif
