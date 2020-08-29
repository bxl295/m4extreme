// -*- C++ -*-

/*!
  \file cuda/bound.h
  \brief Forward declarations for objects that can be bounded.
*/

#if !defined(__cuda_bound_h__)
#define __cuda_bound_h__

#include "../ext/array.h"

namespace geom {

template<typename _T, std::size_t _N>
struct BBox;

// array

template<typename _T, std::size_t _N>
BBox<_T, _N>
bound(const std::tr1::array<_T, _N>& p);

template<typename _T, std::size_t _N>
BBox<_T, _N>
bound(const std::tr1::array<_T, _N>& p, const std::tr1::array<_T, _N>& q);

template<typename _T, std::size_t _N>
inline
BBox<_T, _N>
bound(const std::tr1::array<_T, _N>& p, const std::tr1::array<_T, _N>& q,
      const std::tr1::array<_T, _N>& r);

// BBox

template<typename _T, std::size_t _N>
const BBox<_T, _N>&
bound(const BBox<_T, _N>& x);


// Ball

template<typename _T, std::size_t _N>
class Ball;

template<typename _T, std::size_t _N>
BBox<_T, _N>
bound(const Ball<_T, _N>& ball);

} // namespace geom

#endif
