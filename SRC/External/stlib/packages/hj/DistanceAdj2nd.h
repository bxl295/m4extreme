// -*- C++ -*-

/*!
  \file DistanceAdj2nd.h
  \brief Distance equation.  Second-order, adjacent scheme.
*/

#if !defined(__hj_DistanceAdj2nd_h__)
#define __hj_DistanceAdj2nd_h__

#include "Distance.h"
#include "DistanceScheme.h"

#ifdef DEBUG_stlib
// Include the debugging code.
#include "debug.h"
#endif

namespace hj {

//! Distance equation.  Adjacent difference scheme.  2nd order.
/*!
  \param N is the space dimension.
  \param T is the number type.
*/
template<std::size_t N, typename T>
class DistanceAdj2nd;

} // namespace hj

#define __hj_DistanceAdj2nd2_h__
#include "DistanceAdj2nd2.h"
#undef __hj_DistanceAdj2nd2_h__

//#define __hj_DistanceAdj2nd3_h__
//#include "DistanceAdj2nd3.h"
//#undef __hj_DistanceAdj2nd3_h__

#endif
