// -*- C++ -*-

/*!
  \file EikonalAdj2nd.h
  \brief Eikonal equation.  Second-order, adjacent scheme.
*/

#if !defined(__hj_EikonalAdj2nd_h__)
#define __hj_EikonalAdj2nd_h__

#include "Eikonal.h"
#include "EikonalScheme.h"

#ifdef DEBUG_stlib
// Include the debugging code.
#include "debug.h"
#endif

namespace hj {

//! Eikonal equation.  Adjacent difference scheme.  2nd order.
/*!
  \param N is the space dimension.
  \param T is the number type.
*/
template<std::size_t N, typename T>
class EikonalAdj2nd;

} // namespace hj

#define __hj_EikonalAdj2nd2_h__
#include "EikonalAdj2nd2.h"
#undef __hj_EikonalAdj2nd2_h__

//#define __hj_EikonalAdj2nd3_h__
//#include "EikonalAdj2nd3.h"
//#undef __hj_EikonalAdj2nd3_h__

#endif
