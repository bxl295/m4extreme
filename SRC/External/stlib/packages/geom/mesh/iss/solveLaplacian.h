// -*- C++ -*-

/*!
  \file geom/mesh/iss/solveLaplacian.h
  \brief Implements Laplacian solve.
*/

#if !defined(__geom_mesh_iss_solveLaplacian_h__)
#define __geom_mesh_iss_solveLaplacian_h__

#include "laplacian.h"

#include "../../../third-party/tnt/tnt.h"
#include "../../../third-party/jama/jama_lu.h"

namespace geom {

//-----------------------------------------------------------------------------
/*! \defgroup iss_solveLaplacian Laplacian Smoothing with a Solve
*/
//@{


//! Perform Laplacian smoothing on the interior vertices.
/*!
  \relates IndSimpSetIncAdj

  \param mesh Pointer to the simplicial mesh.

  \note This function is for testing purposes.  Do not use this function
  for large meshes.  This makes an n x n matrix, where n is the number
  of interior nodes.
*/
template<std::size_t N, std::size_t M, typename T>
void
solveLaplacian(IndSimpSetIncAdj<N, M, T>* mesh);


//@}

} // namespace geom

#define __geom_mesh_iss_solveLaplacian_ipp__
#include "solveLaplacian.ipp"
#undef __geom_mesh_iss_solveLaplacian_ipp__

#endif
