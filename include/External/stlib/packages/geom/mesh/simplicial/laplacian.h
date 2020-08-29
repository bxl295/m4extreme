// -*- C++ -*-

/*!
  \file geom/mesh/simplicial/laplacian.h
  \brief Implements Laplacian smoothing.
*/

#if !defined(__geom_mesh_simplicial_laplacian_h__)
#define __geom_mesh_simplicial_laplacian_h__

#include "SimpMeshRed.h"
#include "geometry.h"

#include "../iss/ISS_SignedDistance.h"

#include <set>

namespace geom {

//-----------------------------------------------------------------------------
/*! \defgroup simplicial_laplacian Laplacian Smoothing
*/
//@{

//! Perform Laplacian smoothing on the specified interior node.
/*!
  \relates SimpMeshRed
*/
template < std::size_t N, std::size_t M, typename T,
         template<class> class Node,
         template<class> class Cell,
         template<class, class> class Cont >
void
applyLaplacianAtNode(SimpMeshRed<N, M, T, Node, Cell, Cont>* mesh,
                     typename SimpMeshRed<N, M, T, Node, Cell, Cont>::
                     NodeIterator node);


//! Perform Laplacian smoothing on the specified interior nodes.
/*!
  \relates SimpMeshRed
*/
template < std::size_t N, std::size_t M, typename T,
         template<class> class Node,
         template<class> class Cell,
         template<class, class> class Cont,
         class NodeIterInIter >
void
applyLaplacian(SimpMeshRed<N, M, T, Node, Cell, Cont>* mesh,
               NodeIterInIter begin, NodeIterInIter end,
               std::size_t numSweeps = 1);


//! Perform a sweep of Laplacian smoothing on the interior nodes.
/*!
  \relates SimpMeshRed
*/
template < std::size_t N, std::size_t M, typename T,
         template<class> class Node,
         template<class> class Cell,
         template<class, class> class Cont >
void
applyLaplacian(SimpMeshRed<N, M, T, Node, Cell, Cont>* mesh,
               std::size_t numSweeps = 1);


//! Perform a sweep of Laplacian smoothing on the boundary nodes.
/*!
  \relates SimpMeshRed
*/
template < std::size_t N, typename T,
         template<class> class Node,
         template<class> class Cell,
         template<class, class> class Cont,
         class BoundaryCondition >
void
applyLaplacian(SimpMeshRed<N, N, T, Node, Cell, Cont>* mesh,
               const BoundaryCondition& condition,
               T minAngle, std::size_t numSweeps);


//! Perform a sweep of Laplacian smoothing on the specified interior nodes.
/*!
  \relates SimpMeshRed
*/
template < std::size_t M, typename T,
         template<class> class Node,
         template<class> class Cell,
         template<class, class> class Cont,
         class LevelSet, class NodeIterInIter >
void
applyLaplacian(SimpMeshRed < M + 1, M, T, Node, Cell, Cont > * mesh,
               const LevelSet& levelSet,
               NodeIterInIter begin, NodeIterInIter end,
               std::size_t numSweeps = 1);

//@}

} // namespace geom

#define __geom_mesh_simplicial_laplacian_ipp__
#include "laplacian.ipp"
#undef __geom_mesh_simplicial_laplacian_ipp__

#endif
