// -*- C++ -*-

#if !defined(__geom_mesh_simplicial_laplacian_ipp__)
#error This file is an implementation detail.
#endif

namespace geom {


template < std::size_t N, std::size_t M, typename T,
         template<class> class Node,
         template<class> class Cell,
         template<class, class> class Cont >
inline
void
applyLaplacianAtNode(SimpMeshRed<N, M, T, Node, Cell, Cont>* mesh,
                     typename SimpMeshRed<N, M, T, Node, Cell, Cont>::
                     NodeIterator node) {
   typedef SimpMeshRed<N, M, T, Node, Cell, Cont> SMR;
   typedef typename SMR::Vertex Vertex;
   typedef typename SMR::NodeIteratorSet NodeIteratorSet;

   NodeIteratorSet neighbors;
   // Initialize to zero.
   Vertex pt = {};

   if (! node->isOnBoundary()) {
      // Get the neighboring nodes.
      determineNeighbors(*mesh, node, &neighbors);
      assert(neighbors.size() != 0);
      // For each neighbor.
      for (typename NodeIteratorSet::const_iterator iter = neighbors.begin();
            iter != neighbors.end(); ++iter) {
         pt += (*iter)->getVertex();
      }
      pt /= T(neighbors.size());
      node->setVertex(pt);
   }
}


// Apply Laplacian smoothing to a boundary node.
template < std::size_t N, typename T,
         template<class> class Node,
         template<class> class Cell,
         template<class, class> class Cont >
inline
void
applyLaplacianOnBoundary(SimpMeshRed<N, N, T, Node, Cell, Cont>* mesh,
                         typename SimpMeshRed<N, N, T, Node, Cell, Cont>::
                         NodeIterator node) {
   typedef SimpMeshRed<N, N, T, Node, Cell, Cont> SMR;
   typedef typename SMR::Vertex Vertex;
   typedef typename SMR::NodeIteratorSet NodeIteratorSet;

   NodeIteratorSet neighbors;
   // Initialize to zero.
   Vertex pt = {};

#ifdef DEBUG_stlib
   assert(node->isOnBoundary());
#endif

   // Get the neighboring boundary nodes.
   determineBoundaryNeighbors(*mesh, node, &neighbors);
   assert(neighbors.size() != 0);
   if (N == 2) {
      // CONTINUE
      //std::cerr << "Number of neighbors = " << neighbors.size() << "\n";
      /*
      if (neighbors.size() != 2) {
        printQualityStatistics(std::cerr, *mesh);
        std::ofstream file("problem.vtu");
        writeVtkXml(file, *mesh);
      }
      */
      assert(neighbors.size() == 2);
   }
   // For each boundary neighbor.
   for (typename NodeIteratorSet::const_iterator iter = neighbors.begin();
         iter != neighbors.end(); ++iter) {
      pt += (*iter)->getVertex();
   }
   pt /= T(neighbors.size());
   // Set the location of the node.
   node->setVertex(pt);
}



// Apply Laplacian smoothing subject to a condition to a node.
template < std::size_t N, std::size_t M, typename T,
         template<class> class Node,
         template<class> class Cell,
         template<class, class> class Cont,
         class ISS >
inline
void
applyLaplacianOnBoundary(SimpMeshRed<N, M, T, Node, Cell, Cont>* mesh,
                         const ISS_SD_ClosestPoint<ISS>& condition,
                         typename SimpMeshRed<N, M, T, Node, Cell, Cont>::
                         NodeIterator node) {
   // Laplacian smoothing.
   applyLaplacianOnBoundary(mesh, node);
   // Apply the boundary condition.
   node->setVertex(condition(node->getVertex()));
}



// Apply Laplacian smoothing subject to a condition to a node.
template < std::size_t N, std::size_t M, typename T,
         template<class> class Node,
         template<class> class Cell,
         template<class, class> class Cont,
         class ISS >
inline
void
applyLaplacianOnBoundary(SimpMeshRed<N, M, T, Node, Cell, Cont>* mesh,
                         const ISS_SD_ClosestPointDirection<ISS>& condition,
                         typename SimpMeshRed<N, M, T, Node, Cell, Cont>::
                         NodeIterator node) {
   typedef SimpMeshRed<N, M, T, Node, Cell, Cont> SMR;

   // Laplacian smoothing.
   applyLaplacianOnBoundary(mesh, node);
   // Apply the boundary condition.
   node->setVertex(condition(node->getVertex(), computeNodeNormal<SMR>(node)));
}



// Apply Laplacian smoothing subject to a condition to a node.
template < std::size_t N, std::size_t M, typename T,
         template<class> class Node,
         template<class> class Cell,
         template<class, class> class Cont,
         class ISS >
inline
void
applyLaplacianAtNode(SimpMeshRed<N, M, T, Node, Cell, Cont>* mesh,
                     const ISS_SD_ClosestPoint<ISS>& condition,
                     typename SimpMeshRed<N, M, T, Node, Cell, Cont>::
                     NodeIterator node) {
   // Laplacian smoothing.
   applyLaplacianAtNode(mesh, node);
   // Apply the condition.
   node->setVertex(condition(node->getVertex()));
}



// Apply Laplacian smoothing subject to a condition to a node.
template < std::size_t N, std::size_t M, typename T,
         template<class> class Node,
         template<class> class Cell,
         template<class, class> class Cont,
         class ISS >
inline
void
applyLaplacianAtNode(SimpMeshRed<N, M, T, Node, Cell, Cont>* mesh,
                     const ISS_SD_ClosestPointDirection<ISS>& condition,
                     typename SimpMeshRed<N, M, T, Node, Cell, Cont>::
                     NodeIterator node) {
   typedef SimpMeshRed<N, M, T, Node, Cell, Cont> SMR;

   // Laplacian smoothing.
   applyLaplacianAtNode(mesh, node);
   // Apply the condition.
   node->setVertex(condition(node->getVertex(), computeNodeNormal<SMR>(node)));
}



template < std::size_t N, std::size_t M, typename T,
         template<class> class Node,
         template<class> class Cell,
         template<class, class> class Cont,
         class NodeIterInIter >
inline
void
applyLaplacian(SimpMeshRed<N, M, T, Node, Cell, Cont>* mesh,
               NodeIterInIter begin, NodeIterInIter end,
               std::size_t numSweeps) {
   assert(numSweeps >= 0);

   while (numSweeps-- != 0) {
      for (; begin != end; ++begin) {
         applyLaplacianAtNode(mesh, *begin);
      }
   }
}



template < std::size_t N, std::size_t M, typename T,
         template<class> class Node,
         template<class> class Cell,
         template<class, class> class Cont >
inline
void
applyLaplacian(SimpMeshRed<N, M, T, Node, Cell, Cont>* mesh, std::size_t numSweeps) {
   typedef SimpMeshRed<N, M, T, Node, Cell, Cont> SMR;
   typedef typename SMR::NodeIterator NodeIterator;

   assert(numSweeps >= 0);

   while (numSweeps-- != 0) {
      for (NodeIterator node = mesh->getNodesBeginning();
            node != mesh->getNodesEnd(); ++node) {
         applyLaplacianAtNode(mesh, node);
      }
   }
}


// Perform a sweep of Laplacian smoothing on the boundary nodes.
template < std::size_t N, typename T,
         template<class> class Node,
         template<class> class Cell,
         template<class, class> class Cont,
         class BoundaryCondition >
inline
void
applyLaplacian(SimpMeshRed<N, N, T, Node, Cell, Cont>* mesh,
               const BoundaryCondition& condition,
               const T minAngle, std::size_t numSweeps) {
   typedef SimpMeshRed<N, N, T, Node, Cell, Cont> SMR;
   typedef typename SMR::NodeIterator NodeIterator;

   assert(numSweeps >= 0);

   // Make the functor for determining which nodes can be moved.
   geom::IsNotSharpAngle<SMR> isMovable(minAngle);

   // Determine which nodes are on the boundary and can be moved.
   std::vector<NodeIterator> boundary;
   {
      // Determine which nodes are on the boundary.
      std::vector<NodeIterator> wholeBoundary;
      determineBoundaryNodes(*mesh, std::back_inserter(wholeBoundary));
      // Select the movable boundary nodes.
      for (typename std::vector<NodeIterator>::const_iterator
            i = wholeBoundary.begin(); i != wholeBoundary.end(); ++i) {
         if (isMovable(*i)) {
            boundary.push_back(*i);
         }
      }
   }

   while (numSweeps-- != 0) {
      // For each boundary node.
      for (typename std::vector<NodeIterator>::const_iterator
            i = boundary.begin(); i != boundary.end(); ++i) {
         // Laplacian smoothing subject to the boundary condition.
         applyLaplacianOnBoundary(mesh, condition, *i);
      }
   }
}



template < std::size_t M, typename T,
         template<class> class Node,
         template<class> class Cell,
         template<class, class> class Cont,
         class LevelSet, class NodeIterInIter >
inline
void
applyLaplacian(SimpMeshRed < M + 1, M, T, Node, Cell, Cont > * mesh,
               const LevelSet& levelSet,
               NodeIterInIter begin, NodeIterInIter end, std::size_t numSweeps) {
   assert(numSweeps >= 0);

   while (numSweeps-- != 0) {
      for (; begin != end; ++begin) {
         applyLaplacianAtNode(mesh, levelSet, *begin);
      }
   }
}

} // namespace geom
