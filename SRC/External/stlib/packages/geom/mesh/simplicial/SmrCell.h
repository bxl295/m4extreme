// -*- C++ -*-

/*!
  \file SmrCell.h
  \brief Class for a cell in a simplicial mesh that stores the handles to the adjacent cells.
*/

#if !defined(__geom_mesh_simplicial_SmrCell_h__)
#define __geom_mesh_simplicial_SmrCell_h__

#include "../simplex/topology.h"
#include "../../../loki/static_check.h"

namespace geom {

//! Cell in a simplicial mesh that stores handles to the adjacent cells.
/*!
  \param SMR is the simplicial mesh data structure.
*/
template<class SMR>
class SmrCell {
   //
   // Enumerations.
   //

public:

   //! The simplex dimension.
   enum {M = SMR::M};

   //
   // Public types.
   //

public:

   //! The simplicial mesh.
   typedef SMR Mesh;

   //! An iterator to a cell.
   typedef typename Mesh::CellIterator CellIterator;
   //! An iterator to a const cell.
   typedef typename Mesh::CellConstIterator CellConstIterator;

   //! An iterator to a node.
   typedef typename Mesh::NodeIterator NodeIterator;
   //! An iterator to a const node.
   typedef typename Mesh::NodeConstIterator NodeConstIterator;
   //! The vertex (a Cartesian point).
   typedef typename Mesh::Vertex Vertex;
   //! The number type.
   typedef typename Mesh::Number Number;

   //! The simplex of node iterators.
   typedef std::tr1::array < NodeIterator, M + 1 > NodeIteratorSimplex;
   //! The simplex of cell iterators.
   typedef std::tr1::array < CellIterator, M + 1 > CellIteratorSimplex;

   //
   // Private types.
   //

private:

   //
   // Public types.
   //

public:

   //! A face of the cell is a simplex of \c M vertex iterators.
   typedef std::tr1::array<NodeIterator, M> Face;

   //
   // Data
   //

private:

   //! The incident nodes.
   NodeIteratorSimplex _nodes;
   //! The adjacent cells.
   CellIteratorSimplex _neighbors;
   //! The identifier of this simplex.
   mutable std::size_t _identifier;
   //! An iterator to this cell.
   CellIterator _self;

public:

   //--------------------------------------------------------------------------
   //! \name Constructors and Destructor.
   //! @{

   //! Default constructor.  Null iterators.
   SmrCell() :
      _nodes(),
      _neighbors(),
      _identifier(-1),
      _self(0) {
      std::fill(_nodes.begin(), _nodes.end(), NodeIterator(0));
      std::fill(_neighbors.begin(), _neighbors.end(), CellIterator(0));
   }

   //! Construct a 1-D cell from the vertices and neighbors.
   SmrCell(const NodeIterator v0, const NodeIterator v1,
           const CellIterator c0 = CellIterator(0),
           const CellIterator c1 = CellIterator(0),
           const std::size_t identifier = -1) :
      _nodes(ext::make_array(v0, v1)),
      _neighbors(ext::make_array(c0, c1)),
      _identifier(identifier),
      _self(0) {
      LOKI_STATIC_CHECK(M == 1, TheSimplexDimensionMustBe1);
   }

   //! Construct a 2-D cell from the vertices and neighbors.
   SmrCell(const NodeIterator v0, const NodeIterator v1,
           const NodeIterator v2,
           const CellIterator c0 = CellIterator(0),
           const CellIterator c1 = CellIterator(0),
           const CellIterator c2 = CellIterator(0),
           const std::size_t identifier = -1) :
      _nodes(ext::make_array(v0, v1, v2)),
      _neighbors(ext::make_array(c0, c1, c2)),
      _identifier(identifier),
      _self(0) {
      LOKI_STATIC_CHECK(M == 2, TheSimplexDimensionMustBe2);
   }

   //! Construct a 3-D cell from the vertices and neighbors.
   SmrCell(const NodeIterator v0, const NodeIterator v1,
           const NodeIterator v2, const NodeIterator v3,
           const CellIterator c0 = CellIterator(0),
           const CellIterator c1 = CellIterator(0),
           const CellIterator c2 = CellIterator(0),
           const CellIterator c3 = CellIterator(0),
           const std::size_t identifier = -1) :
      _nodes(ext::make_array(v0, v1, v2, v3)),
      _neighbors(ext::make_array(c0, c1, c2, c3)),
      _identifier(identifier),
      _self(0) {
      LOKI_STATIC_CHECK(M == 3, TheSimplexDimensionMustBe3);
   }

   //! Copy constructor.
   SmrCell(const SmrCell& other) :
      _nodes(other._nodes),
      _neighbors(other._neighbors),
      _identifier(other._identifier),
      _self(other._self) {}

   //! Trivial destructor.
   ~SmrCell() {}

   //! @}
   //--------------------------------------------------------------------------
   //! \name Assignment operators.
   //! @{

   //! Assignment operator.
   SmrCell&
   operator=(const SmrCell& other) {
      if (&other != this) {
         _nodes = other._nodes;
         _neighbors = other._neighbors;
         _identifier = other._identifier;
         _self = other._self;
      }
      return *this;
   }

   //! @}
   //--------------------------------------------------------------------------
   //! \name Accessors.
   //! @{

   //! Return the identifier of this simplex.
   /*!
     Typically, the identifier is in the range [0...num_simplices).
     and a value of -1 indicates that the identifier has not been calculated.
   */
   std::size_t
   getIdentifier() const {
      return _identifier;
   }

   //! Return a const iterator to this cell.
   CellConstIterator
   getSelf() const {
      return _self;
   }

   //! Return a const iterator to the m_th node.
   NodeConstIterator
   getNode(const std::size_t m) const {
      return _nodes[m];
   }

   //! Return the simplex of node iterators.
   const NodeIteratorSimplex&
   getNodes() const {
      return _nodes;
   }

   //! Return the index of the specified node.
   std::size_t
   getIndex(NodeConstIterator node) const {
      return index(_nodes, node);
   }

   //! Return the index of the specified node.
   std::size_t
   getIndex(NodeIterator node) const {
      return index(_nodes, node);
   }

   //! Return the index of the specified vertex.
   /*!
     The vertex is specified by a face of the cell.
   */
   std::size_t
   getIndex(const Face& f) const {
      for (std::size_t m = 0; m != M + 1; ++m) {
         if (! hasElement(f, _nodes[m])) {
            return m;
         }
      }
      assert(false);
      return -1;
   }

   //! Return true if the cell has the specified node.
   bool
   hasNode(NodeConstIterator node) const {
      return hasElement(_nodes, node);
   }

   //! Return true if the cell has the specified node.
   /*!
     If true, compute the index of the node.
   */
   bool
   hasNode(NodeConstIterator node, std::size_t* m) const {
      return hasElement(_nodes, node, m);
   }

   //! Return true if this cell has a boundary face that is incident to the specified node.
   bool
   hasIncidentBoundaryFace(NodeConstIterator node) const {
#ifdef DEBUG_stlib
      assert(hasNode(node));
#endif
      return hasIncidentBoundaryFace(getIndex(node));
   }

   //! Return true if this cell has a boundary face that is incident to the specified node.
   bool
   hasIncidentBoundaryFace(const std::size_t n) const {
#ifdef DEBUG_stlib
      assert(n < M + 1);
#endif
      // For each face incident to the node.
      for (std::size_t m = 0; m != M + 1; ++m) {
         // Exclude the face opposite this node.
         if (m != n) {
            // If the face is on the boundary.
            if (isFaceOnBoundary(m)) {
               return true;
            }
         }
      }
      // We did not encounter any incident boundary faces.
      return false;
   }

   //! Calculate the centroid.
   void
   getCentroid(Vertex* centroid) const {
      // CONTINUE: This is the arithmetic mean.  Look up the formula for
      // the centroid.
      *centroid = _nodes[0]->getVertex();
      for (std::size_t m = 1; m != M + 1; ++m) {
         *centroid += _nodes[m]->getVertex();
      }
      *centroid /= Number(M + 1);
   }

   //! Return a const iterator to the n_th neighboring cell.
   CellConstIterator
   getNeighbor(const std::size_t m) const {
      return _neighbors[m];
   }

   //! Return true if the specified face is on the boundary.
   bool
   isFaceOnBoundary(const std::size_t m) const {
      return getNeighbor(m) == CellConstIterator(0);
   }

   //! Return the number of neighbors.
   std::size_t
   getNumberOfNeighbors() const {
      return M + 1 - std::count(_neighbors.begin(), _neighbors.end(),
                                CellConstIterator(0));
   }

   //! Return the index of the specified neighbor.
   std::size_t
   getIndex(CellConstIterator c) const {
      return index(_neighbors, c);
   }

   //! Return true if the cell has the specified neighbor.
   bool
   hasNeighbor(CellConstIterator c) const {
      return hasElement(_neighbors, c);
   }

   //! Return true if the cell has the specified neighbor.
   /*!
     If true, compute the index of the neighbor.
   */
   bool
   hasNeighbor(CellConstIterator c, std::size_t* m) const {
      return hasElement(_neighbors, c, m);
   }

   //! Return the index of this cell in the m_th neighbor.
   std::size_t
   getMirrorIndex(const std::size_t m) const {
      if (isFaceOnBoundary(m)) {
         return -1;
      }
      return _neighbors[m]->getIndex(_self);
   }

   //! Return the minimum edge length.
   /*!
     Set \c a and \c b to the vertex indices that define the minimum edge.
   */
   Number
   computeMinimumEdgeLength(std::size_t* a, std::size_t* b) const {
      Number x = std::numeric_limits<Number>::max();
      Number d;
      std::size_t j;
      // For each edge.
      for (std::size_t i = 0; i != M; ++i) {
         for (j = i + 1; j != M + 1; ++j) {
            d = euclideanDistance(_nodes[i]->getVertex(), _nodes[j]->getVertex());
            if (d < x) {
               x = d;
               *a = i;
               *b = j;
            }
         }
      }
      return x;
   }

   //! Return the maximum edge length.
   /*!
     Set \c a and \c b to the vertex indices that define the maximum edge.
   */
   Number
   computeMaximumEdgeLength(std::size_t* a, std::size_t* b) const {
      Number x = 0;
      Number d;
      std::size_t j;
      // For each edge.
      for (std::size_t i = 0; i != M; ++i) {
         for (j = i + 1; j != M + 1; ++j) {
            d = euclideanDistance(_nodes[i]->getVertex(), _nodes[j]->getVertex());
            if (d > x) {
               x = d;
               *a = i;
               *b = j;
            }
         }
      }
      return x;
   }

   //! Return the minimum edge length.
   Number
   computeMinimumEdgeLength() const {
      std::size_t a, b;
      return computeMinimumEdgeLength(&a, &b);
   }

   //! Return the maximum edge length.
   Number
   computeMaximumEdgeLength() const {
      std::size_t a, b;
      return computeMaximumEdgeLength(&a, &b);
   }

   //! Get the (vertex) simplex that comprises the cell.
   void
   getSimplex(std::tr1::array < Vertex, M + 1 > * simplex) {
      for (std::size_t m = 0; m != M + 1; ++m) {
         (*simplex)[m] = getNode(m)->getVertex();
      }
   }

   //! @}
   //--------------------------------------------------------------------------
   //! \name Manipulators.
   //! @{

   //! Set the identifier.
   /*!
     \note This member function is const because the identifier is mutable.
     It is intended to be modified as the mesh changes.
   */
   void
   setIdentifier(const std::size_t identifier) const {
      _identifier = identifier;
   }

   //! Return an iterator to this cell.
   CellIterator
   getSelf() {
      return _self;
   }

   //! Set the self iterator for this cell.
   void
   setSelf(const CellIterator self) {
      _self = self;
   }

   //! Return an iterator to the m_th node.
   NodeIterator
   getNode(const std::size_t m) {
      return _nodes[m];
   }

   //! Set the m_th node.
   void
   setNode(const std::size_t m, const NodeIterator node) {
      _nodes[m] = node;
   }

   //! Return an iterator to the m_th neighboring cell.
   CellIterator
   getNeighbor(const std::size_t m) {
      return _neighbors[m];
   }

   //! Set the m_th neighbor.
   void
   setNeighbor(const std::size_t m, const CellIterator c) {
      _neighbors[m] = c;
   }

   //! Return the vertex of the m_th neighbor that is opposite this cell.
   NodeIterator
   getMirrorVertex(const std::size_t m) {
      if (_neighbors[m] == CellIterator(0)) {
         return NodeIterator(0);
      }
      return _neighbors[m]->getNode(getMirrorIndex(m));
   }

   //! Get the face opposite the m_th vertex.
   void
   getFace(const std::size_t m, Face* f) {
      geom::getFace(_nodes, m, f);
   }

   //! Unlink this cell from the mesh.
   void
   unlink() {
      // For each vertex.
      for (typename NodeIteratorSimplex::iterator i = _nodes.begin();
            i != _nodes.end(); ++i) {
         // Remove the vertices' link to this cell.
         (*i)->removeCell(getSelf());
         // Remove the link to the vertex.
         *i = NodeIterator(0);
      }
      // For each neigbor.
      for (typename CellIteratorSimplex::iterator i = _neighbors.begin();
            i != _neighbors.end(); ++i) {
         // If there is a neighbor.
         if (*i != CellIterator(0)) {
            // Remove the neighbor's link to this cell.
            (*i)->removeNeighbor(getSelf());
            // Remove the link to the neighbor.
            *i = CellIterator(0);
         }
      }
   }

   //! Remove the link to the specified neighbor.
   /*!
     The shared face becomes a boundary face.
   */
   void
   removeNeighbor(const CellIterator c) {
      // For each neigbor.
      typename CellIteratorSimplex::iterator i = _neighbors.begin();
      for (; i != _neighbors.end(); ++i) {
         // If this is the specified neighbor.
         if (*i == c) {
            // Remove the link to the neighbor.
            *i = CellIterator(0);
            break;
         }
      }
      // The specified cell must be a neighbor.
      assert(i != _neighbors.end());
   }

   //! Reverse the orientation of this cell.
   void
   reverseOrientation() {
      geom::reverseOrientation(&_nodes);
      geom::reverseOrientation(&_neighbors);
   }

   //! @}
   //--------------------------------------------------------------------------
   //! \name File I/O.
   //! @{

   //! Write the identifier, vertex identifiers and neighbor simplex identifiers.
   void
   put(std::ostream& out) const {
      // Write the simplex identifier.
      out << "Id = " << getIdentifier() << " Vertices = ";
      // Write the incident nodes.
      for (std::size_t m = 0; m != M + 1; ++m) {
         assert(_nodes[m] != NodeIterator(0));
         out << _nodes[m]->getIdentifier() << " ";
      }
      out << " Neighbors = ";
      // Write the neighbors.
      for (std::size_t m = 0; m != M + 1; ++m) {
         if (_neighbors[m] != CellIterator(0)) {
            out << _neighbors[m]->getIdentifier() << " ";
         }
         else {
            out << "-1 ";
         }
      }
      out << "\n";
   }

   //! @}
};


//! Return true if the cell has a boundary face incident to the edge.
/*!
  \pre SMR::M == 3.

  \relates SmrCell
*/
template<class SMR>
bool
doesCellHaveIncidentFaceOnBoundary(const typename SMR::CellConstIterator& c,
                                   std::size_t i, std::size_t j);


//! Return true if the edge is on the boundary.
/*!
  \pre SMR::M == 3.

  \c i and \c j are the indices of the edge in the cell.

  An edge is on the boundary iff an incident face is on the boundary.

  \relates SmrCell
*/
template<class SMR>
bool
isOnBoundary(const typename SMR::CellConstIterator& c, std::size_t i,
             std::size_t j);


//! Return true if the edge is on the boundary.
/*!
  \pre SMR::M == 3.

  \c i and \c j are the indices of the edge in the cell.

  An edge is on the boundary iff an incident face is on the boundary.

  \relates SmrCell
*/
template<class SMR>
inline
bool
isOnBoundary(const typename SMR::Edge& edge) {
   return isOnBoundary<SMR>(edge.first, edge.second, edge.third);
}


//! Return true if the cell has the specified face.  Set the face index.
/*!
  \relates SmrCell
*/
template<class SMR>
inline
bool
hasFace(typename SMR::CellConstIterator cell,
        typename SMR::Cell::Face& face,
        std::size_t* faceIndex) {
   return hasFace(cell->getNodes(), face, faceIndex);
}


//! Return true if the cell has the specified face.  Set the face index.
/*!
  \relates SmrCell
*/
template<class SMR>
inline
bool
hasFace(typename SMR::CellConstIterator cell,
        typename SMR::Cell::Face& face) {
   std::size_t faceIndex;
   return hasFace<SMR>(cell, face, &faceIndex);
}


//! Return true if the cell has the specified face.  Set the face index.
/*!
  \pre SMR::M == 3.

  \relates SmrCell
*/
template<class SMR>
inline
bool
hasFace(typename SMR::CellConstIterator cell,
        typename SMR::NodeConstIterator a,
        typename SMR::NodeConstIterator b,
        typename SMR::NodeConstIterator c,
        std::size_t* faceIndex) {
   return hasFace(cell->getNodes(), a, b, c, faceIndex);
}


//! Return true if the cell has the specified face.
/*!
  \pre SMR::M == 3.

  \relates SmrCell
*/
template<class SMR>
inline
bool
hasFace(typename SMR::CellConstIterator cell,
        typename SMR::NodeConstIterator a,
        typename SMR::NodeConstIterator b,
        typename SMR::NodeConstIterator c) {
   std::size_t dummy;
   return hasFace<SMR>(cell, a, b, c, &dummy);
}


//! Return true if the cell has the specified face.  Set the face index.
/*!
  \pre SMR::M == 3.

  \relates SmrCell
*/
template<class SMR>
inline
bool
hasFace(typename SMR::CellIterator cell,
        typename SMR::NodeIterator a,
        typename SMR::NodeIterator b,
        typename SMR::NodeIterator c,
        std::size_t* faceIndex) {
   return hasFace(cell->getNodes(), a, b, c, faceIndex);
}


//! Return true if the cell has the specified face.
/*!
  \pre SMR::M == 3.

  \relates SmrCell
*/
template<class SMR>
inline
bool
hasFace(typename SMR::CellIterator cell,
        typename SMR::NodeIterator a,
        typename SMR::NodeIterator b,
        typename SMR::NodeIterator c) {
   std::size_t dummy;
   return hasFace<SMR>(cell, a, b, c, &dummy);
}


// Below are the even and odd permutations of the nodes.  We use even
// permutations to get the next node and odd permutations to get the
// previous node.
//
// Even Odd
// 0123 0132
// 0231 0213
// 0312 0321
// 1032 1032
// 1203 1230
// 1320 1302
// 2013 2031
// 2130 2103
// 2301 2310
// 3021 3012
// 3102 3120
// 3210 3201


// CONTINUE: Make this a singleton.  Static data may be a bad idea.
// Also, I should move this to Simplex.
//! Get the next node index.
inline
std::size_t
getNextNodeIndex(const std::size_t i, const std::size_t j) {
   const std::size_t Size = 12;
   static std::size_t even[Size][3] = {{0, 1, 2},
      {0, 2, 3},
      {0, 3, 1},
      {1, 0, 3},
      {1, 2, 0},
      {1, 3, 2},
      {2, 0, 1},
      {2, 1, 3},
      {2, 3, 0},
      {3, 0, 2},
      {3, 1, 0},
      {3, 2, 1}
   };
   for (std::size_t n = 0; n != Size; ++n) {
      if (even[n][0] == i && even[n][1] == j) {
         return even[n][2];
      }
   }
   // Here we check that i and j have sensible values.
   assert(false);
   return std::numeric_limits<std::size_t>::max();
}


// CONTINUE: Make this a singleton.  Static data may be a bad idea.
//! Get the previous node index.
inline
std::size_t
getPreviousNodeIndex(const std::size_t i, const std::size_t j) {
   const std::size_t Size = 12;
   static std::size_t odd[Size][3] = {{0, 1, 3},
      {0, 2, 1},
      {0, 3, 2},
      {1, 0, 2},
      {1, 2, 3},
      {1, 3, 0},
      {2, 0, 3},
      {2, 1, 0},
      {2, 3, 1},
      {3, 0, 1},
      {3, 1, 2},
      {3, 2, 0}
   };
   for (std::size_t n = 0; n != Size; ++n) {
      if (odd[n][0] == i && odd[n][1] == j) {
         return odd[n][2];
      }
   }
   // Here we check that i and j have sensible values.
   assert(false);
   return std::numeric_limits<std::size_t>::max();
}


//! Get the next node index.
template <typename SMR>
inline
std::size_t
getNextNodeIndex(const typename SMR::CellConstIterator cell,
                 const typename SMR::NodeConstIterator a,
                 const typename SMR::NodeConstIterator b) {
   assert(cell->hasNode(a) && cell->hasNode(b));
   return getNextNodeIndex(cell->getIndex(a), cell->getIndex(b));
}


//! Get the previous node index.
template <typename SMR>
inline
std::size_t
getPreviousNodeIndex(const typename SMR::CellConstIterator cell,
                     const typename SMR::NodeConstIterator a,
                     const typename SMR::NodeConstIterator b) {
   assert(cell->hasNode(a) && cell->hasNode(b));
   return getPreviousNodeIndex(cell->getIndex(a), cell->getIndex(b));
}


} // namespace geom

#define __geom_mesh_simplicial_SmrCell_ipp__
#include "SmrCell.ipp"
#undef __geom_mesh_simplicial_SmrCell_ipp__

#endif
