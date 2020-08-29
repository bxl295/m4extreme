// -*- C++ -*-

/*!
  \file IndSimpSet.h
  \brief Implements a mesh that stores vertices and indexed simplices.
*/

#if !defined(__geom_IndSimpSet_h__)
#define __geom_IndSimpSet_h__

#include "SimplexIterator.h"

#if _MSC_VER >= 1600
#include <array>
#include <unordered_map>
#else
#include <tr1/array>
#include <tr1/unordered_map>
#endif

#include <cassert>

#include <vector>

namespace geom {

//! Class for a mesh that stores vertices and indexed simplices.
/*!
  \param N is the space dimension.
  \param M is the simplex dimension  By default it is N.
  \param T is the number type.  By default it is double.

  Note that the indices for indexed simplices follow the C convention of
  starting at 0.

  The free functions that operate on this class are grouped into the
  following categories:
  - \ref iss_accessors
  - \ref iss_build
  - \ref iss_equality
  - \ref iss_file_io
  - \ref iss_optimize
  - \ref iss_quality
  - \ref iss_set
  - \ref iss_tile
  - \ref iss_transfer
  - \ref iss_transform
*/
template < std::size_t _N, std::size_t _M = _N, typename T = double >
class IndSimpSet {
   //
   // Enumerations.
   //

public:

   //! The space dimension.
   static const std::size_t N = _N;
   //! The simplex dimension.
   static const std::size_t M = _M;

   //
   // Public types.
   //

public:

   //! The number type.
   typedef T Number;
   //! A vertex.
   typedef std::tr1::array<T, N> Vertex;
   //! A simplex of vertices.
   typedef std::tr1::array < Vertex, M + 1 > Simplex;
   //! The face of a simplex of vertices.
   typedef std::tr1::array<Vertex, M> SimplexFace;

   //! An indexed simplex.  (A simplex of indices.)
   typedef std::tr1::array < std::size_t, M + 1 > IndexedSimplex;
   //! The face of an indexed simplex.
   typedef std::tr1::array<std::size_t, M> IndexedSimplexFace;

   //! The vertex container.
   typedef std::vector<Vertex> VertexContainer;
   //! A vertex const iterator.
   typedef typename VertexContainer::const_iterator VertexConstIterator;
   //! A vertex iterator.
   typedef typename VertexContainer::iterator VertexIterator;

   //! The indexed simplex container.
   typedef std::vector<IndexedSimplex> IndexedSimplexContainer;
   //! An indexed simplex const iterator.
   typedef typename IndexedSimplexContainer::const_iterator
   IndexedSimplexConstIterator;
   //! An indexed simplex iterator.
   typedef typename IndexedSimplexContainer::iterator
   IndexedSimplexIterator;

   //! A simplex const iterator.
   typedef SimplexIterator<IndSimpSet> SimplexConstIterator;

   //! The size type.
   typedef typename VertexContainer::size_type SizeType;

private:

   //
   // Data
   //

   //! The vertices.
   VertexContainer _vertices;

   //! An indexed simplex is determined by the indices of M+1 vertices.
   IndexedSimplexContainer _indexedSimplices;

   //--------------------------------------------------------------------------
   /*! \name Constructors etc.
     Suppose that we are dealing with a tetrahedron mesh in 3-D.  Below
     we instantiate a mesh that allocates its own memory for the vertices
     and indexed simplices.
     \code
     geom::IndSimpSet<3,3> mesh;
     \endcode
     We can construct the mesh from vertices and indexed simplices stored in
     vectors:
     \code
     typedef geom::IndSimpSet<3,3> ISS;
     typedef typename ISS:Vertex Vertex;
     typedef typename ISS:IndexedSimplex IndexedSimplex;
     std::vector<Vertex> vertices(numberOfVertices);
     std::vector<IndexedSimplex> indexedSimplices(numberOfSimplices);
     ...
     geom::IndSimpSet<3,3> mesh(vertices, indexedSimplices);
     \endcode
     or use C arrays:
     \code
     double* vertices = new double[3 * numberOfVertices]
     std::size_t* indexedSimplices = new std::size_t[4 * numberOfSimplices];
     ...
     geom::IndSimpSet<3,3> mesh(numberOfVertices, vertices, numberOfSimplices, simplices);
     \endcode
   */
   //! @{

public:

   //! Default constructor.  Empty simplex set.
   IndSimpSet() :
      _vertices(),
      _indexedSimplices() {}

   //! Construct from arrays of vertices and indexed simplices.
   /*!
     \param vertices is the array of vertices.
     \param indexedSimplices is the array of indexed simplices.
   */
   IndSimpSet(const std::vector<Vertex>& vertices,
              const std::vector<IndexedSimplex>& indexedSimplices) :
      _vertices(vertices),
      _indexedSimplices(indexedSimplices) {}

   //! Build from arrays of vertices and indexed simplices.
   /*!
     Performs same actions as the constructor.

     \param vertices is the array of vertices.
     \param indexedSimplices is the array of indexed simplices.
   */
   void
   build(const std::vector<Vertex>& vertices,
         const std::vector<IndexedSimplex>& indexedSimplices) {
      _vertices = vertices;
      _indexedSimplices = indexedSimplices;
      updateTopology();
   }

   //! Construct from arrays of vertices, vertex identifiers and indexed simplices.
   /*!
     \param vertices is the array of vertices.
     \param vertexIdentifiers is the array of vertex identifiers.
     \param simplices is the array of indexed simplices, which are
     tuples of the vertex identifiers.
   */
   IndSimpSet(const std::vector<Vertex>& vertices,
              const std::vector<std::size_t>& vertexIdentifiers,
              const std::vector<IndexedSimplex>& simplices) :
      _vertices(vertices),
      // Start with simplices of vertex identifiers.
      _indexedSimplices(simplices) {
      convertFromIdentifiersToIndices(vertexIdentifiers);
   }

   //! Build from arrays of vertices, vertex identifiers and indexed simplices.
   /*!
     \param vertices is the array of vertices.
     \param vertexIdentifiers is the array of vertex identifiers.
     \param simplices is the array of indexed simplices, which are
     tuples of the vertex identifiers.
   */
   void
   build(const std::vector<Vertex>& vertices,
         const std::vector<std::size_t>& vertexIdentifiers,
         const std::vector<IndexedSimplex>& simplices) {
      // Start with simplices of vertex identifiers.
      _vertices = vertices;
      _indexedSimplices = simplices;
      convertFromIdentifiersToIndices(vertexIdentifiers);
      updateTopology();
   }

   //! Construct from pointers to the vertices and indexed simplices.
   /*!
     \param numVertices is the number of vertices.
     \param vertices points to the data for the vertices.
     \param numSimplices is the number of simplices.
     \param indexedSimplices points to the data for the indexed simplices.
   */
   IndSimpSet(const SizeType numVertices,
              const Vertex* vertices,
              const SizeType numSimplices,
              const IndexedSimplex* indexedSimplices) :
      _vertices(vertices, vertices + numVertices),
      _indexedSimplices(indexedSimplices, indexedSimplices + numSimplices) {}

   //! Build from pointers to the vertices and indexed simplices.
   /*!
     Performs same actions as the constructor.

     \param numVertices is the number of vertices.
     \param vertices points to the data for the vertices.
     \param numSimplices is the number of simplices.
     \param indexedSimplices points to the data for the indexed simplices.
   */
   void
   build(const SizeType numVertices,
         const Vertex* vertices,
         const SizeType numSimplices,
         const IndexedSimplex* indexedSimplices) {
      std::vector<Vertex> v(vertices, vertices + numVertices);
      _vertices.swap(v);
      std::vector<IndexedSimplex>
      is(indexedSimplices, indexedSimplices + numSimplices);
      _indexedSimplices.swap(is);
      updateTopology();
   }

   //! Construct from pointers to the vertices and indexed simplices.
   /*!
     \param numVertices is the number of vertices.
     \param coordinates points to the data for the vertex coordinates.
     \param numSimplices is the number of simplices.
     \param indices points to the data for the simplex indices.
   */
   template<typename _Index>
   IndSimpSet(const SizeType numVertices, const Number* coordinates,
              const SizeType numSimplices, const _Index* indices) :
      _vertices(numVertices),
      _indexedSimplices(numSimplices) {
      for (SizeType i = 0; i != numVertices; ++i) {
         for (SizeType n = 0; n != N; ++n) {
            _vertices[i][n] = *coordinates++;
         }
      }
      for (SizeType i = 0; i != numSimplices; ++i) {
         for (SizeType m = 0; m != M + 1; ++m) {
            _indexedSimplices[i][m] = *indices++;
         }
      }
   }

   //! Build from pointers to the vertices and indexed simplices.
   /*!
     Performs same actions as the constructor.

     \param numVertices is the number of vertices.
     \param coordinates points to the data for the vertex coordinates.
     \param numSimplices is the number of simplices.
     \param indices points to the data for the simplex indices.
   */
   template<typename _Index>
   void
   build(const SizeType numVertices, const Number* coordinates,
         const SizeType numSimplices, const _Index* indices) {
      _vertices.resize(numVertices);
      _indexedSimplices.resize(numSimplices);
      for (SizeType i = 0; i != numVertices; ++i) {
         for (SizeType n = 0; n != N; ++n) {
            _vertices[i][n] = *coordinates++;
         }
      }
      for (SizeType i = 0; i != numSimplices; ++i) {
         for (SizeType m = 0; m != M + 1; ++m) {
            _indexedSimplices[i][m] = *indices++;
         }
      }
      updateTopology();
   }

   //! Construct from the number of vertices and simplices.
   /*!
     The vertices and indexed simplices are left uninitialized.

     \param numVertices is the number of vertices.
     \param numSimplices is the number of simplices.
   */
   IndSimpSet(const SizeType numVertices, const SizeType numSimplices) :
      _vertices(numVertices),
      _indexedSimplices(numSimplices) {}

   //! Build from the number of vertices and simplices.
   /*!
     Performs same actions as the constructor.

     \param numVertices is the number of vertices.
     \param numSimplices is the number of simplices.
   */
   void
   build(const SizeType numVertices, const SizeType numSimplices) {
      _vertices.resize(numVertices);
      _indexedSimplices.resize(numSimplices);
   }

   //! Swap data with another mesh.
   void
   swap(IndSimpSet& x) {
      _vertices.swap(x._vertices);
      _indexedSimplices.swap(x._indexedSimplices);
   }

   //! Copy constructor.
   IndSimpSet(const IndSimpSet& other) :
      _vertices(other._vertices),
      _indexedSimplices(other._indexedSimplices) {}

   //! Assignment operator.
   IndSimpSet&
   operator=(const IndSimpSet& other);

   //! Destructor.  Deletes memory only if it was allocated internally.
   virtual
   ~IndSimpSet() {}

private:

   // Convert identifier simplices to index simplices.
   void
   convertFromIdentifiersToIndices
   (const std::vector<std::size_t>& vertexIdentifiers);

   //! @}
   //--------------------------------------------------------------------------
   //! \name Vertex Accessors
   //! @{

public:

   //! Return the dimension of the space.
   std::size_t
   getSpaceDimension() const {
      return N;
   }

   //! Return the number of vertices.
   SizeType
   getVerticesSize() const {
      return _vertices.size();
   }

   //! Return true if there are no vertices.
   SizeType
   areVerticesEmpty() const {
      return getVerticesSize() == 0;
   }

   //! Return a const iterator to the beginning of the vertices.
   VertexConstIterator
   getVerticesBeginning() const {
      return _vertices.begin();
   }

   //! Return a const iterator to the end of the vertices.
   VertexConstIterator
   getVerticesEnd() const {
      return _vertices.end();
   }

   //! Return a const reference to the n_th vertex.
   const Vertex&
   getVertex(const std::size_t n) const {
      return _vertices[n];
   }

   //! Return a const reference to the vertex container.
   const VertexContainer&
   getVertices() const {
      return _vertices;
   }

   //! @}
   //--------------------------------------------------------------------------
   //! \name Simplex Accessors
   //! @{

public:

   //! Return the dimension of the simplices.
   std::size_t
   getSimplexDimension() const {
      return M;
   }

   //! Return the number of simplices.
   SizeType
   getSimplicesSize() const {
      return _indexedSimplices.size();
   }

   //! Return true if there are no simplices.
   SizeType
   areSimplicesEmpty() const {
      return getSimplicesSize() == 0;
   }

   //! Return true if there are no vertices or simplices.
   SizeType
   isEmpty() const {
      return areVerticesEmpty() && areSimplicesEmpty();
   }

   //! Return a const iterator to the beginning of the indexed simplices.
   IndexedSimplexConstIterator
   getIndexedSimplicesBeginning() const {
      return _indexedSimplices.begin();
   }

   //! Return a const iterator to the end of the indexed simplices.
   IndexedSimplexConstIterator
   getIndexedSimplicesEnd() const {
      return _indexedSimplices.end();
   }

   //! Return a const reference to the n_th indexed simplex.
   const IndexedSimplex&
   getIndexedSimplex(const std::size_t n) const {
      return _indexedSimplices[n];
   }

   //! Return a const reference to the indexed simplex container.
   const IndexedSimplexContainer&
   getIndexedSimplices() const {
      return _indexedSimplices;
   }

   //! Return a const iterator to the beginning of the simplices.
   SimplexConstIterator
   getSimplicesBeginning() const {
      return SimplexConstIterator(*this);
   }

   //! Return a const iterator to the end of the simplices.
   SimplexConstIterator
   getSimplicesEnd() const {
      SimplexConstIterator i(*this);
      i += getSimplicesSize();
      return i;
   }

   //! Return a const reference to the m_th vertex of the n_th simplex.
   const Vertex&
   getSimplexVertex(const std::size_t n, const std::size_t m) const {
      return _vertices[_indexedSimplices[n][m]];
   }

   //! Get the n_th simplex.
   void
   getSimplex(const std::size_t n, Simplex* s) const {
      getSimplex(_indexedSimplices.begin() + n, s);
   }

   //! Get the simplex given an iterator to the indexed simplex.
   void
   getSimplex(IndexedSimplexConstIterator i, Simplex* s) const {
      for (std::size_t m = 0; m != M + 1; ++m) {
         (*s)[m] = _vertices[(*i)[m]];
      }
   }

   //! @}
   //--------------------------------------------------------------------------
   //! \name Vertex Manipulators
   //! @{

public:

   //! Return an iterator to the beginning of the vertices.
   VertexIterator
   getVerticesBeginning() {
      return _vertices.begin();
   }

   //! Return an iterator to the end of the vertices.
   VertexIterator
   getVerticesEnd() {
      return _vertices.end();
   }

   //! Set the specified vertex.
   void
   setVertex(const std::size_t n, const Vertex& vertex) {
      _vertices[n] = vertex;
   }

   //! Return a reference to the vertex container.
   VertexContainer&
   getVertices() {
      return _vertices;
   }

   //! @}
   //--------------------------------------------------------------------------
   //! \name Simplex Manipulators
   //! @{

public:

   //! Return an iterator to the beginning of the indexed simplices.
   IndexedSimplexIterator
   getIndexedSimplicesBeginning() {
      return _indexedSimplices.begin();
   }

   //! Return an iterator to the end of the indexed simplices.
   IndexedSimplexIterator
   getIndexedSimplicesEnd() {
      return _indexedSimplices.end();
   }

   //! Return a reference to the indexed simplex container.
   IndexedSimplexContainer&
   getIndexedSimplices() {
      return _indexedSimplices;
   }

   //! @}
   //--------------------------------------------------------------------------
   //! \name Update the topology.
   //! @{

public:

   //! Update the data structure following a change in the topology.
   /*!
     For this class, this function does nothing.  For derived classes,
     it updates data structures that hold auxillary topological information.
   */
   virtual
   void
   updateTopology() {}

   //! @}
};

} // namespace geom

#define __geom_IndSimpSet_ipp__
#include "IndSimpSet.ipp"
#undef __geom_IndSimpSet_ipp__

#endif
