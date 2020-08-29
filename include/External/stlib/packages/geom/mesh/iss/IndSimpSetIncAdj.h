// -*- C++ -*-

/*!
  \file IndSimpSetIncAdj.h
  \brief An indexed simplex set in N-D that optimizes simplex quality.
*/

#if !defined(__geom_IndSimpSetIncAdj_h__)
#define __geom_IndSimpSetIncAdj_h__

#include "IndSimpSet.h"
#include "SimplexAdj.h"
#include "VertexSimplexInc.h"
#include "IssiaFaceIterator.h"

namespace geom {

//! An indexed simplex set that stores vertex-simplex incidences and simplex adjacencies.
/*!
  \param N is the space dimension.
  \param M is the simplex dimension  By default it is N.
  \param T is the number type.  By default it is double.
  \param V is the vertex type, an N-tuple of the number type.  It must be
  subscriptable.  By default it is std::tr1::array<T,N>.
  \param IS is the Indexed Simplex type, a tuple of M+1 integers.
  It must be subscriptable.  By default it is Simplex<M,std::size_t>.

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

  This class derives from geom::IndSimpSet.  Any function that takes a
  geom::IndSimpSet as an argument may also take this class as an argument.
  This includes function that build the mesh or modify the topology.
  This functionality is made possible with the update_topology() virtual
  function.  Any free function that modifies the topology of the mesh
  calls update_topology().  In the base class, the function has no effect,
  but in this class, it builds/rebuilds the vertex-simplex incidences
  and the simplex adjacencies.
*/
template < std::size_t _N,
         std::size_t _M = _N,
         typename T = double >
class IndSimpSetIncAdj :
   public IndSimpSet<_N, _M, T> {
   //
   // The base type.
   //

public:

   //! The base type.
   typedef IndSimpSet<_N, _M, T> Base;

   //
   // Enumerations.
   //

public:

   //! The space dimension.
   static const std::size_t N = Base::N;
   //! The simplex dimension.
   static const std::size_t M = Base::M;

   //
   // Public types.
   //

public:

   //
   // Inherited from IndSimpSet.
   //

   //! The number type.
   typedef typename Base::Number Number;
   //! A vertex.
   typedef typename Base::Vertex Vertex;
   //! A simplex of vertices.
   typedef typename Base::Simplex Simplex;
   //! The face of a simplex of vertices.
   typedef typename Base::SimplexFace SimplexFace;

   //! An indexed simplex.  (A simplex of indices.)
   typedef typename Base::IndexedSimplex IndexedSimplex;
   //! The face of an indexed simplex.
   typedef typename Base::IndexedSimplexFace IndexedSimplexFace;

   //! The vertex container.
   typedef typename Base::VertexContainer VertexContainer;
   //! A vertex const iterator.
   typedef typename Base::VertexConstIterator VertexConstIterator;
   //! A vertex iterator.
   typedef typename Base::VertexIterator VertexIterator;

   //! The indexed simplex container.
   typedef typename Base::IndexedSimplexContainer IndexedSimplexContainer;
   //! An indexed simplex const iterator.
   typedef typename Base::IndexedSimplexConstIterator
   IndexedSimplexConstIterator;
   //! An indexed simplex iterator.
   typedef typename Base::IndexedSimplexIterator IndexedSimplexIterator;

   //! A simplex const iterator.
   typedef typename Base::SimplexConstIterator SimplexConstIterator;

   //! The size type.
   typedef typename Base::SizeType SizeType;

   //
   // New types.
   //

   //! The vertex-simplex incidences.
   typedef VertexSimplexInc<M> IncidenceContainer;
   //! Iterator over the vertex-simplex incidences.
   typedef typename IncidenceContainer::ConstIterator IncidenceConstIterator;

   //! The simplex adjacencies.
   typedef SimplexAdj<M> AdjacencyContainer;

   // Faces.

   //! A face is determined by a simplex index and an integer in [0..M].
   typedef std::pair<std::size_t, std::size_t> Face;
   //! A bidirectional, iterator on the faces.
   typedef IssiaFaceIterator<IndSimpSetIncAdj> FaceIterator;

   //
   // Member data.
   //

private:

   //! The vertex-simplex incidences.
   IncidenceContainer _vertexSimplexIncidence;
   //! The simplex adjacencies.
   AdjacencyContainer _simplexAdjacencies;

public:

   //--------------------------------------------------------------------------
   //! \name Constructors etc.
   //! @{

   //! Default constructor.
   IndSimpSetIncAdj() :
      Base(),
      _vertexSimplexIncidence(),
      _simplexAdjacencies() {}

   //! Copy constructor.
   IndSimpSetIncAdj(const IndSimpSetIncAdj& other) :
      Base(other),
      _vertexSimplexIncidence(other._vertexSimplexIncidence),
      _simplexAdjacencies(other._simplexAdjacencies) {}

   //! Assignment operator.
   IndSimpSetIncAdj&
   operator=(const IndSimpSetIncAdj& other) {
      if (&other != this) {
         Base::operator=(other);
         _vertexSimplexIncidence = other._vertexSimplexIncidence;
         _simplexAdjacencies = other._simplexAdjacencies;
      }
      return *this;
   }

   //! Construct from arrays of vertices and indexed simplices.
   /*!
     \param vertices is the array of vertices.
     \param indexedSimplices is the array of indexed simplices.
   */
   IndSimpSetIncAdj(const std::vector<Vertex>& vertices,
                    const std::vector<IndexedSimplex>& indexedSimplices) :
      Base(vertices, indexedSimplices),
      _vertexSimplexIncidence(getVerticesSize(), getIndexedSimplices()),
      _simplexAdjacencies(getIndexedSimplices(), _vertexSimplexIncidence) {}

   //! Build from arrays of vertices and indexed simplices.
   /*!
     Performs same actions as the constructor.

     \param vertices is the array of vertices.
     \param indexedSimplices is the array of indexed simplices.
   */
   void
   build(const std::vector<Vertex>& vertices,
         const std::vector<IndexedSimplex>& indexedSimplices) {
      Base::build(vertices, indexedSimplices);
      updateTopology();
   }

   //! Construct from arrays of vertices, vertex identifiers and indexed simplices.
   /*!
     \param vertices is the array of vertices.
     \param vertexIdentifiers is the array of vertex identifiers.
     \param simplices is the array of indexed simplices, which are
     tuples of the vertex identifiers.
   */
   IndSimpSetIncAdj(const std::vector<Vertex>& vertices,
                    const std::vector<std::size_t>& vertexIdentifiers,
                    const std::vector<IndexedSimplex>& simplices) :
      Base(vertices, vertexIdentifiers, simplices),
      _vertexSimplexIncidence(getVerticesSize(), getIndexedSimplices()),
      _simplexAdjacencies(getIndexedSimplices(), _vertexSimplexIncidence) {}

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
      Base::build(vertices, vertexIdentifiers, simplices);
      updateTopology();
   }

   //! Construct from pointers to the vertices and indexed simplices.
   /*!
     \param numVertices is the number of vertices.
     \param vertices points to the data for the vertices.
     \param numSimplices is the number of simplices.
     \param indexedSimplices points to the data for the indexed simplices.
    */
   IndSimpSetIncAdj(const SizeType numVertices,
                    const Vertex* vertices,
                    const SizeType numSimplices,
                    const IndexedSimplex* indexedSimplices) :
      Base(numVertices, vertices, numSimplices, indexedSimplices),
      _vertexSimplexIncidence(getVerticesSize(), getIndexedSimplices()),
      _simplexAdjacencies(getIndexedSimplices(), _vertexSimplexIncidence) {}

   //! Build from pointers to the vertices and indexed simplices.
   /*!
     Performs the same action as the constructor.

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
      Base::build(numVertices, vertices, numSimplices, indexedSimplices);
      updateTopology();
   }

   //! Construct from pointers to the vertices and indexed simplices.
   /*!
     \param numVertices is the number of vertices.
     \param coordinates points to the data for the vertices.
     \param numSimplices is the number of simplices.
     \param indices points to the data for the indexed simplices.
    */
   template<typename _Index>
   IndSimpSetIncAdj(const SizeType numVertices, const Number* coordinates,
                    const SizeType numSimplices, const _Index* indices) :
      Base(numVertices, coordinates, numSimplices, indices),
      _vertexSimplexIncidence(getVerticesSize(), getIndexedSimplices()),
      _simplexAdjacencies(getIndexedSimplices(), _vertexSimplexIncidence) {}

   //! Build from pointers to the vertices and indexed simplices.
   /*!
     Performs the same action as the constructor.

     \param numVertices is the number of vertices.
     \param coordinates points to the data for the vertex coordinates.
     \param numSimplices is the number of simplices.
     \param indices points to the data for the simplex indices.
   */
   template<typename _Index>
   void
   build(const SizeType numVertices, const Number* coordinates,
         const SizeType numSimplices, const _Index* indices) {
      Base::build(numVertices, coordinates, numSimplices, indices);
      updateTopology();
   }

   //! Construct from the number of vertices and simplices.
   /*!
     The vertices and indexed simplices are left uninitialized.  The incidence
     and adjacency relations are not built.

     \param numVertices is the number of vertices.
     \param numSimplices is the number of simplices.
   */
   IndSimpSetIncAdj(const SizeType numVertices, const SizeType numSimplices) :
      Base(numVertices, numSimplices),
      _vertexSimplexIncidence(),
      _simplexAdjacencies() {}

   //! Build from the number of vertices and simplices.
   /*!
     Performs the same action as the constructor.

     \param numVertices is the number of vertices.
     \param numSimplices is the number of simplices.
   */
   void
   build(const SizeType numVertices, const SizeType numSimplices) {
      Base::build(numVertices, numSimplices);
      {
         // Clear the vertex-simplex incidence relations.
         IncidenceContainer empty;
         _vertexSimplexIncidence.swap(empty);
      }
      {
         // Clear the adjacency relations.
         AdjacencyContainer empty;
         _simplexAdjacencies.swap(empty);
      }
   }

   //! Swap data with another mesh.
   void
   swap(IndSimpSetIncAdj& x) {
      Base::swap(x);
      _vertexSimplexIncidence.swap(x._vertexSimplexIncidence);
      _simplexAdjacencies.swap(x._simplexAdjacencies);
   }

   //! Construct from an indexed simplex set.
   /*!
     \param iss is the indexed simplex set.
   */
   IndSimpSetIncAdj(const Base& iss) :
      Base(iss),
      _vertexSimplexIncidence(getVerticesSize(), getIndexedSimplices()),
      _simplexAdjacencies(getIndexedSimplices(), _vertexSimplexIncidence) {}

   //! Destructor.  Free any memory that was allocated.
   virtual
   ~IndSimpSetIncAdj() {}

   //! @}
   //--------------------------------------------------------------------------
   /*! \name Vertex Accessors
     Inherited from IndSimpSet.
   */
   //! @{

   //! Return the dimension of the space.
   using Base::getSpaceDimension;

   //! Return the number of vertices.
   using Base::getVerticesSize;

   //! Return true if there are no vertices.
   using Base::areVerticesEmpty;

   //! Return a const iterator to the beginning of the vertices.
   using Base::getVerticesBeginning;

   //! Return a const iterator to the end of the vertices.
   using Base::getVerticesEnd;

   //! Return a const reference to the n_th vertex.
   using Base::getVertex;

   //! Return a const reference to the vertex container.
   using Base::getVertices;

   //! @}
   //--------------------------------------------------------------------------
   /*! \name Simplex Accessors
     Inherited from IndSimpSet.
    */
   //! @{

   //! Return the dimension of the simplices.
   using Base::getSimplexDimension;

   //! Return the number of simplices.
   using Base::getSimplicesSize;

   //! Return true if there are no simplices.
   using Base::areSimplicesEmpty;

   //! Return true if there are no vertices or simplices.
   using Base::isEmpty;

   //! Return a const iterator to the beginning of the indexed simplices.
   using Base::getIndexedSimplicesBeginning;

   //! Return a const iterator to the end of the indexed simplices.
   using Base::getIndexedSimplicesEnd;

   //! Return a const reference to the n_th indexed simplex.
   using Base::getIndexedSimplex;

   //! Return a const reference to the indexed simplex container.
   using Base::getIndexedSimplices;

   //! Return a const iterator to the beginning of the simplices.
   using Base::getSimplicesBeginning;

   //! Return a const iterator to the end of the simplices.
   using Base::getSimplicesEnd;

   //! Return a const reference to the m_th vertex of the n_th simplex.
   using Base::getSimplexVertex;

   //! Get the n_th simplex.
   using Base::getSimplex;

   //! Get the simplex given an iterator to the indexed simplex.
   //using Base::getSimplex;

   //! @}
   //--------------------------------------------------------------------------
   //! \name Vertex-Simplex Incidence Accessors
   //! @{

   //! A const reference to the vertex-simplex incidence.
   const IncidenceContainer&
   getVertexSimplexIncidence() const {
      return _vertexSimplexIncidence;
   }

   //! Return the number of incident simplices to the n_th vertex.
   SizeType
   getIncidentSize(const std::size_t n) const {
      return _vertexSimplexIncidence.getSize(n);
   }

   //! Return true if the n_th vertex has no incident simplices.
   bool
   isIncidentEmpty(const std::size_t n) const {
      return _vertexSimplexIncidence.isEmpty(n);
   }

   //! Return a const iterator on simplex indices to the first incident simplex index of the n_th vertex.
   IncidenceConstIterator
   getIncidentBeginning(const std::size_t n) const {
      return _vertexSimplexIncidence.getBeginning(n);
   }

   //! Return a const iterator on simplex indices to one past the last incident simplex index of the n_th vertex.
   IncidenceConstIterator
   getIncidentEnd(const std::size_t n) const {
      return _vertexSimplexIncidence.getEnd(n);
   }

   //! Return the m_th incident simplex index of the n_th vertex.
   std::size_t
   getIncident(const std::size_t n, const std::size_t m) const {
      return _vertexSimplexIncidence(n, m);
   }

   //! @}
   //--------------------------------------------------------------------------
   //! \name Simplex Adjacency Accessors
   //! @{

   //! Return a const reference to the simplex adjacencies.
   const AdjacencyContainer&
   getSimplexAdjacencies() const {
      return _simplexAdjacencies;
   }

   //! Return number of simplices adjacent to the n_th simplex.
   std::size_t
   getAdjacentSize(const std::size_t n) const {
      return _simplexAdjacencies.getSize(n);
   }

   //! Return the index of the m_th adjacent simplex to the n_th simplex.
   /*!
     An index of -1 indicates that there is no adjacent simplex.
   */
   std::size_t
   getAdjacent(const std::size_t n, const std::size_t m) const {
      return _simplexAdjacencies(n, m);
   }

   //! Return the index of the n_th simplex in its m_th adjacent neighbor.
   std::size_t
   getMirrorIndex(const std::size_t n, const std::size_t m) const {
      const std::size_t a = getAdjacent(n, m);
      if (a == std::size_t(-1)) {
         return std::numeric_limits<std::size_t>::max();
      }
#ifdef DEBUG_stlib
      const std::size_t mi = std::find(_simplexAdjacencies(a).begin(),
                                       _simplexAdjacencies(a).end(), n) -
                             _simplexAdjacencies(a).begin();
      assert(mi < M + 1);
      return mi;
#else
      // Return the index of the value n.
      return std::find(_simplexAdjacencies(a).begin(),
                       _simplexAdjacencies(a).end(), n) -
             _simplexAdjacencies(a).begin();
#endif
   }

   //! @}
   //--------------------------------------------------------------------------
   //! \name Face accessors.
   //! @{

   //! Return the number of faces.
   /*!
     \note This is a slow function.  It counts the faces.
   */
   SizeType
   computeFacesSize() const {
      return std::distance(getFacesBeginning(), getFacesEnd());
   }

   //! Return the beginning of the faces.
   FaceIterator
   getFacesBeginning() const {
      FaceIterator x(this, 0);
      return x;
   }

   //! Return the end of the faces.
   FaceIterator
   getFacesEnd() const {
      return FaceIterator(this, getSimplicesSize());
   }

   //! Return true if the face is on the boundary.
   bool
   isOnBoundary(const Face& f) const {
      return getAdjacent(f.first, f.second) == std::size_t(-1);
   }

   //! Return true if the face is on the boundary.
   bool
   isOnBoundary(const FaceIterator& f) const {
      return isOnBoundary(*f);
   }

   //! @}
   //--------------------------------------------------------------------------
   //! \name Other Accessors
   //! @{

   //! Return true if the vertex is on the boundary of the mesh.
   /*!
     \param index is the index of a vertex.
   */
   bool
   isVertexOnBoundary(std::size_t index) const;

   //! @}
   //--------------------------------------------------------------------------
   // Vertex Manipulators
   // Inherited from IndSimpSet.
   // These are already listed with the vertex accessors.

   // Return an iterator to the beginning of the vertices.
   //using Base::getVerticesBeginning;

   // Return an iterator to the end of the vertices.
   //using Base::getVerticesEnd;

   // Return a reference to the vertex container.
   //using Base::getVertices;

   //--------------------------------------------------------------------------
   //! \name Simplex Manipulators
   //! @{

   // Return an iterator to the beginning of the indexed simplices.
   //using Base::getIndexedSimplicesBeginning;

   // Return an iterator to the end of the indexed simplices.
   //using Base::getIndexedSimplicesEnd;

   // Return a reference to the indexed simplex container.
   //using Base::getIndexedSimplices;

   //! Reverse the orientation of the n_th simplex.
   void
   reverseOrientation(const std::size_t n) {
      std::swap(getIndexedSimplices()[n][0], getIndexedSimplices()[n][1]);
      std::swap(_simplexAdjacencies(n)[0], _simplexAdjacencies(n)[1]);
   }

   //! @}
   //--------------------------------------------------------------------------
   //! \name Update the topology.
   //! @{

   //! Update the data structure following a change in the topology.
   /*!
     Update the vertex-simplex incidences and simplex adjacencies following
     a change in the topology.
   */
   virtual
   void
   updateTopology() {
      _vertexSimplexIncidence.build(getVerticesSize(), getIndexedSimplices());
      _simplexAdjacencies.build(getIndexedSimplices(), _vertexSimplexIncidence);
   }

   //! @}
};

} // namespace geom

#define __geom_IndSimpSetIncAdj_ipp__
#include "IndSimpSetIncAdj.ipp"
#undef __geom_IndSimpSetIncAdj_ipp__

#endif
