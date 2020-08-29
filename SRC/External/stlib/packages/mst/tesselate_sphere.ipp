// -*- C++ -*-

#if !defined(__tesselate_sphere_ipp__)
#error This file is an implementation detail of tesselate_sphere.
#endif

namespace mst {


template<typename T>
inline
void
makeOctahedron(geom::IndSimpSet<3, 2, T>* mesh) {
   typedef geom::IndSimpSet<3, 2, T> ISS;
   typedef typename ISS::Vertex Vertex;
   typedef typename ISS::IndexedSimplex IndexedSimplex;

   std::vector<Vertex> vertices(6);
   vertices[0] = ext::make_array<T>(1, 0, 0);
   vertices[1] = ext::make_array<T>(-1, 0, 0);
   vertices[2] = ext::make_array<T>(0, 1, 0);
   vertices[3] = ext::make_array<T>(0, -1, 0);
   vertices[4] = ext::make_array<T>(0, 0, 1);
   vertices[5] = ext::make_array<T>(0, 0, -1);

   std::vector<IndexedSimplex> indexedSimplices(8);
   indexedSimplices[0] = ext::make_array<std::size_t>(0, 2, 4);
   indexedSimplices[1] = ext::make_array<std::size_t>(2, 1, 4);
   indexedSimplices[2] = ext::make_array<std::size_t>(1, 3, 4);
   indexedSimplices[3] = ext::make_array<std::size_t>(3, 0, 4);
   indexedSimplices[4] = ext::make_array<std::size_t>(2, 0, 5);
   indexedSimplices[5] = ext::make_array<std::size_t>(1, 2, 5);
   indexedSimplices[6] = ext::make_array<std::size_t>(3, 1, 5);
   indexedSimplices[7] = ext::make_array<std::size_t>(0, 3, 5);

   mesh->build(vertices, indexedSimplices);
}



template<typename T>
inline
void
makeIcosahedron(geom::IndSimpSet<3, 2, T>* mesh) {
   typedef geom::IndSimpSet<3, 2, T> ISS;
   typedef typename ISS::Vertex Vertex;
   typedef typename ISS::IndexedSimplex IndexedSimplex;

   const T t = (1.0 + std::sqrt(5.0)) / 2.0;
   const T d = std::sqrt(1.0 + t * t);

   std::vector<Vertex> vertices(12);
   vertices[ 0] = ext::make_array<T>(t,  1,  0) / d;
   vertices[ 1] = ext::make_array<T>(-t,  1,  0) / d;
   vertices[ 2] = ext::make_array<T>(t, -1,  0) / d;
   vertices[ 3] = ext::make_array<T>(-t, -1,  0) / d;
   vertices[ 4] = ext::make_array<T>(1,  0,  t) / d;
   vertices[ 5] = ext::make_array<T>(1,  0, -t) / d;
   vertices[ 6] = ext::make_array<T>(-1,  0,  t) / d;
   vertices[ 7] = ext::make_array<T>(-1,  0, -t) / d;
   vertices[ 8] = ext::make_array<T>(0,  t,  1) / d;
   vertices[ 9] = ext::make_array<T>(0, -t,  1) / d;
   vertices[10] = ext::make_array<T>(0,  t, -1) / d;
   vertices[11] = ext::make_array<T>(0, -t, -1) / d;

   std::vector<IndexedSimplex> indexedSimplices(20);
   indexedSimplices[ 0] = ext::make_array<std::size_t>(0,  8,  4);
   indexedSimplices[ 1] = ext::make_array<std::size_t>(0,  5, 10);
   indexedSimplices[ 2] = ext::make_array<std::size_t>(2,  4,  9);
   indexedSimplices[ 3] = ext::make_array<std::size_t>(2, 11,  5);
   indexedSimplices[ 4] = ext::make_array<std::size_t>(1,  6,  8);
   indexedSimplices[ 5] = ext::make_array<std::size_t>(1, 10,  7);
   indexedSimplices[ 6] = ext::make_array<std::size_t>(3,  9,  6);
   indexedSimplices[ 7] = ext::make_array<std::size_t>(3,  7, 11);
   indexedSimplices[ 8] = ext::make_array<std::size_t>(0, 10,  8);
   indexedSimplices[ 9] = ext::make_array<std::size_t>(1,  8, 10);
   indexedSimplices[10] = ext::make_array<std::size_t>(2,  9, 11);
   indexedSimplices[11] = ext::make_array<std::size_t>(3, 11,  9);
   indexedSimplices[12] = ext::make_array<std::size_t>(4,  2,  0);
   indexedSimplices[13] = ext::make_array<std::size_t>(5,  0,  2);
   indexedSimplices[14] = ext::make_array<std::size_t>(6,  1,  3);
   indexedSimplices[15] = ext::make_array<std::size_t>(7,  3,  1);
   indexedSimplices[16] = ext::make_array<std::size_t>(8,  6,  4);
   indexedSimplices[17] = ext::make_array<std::size_t>(9,  4,  6);
   indexedSimplices[18] = ext::make_array<std::size_t>(10,  5,  7);
   indexedSimplices[19] = ext::make_array<std::size_t>(11,  7,  5);

   mesh->build(vertices, indexedSimplices);
}



template<typename T>
inline
void
tesselateUnitSphere(const T maxEdgeLength,
                    geom::IndSimpSetIncAdj<3, 2, T>* mesh) {
   typedef geom::IndSimpSetIncAdj<3, 2, T> Mesh;
   typedef typename Mesh::VertexIterator VertexIterator;
   typedef T Number;

   assert(maxEdgeLength > 0);

   // The vector of meshes that we have obtained through subdivision.
   // We make this static to avoid recomputing the meshes.
   static std::vector<Mesh> meshes;
   static std::vector<Number> maxEdgeLengths;

   // If we have not yet computed the coarsest mesh.
   if (meshes.size() == 0) {
      // Start with an octahedron and an icosahedron.
      Mesh m;
      makeOctahedron(&m);
      meshes.push_back(m);
      maxEdgeLengths.push_back(geom::computeMaximumEdgeLength(m));
      makeIcosahedron(&m);
      meshes.push_back(m);
      maxEdgeLengths.push_back(geom::computeMaximumEdgeLength(m));
   }

   // If the current meshes are too coarse, compute finer ones.
   while (maxEdgeLength < maxEdgeLengths.back()) {
      Mesh subdividedMesh;
      // Subdivide the second to last mesh.
      geom::subdivide(*(meshes.end() - 2), &subdividedMesh);
      // Normalize the vertices.
      for (VertexIterator i = subdividedMesh.getVerticesBeginning();
            i != subdividedMesh.getVerticesEnd(); ++i) {
         normalize(&*i);
      }
      // Copy the current finest mesh.
      meshes.push_back(subdividedMesh);
      // Compute the maximum edge length for the refined mesh.
      maxEdgeLengths.push_back(computeMaximumEdgeLength(subdividedMesh));
   }

   // Determine the appropriate mesh.
   std::size_t index = 0;
   while (maxEdgeLengths[index] > maxEdgeLength) {
      ++index;
   }
   *mesh = meshes[index];
}



template<typename T>
inline
void
tesselateUnitSphere(const int refinementLevel,
                    geom::IndSimpSetIncAdj<3, 2, T>* mesh) {
   typedef geom::IndSimpSetIncAdj<3, 2, T> Mesh;
   typedef typename Mesh::VertexIterator VertexIterator;
   typedef T Number;

   assert(refinementLevel >= 0);

   // The vector of meshes that we have obtained through subdivision.
   // We make this static to avoid recomputing the meshes.
   static std::vector<Mesh> meshes;

   // If we have not yet computed the coarsest mesh.
   if (meshes.size() == 0) {
      // Start with an octahedron and an icosahedron.
      Mesh m;
      makeOctahedron(&m);
      meshes.push_back(m);
      makeIcosahedron(&m);
      meshes.push_back(m);
   }

   // If the current meshes are too coarse, compute finer ones.
   while (int(meshes.size()) < refinementLevel + 1) {
      Mesh subdividedMesh;
      // Subdivide the second to last mesh.
      geom::subdivide(*(meshes.end() - 2), &subdividedMesh);
      // Normalize the vertices.
      for (VertexIterator i = subdividedMesh.getVerticesBeginning();
            i != subdividedMesh.getVerticesEnd(); ++i) {
         normalize(&*i);
      }
      // Copy the current finest mesh.
      meshes.push_back(subdividedMesh);
   }

   // Copy the appropriate mesh.
   *mesh = meshes[refinementLevel];
}



template<typename T, typename PointInIter, typename NumberInIter >
inline
void
tesselateAllSurfaces(PointInIter centersBegin, PointInIter centersEnd,
                     NumberInIter radiiBegin, NumberInIter radiiEnd,
                     const T maxEdgeLength,
                     geom::IndSimpSetIncAdj<3, 2, T>* mesh) {
   typedef geom::IndSimpSetIncAdj<3, 2, T> Mesh;
   typedef typename Mesh::VertexIterator VertexIterator;
   typedef T Number;
   typedef typename std::iterator_traits<PointInIter>::value_type Point;

   assert(maxEdgeLength > 0);

   // The list of sphere meshes.
   std::list<Mesh> meshes;

   // For each sphere.
   for (; centersBegin != centersEnd; ++centersBegin, ++radiiBegin) {
      // Check that the two ranges are the same size.
      assert(radiiBegin != radiiEnd);
      // The center.
      const Point center = *centersBegin;
      // The radius.
      const Number radius = *radiiBegin;
      // Make an empty mesh.
      Mesh m;
      // Compute the triangulation of the unit sphere.
      tesselateUnitSphere(maxEdgeLength / radius, &m);
      // Transform the unit sphere to match this sphere.
      for (VertexIterator i = m.getVerticesBeginning(); i != m.getVerticesEnd();
            ++i) {
         *i *= radius;
         *i += center;
      }
      // Add the new mesh to the list.
      meshes.push_back(m);
   }
   // Check that the two ranges are the same size.
   assert(radiiBegin == radiiEnd);

   // Merge the meshes into a single mesh.
   geom::merge(meshes.begin(), meshes.end(), mesh);
}


} // namespace mst
