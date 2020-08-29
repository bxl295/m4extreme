// -*- C++ -*-

#if !defined(__geom_mesh_iss_fit_ipp__)
#error This file is an implementation detail.
#endif

namespace geom {


template<typename T, class ISS>
inline
void
fit(IndSimpSetIncAdj<2, 1, T>* mesh,
    const ISS_SignedDistance<ISS, 2>& signedDistance,
    const T deviationTangent) {
   typedef IndSimpSetIncAdj<2, 1, T> Mesh;
   typedef typename Mesh::Vertex Vertex;

   std::size_t i, j;
   std::size_t previousEdgeIndex, previousVertexIndex, nextEdgeIndex, nextVertexIndex;
   Vertex point;
   Vertex tangent;
   T previousDistance, nextDistance, previousDeviation, nextDeviation;

   // Loop over the vertices.
   const std::size_t size = mesh->getVerticesSize();
   for (std::size_t n = 0; n != size; ++n) {
      // Don't alter boundary vertices.
      if (mesh->isVertexOnBoundary(n)) {
         continue;
      }

      //
      // Determine the deviation tangents.
      //

      // The vertex should have two incident faces.
      assert(mesh->getIncidentSize(n) == 2);
      // The two incident edges.
      i = mesh->getIncident(n, 0);
      j = mesh->getIncident(n, 1);
      if (mesh->getIndexedSimplex(i)[0] == n) {
         std::swap(i, j);
      }
      assert(mesh->getIndexedSimplex(i)[1] == n &&
             mesh->getIndexedSimplex(j)[0] == n);
      // Now i is the previous edge and j is the next edge.
      previousEdgeIndex = i;
      previousVertexIndex = mesh->getIndexedSimplex(i)[0];
      nextEdgeIndex = j;
      nextVertexIndex = mesh->getIndexedSimplex(j)[1];

      // The distance at the mid-points of the edges.
      point = mesh->getVertex(n);
      point += mesh->getVertex(previousVertexIndex);
      point *= 0.5;
      previousDistance = signedDistance(point);

      point = mesh->getVertex(n);
      point += mesh->getVertex(nextVertexIndex);
      point *= 0.5;
      nextDistance = signedDistance(point);

      /* REMOVE
      std::cout << "previousDistance = " << previousDistance
            << ", nextDistance = " << nextDistance << "\n";
      */

      // If the deviation is small, don't do anything with this vertex.
      previousDeviation = 2.0 * std::abs(previousDistance) /
                          euclideanDistance(mesh->getVertex(previousVertexIndex),
                                            mesh->getVertex(n));
      nextDeviation = 2.0 * std::abs(nextDistance) /
                      euclideanDistance(mesh->getVertex(nextVertexIndex), mesh->getVertex(n));

      if (previousDeviation < deviationTangent &&
            nextDeviation < deviationTangent) {
         continue;
      }

      // The offset distance.
#if 0
      if (previousDistance > 0 && nextDistance > 0) {
         //offset = - 2.0 * std::abs(previousDistance - nextDistance);
         offset = - 2.0 * std::max(std::abs(previousDistance),
                                   std::abs(nextDistance));
      }
      else if (previousDistance < 0 && nextDistance < 0) {
         //offset = 2.0 * std::abs(previousDistance - nextDistance);
         offset = 2.0 * std::max(std::abs(previousDistance),
                                 std::abs(nextDistance));
      }
      else {
         offset = - 2.0 * (previousDistance + nextDistance);
      }
      offset = - 2.0 * (previousDistance + nextDistance);
#endif

#if 0
      //
      // Compute the outward normal direction.
      //
      normal = 0.0;

      point = mesh->getVertex(previousVertexIndex);
      point -= mesh->getVertex(n);
      rotatePiOver2(&point);
      normalize(&point);
      normal += point;

      point = mesh->getVertex(nextVertexIndex);
      point -= mesh->getVertex(n);
      rotateMinusPiOver2(&point);
      normalize(&point);
      normal += point;

      normalize(&normal);
#endif

      //
      // Find the tangent to the surface.
      //

      tangent = signedDistance.computeNormal(mesh->getVertex(n));
      rotatePiOver2(&tangent);
      // The offset.
      tangent *= 2.0 * (std::abs(nextDistance) - std::abs(previousDistance));

      //
      // Offset the vertex and find the closest point.
      //
      point = mesh->getVertex(n);
      point += tangent;
      mesh->setVertex(n, signedDistance.computeClosestPoint(point));
      /* REMOVE
      std::cout << tangent << " "
            << mesh->vertices()[n] << "\n";
      */
   }
}



template<typename T, class ISS>
inline
void
fit(IndSimpSetIncAdj<2, 1, T>* mesh,
    const ISS_SignedDistance<ISS, 2>& signedDistance,
    const T deviationTangent, std::size_t numSweeps) {
   while (numSweeps-- > 0) {
      fit(mesh, signedDistance, deviationTangent);
   }
}



// Fit the boundary of a mesh to a level-set description.
template<typename T, class ISS>
inline
void
fit(IndSimpSetIncAdj<2, 2, T>* mesh,
    const ISS_SignedDistance<ISS, 2>& signedDistance,
    const T deviationTangent,
    std::size_t numSweeps) {
   typedef IndSimpSetIncAdj<2, 1, T> Boundary;

   Boundary boundary;
   std::vector<std::size_t> boundaryVertexIndices;

   // Get the boundary.
   buildBoundary(*mesh, &boundary,
                 std::back_inserter(boundaryVertexIndices));

   // Fit the boundary to the level-set.
   fit(&boundary, signedDistance, deviationTangent, numSweeps);

   // Update the vertices in the mesh->
   const std::size_t size = boundaryVertexIndices.size();
   assert(size == boundary.getVerticesSize());
   for (std::size_t n = 0; n != size; ++n) {
      mesh->getVertices()[boundaryVertexIndices[n]] = boundary.getVertex(n);
   }
}

} // namespace geom
