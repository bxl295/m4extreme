// -*- C++ -*-

#if !defined(__geom_IndSimpSetIncAdj_ipp__)
#error This file is an implementation detail of the class IndSimpSetIncAdj.
#endif

namespace geom {

template<std::size_t _N, std::size_t _M, typename _T>
inline
bool
IndSimpSetIncAdj<_N, _M, _T>::
isVertexOnBoundary(const std::size_t vertexIndex) const {
   std::size_t si, m, n;
   // For each simplex incident to this vertex.
   for (IncidenceConstIterator sii =
            _vertexSimplexIncidence.getBeginning(vertexIndex);
         sii != _vertexSimplexIncidence.getEnd(vertexIndex); ++sii) {
      // The simplex index.
      si = *sii;
      // The number of the vertex in the simplex.
      n = index(getIndexedSimplices()[si], vertexIndex);
      // For each face incident to the vertex.
      for (m = 0; m != M + 1; ++m) {
         if (m != n) {
            // If the face is on the boundary.
            if (_simplexAdjacencies(si, m) == std::size_t(-1)) {
               return true;
            }
         }
      }
   }
   return false;
}

} // namespace geom
