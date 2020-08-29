// -*- C++ -*-

#if !defined(__geom_IndSimpSet_ipp__)
#error This file is an implementation detail of the class IndSimpSet.
#endif

namespace geom {

// Assignment operator.
template<std::size_t _N, std::size_t _M, typename _T>
inline
IndSimpSet<_N, _M, _T>&
IndSimpSet<_N, _M, _T>::
operator=(const IndSimpSet& other) {
   if (this != &other) {
      _vertices = other._vertices;
      _indexedSimplices = other._indexedSimplices;
   }
   return *this;
}

template<std::size_t _N, std::size_t _M, typename _T>
inline
void
IndSimpSet<_N, _M, _T>::
convertFromIdentifiersToIndices
(const std::vector<std::size_t>& vertexIdentifiers) {
   assert(_vertices.size() == vertexIdentifiers.size());

   //
   // Make the array of indexed simplices.
   //
   // Mapping from vertex identifiers to vertex indices.
   std::tr1::unordered_map<std::size_t, std::size_t> identifierToIndex;
   for (std::size_t i = 0; i != vertexIdentifiers.size(); ++i) {
      identifierToIndex[vertexIdentifiers[i]] = i;
   }
   // Convert to simplices of vertex indices.
   for (std::size_t i = 0; i != _indexedSimplices.size(); ++i) {
      for (std::size_t m = 0; m != M + 1; ++m) {
         _indexedSimplices[i][m] = identifierToIndex[_indexedSimplices[i][m]];
      }
   }
}

} // namespace geom
