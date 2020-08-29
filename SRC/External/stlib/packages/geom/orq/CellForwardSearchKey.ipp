// -*- C++ -*-

#if !defined(__geom_CellForwardSearchKey_ipp__)
#error This file is an implementation detail of the class CellForwardSearchKey.
#endif

namespace geom {


//
// Mathematical member functions
//


template < std::size_t N, typename _Record, typename _MultiKey, typename _Key,
         typename _MultiKeyAccessor >
template<typename OutputIterator>
inline
std::size_t
CellForwardSearchKey<N, _Record, _MultiKey, _Key, _MultiKeyAccessor>::
computeWindowQuery(OutputIterator iter, const BBox& window) const {
   //
   // Convert the multi-key to array indices.
   //
   IndexList mi, lo, hi;
   convertMultiKeyToIndices(window.lower, &lo);
   convertMultiKeyToIndices(window.upper, &hi);

   //
   // Truncate the index window to lie within the cell array.
   //
   for (std::size_t n = 0; n != N - 1; ++n) {
      lo[n] = std::max(Index(0), lo[n]);
      hi[n] = std::min(Index(getCellArrayExtents()[n] - 1), hi[n]);
   }

   // The interior portion of the index window.
   geom::BBox<Index, N - 1> interior = {lo + Index(1), hi - Index(1)};

   // The number of records in the window.
   std::size_t count = 0;
   typename Cell::const_iterator recordIterator, recordIteratorEnd;
   typename Cell::KeyConstIterator keyIterator, keyIteratorEnd;

   //
   // Iterate over the cells in the index window.
   //

   std::tr1::array<Index, N-1> indexPoint;
   const Key coord_min = window.lower[N - 1];
   const Key coord_max = window.upper[N - 1];
   std::size_t i;
   std::size_t n = N - 2;
   mi = lo;
   while (mi[N-2] <= hi[N-2]) {
      if (n == 0) {
         for (mi[0] = lo[0]; mi[0] <= hi[0]; ++mi[0]) {
            // Iterate over the records in the cell.
            const Cell& cell = getCell(mi);

            // Forward search to find the beginning of the records in window.
            i = cell.search(coord_min);
            recordIterator = cell.begin() + i;
            keyIterator = cell.getLastKeys().begin() + i;

            // Project to a lower dimension.
            for (std::size_t i = 0; i != indexPoint.size(); ++i) {
               indexPoint[i] = mi[i];
            }
            // If this is an interior cell.
            if (interior.isIn(indexPoint)) {
               keyIteratorEnd = cell.getLastKeys().end();
               for (; keyIterator != keyIteratorEnd && *keyIterator <= coord_max;
                     ++keyIterator, ++recordIterator) {
                  // There is no need to check if the record is in the window.
                  *iter = *recordIterator;
                  ++iter;
                  ++count;
               }
            }
            else { // This is a boundary cell.
               recordIteratorEnd = cell.end();
               for (; recordIterator != recordIteratorEnd &&
                     _multiKeyAccessor(*recordIterator)[N-1] <= coord_max;
                     ++recordIterator) {
                  if (window.isIn(_multiKeyAccessor(*recordIterator))) {
                     *iter = *recordIterator;
                     ++iter;
                     ++count;
                  }
               }
            }
         }
         ++n;
      }
      else if (mi[n-1] > hi[n-1]) {
         mi[n-1] = lo[n-1];
         ++mi[n];
         ++n;
      }
      else {
         --n;
      }
   }

   return count;
}

} // namespace geom
