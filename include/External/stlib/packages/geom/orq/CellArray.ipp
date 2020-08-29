// -*- C++ -*-

#if !defined(__geom_CellArray_ipp__)
#error This file is an implementation detail of the class CellArray.
#endif

namespace geom {


//
// Memory usage.
//


template < std::size_t N, typename _Record, typename _MultiKey, typename _Key,
         typename _MultiKeyAccessor >
inline
std::size_t
CellArray<N, _Record, _MultiKey, _Key, _MultiKeyAccessor>::
getMemoryUsage() const {
   std::size_t usage = 0;
   for (typename DenseArray::const_iterator i = _cellArray.begin();
         i != _cellArray.end(); ++i) {
      usage += i->size() * sizeof(Record);
   }
   usage += _cellArray.size() * sizeof(Cell);
   return usage;
}


//
// Mathematical member functions
//


template < std::size_t N, typename _Record, typename _MultiKey, typename _Key,
         typename _MultiKeyAccessor >
template<class OutputIterator>
inline
std::size_t
CellArray<N, _Record, _MultiKey, _Key, _MultiKeyAccessor>::
computeWindowQuery(OutputIterator iter, const BBox& window) const {
   //
   // Convert the multi-key to array index coordinates.
   //

   // A multi-index and the lower and upper bounds of the index window.
   IndexList mi, lo, hi;
   convertMultiKeyToIndices(window.lower, &lo);
   convertMultiKeyToIndices(window.upper, &hi);

   //
   // Truncate the index window to lie within the cell array.
   //

   {
      for (std::size_t n = 0; n != N; ++n) {
         lo[n] = std::max(Index(0), lo[n]);
         hi[n] = std::min(Index(getExtents()[n] - 1), hi[n]);
      }
   }

   // The interior portion of the index window.
   Point lower, upper;
   for (std::size_t n = 0; n != N; ++n) {
      lower[n] = lo[n] + 1;
      upper[n] = hi[n] - 1;
   }
   BBox interior = {lower, upper};

   // The number of records in the window.
   std::size_t count = 0;
   typename Cell::const_iterator recordIter;
   typename Cell::const_iterator recordIterEnd;

   //
   // Iterate over the cells in the index window.
   //

   Point p;
   std::size_t n = N - 1;
   mi = lo;
   while (mi[N-1] <= hi[N-1]) {
      if (n == 0) {
         for (mi[0] = lo[0]; mi[0] <= hi[0]; ++mi[0]) {
            // Iterate over the records in the cell.
            const Cell& cell = _cellArray(mi);
            recordIter = cell.begin();
            recordIterEnd = cell.end();

            // Convert the index to a Cartesian point.
            for (std::size_t i = 0; i != p.size(); ++i) {
               p[i] = mi[i];
            }
            // If this is an interior cell.
            if (interior.isIn(p)) {
               for (; recordIter != recordIterEnd; ++recordIter) {
                  // No need to check if it is in the window.
                  *iter = *recordIter;
                  ++iter;
               }
               count += std::size_t(cell.size());
            }
            else {
               for (; recordIter != recordIterEnd; ++recordIter) {
                  // If the record is in the window.
                  if (window.isIn(_multiKeyAccessor(*recordIter))) {
                     *iter = *recordIter;
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


//
// File I/O
//


template < std::size_t N, typename _Record, typename _MultiKey, typename _Key,
         typename _MultiKeyAccessor >
inline
void
CellArray<N, _Record, _MultiKey, _Key, _MultiKeyAccessor>::
put(std::ostream& out) const {
   Base::put(out);

   for (typename DenseArray::const_iterator i = _cellArray.begin();
         i != _cellArray.end(); ++i) {
      const Cell& b = *i;
      typename Cell::const_iterator iter(b.begin());
      while (iter != b.end()) {
         out << _multiKeyAccessor(*(iter++)) << '\n';
      }
   }
}

} // namespace geom
