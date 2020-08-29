// -*- C++ -*-

#if !defined(__geom_SparseCellArray_ipp__)
#error This file is an implementation detail of the class SparseCellArray.
#endif

namespace geom {


//
// Constructors
//


template < std::size_t N, typename _Record, typename _MultiKey, typename _Key,
         typename _MultiKeyAccessor >
inline
void
SparseCellArray<N, _Record, _MultiKey, _Key, _MultiKeyAccessor>::
build() {
   typename VectorArray::SizeList arrayExtents;
   for (std::size_t n = 0; n != N - 1; ++n) {
      arrayExtents[n] = getExtents()[n];
   }
   _vectorArray.rebuild(arrayExtents);
}


//
// Memory usage.
//


template < std::size_t N, typename _Record, typename _MultiKey, typename _Key,
         typename _MultiKeyAccessor >
inline
std::size_t
SparseCellArray<N, _Record, _MultiKey, _Key, _MultiKeyAccessor>::
getMemoryUsage() const {
   std::size_t usage = 0;
   // The (N-1)-D array of sparse vectors.
   usage += _vectorArray.size()
            * sizeof(SparseCellVector<Cell>);
   // For each sparse vector.
   for (typename VectorArray::const_iterator i = _vectorArray.begin();
         i != _vectorArray.end(); ++i) {
      // The memory for the structure of the cell.
      usage += i->size() *
               sizeof(IndexAndCell<Cell>);
      for (typename SparseCellVector<Cell>::
            const_iterator j = i->begin();
            j != i->end();
            ++j) {
         // The memory for the contents of the cell.
         usage += j->cell.size() * sizeof(Record);
      }
   }
   return usage;
}


//
// Mathematical member functions
//


template < std::size_t N, typename _Record, typename _MultiKey, typename _Key,
         typename _MultiKeyAccessor >
template< class OutputIterator >
inline
std::size_t
SparseCellArray<N, _Record, _MultiKey, _Key, _MultiKeyAccessor>::
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

   for (std::size_t n = 0; n != N; ++n) {
      lo[n] = std::max(Index(0), lo[n]);
      hi[n] = std::min(Index(getExtents()[n] - 1), hi[n]);
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
   typename Cell::const_iterator recordIterator, recordIteratorEnd;
   typename SparseCellVector<Cell>::const_iterator cellIterator,
            cellIteratorEnd;

   //
   // Iterate over the cells in the index window.
   //

   Point p;
   // Array Index.
   std::tr1::array < typename IndexList::value_type, N - 1 > ai;
   std::size_t n = N - 2;
   mi = lo;
   for (std::size_t i = 0; i != N - 1; ++i) {
      ai[i] = mi[i];
   }
   while (mi[N-2] <= hi[N-2]) {
      if (n == 0) {
         for (mi[0] = lo[0]; mi[0] <= hi[0]; ++mi[0]) {
            ai[0] = mi[0];
            const SparseCellVector<Cell>& cell_vector = _vectorArray(ai);
            cellIterator = cell_vector.lower_bound(lo[N-1]);
            cellIteratorEnd = cell_vector.end();
            for (; cellIterator != cellIteratorEnd &&
                  Index(cellIterator->index) <= hi[N-1];
                  ++cellIterator) {
               const Cell& cell = cellIterator->cell;
               recordIterator = cell.begin();
               recordIteratorEnd = cell.end();

               mi[N-1] = cellIterator->index;
               // Convert the index to a Cartesian point.
               for (std::size_t i = 0; i != p.size(); ++i) {
                  p[i] = mi[i];
               }
               // If this is an interior cell.
               if (interior.isIn(p)) {
                  for (; recordIterator != recordIteratorEnd; ++recordIterator) {
                     // There is no need to check if the record is in the window.
                     *iter = *recordIterator;
                     ++iter;
                  }
                  count += cell.size();
               }
               else { // This is a boundary cell.
                  for (; recordIterator != recordIteratorEnd; ++recordIterator) {
                     if (window.isIn(_multiKeyAccessor(*recordIterator))) {
                        *iter = *recordIterator;
                        ++iter;
                        ++count;
                     }
                  }
               }
            }
         }
         ++n;
      }
      else if (mi[n-1] > hi[n-1]) {
         ai[n-1] = mi[n-1] = lo[n-1];
         ++mi[n];
         ++ai[n];
         ++n;
      }
      else {
         --n;
      }
   }

   return count;
}


// Indexing by multi-key.  Return a reference to a cell.
template < std::size_t N, typename _Record, typename _MultiKey, typename _Key,
         typename _MultiKeyAccessor >
inline
typename SparseCellArray<N, _Record, _MultiKey, _Key, _MultiKeyAccessor>::Cell&
SparseCellArray<N, _Record, _MultiKey, _Key, _MultiKeyAccessor>::
operator()(const MultiKey& multiKey) {
   // CONTINUE: this is not efficient.
   IndexList mi;
   convertMultiKeyToIndices(multiKey, &mi);
   std::tr1::array < typename IndexList::value_type, N - 1 > ai;
   for (std::size_t n = 0; n != N - 1; ++n) {
      ai[n] = mi[n];;
   }
   return _vectorArray(ai).find(mi[N-1]);
}


//
// File IO
//


template < std::size_t N, typename _Record, typename _MultiKey, typename _Key,
         typename _MultiKeyAccessor >
inline
void
SparseCellArray<N, _Record, _MultiKey, _Key, _MultiKeyAccessor>::
put(std::ostream& out) const {
   Base::put(out);

   typename Cell::const_iterator recordIterator;
   typename SparseCellVector<Cell>::const_iterator vectorIterator;

   typename VectorArray::const_iterator
   i = _vectorArray.begin(),
   iEnd = _vectorArray.end();
   for (; i != iEnd; ++i) {
      const SparseCellVector<Cell>& cellVector = *i;
      vectorIterator = cellVector.begin();
      for (; vectorIterator != cellVector.end(); ++vectorIterator) {
         const Cell& b = vectorIterator->cell;
         recordIterator = b.begin();
         while (recordIterator != b.end()) {
            out << _multiKeyAccessor(*(recordIterator++)) << '\n';
         }
      }
   }
}

} // namespace geom
