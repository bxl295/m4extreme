// -*- C++ -*-

#if !defined(__geom_SortRankProject_ipp__)
#error This file is an implementation detail of the class SortRankProject.
#endif

namespace geom {

//
// Constructors
//

template < std::size_t N, typename _Record, typename _MultiKey, typename _Key,
         typename _MultiKeyAccessor >
inline
SortRankProject<N, _Record, _MultiKey, _Key, _MultiKeyAccessor>::
SortRankProject() :
   Base(),
   _records(),
   _sorted(),
   _rank() {}


template < std::size_t N, typename _Record, typename _MultiKey, typename _Key,
         typename _MultiKeyAccessor >
inline
SortRankProject<N, _Record, _MultiKey, _Key, _MultiKeyAccessor>::
SortRankProject(const std::size_t size) :
   Base(),
   _records(),
   _sorted(),
   _rank() {
   _records.reserve(size);
   for (std::size_t n = 0; n != N; ++n) {
      _sorted[n].reserve(size);
      _rank[n].reserve(size);
   }
}


template < std::size_t N, typename _Record, typename _MultiKey, typename _Key,
         typename _MultiKeyAccessor >
template<typename InputIterator>
inline
SortRankProject<N, _Record, _MultiKey, _Key, _MultiKeyAccessor>::
SortRankProject(InputIterator first, InputIterator last) :
   Base(),
   _records(),
   _sorted(),
   _rank() {
   insert(first, last);
   sortAndRank();
}


//
// Memory usage.
//


template < std::size_t N, typename _Record, typename _MultiKey, typename _Key,
         typename _MultiKeyAccessor >
inline
std::size_t
SortRankProject<N, _Record, _MultiKey, _Key, _MultiKeyAccessor>::
getMemoryUsage() const {
   return (sizeof(SortRankProject) +
           _records.size() *(sizeof(Record) + N * sizeof(const_pointer) +
                             N * sizeof(int)));
}


//
// Mathematical member functions
//


template < std::size_t N, typename _Record, typename _MultiKey, typename _Key,
         typename _MultiKeyAccessor >
inline
void
SortRankProject<N, _Record, _MultiKey, _Key, _MultiKeyAccessor>::
sortAndRank() {
   //
   // Make the sorted pointer and rank vectors the right size.
   //
   for (std::size_t n = 0; n != N; ++n) {
      _sorted[n].resize(_records.size());
      _rank[n].resize(_records.size());
   }

   //
   // Initialize the vectors of sorted pointers.
   //
   const_pointer recordsBeginning = &*(_records.begin());
   for (std::size_t i = 0; i != getSize(); ++i) {
      for (std::size_t n = 0; n != N; ++n) {
         _sorted[n][i] = recordsBeginning + i;
      }
   }

   //
   // Sort in each direction.
   //
   LessThanCompare compare;
   for (std::size_t n = 0; n != N; ++n) {
      compare.set(n);
      std::sort(_sorted[n].begin(), _sorted[n].end(), compare);
   }

   //
   // Make the rank vectors.
   //
   for (std::size_t n = 0; n != N; ++n) {
      for (std::size_t i = 0; i != getSize(); ++i) {
         _rank[n][_sorted[n][i] - recordsBeginning] = i;
      }
   }
}


template < std::size_t N, typename _Record, typename _MultiKey, typename _Key,
         typename _MultiKeyAccessor >
template<typename OutputIterator>
inline
std::size_t
SortRankProject<N, _Record, _MultiKey, _Key, _MultiKeyAccessor>::
computeWindowQuery(OutputIterator iter, const BBox& window) const {
   //
   // Get each slice.
   //

   std::tr1::array<typename std::vector<const_pointer>::const_iterator, N>
   begin, end;
   std::tr1::array<std::size_t, N> size;
   LessThanComparePointerAndMultiKey comparePointerAndMultiKey;
   LessThanCompareMultiKeyAndPointer compareMultiKeyAndPointer;
   for (std::size_t n = 0; n != N; ++n) {
      comparePointerAndMultiKey.set(n);
      compareMultiKeyAndPointer.set(n);
      begin[n] = std::lower_bound(_sorted[n].begin(), _sorted[n].end(),
                                  window.lower,
                                  comparePointerAndMultiKey);
      end[n] = std::upper_bound(_sorted[n].begin(), _sorted[n].end(),
                                window.upper,
                                compareMultiKeyAndPointer);
      size[n] = end[n] - begin[n];
   }

   //
   // Get the intersection of the N slices.
   //

   // The ranks corresponding to the query window.
   geom::BBox<int, N> rankWindow;
   for (std::size_t n = 0; n != N; ++n) {
      rankWindow.lower[n] = begin[n] - _sorted[n].begin();
      rankWindow.upper[n] = end[n] - _sorted[n].begin() - 1;
   }

   // The smallest slice.
   const std::size_t smallestSliceIndex =
      std::min_element(size.begin(), size.end()) - size.begin();

   std::size_t count = 0;
   std::tr1::array<int, N> r;
   typename std::vector<const_pointer>::const_iterator i;
   const const_pointer recordsBeginning
   = &*(_records.begin());
   std::size_t index;
   // Iterate over the smallest slice.
   for (i = begin[smallestSliceIndex]; i != end[smallestSliceIndex]; ++i) {
      // The index of the record.
      index = *i - recordsBeginning;
      // Get the ranks of the record.
      for (std::size_t n = 0; n != N; ++n) {
         r[n] = _rank[n][index];
      }
      // If the ranks are in the query window.
      if (rankWindow.isIn(r)) {
         // Add the record.
         *(iter++) = **i;
         ++count;
      }
   }

   return count;
}


//
// File IO
//


template < std::size_t N, typename _Record, typename _MultiKey, typename _Key,
         typename _MultiKeyAccessor >
inline
void
SortRankProject<N, _Record, _MultiKey, _Key, _MultiKeyAccessor>::
put(std::ostream& out) const {
   for (typename std::vector<Record>::const_iterator i = _records.begin();
         i != _records.end(); ++i) {
      out << _multiKeyAccessor(*i) << '\n';
   }
}


//
// Validity check.
//


template < std::size_t N, typename _Record, typename _MultiKey, typename _Key,
         typename _MultiKeyAccessor >
inline
bool
SortRankProject<N, _Record, _MultiKey, _Key, _MultiKeyAccessor>::
isValid() const {
   LessThanCompare compare;
   for (std::size_t n = 0; n != N; ++n) {
      compare.set(n);
      if (! ads::is_sorted(_sorted[n].begin(), _sorted[n].end(), compare)) {
         return false;
      }
   }
   return true;
}

} // namespace geom
