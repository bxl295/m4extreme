// -*- C++ -*-

#if !defined(__geom_SortFirstDynamic_ipp__)
#error This file is an implementation detail of the class SortFirstDynamic.
#endif

namespace geom {


//
// Window queries.
//


template < std::size_t N, typename _Record, typename _MultiKey, typename _Key,
         typename _MultiKeyAccessor >
template<typename OutputIterator>
inline
std::size_t
SortFirstDynamic<N, _Record, _MultiKey, _Key, _MultiKeyAccessor>::
computeWindowQuery(OutputIterator output, const BBox& window) const {
   ConstIterator i = _records.lower_bound(window.lower[0]);
   const ConstIterator iEnd = _records.end();

   const Key upperBound = window.upper[0];
   std::size_t count = 0;
   for (; i != iEnd && i->first <= upperBound; ++i) {
      if (window.isIn(_multiKeyAccessor(i->second))) {
         *output++ = i->second;
         ++count;
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
SortFirstDynamic<N, _Record, _MultiKey, _Key, _MultiKeyAccessor>::
put(std::ostream& out) const {
   for (ConstIterator i = _records.begin(); i != _records.end(); ++i) {
      out << _multiKeyAccessor(i->second) << '\n';
   }
}

} // namespace geom
