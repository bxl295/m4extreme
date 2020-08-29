// -*- C++ -*-

#if !defined(__geom_SequentialScan_ipp__)
#error This file is an implementation detail of the class SequentialScan.
#endif

namespace geom {


//
// Mathematical member functions
//


template < std::size_t N, typename _Record, typename _MultiKey, typename _Key,
         typename _MultiKeyAccessor >
template<class OutputIterator>
inline
std::size_t
SequentialScan<N, _Record, _MultiKey, _Key, _MultiKeyAccessor>::
computeWindowQuery(OutputIterator output, const BBox& window) const {
   std::size_t count = 0;
   ConstIterator iter = _recordIterators.begin();
   const ConstIterator iterEnd = _recordIterators.end();
   for (; iter != iterEnd; ++iter) {
      if (window.isIn(_multiKeyAccessor(*iter))) {
         *(output++) = *iter;
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
SequentialScan<N, _Record, _MultiKey, _Key, _MultiKeyAccessor>::
put(std::ostream& out) const {
   for (ConstIterator i = _recordIterators.begin();
         i != _recordIterators.end(); ++i) {
      out << **i << '\n';
   }
}

} // namespace geom
