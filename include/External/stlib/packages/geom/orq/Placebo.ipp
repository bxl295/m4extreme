// -*- C++ -*-

#if !defined(__geom_Placebo_ipp__)
#error This file is an implementation detail of the class Placebo.
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
Placebo<N, _Record, _MultiKey, _Key, _MultiKeyAccessor>::
computeWindowQuery(OutputIterator iter, const BBox& /*window*/) const {
   typename std::vector<Record>::const_iterator recordIterator
   = _records.begin() + getStartingPoint();
   const typename std::vector<Record>::const_iterator recordEnd
   = recordIterator + _querySize;
   while (recordIterator != recordEnd) {
      *(iter++) = *(recordIterator++);
   }
   return _querySize;
}


//
// File I/O
//


template < std::size_t N, typename _Record, typename _MultiKey, typename _Key,
         typename _MultiKeyAccessor >
inline
void
Placebo<N, _Record, _MultiKey, _Key, _MultiKeyAccessor>::
put(std::ostream& out) const {
   for (typename std::vector<Record>::const_iterator i = _records.begin();
         i != _records.end(); ++i) {
      out << _multiKeyAccessor(*i) << '\n';
   }
}

} // namespace geom
