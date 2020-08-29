// -*- C++ -*-

#if !defined(__geom_FixedRadiusNeighborSearch_ipp__)
#error This file is an implementation detail of the class FixedRadiusNeighborSearch.
#endif

namespace geom {

//
// Mathematical member functions
//

template<std::size_t N, typename _Record, typename _MultiKey, typename _Key,
         typename _MultiKeyAccessor>
template<class _IndexOutputIterator>
inline
void
FixedRadiusNeighborSearch<N, _Record, _MultiKey, _Key, _MultiKeyAccessor>::
findNeighbors(_IndexOutputIterator iter, const std::size_t recordIndex) {
   // Make a bounding box that contains the sphere.
   const Point center = _multiKeyAccessor(_recordsBegin + recordIndex);
   Point p = center;
   p -= _radius;
   _boundingBox.lower = p;
   p = center;
   p += _radius;
   _boundingBox.upper = p;
   // Get the records that are in the box.
   _recordsInBox.clear();
   Base::computeWindowQuery(std::back_inserter(_recordsInBox), _boundingBox);
   // Select the records that are within the search radius (and are not the
   // specified record).
   std::size_t index;
   for (std::size_t i = 0; i != _recordsInBox.size(); ++i) {
      const Record& r = _recordsInBox[i];
      index = std::distance(_recordsBegin, r);
      if (index != recordIndex &&
            squaredDistance(center, _multiKeyAccessor(r)) <= _squaredRadius) {
         *iter++ = index;
      }
   }
}

} // namespace geom
