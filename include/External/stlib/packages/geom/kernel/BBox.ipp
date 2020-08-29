// -*- C++ -*-

#if !defined(__geom_BBox_ipp__)
#error This file is an implementation detail of the class BBox.
#endif

namespace geom {

//
// Member functions
//

template<typename _T, std::size_t _Dim>
inline
bool
BBox<_T, _Dim>::
isEmpty() const {
   for (std::size_t i = 0; i != Dimension; ++i) {
      if (lower[i] > upper[i]) {
         return true;
      }
   }
   return false;
}


template<typename _T, std::size_t _Dim>
inline
bool
BBox<_T, _Dim>::
isIn(const Point& p) const {
   for (std::size_t i = 0; i != Dimension; ++i) {
      if (p[i] < lower[i] || upper[i] < p[i]) {
         return false;
      }
   }
   return true;
}


template<typename _T, std::size_t _Dim>
inline
void
BBox<_T, _Dim>::
add(const Point& p) {
   for (std::size_t i = 0; i != Dimension; ++i) {
      if (p[i] < lower[i]) {
         lower[i] = p[i];
      }
      if (p[i] > upper[i]) {
         upper[i] = p[i];
      }
   }
}


template<typename _T, std::size_t _Dim>
template<typename _InputIterator>
inline
void
BBox<_T, _Dim>::
bound(_InputIterator begin, _InputIterator end) {
   // If there are no objects, make the box invalid.
   if (begin == end) {
      invalidate();
      return;
   }
   // Bound the first object.
   *this = geom::bound(*begin++);
   // Add the rest of the objects.
   while (begin != end) {
      add(*begin++);
   }
}


//
// File I/O
//


template<typename _T, std::size_t _Dim>
inline
void
printFormatted(std::ostream& out, const BBox<_T, _Dim>& x) {
   for (std::size_t i = 0; i != _Dim; ++i) {
      if (i != 0) {
         out << " x ";
      }
      out << "[ " << x.lower[i] << " ... " << x.upper[i] << " ]";
   }
}

// CONTINUE: REMOVE
#if 0
// Compute a bounding box around the sequence of objects.
template<typename _InputIterator, std::size_t _Dim, typename _T>
inline
void
bound(_InputIterator begin, _InputIterator end, BBox<_T, _Dim>* box) {
   assert(begin != end);
   // Make a bounding box around the first object.
   bound(*begin++, box);
   // Add the rest of the objects.
   BBox<_T, _Dim> b;
   while (begin != end) {
      bound(*begin++, &b);
      box->add(b);
   }
}
#endif


//
// Mathematical free functions
//

template<typename _T, std::size_t _Dim>
inline
bool
doOverlap(const BBox<_T, _Dim>& a, const BBox<_T, _Dim>& b) {
   for (std::size_t i = 0; i != _Dim; ++i) {
      if (std::max(a.lower[i], b.lower[i]) >
          std::min(a.upper[i], b.upper[i])) {
         return false;
      }
   }
   return true;
}


// CONTINUE
#if 0
template<typename _T, std::size_t _Dim>
inline
bool
doOverlap(const BBox<_T, _Dim>& a, const SemiOpenInterval<N, T>& b) {
   for (std::size_t i = 0; i != N; ++i) {
      if (a.lower[i] > a.upper[i] ||
          b.lower[i] >= b.upper[i] ||
          a.lower[i] >= b.upper[i] ||
          b.lower[i] > a.upper[i]) {
         return false;
      }
   }
   return true;
}
#endif


// CONTINUE: REMOVE
#if 0
template<typename T>
inline
bool
isIn(const BBox<T, 3>& bb, typename BBox<T, 3>::point_type& p) {
   if (p[0] >= bb.lower[0] &&
       p[1] >= bb.lower[1] &&
       p[2] >= bb.lower[2] &&
       p[0] <= bb.upper[0] &&
       p[1] <= bb.upper[1] &&
       p[2] <= bb.upper[2]) {
      return true;
   }
   return false;
}
#endif


// CONTINUE: REMOVE
#if 0
template<typename T>
inline
bool
isIn(const Interval<3, T>& a, const Interval<3, T>& b) {
   if (b.lower[0] < a.lower[0] ||
       a.upper[0] < b.upper[0] ||
       b.lower[1] < a.lower[1] ||
       a.upper[1] < b.upper[1] ||
       b.lower[2] < a.lower[2] ||
       a.upper[2] < b.upper[2]) {
      return false;
   }
   return true;
}
#endif

// Return the squared distance between two 1-D intervals.
template<typename _T>
inline
_T
squaredDistanceBetweenIntervals(const _T lower1, const _T upper1,
                                const _T lower2, const _T upper2) {
   // Consider the intervals [l1..u1] and [l2..u2].
   // l1 u1 l2 u2
   if (upper1 < lower2) {
      return (upper1 - lower2) *(upper1 - lower2);
   }
   // l2 u2 l1 u2
   if (upper2 < lower1) {
      return (upper2 - lower1) *(upper2 - lower1);
   }
   return 0;
}


// Return the squared distance between two bounding boxes.
template<typename _T, std::size_t _Dim>
inline
_T
squaredDistance(const BBox<_T, _Dim>& x, const BBox<_T, _Dim>& y) {
   _T d2 = 0;
   for (std::size_t i = 0; i != _Dim; ++i) {
      d2 += squaredDistanceBetweenIntervals(x.lower[i], x.upper[i],
                                            y.lower[i], y.upper[i]);
   }
   return d2;
}


template<typename _T, std::size_t _Dim>
inline
BBox<_T, _Dim>
intersection(const BBox<_T, _Dim>& a, const BBox<_T, _Dim>& b) {
   BBox<_T, _Dim> box = {max(a.lower, b.lower), min(a.upper, b.upper)};
   return box;
}


// CONTINUE: Try to get rid of floor and ceil.
template<typename Index, typename MultiIndexOutputIterator>
inline
void
scanConvert(MultiIndexOutputIterator indices, const BBox<double, 3>& box) {
   // Make the index bounding box.
   BBox<Index, 3> ib = {{{Index(std::ceil(box.lower[0])),
                          Index(std::ceil(box.lower[1])),
                          Index(std::ceil(box.lower[2]))}},
                        {{Index(std::floor(box.upper[0])),
                          Index(std::floor(box.upper[1])),
                          Index(std::floor(box.upper[2]))}}};

   // Scan convert the integer bounding box.
   scanConvertIndex(indices, ib);
}


template<typename MultiIndexOutputIterator, typename Index>
inline
void
scanConvertIndex(MultiIndexOutputIterator indices, const BBox<Index, 3>& box) {

   const Index iStart = box.lower[0];
   const Index iEnd = box.upper[0];
   const Index jStart = box.lower[1];
   const Index jEnd = box.upper[1];
   const Index kStart = box.lower[2];
   const Index kEnd = box.upper[2];

   std::tr1::array<Index, 3> index;
   for (index[2] = kStart; index[2] <= kEnd; ++index[2]) {
      for (index[1] = jStart; index[1] <= jEnd; ++index[1]) {
         for (index[0] = iStart; index[0] <= iEnd; ++index[0]) {
            *indices++ = index;
         }
      }
   }
}


template<typename MultiIndexOutputIterator, typename Index>
inline
void
scanConvert(MultiIndexOutputIterator indices, const BBox<double, 3>& box,
            const BBox<Index, 3>& domain) {
   // Make the integer bounding box.
   BBox<Index, 3> ib = {{{Index(std::ceil(box.lower[0])),
                          Index(std::ceil(box.lower[1])),
                          Index(std::ceil(box.lower[2]))}},
                        {{Index(std::floor(box.upper[0])),
                          Index(std::floor(box.upper[1])),
                          Index(std::floor(box.upper[2]))}}};

   // Scan convert the integer bounding box on the specified domain.
   scanConvertIndex(indices, ib, domain);
}


template<typename MultiIndexOutputIterator, typename Index>
inline
void
scanConvertIndex(MultiIndexOutputIterator indices, const BBox<Index, 3>& box,
                 const BBox<Index, 3>& domain) {
   BBox<Index, 3> inter = intersection(box, domain);
   if (! inter.isEmpty()) {
      scanConvertIndex(indices, inter);
   }
}

// CONTINUE: Remove this.
#if 0
template<typename T, typename _IndexType>
inline
int
scanConvert(std::vector<_IndexType>& indexSet, const BBox<T, 3>& box) {
   const int old_size = indexSet.size();

   const int iStart = int(std::ceil(box.lower[0]));
   const int iEnd = int(std::floor(box.upper[0]));
   const int jStart = int(std::ceil(box.lower[1]));
   const int jEnd = int(std::floor(box.upper[1]));
   const int kStart = int(std::ceil(box.lower[2]));
   const int kEnd = int(std::floor(box.upper[2]));

   _IndexType index;
   for (index[2] = kStart; index[2] <= kEnd; ++index[2]) {
      for (index[1] = jStart; index[1] <= jEnd; ++index[1]) {
         for (index[0] = iStart; index[0] <= iEnd; ++index[0]) {
            indexSet.push_back(index);
         }
      }
   }
   return indexSet.size() - old_size;
}
#endif


// CONTINUE: Remove this.
#if 0
template<typename T, typename _IndexType>
inline
int
scanConvert(std::vector<_IndexType>& indexSet,
            const BBox<T, 3>& box, const BBox<T, 3>& domain) {
   BBox<T, 3> inter = intersection(box, domain);
   if (!inter.isEmpty()) {
      return scanConvert(indexSet, inter);
   }
   return 0;
}
#endif


} // namespace geom
