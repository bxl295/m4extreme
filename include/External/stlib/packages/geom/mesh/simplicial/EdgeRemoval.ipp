// -*- C++ -*-

#if !defined(__geom_mesh_simplicial_EdgeRemoval_ipp__)
#error This file is an implementation detail of the class EdgeRemoval.
#endif

namespace geom {

//
// Tetrahedralization.
//

template<class _QualityMetric, class _Point, typename _Number>
inline
bool
EdgeRemoval<_QualityMetric, _Point, _Number>::
solve() {
   // Calculate the best quality if we remove the edge.
   fillTables();
   // If the new quality is no better than the old quality with the edge.
   // CONTINUE
   typename array::MultiArray<Number, 2>::IndexList i = {{0, _ring.size() - 1}};
   if (_quality(i) <= computeQualityWithEdge() *
   (1.0 + std::sqrt(std::numeric_limits<Number>::epsilon()))) {
      return false;
   }

   // Otherwise, removing the edge improves the quality.
   // Build the triangulation of the ring.
   buildTriangles();
   return true;
}


//
// Private member functions.
//


template<class _QualityMetric, class _Point, typename _Number>
inline
void
EdgeRemoval<_QualityMetric, _Point, _Number>::
fillTables() {
   const std::size_t N = _ring.size();
   Number q;
   for (std::ptrdiff_t i = N - 3; i >= 0; --i) {
      for (std::ptrdiff_t j = i + 2; j != std::ptrdiff_t(N); ++j) {
         for (std::ptrdiff_t k = i + 1; k <= j - 1; ++k) {
            q = computeQuality(i, k, j);
            if (k < j - 1) {
               q = std::min(q, _quality(ext::make_array(k, j)));
            }
            if (k > i + 1) {
               q = std::min(q, _quality(ext::make_array(i, k)));
            }
            if (k == i + 1 || q > _quality(ext::make_array(i, j))) {
               _quality(ext::make_array(i, j)) = q;
               _index(ext::make_array(i, j)) = k;
            }
         }
      }
   }
}


template<class _QualityMetric, class _Point, typename _Number>
inline
void
EdgeRemoval<_QualityMetric, _Point, _Number>::
buildTriangles() {
   _triangles.clear();
   buildTrianglesRecurse(0, _ring.size() - 1);
#ifdef DEBUG_stlib
   assert(_triangles.size() == _ring.size() - 2);
#endif
}


template<class _QualityMetric, class _Point, typename _Number>
inline
void
EdgeRemoval<_QualityMetric, _Point, _Number>::
buildTrianglesRecurse(const std::size_t i, const std::size_t j) {
#ifdef DEBUG_stlib
   assert(0 <= i && i < j && j < _ring.size());
#endif

   if (i + 1 < j) {
      // Get the third index of the triangle.
      const std::size_t k = _index(ext::make_array<Index>(i, j));
#ifdef DEBUG_stlib
      assert(i < k && k < j);
#endif
      // Add the triangle.
      _triangles.push_back(ext::make_array(i, k, j));
      // Build the rest of the triangles.
      buildTrianglesRecurse(i, k);
      buildTrianglesRecurse(k, j);
   }
}


// Return the worse quality of the tetrahedra:
// _ring[i], _ring[k], _ring[j], _target
// and
// _source, _ring[i], _ring[k], _ring[j]
template<class _QualityMetric, class _Point, typename _Number>
inline
typename EdgeRemoval<_QualityMetric, _Point, _Number>::Number
EdgeRemoval<_QualityMetric, _Point, _Number>::
computeQuality(const std::size_t i, const std::size_t k, const std::size_t j)
const {
#ifdef DEBUG_stlib
   assert(0 <= i && i <= _ring.size() - 3);
   assert(j >= i + 2 && j < _ring.size());
   assert(k >= i + 1 && k <= j - 1);
#endif
   static Simplex tet;

   // Make a simplex from the first tetrahedra.
   tet = ext::make_array(_ring[i], _ring[k], _ring[j], _target);
   // Calculate the Jacobian.
   _qualityFunction.setFunction(tet);
   // Calculate the mean ratio quality function.
   const Number q1 = 1.0 / _qualityFunction();

   // Make a simplex from the first tetrahedra.
   tet = ext::make_array(_source, _ring[i], _ring[k], _ring[j]);
   // Calculate the Jacobian.
   _qualityFunction.setFunction(tet);
   // Calculate the mean ratio quality function.
   const Number q2 = 1.0 / _qualityFunction();

   return std::min(q1, q2);
}


// Return the quality of the complex with the center edge.
// The quality of the complex is the quality of the worst tetrahedron.
template<class _QualityMetric, class _Point, typename _Number>
inline
typename EdgeRemoval<_QualityMetric, _Point, _Number>::Number
EdgeRemoval<_QualityMetric, _Point, _Number>::
computeQualityWithEdge() const {
#ifdef DEBUG_stlib
   assert(_ring.size() >= 3);
#endif

   Number qualityNew;
   Number q = computeQuality(0);
   const std::size_t end = _ring.size() - 1;
   for (std::size_t i = 1; i != end; ++i) {
      qualityNew = computeQuality(i);
      if (qualityNew < q) {
         q = qualityNew;
      }
   }
   return q;
}


// Return the quality of the tetrahedra:
// _source, _target, _ring[i], _ring[i+1]
template<class _QualityMetric, class _Point, typename _Number>
inline
typename EdgeRemoval<_QualityMetric, _Point, _Number>::Number
EdgeRemoval<_QualityMetric, _Point, _Number>::
computeQuality(const std::size_t i) const {
#ifdef DEBUG_stlib
   assert(0 <= i && i <  _ring.size() - 1);
#endif

   static Simplex tet;

   // Make the simplex.
   tet = ext::make_array(_source, _target, _ring[i], _ring[i+1]);
   // Calculate the Jacobian.
   _qualityFunction.setFunction(tet);
   // REMOVE
   /*
   Number q = 1.0 / _qualityFunction();
   if (q < 0.1) {
     std::cout << "Low quality: " << q << '\n';
   }
   */
   // Return the quality.
   return 1.0 / _qualityFunction();
}

} // namespace geom
