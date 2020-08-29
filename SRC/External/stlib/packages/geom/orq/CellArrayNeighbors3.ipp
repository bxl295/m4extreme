// -*- C++ -*-

#if !defined(__geom_CellArrayNeighbors3_ipp__)
#error This file is an implementation detail of the class CellArrayNeighbors.
#endif

#ifdef __SSE4_1__
#include <smmintrin.h>
#elif defined(__SSE3__)
#include <pmmintrin.h>
#else
#error SSE3 or greater must be available.
#endif

namespace geom {

//
// 3-D specialization.
//

template<typename _Record, typename _Location>
class CellArrayNeighbors<3, _Record, _Location> {
   //
   // Constants.
   //
public:

   //! The space dimension.
   static const std::size_t D = 3;

   //
   // Types.
   //
public:

   //! A pointer to the record type.
   typedef _Record Record;
   //! A Cartesian point.
   typedef std::tr1::array<double, D> Point;

protected:

   //! The number type used internally.
   typedef float Number;
   //! The internal representation of a point.
   typedef std::tr1::array<Number, D> Pt;
   //! Bounding box.
   typedef geom::BBox<Number, D> BBox;
   //! A multi-index.
   typedef typename array::SimpleMultiArray<int, D>::IndexList IndexList;
   //! A single index.
   typedef typename array::SimpleMultiArray<int, D>::Index Index;

   //
   // Nested classes.
   //

private:

   struct RecLoc {
      Record record;
      __m128 location;
   };

   typedef typename std::vector<RecLoc>::const_iterator ConstIterator;

   //
   // Data
   //

private:

   //! The functor for computing record locations.
   _Location _location;
   //! The records along with the locations.
   std::vector<RecLoc> _recordData;
   //! Buffer for recording the records in a neighbor query.
   /*! We use this to avoid the possible inefficiencies of output iterators. */
   std::vector<Record> _recordsBuffer;
   //! The array of cells.
   /*! The array is padded by one empty slice in the first dimension. 
     This makes it easier to iterate over a range of cells. */
   array::SimpleMultiArray<ConstIterator, D> _cellArray;
   //! The Cartesian location of the lower corner of the cell array.
   Pt _lowerCorner;
   //! The inverse cell lengths.
   /*! This is used for converting locations to cell multi-indices. */
   Pt _inverseCellLengths;

   // SIMD vectors for locationToIndices().
   const __m128 _zero;
   __m128 _cellArrayUpper;
   __m128 _lc;
   __m128 _icl;

   // Scratch data for cellSort().
   std::vector<Index> _cellIndices;
   std::vector<std::size_t> _cellCounts;
   std::vector<RecLoc> _recordDataCopy;


private:

   //
   // Not implemented
   //

   // Copy constructor not implemented. _cellArray has iterators to _recordData
   // so the synthesized copy constructor would not work.
   CellArrayNeighbors(const CellArrayNeighbors&);

   // Assignment operator not implemented.
   CellArrayNeighbors&
   operator=(const CellArrayNeighbors&);

   //--------------------------------------------------------------------------
   //! \name Constructors.
   // @{
public:

   //! Construct from the search radius.
   CellArrayNeighbors(const _Location& location = _Location()) :
      _location(location),
      _recordData(),
      _cellArray(),
      // Fill with invalid values.
      _lowerCorner(ext::filled_array<Pt>
                          (std::numeric_limits<Number>::quiet_NaN())),
      _inverseCellLengths(ext::filled_array<Pt>
                          (std::numeric_limits<Number>::quiet_NaN())),
      _zero(_mm_set1_ps(0)) {
   }

   // @}
   //--------------------------------------------------------------------------
   //! \name Window Queries.
   // @{
public:

   //! Initialize with the given sequence of records.
   template<typename _InputIterator>
   void
   initialize(_InputIterator begin, _InputIterator end);

   //! Find the records that are in the specified ball.
   template<typename _OutputIterator>
   void
   neighborQuery(const Point& center, double radius, _OutputIterator neighbors);

   //! Find the records that are in the specified ball.
   /*! Store the neighbors in the supplied container, which will first be 
     cleared. */
   void
   neighborQuery(const Point& center, double radius,
                 std::vector<Record>* neighbors) {
      neighbors->clear();
      neighborQuery(center, radius, std::back_inserter(*neighbors));
   }

protected:

   //! Convert a location to a valid cell array multi-index.
   IndexList
   locationToIndices(const Pt& x);

   //! Convert a location to a container index.
   Index
   containerIndex(const Pt& x);

   //! Sort _recordData by cell container indices.
   void
   cellSort();

   //! Compute the array extents and the sizes for the cells.
   IndexList
   computeExtentsAndSizes(std::size_t numberOfCells, const BBox& domain);

   // @}
};


template<typename _Record, typename _Location>
template<typename _InputIterator>
inline
void
CellArrayNeighbors<3, _Record, _Location>::
initialize(_InputIterator begin, _InputIterator end) {
   // Clear results from previous queries.
   _recordData.clear();
   // Dispense with the trivial case.
   if (begin == end) {
      return;
   }

   // Copy the records and compute the locations.
   RecLoc rl;
   Pt p;
   while (begin != end) {
      rl.record = begin++;
      p = ext::convert_array<Number>(_location(rl.record));
      rl.location = _mm_set_ps(0, p[2], p[1], p[0]);
      _recordData.push_back(rl);
   }
   _recordsBuffer.resize(_recordData.size());

   // Bound the locations.
   _mm_store_ps(&p[0], _recordData[0].location);
   BBox box = {p, p};
   for (std::size_t i = 1; i != _recordData.size(); ++i) {
      _mm_store_ps(&p[0], _recordData[i].location);
      box.add(p);
   }

   // Expand to avoid errors in converting locations to cell multi-indices.
   box.offset(std::sqrt(std::numeric_limits<Number>::epsilon()) * 
              (1. + max(box.upper - box.lower)));
   // Set the lower corner and the cell array extents.
   _lowerCorner = box.lower;
   _cellArray.rebuild(computeExtentsAndSizes(_recordData.size(), box));
   // Sort by the cell indices and define the cells.
   cellSort();

   // Set data for locationToIndices().
   // Note that we adjust the first index to map to a non-empty cell.
   _cellArrayUpper = _mm_set_ps(0,
                                _cellArray.extents()[2] - 1,
                                _cellArray.extents()[1] - 1,
                                _cellArray.extents()[0] - 2);
   _lc = _mm_set_ps(0, _lowerCorner[2], _lowerCorner[1], _lowerCorner[0]);
   _icl = _mm_set_ps(0, _inverseCellLengths[2], _inverseCellLengths[1],
                     _inverseCellLengths[0]);
}


// 3-D specialization has little affect on performance.
template<typename _Record, typename _Location>
template<typename _OutputIterator>
inline
void
CellArrayNeighbors<3, _Record, _Location>::
neighborQuery(const Point& center_, const double radius_,
              _OutputIterator neighbors) {
   // Check trivial case.
   if (_recordData.empty()) {
      return;
   }
#ifdef DEBUG_stlib
   // Check that the records were initialized.
   assert(_lowerCorner[0] == _lowerCorner[0]);
#endif

   // The window for the ORQ has corners at center - radius and
   // center + radius. Convert the corners to cell array indices.
   const Pt center = ext::convert_array<Number>(center_);
   const Number radius = radius_;
   IndexList lo = locationToIndices(center - radius);
   IndexList hi = locationToIndices(center + radius);
   const Number squaredRadius = radius * radius;
   // Iterate over cells in all except the first dimension.
   IndexList index;
   const __m128 c = _mm_set_ps(0, center[2], center[1], center[0]);
   __m128 t;
   Record* r = &_recordsBuffer[0];
   const __m128 sr = _mm_set1_ps(squaredRadius);
   for (index[2] = lo[2]; index[2] <= hi[2]; ++index[2]) {
      for (index[1] = lo[1]; index[1] <= hi[1]; ++index[1]) {
         index[0] = hi[0] + 1;
         const ConstIterator jEnd = _cellArray(index);
         index[0] = lo[0];
         // Iterate over the records in the row.
         for (ConstIterator j = _cellArray(index); j != jEnd; ++j) {
            t = c - j->location;
#ifdef __SSE4_1__
            t = _mm_dp_ps(t, t, 0x71);
#else
            // If the dot product is not available, use two horizontal
            // additions.
            t = _mm_hadd_ps(t, t);
            t = _mm_hadd_ps(t, t);
#endif
            t = _mm_cmplt_ss(t, sr);
            *r = j->record;
            r += *reinterpret_cast<const unsigned*>(&t) & 0x1;
         }
      }
   }
   for (std::size_t i = 0; i != std::size_t(r - &_recordsBuffer[0]); ++i) {
      *neighbors++ = _recordsBuffer[i];
   }
}


template<typename _Record, typename _Location>
inline
typename CellArrayNeighbors<3, _Record, _Location>::IndexList
CellArrayNeighbors<3, _Record, _Location>::
locationToIndices(const Pt& x) {
   IndexList index;
#if 0
   for (std::size_t i = 0; i != D; ++i) {
      index[i] = std::min(_cellArray.extents()[i] - 1,
                          Index(std::max(Number(0), x[i] - _lowerCorner[i]) *
                                _inverseCellLengths[i]));
   }
   // Adjust to map to a non-empty cell.
   if (index[0] == _cellArray.extents()[0] - 1) {
      --index[0];
   }
#else
   // This is slightly faster.
   __m128 idx = _mm_set_ps(0, x[2], x[1], x[0]);
   idx = _mm_min_ps(_cellArrayUpper, _mm_max_ps(_zero, idx - _lc) * _icl);
   // This is a little faster than using _mm_cvttps_epi32().
   const float* i = reinterpret_cast<const float*>(&idx);
   index[0] = int(i[0]);
   index[1] = int(i[1]);
   index[2] = int(i[2]);
#endif
   return index;
}


template<typename _Record, typename _Location>
inline
typename CellArrayNeighbors<3, _Record, _Location>::Index
CellArrayNeighbors<3, _Record, _Location>::
containerIndex(const Pt& x) {
   IndexList index;
   for (std::size_t i = 0; i != D; ++i) {
      index[i] = Index((x[i] - _lowerCorner[i]) * _inverseCellLengths[i]);
   }
   return _cellArray.arrayIndex(index);
}


template<typename _Record, typename _Location>
inline
void
CellArrayNeighbors<3, _Record, _Location>::
cellSort() {
   // Calculate the cell container index for each record.
   _cellIndices.resize(_recordData.size());
   // Count the number of records in each cell.
   _cellCounts.resize(_cellArray.size());
   // Initialize to zero.
   Pt p;
   memset(&_cellCounts[0], 0, _cellCounts.size() * sizeof(std::size_t));
   for (std::size_t i = 0; i != _recordData.size(); ++i) {
      _mm_store_ps(&p[0], _recordData[i].location);
      _cellIndices[i] = containerIndex(p);
      ++_cellCounts[_cellIndices[i]];
   }
   // Turn the counts into end delimiters.
   std::partial_sum(_cellCounts.begin(), _cellCounts.end(),
                    _cellCounts.begin());
   // Copy the record data.
   _recordDataCopy.swap(_recordData);
   _recordData.resize(_recordDataCopy.size());
   // Define the cells.
   _cellArray[0] = _recordData.begin();
   for (std::size_t i = 1; i != _cellArray.size(); ++i) {
      _cellArray[i] = _recordData.begin() + _cellCounts[i-1];
   }
   // Sort by the cell indices.
   for (std::size_t i = 0; i != _recordData.size(); ++i) {
      _recordData[--_cellCounts[_cellIndices[i]]] = _recordDataCopy[i];
   }
}


// Compute the array extents and the sizes for the cells.
template<typename _Record, typename _Location>
inline
typename CellArrayNeighbors<3, _Record, _Location>::IndexList
CellArrayNeighbors<3, _Record, _Location>::
computeExtentsAndSizes(const std::size_t numberOfCells, const BBox& domain) {
   assert(numberOfCells > 0);

   // Work from the the least to greatest Cartesian extent to compute the
   // grid extents.
   IndexList extents;
   IndexList order;
   Pt ext = domain.upper - domain.lower;
   ads::computeOrder(ext.begin(), ext.end(), order.begin());
   for (std::size_t i = 0; i != D; ++i) {
      // Normalize the domain to numberOfCells content.
      Number content = product(ext);
      assert(content != 0);
      const Number factor = std::pow(numberOfCells / content,
                                     Number(1.0 / (D - i)));
      for (std::size_t j = i; j != D; ++j) {
         ext[order[j]] *= factor;
      }
      // The current index;
      const std::size_t n = order[i];
      // Add 0.5 and truncate to round to the nearest integer.
      ext[n] = extents[n] = std::max(std::size_t(1),
                                     std::size_t(ext[n] + 0.5));
   }

   // From the domain and the cell array extents, compute the cell size.
   for (std::size_t n = 0; n != D; ++n) {
      const Number d = domain.upper[n] - domain.lower[n];
      if (d == 0 || extents[n] == 1) {
         // The cell covers the entire domain.
         _inverseCellLengths[n] = 0;
      }
      else {
         _inverseCellLengths[n] = extents[n] / d;
      }
   }

   // Add one cell in the first dimension so that we will have convenient 
   // end iterators when doing the window queries.
   ++extents[0];

   return extents;
}


} // namespace geom
