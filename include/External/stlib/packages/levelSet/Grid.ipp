// -*- C++ -*-

#if !defined(__levelSet_Grid_ipp__)
#error This file is an implementation detail of Grid.
#endif

namespace levelSet {


// Construct from the Cartesian domain and the suggested grid patch spacing.
template<typename _T, std::size_t _D, std::size_t _N, typename _R>
inline
Grid<_T, _D, _N, _R>::
Grid(const BBox& domain, const _R targetSpacing) :
   Base(calculateExtents(domain, targetSpacing)),
   // Don't allocate memory for the patch arrays until the grid is refined.
   lowerCorner(domain.lower),
   // Choose the maximum over the dimensions in order to cover the domain.
   spacing(max((domain.upper - domain.lower) / (Base::_extents * _N - 1))),
   _data() {
   // Ensure that the patch extent is a power of 2.
   // http://graphics.stanford.edu/~seander/bithacks.html#DetermineIfPowerOf2
   LOKI_STATIC_CHECK(_N && !(_N & (_N - 1)), PatchExtentMustBeAPowerOf2);
   assert(spacing <= targetSpacing);
}


// Construct from the Cartesian domain and the grid extents.
template<typename _T, std::size_t _D, std::size_t _N, typename _R>
inline
Grid<_T, _D, _N, _R>::
Grid(const BBox& domain, const typename Base::IndexList& extents) :
   Base(extents),
   // Don't allocate memory for the patch arrays until the grid is refined.
   lowerCorner(domain.lower),
   // Choose the maximum over the dimensions in order to cover the domain.
   spacing(max((domain.upper - domain.lower) / (Base::_extents * _N - 1))),
   _data() {
   // Ensure that the patch extent is a power of 2.
   // http://graphics.stanford.edu/~seander/bithacks.html#DetermineIfPowerOf2
   LOKI_STATIC_CHECK(_N && !(_N & (_N - 1)), PatchExtentMustBeAPowerOf2);
}


template<typename _T, std::size_t _D, std::size_t _N, typename _R>
inline
typename Grid<_T, _D, _N, _R>::Base::IndexList
Grid<_T, _D, _N, _R>::
calculateExtents(const BBox& domain, const _R targetSpacing) {
   // Number of cells = ext * _N. Note that we don't count the layer of
   // boundary vertices along the upper sides.
   typename Base::IndexList extents;
   for (std::size_t i = 0; i != extents.size(); ++i) {
      // Include a fudge factor for the length.
      const _R length = (domain.upper[i] - domain.lower[i]) * 
         (1 + std::numeric_limits<_R>::epsilon());
      // length = dx * (extents * _N - 1)
      // extents = (length / dx + 1) / _N
      extents[i] = std::size_t(std::ceil((length / targetSpacing + 1)
                                         / _N));
      assert(extents[i] != 0);
      // Ensure that there are at least two grid points.
      if (_N == 1 && extents[i] == 1) {
         extents[i] = 2;
      }
   }
   return extents;
}


// Return the value at the specified grid point.
template<typename _T, std::size_t _D, std::size_t _N, typename _R>
inline
_T
Grid<_T, _D, _N, _R>::
operator()(const typename Base::IndexList& patch,
           const typename Base::IndexList& index) const {
   if (Base::operator()(patch).isRefined()) {
      return Base::operator()(patch)(index);
   }
   return Base::operator()(patch).fillValue;
}


// Return the Cartesian domain spanned by the grid.
template<typename _T, std::size_t _D, std::size_t _N, typename _R>
inline
typename Grid<_T, _D, _N, _R>::BBox
Grid<_T, _D, _N, _R>::
domain() const {
   BBox domain = {lowerCorner,
                  lowerCorner + (Base::extents() * _N - 1) * spacing};
   return domain;
}


// Return the domain for the specified vertex patch.
template<typename _T, std::size_t _D, std::size_t _N, typename _R>
inline
typename Grid<_T, _D, _N, _R>::BBox
Grid<_T, _D, _N, _R>::
getVertexPatchDomain(const typename Base::IndexList& i) const {
   const Point lower = lowerCorner + getVoxelPatchLength() * i;
   const BBox domain = {lower, lower + spacing * (_N - 1)};
   return domain;
}


// Return the lower corner of the Cartesian domain for the specified patch.
template<typename _T, std::size_t _D, std::size_t _N, typename _R>
inline
typename Grid<_T, _D, _N, _R>::Point
Grid<_T, _D, _N, _R>::
getPatchLowerCorner(const typename Base::IndexList& i) const {
   return lowerCorner + getVoxelPatchLength() * i;
}

template<typename _T, std::size_t _N, typename _R>
inline
void
_getVoxelPatch(const Grid<_T, 1, _N, _R>& grid,
              typename Grid<_T, 1, _N, _R>::IndexList index,
              typename Grid<_T, 1, _N, _R>::VoxelPatch* patch) {
   typedef Grid<_T, 1, _N, _R> Grid;
   typedef typename Grid::VertexPatch VertexPatch;

   const VertexPatch& p = grid(index);
   assert(p.isRefined());
   for (std::size_t i = 0; i != _N; ++i) {
      (*patch)[i] = p[i];
   }
   (*patch)[_N] = (*patch)[_N-1];
   ++index[0];
   if (index[0] != grid.extents()[0] && grid(index).isRefined()) {
      (*patch)[_N] = grid(index)[0];
   }
}


template<typename _T, std::size_t _N, typename _R>
inline
void
_getVoxelPatch(const Grid<_T, 2, _N, _R>& grid,
              typename Grid<_T, 2, _N, _R>::IndexList index,
              typename Grid<_T, 2, _N, _R>::VoxelPatch* patch) {
   const std::size_t Dimension = 2;
   typedef Grid<_T, Dimension, _N, _R> Grid;
   typedef typename Grid::IndexList IndexList;
   typedef typename Grid::VertexPatch VertexPatch;

   IndexList i, j;
   {
      const VertexPatch& p = grid(index);
      assert(p.isRefined());
      // Body.
      for (i[0] = 0; i[0] != _N; ++i[0]) {
         for (i[1] = 0; i[1] != _N; ++i[1]) {
            (*patch)(i) = p(i);
         }
      }
      // Sides.
      for (std::size_t a = 0; a != Dimension; ++a) {
         const std::size_t b = (a + 1) % Dimension;
         i[a] = _N;
         j[a] = _N - 1;
         for (i[b] = 0; i[b] != _N; ++i[b]) {
            j[b] = i[b];
            (*patch)(i) = p(j);
         }
      }
      // Corner.
      (*patch)(_N, _N) = p(_N-1, _N-1);
   }

   // Sides.
   for (std::size_t a = 0; a != Dimension; ++a) {
      ++index[a];
      if (index[a] != grid.extents()[a] && grid(index).isRefined()) {
         const VertexPatch& p = grid(index);
         const std::size_t b = (a + 1) % Dimension;
         i[a] = _N;
         j[a] = 0;
         for (i[b] = 0; i[b] != _N; ++i[b]) {
            j[b] = i[b];
            (*patch)(i) = p(j);
         }
      }
      --index[a];
   }

   // Corner.
   ++index[0];
   ++index[1];
   if (index[0] != grid.extents()[0] && index[1] != grid.extents()[1] &&
       grid(index).isRefined()) {
      const VertexPatch& p = grid(index);
      (*patch)(_N, _N) = p(0, 0);
   }
}


template<typename _T, std::size_t _N, typename _R>
inline
void
_getVoxelPatch(const Grid<_T, 3, _N, _R>& grid,
              typename Grid<_T, 3, _N, _R>::IndexList index,
              typename Grid<_T, 3, _N, _R>::VoxelPatch* patch) {
   const std::size_t Dimension = 3;
   typedef Grid<_T, Dimension, _N, _R> Grid;
   typedef typename Grid::IndexList IndexList;
   typedef typename Grid::VertexPatch VertexPatch;

   IndexList i, j;
   {
      const VertexPatch& p = grid(index);
      assert(p.isRefined());
      // Body.
      for (i[0] = 0; i[0] != _N; ++i[0]) {
         for (i[1] = 0; i[1] != _N; ++i[1]) {
            for (i[2] = 0; i[2] != _N; ++i[2]) {
               (*patch)(i) = p(i);
            }
         }
      }
      // Sides.
      for (std::size_t a = 0; a != Dimension; ++a) {
         const std::size_t b = (a + 1) % Dimension;
         const std::size_t c = (a + 2) % Dimension;
         i[a] = _N;
         j[a] = _N - 1;
         for (i[b] = 0; i[b] != _N; ++i[b]) {
            for (i[c] = 0; i[c] != _N; ++i[c]) {
               j[b] = i[b];
               j[c] = i[c];
               (*patch)(i) = p(j);
            }
         }
      }
      // Edges.
      for (std::size_t a = 0; a != Dimension; ++a) {
         const std::size_t b = (a + 1) % Dimension;
         const std::size_t c = (a + 2) % Dimension;
         i[b] = i[c] = _N;
         j[b] = j[c] = _N - 1;
         for (i[a] = 0; i[a] != _N; ++i[a]) {
            j[a] = i[a];
            (*patch)(i) = p(j);
         }
      }
      // Corner.
      (*patch)(_N, _N, _N) = p(_N-1, _N-1, _N-1);
   }

   // Sides.
   for (std::size_t a = 0; a != Dimension; ++a) {
      ++index[a];
      if (index[a] != grid.extents()[a] && grid(index).isRefined()) {
         const VertexPatch& p = grid(index);
         const std::size_t b = (a + 1) % Dimension;
         const std::size_t c = (a + 2) % Dimension;
         i[a] = _N;
         j[a] = 0;
         for (i[b] = 0; i[b] != _N; ++i[b]) {
            for (i[c] = 0; i[c] != _N; ++i[c]) {
               j[b] = i[b];
               j[c] = i[c];
               (*patch)(i) = p(j);
            }
         }
      }
      --index[a];
   }
   // Edges.
   for (std::size_t a = 0; a != Dimension; ++a) {
      const std::size_t b = (a + 1) % Dimension;
      const std::size_t c = (a + 2) % Dimension;
      ++index[b];
      ++index[c];
      if (index[b] != grid.extents()[b] && index[c] != grid.extents()[c] &&
          grid(index).isRefined()) {
         const VertexPatch& p = grid(index);
         i[b] = i[c] = _N;
         j[b] = j[c] = 0;
         for (i[a] = 0; i[a] != _N; ++i[a]) {
            j[a] = i[a];
            (*patch)(i) = p(j);
         }
      }
      --index[b];
      --index[c];
   }
   // Corner.
   ++index[0];
   ++index[1];
   ++index[2];
   if (index[0] != grid.extents()[0] && index[1] != grid.extents()[1] &&
       index[2] != grid.extents()[2] && grid(index).isRefined()) {
      const VertexPatch& p = grid(index);
      (*patch)(_N, _N, _N) = p(0, 0, 0);
   }
}


// Get the specified voxel patch.
template<typename _T, std::size_t _D, std::size_t _N, typename _R>
inline
void
Grid<_T, _D, _N, _R>::
getVoxelPatch(const typename Base::IndexList& i, VoxelPatch* patch) const {
   // A level of indirection for the space dimension.
   _getVoxelPatch(*this, i, patch);
}


// Return the Cartesian position of the specified vertex.
template<typename _T, std::size_t _D, std::size_t _N, typename _R>
inline
typename Grid<_T, _D, _N, _R>::Point
Grid<_T, _D, _N, _R>::
indexToLocation(const typename Base::IndexList& patch,
                const typename Base::IndexList& index) const {
   Point x;
   for (std::size_t d = 0; d != Dimension; ++d) {
      x[d] = lowerCorner[d] + patch[d] * getVoxelPatchLength() +
         index[d] * spacing;
   }
   return x;
}


// Get the number of refined patches.
template<typename _T, std::size_t _D, std::size_t _N, typename _R>
inline
std::size_t
Grid<_T, _D, _N, _R>::
numRefined() const {
   std::size_t count = 0;
   for (std::size_t i = 0; i != Base::size(); ++i) {
      count += Base::operator[](i).isRefined();
   }
   return count;
}


// Report the specified set of grid points as patch/grid multi-index pairs.
template<typename _T, std::size_t _D, std::size_t _N, typename _R>
template<typename _OutputIterator>
inline
void
Grid<_T, _D, _N, _R>::
report(const typename Base::IndexList& patch,
       const typename Base::Range& range,
       _OutputIterator neighbors) const {
   typedef typename Base::IndexList IndexList;
   typedef array::SimpleMultiIndexRangeIterator<_D> Iterator;
   std::pair<IndexList, IndexList> dual;
   dual.first = patch;
   const Iterator end = Iterator::end(range);
   for (Iterator i = Iterator::begin(range); i != end; ++i) {
      dual.second = *i;
      *neighbors++ = dual;
   }
}


// Report the set of adjacent neighbors in the specified direction.
template<typename _T, std::size_t _D, std::size_t _N, typename _R>
template<typename _OutputIterator>
inline
void
Grid<_T, _D, _N, _R>::
adjacentNeighbors(const typename Base::IndexList& patch,
                  const typename Base::IndexList& index,
                  const std::size_t direction,
                  _OutputIterator neighbors) const {
#ifdef DEBUG_stlib
   assert(direction < 2 * _D);
#endif
   typedef typename Base::IndexList IndexList;
   typedef typename Base::Range Range;
   const IndexList NullIndex = 
      ext::filled_array<IndexList>(std::numeric_limits<std::size_t>::max());
   const std::size_t coordinate = direction / 2;
   IndexList p = patch;
   IndexList i = index;
   // Refined.
   if (Base::operator()(patch).isRefined()) {
      // Negative direction.
      if (direction % 2 == 0) {
         // On lower face.
         if (index[coordinate] == 0) {
            if (patch[coordinate] != 0) {
               --p[coordinate];
               if (Base::operator()(p).isRefined()) {
                  i[coordinate] = _N - 1;
                  *neighbors++ = std::make_pair(p, i);
               }
               else {
                  *neighbors++ = std::make_pair(p, NullIndex);
               }
            }
         }
         else {
            --i[coordinate];
            *neighbors++ = std::make_pair(p, i);
         }
      }
      // Positive direction.
      else {
         // On upper face.
         if (index[coordinate] == _N - 1) {
            if (patch[coordinate] != Base::extents()[coordinate] - 1) {
               ++p[coordinate];
               if (Base::operator()(p).isRefined()) {
                  i[coordinate] = 0;
                  *neighbors++ = std::make_pair(p, i);
               }
               else {
                  *neighbors++ = std::make_pair(p, NullIndex);
               }
            }
         }
         else {
            ++i[coordinate];
            *neighbors++ = std::make_pair(p, i);
         }
      }
   }
   // Unrefined.
   else {
      // Negative direction.
      if (direction % 2 == 0) {
         if (patch[coordinate] != 0) {
            --p[coordinate];
            // Neighbor refined.
            if (Base::operator()(p).isRefined()) {
               IndexList extents = ext::filled_array<IndexList>(_N);
               extents[coordinate] = 1;
               IndexList bases = ext::filled_array<IndexList>(0);
               bases[coordinate] = _N - 1;
               Range range = {extents, bases};
               report(p, range, neighbors);
            }
            // Neighbor unrefined.
            else {
               *neighbors++ = std::make_pair(p, NullIndex);
            }
         }
      }
      // Positive direction.
      else {
         if (patch[coordinate] != Base::extents()[coordinate] - 1) {
            ++p[coordinate];
            // Neighbor refined.
            if (Base::operator()(p).isRefined()) {
               IndexList extents = ext::filled_array<IndexList>(_N);
               extents[coordinate] = 1;
               IndexList bases = ext::filled_array<IndexList>(0);
               Range range = {extents, bases};
               report(p, range, neighbors);
            }
            // Neighbor unrefined.
            else {
               *neighbors++ = std::make_pair(p, NullIndex);
            }
         }
      }
   }
}


// Report the set of all adjacent neighbors.
template<typename _T, std::size_t _D, std::size_t _N, typename _R>
template<typename _OutputIterator>
inline
void
Grid<_T, _D, _N, _R>::
adjacentNeighbors(const DualIndices& pair, _OutputIterator neighbors) const {
   for (std::size_t direction = 0; direction != 2 * _D; ++direction) {
      adjacentNeighbors(pair.first, pair.second, direction, neighbors);
   }
}


// Report the set of all neighbors.
template<typename _T, std::size_t _D, std::size_t _N, typename _R>
template<typename _OutputIterator>
inline
void
Grid<_T, _D, _N, _R>::
allNeighbors(const DualIndices& pair, _OutputIterator neighbors) const {
   typedef typename Base::IndexList IndexList;
   typedef typename Base::Range Range;
   typedef array::SimpleMultiIndexRangeIterator<_D> Iterator;
   const IndexList NullIndex = 
      ext::filled_array<IndexList>(std::numeric_limits<std::size_t>::max());

   const IndexList& patch = pair.first;
   const IndexList& index = pair.second;

   // Note that in this function we will work with global indices.
   // g = p * _N + i.
   IndexList lower, upper;
   DualIndices di;
   // Refined.
   if (Base::operator()(patch).isRefined()) {
      // Define a range that contains the neighbors of the point.
      for (std::size_t d = 0; d != _D; ++d) {
         lower[d] = patch[d] * _N + index[d];
         if (lower[d] != 0) {
            --lower[d];
         }
         if (patch[d] == Base::extents()[d] - 1 && index[d] == _N - 1) {
            upper[d] = Base::extents()[d] * _N;
         }
         else {
            upper[d] = patch[d] * _N + index[d] + 2;
         }
      }
      // Loop over the range. Note that unrefined patches may be reported 
      // multiple times.
      const Range range = {upper - lower, lower};
      const Iterator end = Iterator::end(range);
      for (Iterator i = Iterator::begin(range); i != end; ++i) {
         di.first = *i / _N;
         di.second = *i % _N;
         // Skip the grid point itself.
         if (di == pair) {
            continue;
         }
         *neighbors++ = di;
      }
   }
   // Unrefined.
   else {
      // Define a range that contains the neighbors of the patch.
      for (std::size_t d = 0; d != _D; ++d) {
         if (patch[d] == 0) {
            lower[d] = 0;
         }
         else {
            lower[d] = patch[d] * _N - 1;
         }
         if (patch[d] == Base::extents()[d] - 1) {
            upper[d] = Base::extents()[d] * _N;
         }
         else {
            upper[d] = (patch[d] + 1) * _N + 1;
         }
      }
      // Loop over the range and skip the patch itself.
      const Range range = {upper - lower, lower};
      const Iterator end = Iterator::end(range);
      for (Iterator i = Iterator::begin(range); i != end; ++i) {
         di.first = *i / _N;
         if (di.first == patch) {
            continue;
         }
         di.second = *i % _N;
         *neighbors++ = di;
      }
   }
}


// Return true if the grid is valid.
template<typename _T, std::size_t _D, std::size_t _N, typename _R>
inline
bool
Grid<_T, _D, _N, _R>::
isValid() const {
   // Check each patch.
   for (std::size_t i = 0; i != Base::size(); ++i) {
      if (! Base::operator[](i).isValid()) {
         return false;
      }
   }
   if (spacing <= 0) {
      return false;
   }
   // Check the memory.
   typename std::vector<_T>::const_pointer p = &_data[0];
   for (std::size_t i = 0; i != Base::size(); ++i) {
      if (Base::operator[](i).isRefined()) {
         if (Base::operator[](i).data() != p) {
            return false;
         }
         p += NumVerticesPerPatch;
      }
   }
   if (p != &_data[0] + _data.size()) {
      return false;
   }
   return true;
}


// Return a reference to the value at the specified grid point.
template<typename _T, std::size_t _D, std::size_t _N, typename _R>
inline
_T&
Grid<_T, _D, _N, _R>::
operator()(const typename Base::IndexList& patch,
           const typename Base::IndexList& index) {
   if (Base::operator()(patch).isRefined()) {
      return Base::operator()(patch)(index);
   }
   return Base::operator()(patch).fillValue;
}


// Make all of the grids unrefined.
template<typename _T, std::size_t _D, std::size_t _N, typename _R>
inline
void
Grid<_T, _D, _N, _R>::
clear() {
   for (std::size_t i = 0; i != Base::size(); ++i) {
      Base::operator[](i).clear();
   }
   _data.clear();
}


// Refine the specified grids.
template<typename _T, std::size_t _D, std::size_t _N, typename _R>
inline
void
Grid<_T, _D, _N, _R>::
refine(const std::vector<std::size_t>& indices) {
   // First clear the grid.
   clear();
   // Then refine the specified patches.
   _data.resize(indices.size() * NumVerticesPerPatch);
   // CONTINUE: Fix for integer value types.
   std::fill(_data.begin(), _data.end(),
             std::numeric_limits<_T>::quiet_NaN());
   std::size_t dataIndex = 0;
   for (std::size_t i = 0; i != indices.size(); ++i) {
      const std::size_t index = indices[i];
      Base::operator[](index).refine(&_data[dataIndex]);
      dataIndex += NumVerticesPerPatch;
   }
   // Check that we used all of the data.
   assert(dataIndex == _data.size());
}


// Refine the patches that have one or more dependencies.
template<typename _T, std::size_t _D, std::size_t _N, typename _R>
inline
void
Grid<_T, _D, _N, _R>::
refine(const array::StaticArrayOfArrays<std::size_t>& dependencies) {
   std::vector<std::size_t> patchIndices;
   for (std::size_t i = 0; i != dependencies.getNumberOfArrays(); ++i) {
      if (! dependencies.empty(i)) {
         patchIndices.push_back(i);
      }
   }
   // Refine the appropriate patches and set the rest to have an unknown
   // distance.
   refine(patchIndices);
}


template<typename _T, std::size_t _D, std::size_t _N, typename _R>
inline
void
Grid<_T, _D, _N, _R>::
coarsen() {
   // Pack the array.
   std::size_t target = 0;
   // For each refined patch.
   for (std::size_t i = 0; i != Base::size(); ++i) {
      VertexPatch& patch = Base::operator[](i);
      if (! patch.isRefined()) {
         continue;
      }
      if (patch.shouldBeCoarsened()) {
         patch.coarsen();
      }
      else {
         // If the patch data needs to be moved.
         if (patch.data() != &_data[target]) {
            std::copy(patch.begin(), patch.end(), &_data[target]);
         }
         patch.refine(&_data[target]);
         target += NumVerticesPerPatch;
      }
   }

   // Resize.
   _data.resize(target);
}


// Add a constant to all vertices and fill values.
template<typename _T, std::size_t _D, std::size_t _N, typename _R>
inline
Grid<_T, _D, _N, _R>&
operator+=(Grid<_T, _D, _N, _R>& grid, const _T x) {
   for (std::size_t i = 0; i != grid.size(); ++i) {
      if (grid[i].isRefined()) {
         grid[i] += x;
      }
      else {
         grid[i].fillValue += x;
      }
   }
   return grid;
}


// Subtract a constant from all vertices and fill values.
template<typename _T, std::size_t _D, std::size_t _N, typename _R>
inline
Grid<_T, _D, _N, _R>&
operator-=(Grid<_T, _D, _N, _R>& grid, const _T x) {
   for (std::size_t i = 0; i != grid.size(); ++i) {
      if (grid[i].isRefined()) {
         grid[i] -= x;
      }
      else {
         grid[i].fillValue -= x;
      }
   }
   return grid;
}


// Determine the objects whose bounding boxes overlap each patch.
template<typename _T, std::size_t _D, std::size_t _N, typename _R,
         typename _InputIterator>
inline
void
patchDependencies(const Grid<_T, _D, _N, _R>& grid, _InputIterator begin,
                  _InputIterator end, 
                  array::StaticArrayOfArrays<std::size_t>* dependencies) {
   typedef Grid<_T, _D, _N, _R> Grid;
   typedef typename Grid::Index Index;
   typedef typename Grid::IndexList IndexList;
   typedef array::SimpleMultiIndexRange<_D> Range;
   typedef array::SimpleMultiIndexRangeIterator<_D> Iterator;
   typedef std::pair<std::size_t, std::size_t> PatchObjectPair;

   // The inverse length of a voxel patch.
   const _R inverseLength = 1. / grid.getVoxelPatchLength();

   geom::BBox<_R, _D> box;
   IndexList lower, upper;
   Range range;
   PatchObjectPair p;
   // The patch/object dependency pairs.
   std::vector<PatchObjectPair> pairs;
   // For each object, record the patches that depend on it.
   for (std::size_t objectIndex = 0; begin != end; ++begin, ++objectIndex) {
      // Make a bounding box around the object.
      box = bound(*begin);
      //
      // Convert the Cartesian coordinates to index coordinates.
      //
      box.lower -= grid.lowerCorner;
      box.lower *= inverseLength;
      for (std::size_t i = 0; i != _D; ++i) {
         // Round down. Use only voxels in the grid range.
         lower[i] = Index(std::max(std::floor(box.lower[i]), _R(0.)));
      }
      box.upper -= grid.lowerCorner;
      box.upper *= inverseLength;
      for (std::size_t i = 0; i != _D; ++i) {
         // Round up for an open upper bound. Use only voxels in the grid range.
         upper[i] = Index(std::max(std::min(std::ceil(box.upper[i]),
                                            _R(grid.extents()[i])),
                                   _R(lower[i])));
      }
      range.extents = upper - lower;
      range.bases = lower;

      // Loop over the range.
      p.second = objectIndex;
      const Iterator end = Iterator::end(range);
      for (Iterator i = Iterator::begin(range); i != end; ++i) {
         p.first = grid.arrayIndex(*i);
         pairs.push_back(p);
      }
   }

   // Determine the number of dependencies for each patch.
   std::vector<std::size_t> sizes(grid.size(), 0);
   for (std::size_t i = 0; i != pairs.size(); ++i) {
      ++sizes[pairs[i].first];
   }

   // Allocate memory for the dependencies array.
   dependencies->rebuild(sizes.begin(), sizes.end());

   // Record the dependencies.
   std::fill(sizes.begin(), sizes.end(), 0);
   for (std::size_t i = 0; i != pairs.size(); ++i) {
      const std::size_t j = pairs[i].first;
      (*dependencies)(j, sizes[j]++) = pairs[i].second;
   }
}


template<typename _T, std::size_t _D, std::size_t _PatchExtent, typename _R,
         typename _OutputIterator>
inline
void
getIntersectingPatches(const Grid<_T, _D, _PatchExtent, _R>& grid,
                       geom::BBox<_R, _D> box, _OutputIterator indices) {
   typedef std::tr1::array<std::size_t, _D> IndexList;
   typedef array::SimpleMultiIndexRangeIterator<_D> Iterator;

   // Convert the box to index coordinates.
   const _R inverseGridSpacing = 1. / grid.spacing;
   box.lower -= grid.lowerCorner;
   box.lower *= inverseGridSpacing;
   box.upper -= grid.lowerCorner;
   box.upper *= inverseGridSpacing;

   // Determine the index range.
   // Closed lower and open upper bounds.
   IndexList lower, upper;
   for (std::size_t i = 0; i != _D; ++i) {
      // Lower.
      if (box.lower[i] < 0) {
         lower[i] = 0;
      }
      else if (box.lower[i] >= grid.extents()[i]) {
         lower[i] = grid.extents()[i];
      }
      else {
         lower[i] = std::size_t(box.lower[i]);
      }
      // Upper.
      if (box.upper[i] < 0) {
         upper[i] = 0;
      }
      else if (box.upper[i] >= grid.extents()[i]) {
         upper[i] = grid.extents()[i];
      }
      else {
         upper[i] = std::size_t(box.upper[i]) + 1;
      }
   }
   // Construct from the extents and bases.
   const array::SimpleMultiIndexRange<_D> range = {upper - lower, lower};
   // Loop over the index range and report the array indices.
   const Iterator end = Iterator::end(range);
   for (Iterator i = Iterator::begin(range); i != end; ++i) {
      *indices++ = grid.arrayIndex(*i);
   }
}


template<typename _T, std::size_t _N, typename _R>
inline
void
writeVtkXml(const Grid<_T, 3, _N, _R>& grid, std::ostream& out) {
   const std::size_t Dimension = 3;
   typedef array::SimpleMultiIndexRangeIterator<Dimension> Iterator;
   typedef typename Grid<_T, Dimension, _N, _R>::VertexPatch VertexPatch;

   out << "<?xml version=\"1.0\"?>\n"
       << "<VTKFile type=\"ImageData\">\n";

   out << "<ImageData WholeExtent=\""
       << "0 " << grid.extents()[0] * _N - 1 << ' '
       << "0 " << grid.extents()[1] * _N - 1 << ' '
       << "0 " << grid.extents()[2] * _N - 1
       << "\" Origin=\"" << grid.lowerCorner
       << "\" Spacing=\""
       << grid.spacing << ' '
       << grid.spacing << ' '
       << grid.spacing
       << "\">\n";

   // For each refined patch.
   const Iterator end = Iterator::end(grid.extents());
   for (Iterator i = Iterator::begin(grid.extents()); i != end; ++i) {
      out << "<Piece Extent=\""
          << (*i)[0] * _N << ' '
          << ((*i)[0] + 1) * _N - 1 << ' '
          << (*i)[1] * _N << ' '
          << ((*i)[1] + 1) * _N - 1 << ' '
          << (*i)[2] * _N << ' '
          << ((*i)[2] + 1) * _N - 1
          << "\">\n";

      out << "<PointData Scalars=\"Distance\">\n"
          << "<DataArray type=\"Float64\" Name=\"Distance\" "
          << "NumberOfComponents=\"1\" format=\"ascii\">\n";
      // Loop over the vertices in the patch.
      const VertexPatch& patch = grid(*i);
      if (patch.isRefined()) {
         std::copy(patch.begin(), patch.end(),
                   std::ostream_iterator<_T>(out, " "));
      }
      else {
         for (std::size_t i = 0; i != patch.size(); ++i) {
            out << patch.fillValue << ' ';
         }
      }
      out << "</DataArray>\n";
      out << "</PointData>\n";

      out << "<CellData></CellData>\n";
      out << "</Piece>\n";
   }
   out << "</ImageData>\n";
   out << "</VTKFile>\n";
}


template<typename _T, std::size_t _N, typename _R>
inline
void
writeVtkXml(const Grid<_T, 2, _N, _R>& grid, std::ostream& out) {
   typedef array::SimpleMultiIndexRangeIterator<2> Iterator;
   typedef Grid<_T, 2, _N, _R> Grid;
   typedef typename Grid::VertexPatch Patch;

   // Construct a uniform grid.
   GridUniform<_T, 2> uniform(grid.extents() * _N, grid.lowerCorner,
                              grid.spacing);

   //
   // Copy the elements.
   //
   // For each patch.
   const Iterator iEnd = Iterator::end(grid.extents());
   for (Iterator i = Iterator::begin(grid.extents()); i != iEnd; ++i) {
      const Patch& patch = grid(*i);
      // For each vertex.
      if (patch.isRefined()) {
         const Iterator jEnd = Iterator::end(patch.extents());
         for (Iterator j = Iterator::begin(patch.extents()); j != jEnd;
              ++j) {
            uniform(*i * _N + *j) = patch(*j);
         }
      }
      else {
         const Iterator jEnd = Iterator::end(patch.extents());
         for (Iterator j = Iterator::begin(patch.extents()); j != jEnd;
              ++j) {
            uniform(*i * _N + *j) = patch.fillValue;
         }
      }
   }

   // Write the uniform grid.
   writeVtkXml(uniform, out);
}


// CONTINUE Writing in multiple pieces does not work. The result doesn't 
// display correctly in Paraview, even when all the pieces are output.
template<typename _T, std::size_t _N, typename _R>
inline
void
writeVtkXmlOld(const Grid<_T, 2, _N, _R>& grid, std::ostream& out) {
   const std::size_t Dimension = 2;
   typedef array::SimpleMultiIndexRangeIterator<Dimension> Iterator;
   typedef typename Grid<_T, Dimension, _N, _R>::VertexPatch VertexPatch;

   out << "<?xml version=\"1.0\"?>\n"
       << "<VTKFile type=\"ImageData\">\n";

   out << "<ImageData WholeExtent=\""
       << "0 " << grid.extents()[0] * _N - 1 << ' '
       << "0 " << grid.extents()[1] * _N - 1 << ' '
       << "0 0"
       << "\" Origin=\"" << grid.lowerCorner << " 0"
       << "\" Spacing=\""
       << grid.spacing << ' '
       << grid.spacing << ' '
       << grid.spacing
       << "\">\n";

   // For each refined patch.
   const Iterator end = Iterator::end(grid.extents());
   for (Iterator i = Iterator::begin(grid.extents()); i != end; ++i) {
      if (! grid(*i).isRefined()) {
         continue;
      }
      out << "<Piece Extent=\""
          << (*i)[0] * _N << ' '
          << ((*i)[0] + 1) * _N - 1 << ' '
          << (*i)[1] * _N << ' '
          << ((*i)[1] + 1) * _N - 1 << ' '
          << " 0 0\">\n";

      out << "<PointData Scalars=\"Distance\">\n"
          << "<DataArray type=\"Float64\" Name=\"Distance\" "
          << "NumberOfComponents=\"1\" format=\"ascii\">\n";
      // Loop over the vertices in the patch.
      const VertexPatch& patch = grid(*i);
      std::copy(patch.begin(), patch.end(),
                std::ostream_iterator<_T>(out, " "));
      out << "</DataArray>\n";
      out << "</PointData>\n";

      out << "<CellData></CellData>\n";
      out << "</Piece>\n";
   }
   out << "</ImageData>\n";
   out << "</VTKFile>\n";
}


template<typename _T, std::size_t _N, typename _R>
inline
void
writeVtkXml(const Grid<_T, 2, _N, _R>& grid,
            const std::tr1::array<std::size_t, 2>& patchIndices,
            std::ostream& out) {
   const std::size_t Dimension = 2;
   typedef typename Grid<_T, Dimension, _N, _R>::VertexPatch VertexPatch;

   out << "<?xml version=\"1.0\"?>\n"
       << "<VTKFile type=\"ImageData\">\n";

   out << "<ImageData WholeExtent=\""
       << "0 " << grid.extents()[0] * _N - 1 << ' '
       << "0 " << grid.extents()[1] * _N - 1 << ' '
       << "0 0"
       << "\" Origin=\"" << grid.lowerCorner << " 0"
       << "\" Spacing=\""
       << grid.spacing << ' '
       << grid.spacing << ' '
       << grid.spacing
       << "\">\n";

   out << "<Piece Extent=\""
       << patchIndices[0] * _N << ' '
       << (patchIndices[0] + 1) * _N - 1 << ' '
       << patchIndices[1] * _N << ' '
       << (patchIndices[1] + 1) * _N - 1 << ' '
       << " 0 0\">\n";

   out << "<PointData Scalars=\"Distance\">\n"
       << "<DataArray type=\"Float64\" Name=\"Distance\" "
       << "NumberOfComponents=\"1\" format=\"ascii\">\n";
   // Loop over the vertices in the patch.
   const VertexPatch& patch = grid(patchIndices);
   if (patch.isRefined()) {
      std::copy(patch.begin(), patch.end(),
                std::ostream_iterator<_T>(out, " "));
   }
   else {
      for (std::size_t i = 0; i != patch.size(); ++i) {
         out << patch.fillValue << ' ';
      }
   }
   out << "</DataArray>\n";
   out << "</PointData>\n";

   out << "<CellData></CellData>\n";
   out << "</Piece>\n";

   out << "</ImageData>\n";
   out << "</VTKFile>\n";
}


// CONTINUE Writing in multiple pieces does not work. The result doesn't 
// display correctly in Paraview, even when all the pieces are output.
template<typename _T, std::size_t _N, typename _R>
inline
void
writeParaview(const Grid<_T, 2, _N, _R>& grid, const std::string& name) {
   const std::size_t Dimension = 2;
   typedef array::SimpleMultiIndexRangeIterator<Dimension> Iterator;
   typedef typename Grid<_T, Dimension, _N, _R>::VertexPatch VertexPatch;

   // Open the ParaView file.
   std::string paraviewName = name;
   paraviewName += ".pvd";
   std::ofstream paraviewFile(paraviewName.c_str());
   paraviewFile << "<?xml version=\"1.0\"?>\n"
                << "<VTKFile type=\"Collection\">\n"
                << "<Collection>\n";
   // For each refined patch.
   const Iterator end = Iterator::end(grid.extents());
   for (Iterator i = Iterator::begin(grid.extents()); i != end; ++i) {
      std::string vtkName = name;
      std::ostringstream oss;
      oss << std::setfill('0') << std::setw(6) << grid.arrayIndex(*i);
      vtkName += oss.str();
      vtkName += ".vti";

      paraviewFile << "<DataSet part=\"1\" file=\"" << vtkName << "\"/>\n";

      std::ofstream vtkFile(vtkName.c_str());
      writeVtkXml(grid, *i, vtkFile);
   }
   paraviewFile << "</Collection>\n";
   paraviewFile << "</VTKFile>\n";
}


template<typename _T, std::size_t _D, std::size_t _N, typename _R>
inline
void
printInfo(const Grid<_T, _D, _N, _R>& grid, std::ostream& out) {
   const _T Inf = std::numeric_limits<_T>::infinity();

   std::size_t size = 0;
   std::size_t nonNegative = 0;
   std::size_t negative = 0;
   std::size_t unknown = 0;
   std::size_t positiveFar = 0;
   std::size_t negativeFar = 0;
   _T x;
   // For each refined patch.
   for (std::size_t i = 0; i != grid.size(); ++i) {
      if (! grid[i].isRefined()) {
         continue;
      }
      // For each grid point in the patch.
      for (std::size_t j = 0; j != grid[i].size(); ++j) {
         x = grid[i][j];
         ++size;
         if (x >= 0) {
            ++nonNegative;
         }
         else if (x < 0) {
            ++negative;
         }
         else {
            ++unknown;
         }
         if (x == Inf) {
            ++positiveFar;
         }
         else if (x == -Inf) {
            ++negativeFar;
         }
      }
   }
   out << "Domain = " << grid.domain() << '\n'
       << "Extents = " << grid.extents() << '\n'
       << "Patch extent = " << _N << '\n'
       << "Number of patches = " << grid.size() << '\n'
       << "Number of refined patches = " << grid.numRefined() << '\n'
       << "Number of grid points = " << size << '\n'
       << "known/unknown = " << size - unknown << " / " << unknown << '\n'
       << "non-negative/negative = " << nonNegative << " / " << negative
       << '\n'
       << "positive far/negative far = " << positiveFar << " / " << negativeFar
       << '\n';
}


} // namespace levelSet
