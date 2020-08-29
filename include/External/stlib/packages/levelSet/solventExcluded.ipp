// -*- C++ -*-

#if !defined(__levelSet_solventExcluded_ipp__)
#error This file is an implementation detail of solventExcluded.
#endif

namespace levelSet {


template<typename _T, std::size_t _D>
inline
void
solventExcluded(array::SimpleMultiArrayRef<_T, _D>* grid,
                const geom::BBox<_T, _D>& domain,
                const std::vector<geom::Ball<_T, _D> >& balls,
                const _T probeRadius) {
   // Compute the distance to the offset balls.
   {
      std::vector<geom::Ball<_T, _D> > offset = balls;
      for (std::size_t i = 0; i != offset.size(); ++i) {
         offset[i].radius += probeRadius;
      }
      negativeDistance(grid, domain, offset);
   }
   // Offset by the probe radius.
   *grid += probeRadius;
}


template<typename _T, std::size_t _D>
inline
void
solventExcluded(GridUniform<_T, _D>* grid,
                const std::vector<geom::Ball<_T, _D> >& balls,
                const _T probeRadius) {
   solventExcluded(grid, grid->domain(), balls, probeRadius);
}


// Determine the patch/ball dependencies for the solvent-excluded surface.
template<typename _T, std::size_t _D, std::size_t _N>
inline
void
solventExcludedDependencies
(const Grid<_T, _D, _N>& grid,
 const std::vector<geom::Ball<_T, _D> >& balls,
 const _T probeRadius,
 array::StaticArrayOfArrays<std::size_t>* dependencies) {
   // Calculate the largest radius of the solvent-expanded balls.
   _T maxRadius = 0;
   for (std::size_t i = 0; i != balls.size(); ++i) {
      if (balls[i].radius > maxRadius) {
         maxRadius = balls[i].radius;
      }
   }
   maxRadius += probeRadius;
   // First expand the balls by the probe radius, then consider
   // their influence. A ball may influence points up to maxRadius
   // beyond its surface. Consider the case that one ball barely
   // intersects the largest ball. Then points all the way up to
   // the center of the largest ball may be closest to the
   // intersection of the two balls.
   const _T offset = probeRadius + maxRadius;
   std::vector<geom::Ball<_T, _D> > offsetBalls = balls;
   for (std::size_t i = 0; i != offsetBalls.size(); ++i) {
      offsetBalls[i].radius += offset;
   }
   patchDependencies(grid, offsetBalls.begin(), offsetBalls.end(),
                     dependencies);
}


template<typename _T, std::size_t _D, std::size_t _PatchExtent>
inline
void
solventExcluded(Grid<_T, _D, _PatchExtent>* grid,
                const std::vector<geom::Ball<_T, _D> >& balls,
                const _T probeRadius) {
   // Compute the distance to the offset balls.
   {
      std::vector<geom::Ball<_T, _D> > offset = balls;
      for (std::size_t i = 0; i != offset.size(); ++i) {
         offset[i].radius += probeRadius;
      }
      negativeDistance(grid, offset);
   }
   // Offset by the probe radius.
   *grid += probeRadius;
}


template<typename _T, std::size_t _D>
inline
void
solventAccessibleCavities(array::SimpleMultiArrayRef<_T, _D>* grid,
                          const geom::BBox<_T, _D>& domain,
                          const std::vector<geom::Ball<_T, _D> >& balls,
                          const _T probeRadius) {
   typedef array::SimpleMultiArrayRef<_T, _D> SimpleMultiArray;
   typedef typename SimpleMultiArray::IndexList IndexList;
   typedef typename SimpleMultiArray::Range Range;
   typedef array::SimpleMultiIndexRangeIterator<_D> Iterator;

   const _T Inf = std::numeric_limits<_T>::infinity();
   assert(grid->extents()[0] > 1);
   const _T dx = (domain.upper[0] - domain.lower[0]) / (grid->extents()[0] - 1);
   // 1.1 * (diagonal length of a voxel)
   const _T buffer = 1.1 * std::sqrt(double(_D)) * dx;

   // Compute the solvent-center excluded domain.
   {
      // Compute the distance a little bit past the zero iso-surface in order
      // to get an accurate initial condition for the static H-J problem.
      std::vector<geom::Ball<_T, _D> > offset = balls;
      for (std::size_t i = 0; i != offset.size(); ++i) {
         offset[i].radius += probeRadius;
      }
      std::fill(grid->begin(), grid->end(),
                std::numeric_limits<_T>::quiet_NaN());
      positiveDistance(grid, domain, offset, _T(0), buffer);
   }

   // Reverse the sign to get the solvent-center accessible domain.
   for (typename SimpleMultiArray::iterator i = grid->begin(); i != grid->end();
        ++i) {
      *i = -*i;
   }

   // Mark the outside component of the solvent-center-accesible domain as 
   // positive far-away.
   const IndexList closedUpper = grid->extents() - 1;
   std::deque<IndexList> queue;
   IndexList lower, upper;
   // Start with the lower corner.
   IndexList i = ext::filled_array<IndexList>(0);
   assert((*grid)(i) != (*grid)(i) || (*grid)(i) <= 0);
   (*grid)(i) = Inf;
   queue.push_back(i);
   // Set the rest of the points outside to positive far-away.
   while (! queue.empty()) {
      i = queue.front();
      queue.pop_front();
      // Examine all neighboring grid points. There are 3^_D - 1 neighbors.
      // Note that it is not sufficient to just examine the 2*_D adjacent
      // neighbors. Crevices on the surface of the solvent-center-accessible
      // surface can isolate grid points.
      for (std::size_t n = 0; n != _D; ++n) {
         if (i[n] == 0) {
            lower[n] = 0;
         }
         else {
            lower[n] = i[n] - 1;
         }
         upper[n] = std::min(i[n] + 2, grid->extents()[n]);
      }
      const Range range = {upper - lower, lower};
      const Iterator end = Iterator::end(range);
      for (Iterator j = Iterator::begin(range); j != end; ++j) {
         if (*j == i) {
            continue;
         }
         if ((*grid)(*j) != (*grid)(*j) || (*grid)(*j) <= 0) {
            queue.push_back(*j);
            (*grid)(*j) = Inf;
         }
      }
   }

   // The NaN points are inside the solvent-center-accessible cavities.
   // Erase the points with positive distance. The distance from the
   // SCA will be computed with the H-J solver.
   for (typename SimpleMultiArray::iterator i = grid->begin(); i != grid->end();
        ++i) {
      if (*i > 0) {
         // Note that the H-J solver uses max() to denote an unknown point.
         *i = std::numeric_limits<_T>::max();
      }
      else if (*i != *i) {
         *i = -Inf;
      }
   }

   // If there are no non-positive points, there are no solvent-accessible 
   // cavities.
   if (! hasNonPositive(grid->begin(), grid->end())) {
      // Fill the grid with positive far-away values and return.
      std::fill(grid->begin(), grid->end(), Inf);
      return;
   }

   // CONTINUE REMOVE
#if 1
   std::cerr << "Before H-J solver.\n";
   printLevelSetInfo(grid->begin(), grid->end(), std::cerr);
#endif
   // Compute the distance from the SCA interface.
   array::MultiArrayRef<_T, _D> gridRef(grid->data(), grid->extents());
   hj::computeUnsignedDistance(gridRef, dx, probeRadius + buffer);

   // Offset to get the solvent accessible surface.
   *grid -= probeRadius;
   // CONTINUE REMOVE
#if 1
   std::cerr << "After H-J solver.\n";
   printLevelSetInfo(grid->begin(), grid->end(), std::cerr);
#endif
}


template<typename _T, std::size_t _D>
inline
void
solventAccessibleCavitiesOld(array::SimpleMultiArrayRef<_T, _D>* grid,
                             const geom::BBox<_T, _D>& domain,
                             const std::vector<geom::Ball<_T, _D> >& balls,
                             const _T probeRadius) {
   typedef array::SimpleMultiArrayRef<_T, _D> SimpleMultiArray;
   typedef typename SimpleMultiArray::IndexList IndexList;
   typedef typename SimpleMultiArray::Range Range;
   typedef array::SimpleMultiIndexRangeIterator<_D> Iterator;

   const _T Inf = std::numeric_limits<_T>::infinity();

   // Compute the solvent-center excluded domain.
   {
      std::vector<geom::Ball<_T, _D> > offset = balls;
      for (std::size_t i = 0; i != offset.size(); ++i) {
         offset[i].radius += probeRadius;
      }
      negativeDistance(grid, domain, offset);
   }

#if 0
   std::cerr << "Solvent-center excluded domain\n";
   printLevelSetInfo(*grid, std::cerr);
#endif

   // Mark the outside component of the solvent-center-accesible domain as 
   // positive far-away.
   const IndexList closedUpper = grid->extents() - 1;
   std::deque<IndexList> queue;
   IndexList lower, upper;
   // Start with the lower corner.
   IndexList i = ext::filled_array<IndexList>(0);
   assert((*grid)(i) != (*grid)(i));
   (*grid)(i) = Inf;
   queue.push_back(i);
   // Set the rest of the points outside to positive far-away.
   while (! queue.empty()) {
      i = queue.front();
      queue.pop_front();
      // Examine all neighboring grid points. There are 3^_D - 1 neighbors.
      // Note that it is not sufficient to just examine the 2*_D adjacent
      // neighbors. Crevices on the surface of the solvent-center-accessible
      // surface can isolate grid points.
      for (std::size_t n = 0; n != _D; ++n) {
         if (i[n] == 0) {
            lower[n] = 0;
         }
         else {
            lower[n] = i[n] - 1;
         }
         upper[n] = std::min(i[n] + 2, grid->extents()[n]);
      }
      const Range range = {upper - lower, lower};
      const Iterator end = Iterator::end(range);
      for (Iterator j = Iterator::begin(range); j != end; ++j) {
         if (*j == i) {
            continue;
         }
         if ((*grid)(*j) != (*grid)(*j)) {
            queue.push_back(*j);
            (*grid)(*j) = Inf;
         }
      }
   }

#if 0
   std::cerr << "\nOutside far away:\n";
   printLevelSetInfo(*grid, std::cerr);
#endif   

   // If there are no unknown points, there are no solvent-accessible cavities.
   if (! hasUnknown(grid->begin(), grid->end())) {
      // Fill the grid with positive far-away values and return.
      std::fill(grid->begin(), grid->end(), Inf);
      return;
   }

   // Get the boundary of the solvent-center-accessible cavities.
   Iterator end = Iterator::end(grid->extents());
   for (Iterator i = Iterator::begin(grid->extents()); i != end; ++i) {
      if ((*grid)(*i) < 0 && hasUnknownAdjacentNeighbor(*grid, *i)) {
         (*grid)(*i) = -(*grid)(*i);
         queue.push_back(*i);
      }
   }
   
   // Walk downhill and reverse the sign of the distance.
   while (! queue.empty()) {
      i = queue.front();
      queue.pop_front();
      const _T x = (*grid)(i);
      // Examine each adjacent grid point.
      for (std::size_t n = 0; n != _D; ++n) {
         if (i[n] != closedUpper[n]) {
            ++i[n];
            if ((*grid)(i) < 0 && -(*grid)(i) > x) {
               (*grid)(i) = -(*grid)(i);
               queue.push_back(i);
            }
            --i[n];
         }
         if (i[n] != 0) {
            --i[n];
            if ((*grid)(i) < 0 && -(*grid)(i) > x) {
               (*grid)(i) = -(*grid)(i);
               queue.push_back(i);
            }
            ++i[n];
         }
      }
   }

   // Set the unknown distances to negative far-away. These are the 
   // solvent-center-accessible portion of the cavities.
   // Set all of the negative distances to positive far-away. These are not
   // part of the solvent-accessible cavities.
   for (std::size_t i = 0; i != grid->size(); ++i) {
      _T& x = (*grid)[i];
      if (x != x) {
         x = -Inf;
      }
      else if (x < 0) {
         x = Inf;
      }
   }

   // Offset to get the solvent accessible surface.
   *grid -= probeRadius;

#if 0
   std::cerr << "\nSolvent-accessible cavities:\n";
   printLevelSetInfo(*grid, std::cerr);
#endif
}


// This is an old method.
template<typename _T, std::size_t _D, std::size_t _N>
inline
void
solventAccessibleCavitiesWithSeedInterpolation
(Grid<_T, _D, _N>* grid,
 const std::vector<geom::Ball<_T, _D> >& balls,
 const _T probeRadius) {
   typedef Grid<_T, _D, _N> Grid;
   typedef typename Grid::IndexList IndexList;
   typedef typename Grid::DualIndices DualIndices;
   typedef typename Grid::VertexPatch VertexPatch;
   typedef typename Grid::Range Range;
   typedef array::SimpleMultiIndexRangeIterator<_D> Iterator;

   const _T Inf = std::numeric_limits<_T>::infinity();

   // For this method to work, the probe radius must be at least as large
   // as the diagonal length of a patch. Consider the case that a patch
   // has a single seed in its lower corner. The seed sphere must cover 
   // all of the patch grid points in order to ensure that they are 
   // correctly marked as being in the solvent-accessible cavity.
   assert(probeRadius > grid->spacing * _N * std::sqrt(3));

   // Compute the solvent-center excluded domain.
   {
      std::vector<geom::Ball<_T, _D> > offset(balls);
      for (std::size_t i = 0; i != offset.size(); ++i) {
         offset[i].radius += probeRadius;
      }
      negativePowerDistance(grid, offset);
   }

   // Mark the outside component of the solvent-center-accesible domain as 
   // negative far-away.
   std::deque<DualIndices> queue;
   std::vector<DualIndices> neighbors;
   std::back_insert_iterator<std::vector<DualIndices> >
      insertIterator(neighbors);
   // Start with the lower corner.
   DualIndices di;
   di.first = di.second = ext::filled_array<IndexList>(0);
   assert((*grid)(di) != (*grid)(di) || (*grid)(di) > 0);
   (*grid)(di) = -Inf;
   queue.push_back(di);
   // Set the rest of the points outside to negative far-away.
   while (! queue.empty()) {
      di = queue.front();
      queue.pop_front();
      // Examine the set of all neighboring grid points.
      neighbors.clear();
      grid->allNeighbors(di, insertIterator);
      for (std::size_t i = 0; i != neighbors.size(); ++i) {
         const DualIndices& n = neighbors[i];
         if ((*grid)(n) != (*grid)(n) || (*grid)(n) > 0) {
            queue.push_back(n);
            (*grid)(n) = -Inf;
         }
      }
   }

   // Record the points at the boundary of the solvent-center-accessible
   // cavities. These are positive points with negative neighbors, and will
   // be used as seed points to compute the solvent-accessible cavities.
   std::vector<geom::Ball<_T, _D> > seeds;
   {
      geom::Ball<_T, _D> seed;
      seed.radius = probeRadius;
      std::vector<DualIndices> neighbors;
      std::back_insert_iterator<std::vector<DualIndices> >
         insertIterator(neighbors);
      const Iterator pEnd = Iterator::end(grid->extents());
      for (Iterator p = Iterator::begin(grid->extents()); p != pEnd; ++p) {
         VertexPatch& patch = (*grid)(*p);
         if (! patch.isRefined()) {
            continue;
         }
         // Loop over the grid points in the patch.
         const Iterator iEnd = Iterator::end(patch.extents());
         for (Iterator i = Iterator::begin(patch.extents()); i != iEnd; ++i) {
            if (patch(*i) > 0) {
               const _T a = patch(*i);
               // For each adjacent neighbor.
               for (std::size_t direction = 0; direction != 2 * _D;
                    ++direction) {
                  neighbors.clear();
                  grid->adjacentNeighbors(*p, *i, direction, insertIterator);
                  if (neighbors.empty()) {
                     continue;
                  }
                  const _T b = (*grid)(neighbors[0]);
                  if (b <= 0) {
                     seed.center = grid->indexToLocation(*p, *i);
                     // Use linear interpolation to place the seed between
                     // this grid point and its adjacent neighbor.
                     // a + t (b - a) == 0
                     // t = a / (a - b)
                     seed.center[direction / 2] += (2 * (direction % 2) - 1) *
                        a / (a - b) * grid->spacing;
                     seeds.push_back(seed);
                  }
               }
            }
         }
      }
   }

   // If there are no seed points, there are no solvent-accessible cavities.
   if (seeds.empty()) {
      // Fill the grid with positive far-away values and return.
      grid->clear();
      for (std::size_t i = 0; i != grid->size(); ++i) {
         (*grid)[i].fillValue = Inf;
      }
      return;
   }

   // Compute the solvent-accessible cavities from the seed points.
   negativePowerDistance(grid, seeds);
}


template<typename _T, std::size_t _D, std::size_t _N>
inline
void
markOutsideAsNegativeInf(Grid<_T, _D, _N>* grid) {
   typedef Grid<_T, _D, _N> Grid;
   typedef typename Grid::IndexList IndexList;
   typedef typename Grid::DualIndices DualIndices;

   const _T Inf = std::numeric_limits<_T>::infinity();

   std::deque<DualIndices> queue;
   std::vector<DualIndices> neighbors;
   std::back_insert_iterator<std::vector<DualIndices> >
      insertIterator(neighbors);
   // Start with the lower corner.
   DualIndices di;
   di.first = di.second = ext::filled_array<IndexList>(0);
   assert((*grid)(di) != (*grid)(di) || (*grid)(di) > 0);
   (*grid)(di) = -Inf;
   queue.push_back(di);
   // Set the rest of the points outside to negative far-away.
   while (! queue.empty()) {
      di = queue.front();
      queue.pop_front();
      // Examine the set of all neighboring grid points.
      neighbors.clear();
      grid->allNeighbors(di, insertIterator);
      for (std::size_t i = 0; i != neighbors.size(); ++i) {
         const DualIndices& n = neighbors[i];
         if ((*grid)(n) != (*grid)(n) || (*grid)(n) > 0) {
            queue.push_back(n);
            (*grid)(n) = -Inf;
         }
      }
   }
}


template<typename _T, std::size_t _D, std::size_t _N>
inline
void
solventAccessibleCavitySeeds(const Grid<_T, _D, _N>& grid, const _T probeRadius,
                             std::vector<geom::Ball<_T, _D> >* seeds) {
   typedef Grid<_T, _D, _N> Grid;
   typedef typename Grid::DualIndices DualIndices;
   typedef typename Grid::VertexPatch VertexPatch;
   typedef array::SimpleMultiIndexRangeIterator<_D> Iterator;

   const _T PatchDiagonalHalfLength = 0.5 * grid.spacing * _N * std::sqrt(_D);
   // The threshhold for accepting a point as a seed.
   const _T SeedThreshhold = 1.01 * grid.spacing;

   // For this method to work, the probe radius must be at least as large
   // as the diagonal length of a patch. Consider the case that a patch
   // has a single seed in its lower corner. The seed sphere must cover 
   // all of the patch grid points in order to ensure that they are 
   // correctly marked as being in the solvent-accessible cavity.
   assert(probeRadius > grid.spacing * _N * std::sqrt(3));

   geom::Ball<_T, _D> seed;
   std::vector<DualIndices> neighbors;
   std::back_insert_iterator<std::vector<DualIndices> >
      insertIterator(neighbors);
   DualIndices di;
   const Iterator pEnd = Iterator::end(grid.extents());
   for (Iterator p = Iterator::begin(grid.extents()); p != pEnd; ++p) {
      const VertexPatch& patch = grid(*p);
      di.first = *p;
      if (patch.isRefined()) {
         // Loop over the grid points in the patch.
         const Iterator iEnd = Iterator::end(patch.extents());
         for (Iterator i = Iterator::begin(patch.extents()); i != iEnd; ++i) {
            const _T a = patch(*i);
            // A seed must be near the interface.
            if (a > 0 && a < SeedThreshhold) {
               // A seed must have a neighbor with positive distance.
               // This requirement excludes cavities that result from the 
               // solvent being pinched between two atoms.
               // Examine the set of all neighboring grid points.
               neighbors.clear();
               di.second = *i;
               grid.allNeighbors(di, insertIterator);
               for (std::size_t j = 0; j != neighbors.size(); ++j) {
                  const DualIndices& n = neighbors[j];
                  if (grid(n) > 0) {
                     seed.radius = probeRadius + a;
                     seed.center = grid.indexToLocation(*p, *i);
                     seeds->push_back(seed);
                     break;
                  }
               }
            }
         }
      }
      // Unrefined patches.
      else {
         // Check for unrefined patches in the middle of cavities.
         if (patch.fillValue != patch.fillValue) {
            // Add a seed that covers the patch.
            seed.center = grid.getPatchLowerCorner(*p);
            seed.center += PatchDiagonalHalfLength;
            seed.radius = 1.01 * PatchDiagonalHalfLength;
            seeds->push_back(seed);
         }
      }
   }
}


template<typename _T, std::size_t _D, std::size_t _N>
inline
void
solventAccessibleCavities(Grid<_T, _D, _N>* grid,
                          const std::vector<geom::Ball<_T, _D> >& balls,
                          const _T probeRadius) {
   const _T Inf = std::numeric_limits<_T>::infinity();

   // Compute the solvent-center excluded domain.
   positiveDistanceOutside(grid, balls, probeRadius, 2 * grid->spacing);

   // Mark the outside component of the solvent-center-accesible domain as 
   // negative far-away.
   markOutsideAsNegativeInf(grid);

   // Record the points at the boundary of the solvent-center-accessible
   // cavities. These are positive points with negative neighbors, and will
   // be used as seed points to compute the solvent-accessible cavities.
   std::vector<geom::Ball<_T, _D> > seeds;
   solventAccessibleCavitySeeds(*grid, probeRadius, &seeds);

   // If there are no seed points, there are no solvent-accessible cavities.
   if (seeds.empty()) {
      // Fill the grid with positive far-away values and return.
      grid->clear();
      for (std::size_t i = 0; i != grid->size(); ++i) {
         (*grid)[i].fillValue = Inf;
      }
      return;
   }

   // Compute the solvent-accessible cavities from the seed points.
   negativePowerDistance(grid, seeds);
}


template<typename _T, std::size_t _D, std::size_t _N>
inline
void
solventAccessibleCavitiesOld(Grid<_T, _D, _N>* grid,
                          const std::vector<geom::Ball<_T, _D> >& balls,
                          const _T probeRadius) {
   typedef Grid<_T, _D, _N> Grid;
   typedef typename Grid::IndexList IndexList;
   typedef typename Grid::DualIndices DualIndices;
   typedef typename Grid::VertexPatch VertexPatch;
   typedef array::SimpleMultiIndexRangeIterator<_D> Iterator;

   const _T Inf = std::numeric_limits<_T>::infinity();

   // Compute the solvent-center excluded domain.
   {
      std::vector<geom::Ball<_T, _D> > offset = balls;
      for (std::size_t i = 0; i != offset.size(); ++i) {
         offset[i].radius += probeRadius;
      }
      negativeDistance(grid, offset);
   }

   // Mark the outside component of the solvent-center-accesible domain as 
   // positive far-away.
   std::deque<DualIndices> queue;
   std::vector<DualIndices> neighbors;
   std::back_insert_iterator<std::vector<DualIndices> >
      insertIterator(neighbors);
   // Start with the lower corner.
   DualIndices di;
   di.first = di.second = ext::filled_array<IndexList>(0);
   assert((*grid)(di) != (*grid)(di));
   (*grid)(di) = Inf;
   queue.push_back(di);
   // Set the rest of the points outside to positive far-away.
   while (! queue.empty()) {
      di = queue.front();
      queue.pop_front();
      // Examine the set of all neighboring grid points.
      neighbors.clear();
      grid->allNeighbors(di, insertIterator);
      for (std::size_t i = 0; i != neighbors.size(); ++i) {
         const DualIndices& n = neighbors[i];
         if ((*grid)(n) != (*grid)(n)) {
            queue.push_back(n);
            (*grid)(n) = Inf;
         }
      }
   }

#if 1
   // CONTINUE: REMOVE
   printLevelSetInfo(*grid, std::cout);
#endif
#if 0
   {
      // CONTINUE: REMOVE
      // Print the unknown points.
      const Iterator pEnd = Iterator::end(grid->extents());
      for (Iterator p = Iterator::begin(grid->extents()); p != pEnd; ++p) {
         VertexPatch& patch = (*grid)(*p);
         if (! patch.isRefined()) {
            continue;
         }
         // Loop over the grid points in the patch.
         const Iterator iEnd = Iterator::end(patch.extents());
         for (Iterator i = Iterator::begin(patch.extents()); i != iEnd; ++i) {
            if (patch(*i) != patch(*i)) {
               std::cout << "Unknown: " << *p << ", " << *i << '\n';
               // Print the neighboring values.
               neighbors.clear();
               grid->allNeighbors(std::make_pair(*p, *i), insertIterator);
               for (std::size_t j = 0; j != neighbors.size(); ++j) {
                  std::cout << (*grid)(neighbors[j]) << ' ';
               }
               std::cout << '\n';
            }
         }
      }
   }
#endif

   // If there are no unknown points, there are no solvent-accessible cavities.
   if (! hasUnknown(*grid)) {
      // Fill the grid with positive far-away values and return.
      grid->clear();
      for (std::size_t i = 0; i != grid->size(); ++i) {
         (*grid)[i].fillValue = Inf;
      }
      return;
   }

   // Get the boundary of the solvent-center-accessible cavities.
   // Loop over the refined patches.
   const Iterator pEnd = Iterator::end(grid->extents());
   for (Iterator p = Iterator::begin(grid->extents()); p != pEnd; ++p) {
      VertexPatch& patch = (*grid)(*p);
      if (! patch.isRefined()) {
         continue;
      }
      // Loop over the grid points in the patch.
      const Iterator iEnd = Iterator::end(patch.extents());
      for (Iterator i = Iterator::begin(patch.extents()); i != iEnd; ++i) {
         di.first = *p;
         di.second = *i;
         if (patch(*i) < 0 && hasUnknownAdjacentNeighbor(*grid, di)) {
            patch(*i) = - patch(*i);
            queue.push_back(di);
         }
      }
   }

   // Walk downhill and reverse the sign of the distance.
   while (! queue.empty()) {
      di = queue.front();
      queue.pop_front();
      const _T x = (*grid)(di);
      // Examine the set of adjacent grid points.
      neighbors.clear();
      grid->adjacentNeighbors(di, insertIterator);
      for (std::size_t i = 0; i != neighbors.size(); ++i) {
         _T& y = (*grid)(neighbors[i]);
         if (y < 0 && -y > x) {
            y = - y;
            queue.push_back(neighbors[i]);
         }
      }
   }

   // Set the unknown distances to negative far-away. These are the 
   // solvent-center-accessible portion of the cavities.
   // Set all of the negative distances to positive far-away. These are not
   // part of the solvent-accessible cavities.
   for (std::size_t i = 0; i != grid->size(); ++i) {
      VertexPatch& p = (*grid)[i];
      if (p.isRefined()) {
         for (std::size_t j = 0; j != p.size(); ++j) {
            _T& x = p[j];
            if (x != x) {
               x = -Inf;
            }
            else if (x < 0) {
               x = Inf;
            }
         }
      }
      else {
         _T& x = p.fillValue;
         if (x != x) {
            x = -Inf;
         }
         else if (x < 0) {
            x = Inf;
         }
      }
   }

   // Offset to get the solvent accessible surface.
   *grid -= probeRadius;
}


// Subtract the balls from the level set. Compute the distance up to
// maxDistance past the surface of the balls.
template<typename _T, std::size_t _D>
inline
void
subtract(array::SimpleMultiArrayRef<_T, _D>* grid,
         const geom::BBox<_T, _D>& domain,
         const std::vector<geom::Ball<_T, _D> >& balls,
         const _T maxDistance) {
   typedef array::SimpleMultiArrayRef<_T, _D> SimpleMultiArray;
   typedef typename SimpleMultiArray::IndexList IndexList;
   typedef typename SimpleMultiArray::Range Range;
   typedef array::SimpleMultiIndexRangeIterator<_D> Iterator;

   geom::SimpleRegularGrid<_T, _D> regularGrid(grid->extents(), domain);
   geom::BBox<_T, _D> window;
   Range range;
   std::tr1::array<_T, _D> p;
   for (std::size_t i = 0; i != balls.size(); ++i) {
      // Build a bounding box around the ball.
      window = bound(balls[i]);
      // Expand by two grid cells.
      window.offset(maxDistance);
      // Compute the index range for the bounding box.
      range = regularGrid.computeRange(window);
      // For each index in the range.
      const Iterator end = Iterator::end(range);
      for (Iterator index = Iterator::begin(range); index != end; ++index) {
         // Convert the index to a Cartesian location.
         p = regularGrid.indexToLocation(*index);
         // Compute the signed distance to the surface of the ball.
         // Use this to subtract the ball at this point.
         (*grid)(*index) = difference((*grid)(*index), distance(balls[i], p));
      }
   }
}


// Subtract the balls from the level set. 
template<typename _T, std::size_t _D, std::size_t _N, typename _Base>
inline
void
subtract(array::EquilateralArrayImp<_T, _D, _N, _Base>* patch,
         const std::tr1::array<_T, _D>& lowerCorner,
         const _T spacing,
         const std::vector<geom::Ball<_T, _D> >& balls) {
   typedef array::SimpleMultiIndexRangeIterator<_D> Iterator;

   std::tr1::array<_T, _D> p;
   // For each grid point.
   const Iterator end = Iterator::end(patch->extents());
   for (Iterator i = Iterator::begin(patch->extents()); i != end; ++i) {
      _T& g = (*patch)(*i);
      p = lowerCorner + spacing * *i;
      // For each ball.
      for (std::size_t j = 0; j != balls.size(); ++j) {
         // Compute the signed distance to the surface of the ball.
         // Use this to subtract the ball at this point.
         g = difference(g, distance(balls[j], p));
      }
   }
}


template<typename _T, std::size_t _D>
inline
void
solventExcludedCavities(GridUniform<_T, _D>* grid,
                        const std::vector<geom::Ball<_T, _D> >& balls,
                        const _T probeRadius, const _T maxDistance) {
   // Compute the solvent excluded surface.
   solventExcluded(grid, balls, probeRadius);
   // Subtract the balls. 
   subtract(grid, grid->domain(), balls, maxDistance);
}


template<typename _T, std::size_t _D, std::size_t _N>
inline
void
solventExcludedCavities(Grid<_T, _D, _N>* grid,
                        const std::vector<geom::Ball<_T, _D> >& balls,
                        const _T probeRadius, const _T maxDistance) {
   typedef array::SimpleMultiArrayRef<_T, _D> MultiArrayRef;
   typedef typename MultiArrayRef::IndexList IndexList;
   typedef typename MultiArrayRef::Range Range;
   typedef array::SimpleMultiIndexRangeIterator<_D> Iterator;

   // Determine the patch/ball dependencies for the solvent-excluded surface.
   array::StaticArrayOfArrays<std::size_t> solvExcDependencies;
   solventExcludedDependencies(*grid, balls, probeRadius,
                               &solvExcDependencies);

   // Determine the patch/ball dependencies for subtracting the atoms.
   array::StaticArrayOfArrays<std::size_t> differenceDependencies;
   {
      std::vector<geom::Ball<_T, _D> > offsetBalls = balls;
      for (std::size_t i = 0; i != offsetBalls.size(); ++i) {
         offsetBalls[i].radius += maxDistance;
      }
      patchDependencies(*grid, offsetBalls.begin(), offsetBalls.end(),
                        &differenceDependencies);
   }

   // Use the solvent-excluded dependencies to determine the patches to refine.
   // Refine the appropriate patches and set the rest to have an unknown
   // distance.
   grid->refine(solvExcDependencies);

   // Use a multi-array to wrap the patches.
   const IndexList patchExtents = ext::filled_array<IndexList>(_N);
   array::SimpleMultiArrayRef<_T, _D> patch(0, patchExtents);
   std::vector<geom::Ball<_T, _D> > influencingBalls;
   geom::BBox<_T, _D> domain;
   // Loop over the patches.
   const Iterator end = Iterator::end(grid->extents());
   for (Iterator i = Iterator::begin(grid->extents()); i != end; ++i) {
      const std::size_t index = grid->arrayIndex(*i);
      if (! (*grid)[index].isRefined()) {
         continue;
      }
      // Build the parameters for the solvent-excluded surface.
      patch.rebuild((*grid)[index].data(), patchExtents);
      influencingBalls.clear();
      for (std::size_t n = 0; n != solvExcDependencies.size(index);
           ++n) {
         influencingBalls.push_back(balls[solvExcDependencies(index, n)]);
      }
      // Compute the solvent excluded surface.
      domain = grid->getVertexPatchDomain(*i);
      solventExcluded(&patch, domain, influencingBalls, probeRadius);

      // Build the parameters for subtracting the balls.
      influencingBalls.clear();
      for (std::size_t n = 0; n != differenceDependencies.size(index); ++n) {
         influencingBalls.push_back(balls[differenceDependencies(index, n)]);
      }
      // Subtract the balls.
      subtract(&patch, domain, influencingBalls, maxDistance);
   }
}


} // namespace levelSet
