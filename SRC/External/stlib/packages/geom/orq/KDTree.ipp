// -*- C++ -*-

#if !defined(__geom_KDTree_ipp__)
#error This file is an implementation detail of the class KDTree.
#endif

namespace geom {

//-----------------------------KDTreeBranch-------------------------------

//
// Constructors
//

// Construct from a set of grid elements.
template < std::size_t N, typename _Record, typename _MultiKey, typename _Key,
         typename _MultiKeyAccessor, typename _RecordOutputIterator >
inline
KDTreeBranch<N, _Record, _MultiKey, _Key, _MultiKeyAccessor, _RecordOutputIterator>::
KDTreeBranch(const std::tr1::array<std::vector<Record>, N>& sorted,
             const std::size_t leafSize) {
#ifdef DEBUG_stlib
   // There must be more records than in a single leaf.
   assert(std::size_t(sorted[0].size()) > leafSize);
   // Each of the sorted arrays should be of the same size.
   for (std::size_t n = 1; n != N; ++n) {
      assert(sorted[n].size() == sorted[0].size());
   }
#endif
   //
   // Determine the splitting direction.
   //

   std::tr1::array<Key, N> spreads;
   for (std::size_t n = 0; n != N; ++n) {
      spreads[n] = _multiKeyAccessor(sorted[n].back())[n]
                   - _multiKeyAccessor(sorted[n].front())[n];
   }

   // The index of the maximum spread.
   _splitDimension =
      std::max_element(spreads.begin(), spreads.end()) - spreads.begin();

   //
   // Name the input vectors.
   //

   const std::vector<Record>& splitSorted = sorted[_splitDimension];

   //
   // Compute the median.
   //

   const std::size_t medianIndex = sorted[0].size() / 2;
   const std::size_t leftSize = medianIndex;
   const std::size_t rightSize = sorted[0].size() - medianIndex;
   const MultiKey& medianPoint = _multiKeyAccessor(splitSorted[medianIndex]);
   _splitValue = medianPoint[_splitDimension];

   //
   // Vectors for the subtrees.
   //
   std::tr1::array<std::vector<Record>, N> sub;
   for (std::size_t n = 0; n != N; ++n) {
      sub[n].reserve(rightSize);
   }
   typename std::vector<Record>::const_iterator iter;

   //
   // Make the left subtree.
   //

   std::copy(splitSorted.begin(), splitSorted.begin() + medianIndex,
             std::back_inserter(sub[_splitDimension]));

   // For each dimension except the split dimension.
   for (std::size_t n = 0; n != N; ++n) {
      if (n != _splitDimension) {
         for (iter = sorted[n].begin(); iter != sorted[n].end(); ++iter) {
            if (ads::less_composite_fcn<N>(_splitDimension,
                                           _multiKeyAccessor(*iter),
                                           medianPoint)) {
               sub[n].push_back(*iter);
            }
         }
      }
   }

   // If the left subtree is a leaf.
   if (leftSize <= leafSize) {
      std::vector<Record> leftLeaf(splitSorted.begin(),
                                   splitSorted.begin() + leftSize);
      _left = new Leaf(leftLeaf);
   }
   else {
      _left = new KDTreeBranch(sub, leafSize);
   }

   for (std::size_t n = 0; n != N; ++n) {
      sub[n].clear();
   }

   //
   // Make the right subtree.
   //

   std::copy(splitSorted.begin() + medianIndex, splitSorted.end(),
             std::back_inserter(sub[_splitDimension]));

   // For each dimension except the split dimension.
   for (std::size_t n = 0; n != N; ++n) {
      if (n != _splitDimension) {
         for (iter = sorted[n].begin(); iter != sorted[n].end(); ++iter) {
            if (! ads::less_composite_fcn<N>(_splitDimension,
                                             _multiKeyAccessor(*iter),
                                             medianPoint)) {
               sub[n].push_back(*iter);
            }
         }
      }
   }

   // If the right subtree is a leaf.
   if (rightSize <= leafSize) {
      std::vector<Record>
      rightLeaf(splitSorted.begin() + medianIndex, splitSorted.end());
      _right = new Leaf(rightLeaf);
   }
   else {
      _right = new KDTreeBranch(sub, leafSize);
   }
}


//
// Window queries.
//


template < std::size_t N, typename _Record, typename _MultiKey, typename _Key,
         typename _MultiKeyAccessor, typename _RecordOutputIterator >
inline
std::size_t
KDTreeBranch<N, _Record, _MultiKey, _Key, _MultiKeyAccessor, _RecordOutputIterator>::
computeWindowQuery(RecordOutputIterator iter, const BBox& window) const {
   if (_splitValue < window.lower[_splitDimension]) {
      return _right->computeWindowQuery(iter, window);
   }
   else if (_splitValue > window.upper[_splitDimension]) {
      return _left->computeWindowQuery(iter, window);
   }
   return (_left->computeWindowQuery(iter, window) +
           _right->computeWindowQuery(iter, window));
}


template < std::size_t N, typename _Record, typename _MultiKey, typename _Key,
         typename _MultiKeyAccessor, typename _RecordOutputIterator >
inline
std::size_t
KDTreeBranch<N, _Record, _MultiKey, _Key, _MultiKeyAccessor, _RecordOutputIterator>::
computeWindowQuery(RecordOutputIterator iter, BBox* domain,
                   const BBox& window)
const {
   std::size_t count = 0;

   // If the domain of the left sub-tree intersects the window.
   if (_splitValue >= window.lower[_splitDimension]) {
      // Make the domain of the left sub-tree.
      Key max = domain->upper[_splitDimension];
      domain->upper[_splitDimension] = _splitValue;

      // If the domain lies inside the window.
      if (window.isIn(*domain)) {
         // Report the records in the left sub-tree.
         count += _left->report(iter);
      }
      else {
         // Do a window query of the left sub-tree.
         count += _left->computeWindowQuery(iter, domain, window);
      }

      // Reset the domain.
      domain->upper[_splitDimension] = max;
   }

   // If the domain of the right sub-tree intersects the window.
   if (_splitValue <= window.upper[_splitDimension]) {
      // Make the domain of the right sub-tree.
      Key min = domain->lower[_splitDimension];
      domain->lower[_splitDimension] = _splitValue;

      // If the domain lies inside the window.
      if (window.isIn(*domain)) {
         // Report the records in the right sub-tree.
         count += _right->report(iter);
      }
      // If the domain intersects the window.
      else {
         // Do a window query of the right sub-tree.
         count += _right->computeWindowQuery(iter, domain, window);
      }

      // Reset the domain.
      domain->lower[_splitDimension] = min;
   }

   return count;
}


//
// Validity check.
//


template < std::size_t N, typename _Record, typename _MultiKey, typename _Key,
         typename _MultiKeyAccessor, typename _RecordOutputIterator >
inline
bool
KDTreeBranch<N, _Record, _MultiKey, _Key, _MultiKeyAccessor, _RecordOutputIterator>::
isValid(const BBox& window) const {
   if (window.lower[_splitDimension] > _splitValue ||
         _splitValue > window.upper[_splitDimension]) {
      return false;
   }

   BBox win(window);
   win.upper[_splitDimension] = _splitValue;
   if (! _left->isValid(win)) {
      return false;
   }

   win = window;
   win.lower[_splitDimension] = _splitValue;
   if (! _right->isValid(win)) {
      return false;
   }

   return true;
}


//-----------------------------KDTreeLeaf-------------------------------


//
// Mathematical member functions
//


template < std::size_t N, typename _Record, typename _MultiKey, typename _Key,
         typename _MultiKeyAccessor, typename _RecordOutputIterator >
inline
std::size_t
KDTreeLeaf<N, _Record, _MultiKey, _Key, _MultiKeyAccessor, _RecordOutputIterator>::
computeWindowQuery(RecordOutputIterator iter, const BBox& window) const {
   std::size_t count = 0;
   ConstIterator recordsEnd = _records.end();
   for (ConstIterator i = _records.begin(); i != recordsEnd; ++i) {
      if (window.isIn(_multiKeyAccessor(*i))) {
         *(iter++) = (*i);
         ++count;
      }
   }
   return count;
}


template < std::size_t N, typename _Record, typename _MultiKey, typename _Key,
         typename _MultiKeyAccessor, typename _RecordOutputIterator >
inline
std::size_t
KDTreeLeaf<N, _Record, _MultiKey, _Key, _MultiKeyAccessor, _RecordOutputIterator>::
computeWindowQuery(RecordOutputIterator iter, BBox* domain,
                   const BBox& window) const {
   std::size_t count = 0;
   if (window.isIn(*domain)) {
      for (ConstIterator i = _records.begin(); i != _records.end(); ++i) {
         *(iter++) = (*i);
      }
      count += _records.size();
   }
   else {
      ConstIterator recordsEnd = _records.end();
      for (ConstIterator i = _records.begin(); i != recordsEnd; ++i) {
         if (window.isIn(_multiKeyAccessor(*i))) {
            *(iter++) = (*i);
            ++count;
         }
      }
   }
   return count;
}









//-----------------------------KDTree-----------------------------------

//
// Constructors
//

template < std::size_t N, typename _Record, typename _MultiKey, typename _Key,
         typename _MultiKeyAccessor, typename _RecordOutputIterator >
template<class InputIterator>
inline
KDTree<N, _Record, _MultiKey, _Key, _MultiKeyAccessor, _RecordOutputIterator>::
KDTree(InputIterator first, InputIterator last, const std::size_t leafSize) :
   Base(),
   _root(0),
   _domain() {
   if (first == last) {
      std::vector<Record> empty;
      _root = new Leaf(empty);
      _domain.lower = ext::filled_array<Point>(0);
      _domain.upper = ext::filled_array<Point>(-1);
      return;
   }

   // Make N vectors of pointers to the records.
   std::tr1::array<std::vector<Record>, N> sorted;
   while (first != last) {
      sorted[0].push_back(first);
      ++first;
   }
   for (std::size_t n = 1; n != N; ++n) {
      sorted[n] = sorted[0];
   }

   // Sort these vectors in each coordinate.
   LessThanComposite comparison;
   for (std::size_t n = 0; n != N; ++n) {
      comparison.set(n);
      std::sort(sorted[n].begin(), sorted[n].end(), comparison);
   }

   // Determine the domain.
   for (std::size_t n = 0; n != N; ++n) {
      _domain.lower[n] = _multiKeyAccessor(sorted[n].front())[n];
      _domain.upper[n] = _multiKeyAccessor(sorted[n].back())[n];
   }

   // The number of records.
   setSize(sorted[0].size());

   // Make the tree.
   if (getSize() > leafSize) {
      _root = new Branch(sorted, leafSize);
   }
   else {
      _root = new Leaf(sorted[0]);
   }
}


//
// File I/O
//


template < std::size_t N, typename _Record, typename _MultiKey, typename _Key,
         typename _MultiKeyAccessor, typename _RecordOutputIterator >
inline
void
KDTree<N, _Record, _MultiKey, _Key, _MultiKeyAccessor, _RecordOutputIterator>::
put(std::ostream& out) const {
   out << getSize() << " records"
       << '\n'
       << "domain = " << getDomain()
       << '\n';
   _root->put(out);
}

} // namespace geom
