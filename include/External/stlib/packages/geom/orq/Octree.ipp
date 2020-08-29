// -*- C++ -*-

// geom/orq/Octree.ipp
// A class for an Octree in 3-D.

#if !defined(__geom_orq_Octree_ipp__)
#error This file is an implementation detail of the class Octree.
#endif

namespace geom {

//-----------------------------OctreeBranch-------------------------------


//
// OctreeBranch Constructors and Destructors
//


template < typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor,
         typename _RecordOutputIterator >
inline
OctreeBranch<_Record, _MultiKey, _Key, _MultiKeyAccessor, _RecordOutputIterator>::
OctreeBranch(const BBox& domain) :
   _midpoint(domain.lower + (domain.upper - domain.lower) / 2.0) {
   Base::_domain = domain;
   for (std::size_t i = 0; i < 8; ++i) {
      _octant[i] = 0;
   }
}


template < typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor,
         typename _RecordOutputIterator >
inline
OctreeBranch<_Record, _MultiKey, _Key, _MultiKeyAccessor, _RecordOutputIterator>::
~OctreeBranch() {
   for (std::size_t i = 0; i < 8; ++i) {
      if (_octant[i]) {
         delete _octant[i];
      }
   }
}


//
// Accesors
//


template < typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor,
         typename _RecordOutputIterator >
inline
std::size_t
OctreeBranch<_Record, _MultiKey, _Key, _MultiKeyAccessor, _RecordOutputIterator>::
getOctantIndex(Record record) const {
   const MultiKey& mk = _multiKeyAccessor(record);
   if (mk[2] < _midpoint[2]) {
      if (mk[1] < _midpoint[1]) {
         if (mk[0] < _midpoint[0])
            return 0;
         else
            return 1;
      }
      else {
         if (mk[0] < _midpoint[0])
            return 2;
         else
            return 3;
      }
   }
   else {
      if (mk[1] < _midpoint[1]) {
         if (mk[0] < _midpoint[0])
            return 4;
         else
            return 5;
      }
      else {
         if (mk[0] < _midpoint[0])
            return 6;
         else
            return 7;
      }
   }
}


template < typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor,
         typename _RecordOutputIterator >
inline
typename OctreeBranch<_Record, _MultiKey, _Key, _MultiKeyAccessor, _RecordOutputIterator>::BBox
OctreeBranch<_Record, _MultiKey, _Key, _MultiKeyAccessor, _RecordOutputIterator>::
getOctantDomain(const std::size_t index) const {
   BBox x;

   if (index & 1) {
      x.lower[0] = _midpoint[0];
      x.upper[0] = Base::_domain.upper[0];
   }
   else {
      x.lower[0] = Base::_domain.lower[0];
      x.upper[0] = _midpoint[0];
   }

   if (index & 2) {
      x.lower[1] = _midpoint[1];
      x.upper[1] = Base::_domain.upper[1];
   }
   else {
      x.lower[1] = Base::_domain.lower[1];
      x.upper[1] = _midpoint[1];
   }

   if (index & 4) {
      x.lower[2] = _midpoint[2];
      x.upper[2] = Base::_domain.upper[2];
   }
   else {
      x.lower[2] = Base::_domain.lower[2];
      x.upper[2] = _midpoint[2];
   }

   return x;
}


//
// Add grid elements.
//


template < typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor,
         typename _RecordOutputIterator >
inline
OctreeNode<_Record, _MultiKey, _Key, _MultiKeyAccessor, _RecordOutputIterator>*
OctreeBranch<_Record, _MultiKey, _Key, _MultiKeyAccessor, _RecordOutputIterator>::
insert(Record record, const std::size_t leafSize) {
   std::size_t index = getOctantIndex(record);
   if (_octant[index] == 0) {
      _octant[index] = new Leaf(getOctantDomain(index));
   }
   _octant[index] = _octant[index]->insert(record, leafSize);
   return this;
}


//
// Mathematical member functions
//


template < typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor,
         typename _RecordOutputIterator >
inline
std::size_t
OctreeBranch<_Record, _MultiKey, _Key, _MultiKeyAccessor, _RecordOutputIterator>::
report(RecordOutputIterator iter) const {
   std::size_t count = 0;
   for (std::size_t i = 0; i < 8; ++i) {
      if (_octant[i]) {
         count += _octant[i]->report(iter);
      }
   }
   return count;
}


template < typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor,
         typename _RecordOutputIterator >
inline
std::size_t
OctreeBranch<_Record, _MultiKey, _Key, _MultiKeyAccessor, _RecordOutputIterator>::
computeWindowQuery(RecordOutputIterator iter, const BBox& window) const {
   std::size_t count = 0;
   if (doOverlap(window, Base::_domain)) {
      for (std::size_t i = 0; i < 8; ++i) {
         if (_octant[i]) {
            count += _octant[i]->computeWindowQuery(iter, window);
         }
      }
   }
   return count;
}


template < typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor,
         typename _RecordOutputIterator >
inline
std::size_t
OctreeBranch<_Record, _MultiKey, _Key, _MultiKeyAccessor, _RecordOutputIterator>::
computeWindowQueryCheckDomain(RecordOutputIterator iter, const BBox& window)
const {
   std::size_t count = 0;
   if (window.isIn(Base::_domain)) {
      for (std::size_t i = 0; i < 8; ++i) {
         if (_octant[i]) {
            count += _octant[i]->report(iter);
         }
      }
   }
   else if (doOverlap(window, Base::_domain)) {
      for (std::size_t i = 0; i < 8; ++i) {
         if (_octant[i]) {
            count += _octant[i]->computeWindowQueryCheckDomain(iter, window);
         }
      }
   }
   return count;
}


//
// File IO
//


template < typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor,
         typename _RecordOutputIterator >
inline
void
OctreeBranch<_Record, _MultiKey, _Key, _MultiKeyAccessor, _RecordOutputIterator>::
put(std::ostream& out) const {
   for (std::size_t i = 0; i < 8; ++i) {
      if (_octant[i]) {
         _octant[i]->put(out);
      }
   }
}


template < typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor,
         typename _RecordOutputIterator >
inline
void
OctreeBranch<_Record, _MultiKey, _Key, _MultiKeyAccessor, _RecordOutputIterator>::
print(std::ostream& out, std::string tabbing) const {
   out << tabbing << Base::_domain << '\n';
   std::string newTabbing(tabbing);
   newTabbing.append("  ");
   for (std::size_t i = 0; i < 8; ++i) {
      if (_octant[i]) {
         _octant[i]->print(out, newTabbing);
      }
   }
}


//
// Memory usage.
//


template < typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor,
         typename _RecordOutputIterator >
inline
std::size_t
OctreeBranch<_Record, _MultiKey, _Key, _MultiKeyAccessor, _RecordOutputIterator>::
getMemoryUsage() const {
   std::size_t count = sizeof(OctreeBranch);
   for (std::size_t i = 0; i < 8; ++i) {
      if (_octant[i]) {
         count += _octant[i]->getMemoryUsage();
      }
   }
   return count;
}


//
// Validity check.
//


template < typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor,
         typename _RecordOutputIterator >
inline
bool
OctreeBranch<_Record, _MultiKey, _Key, _MultiKeyAccessor, _RecordOutputIterator>::
isValid() const {
   for (std::size_t i = 0; i < 8; ++i) {
      if (_octant[i]) {
         if (! _octant[i]->isValid()) {
            return false;
         }
      }
   }
   return true;
}


//-----------------------------OctreeLeaf-------------------------------


//
// Add grid elements.
//


template < typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor,
         typename _RecordOutputIterator >
inline
OctreeNode<_Record, _MultiKey, _Key, _MultiKeyAccessor, _RecordOutputIterator>*
OctreeLeaf<_Record, _MultiKey, _Key, _MultiKeyAccessor, _RecordOutputIterator>::
insert(Record record, const std::size_t leafSize) {
   // If this leaf is full.
   if (_records.size() == leafSize) {
      // Replace the leaf with a branch.
      Branch* branch = new Branch(Base::_domain);
      for (ConstIterator i = _records.begin(); i != _records.end(); ++i) {
         branch->insert(*i, leafSize);
      }
      branch->insert(record, leafSize);
      delete this;
      return branch;
   }
   _records.push_back(record);
   return this;
}


//
// Mathematical member functions
//


template < typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor,
         typename _RecordOutputIterator >
inline
std::size_t
OctreeLeaf<_Record, _MultiKey, _Key, _MultiKeyAccessor, _RecordOutputIterator>::
report(RecordOutputIterator iter) const {
   for (ConstIterator i = _records.begin(); i != _records.end(); ++i) {
      *(iter++) = (*i);
   }
   return _records.size();
}


template < typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor,
         typename _RecordOutputIterator >
inline
std::size_t
OctreeLeaf<_Record, _MultiKey, _Key, _MultiKeyAccessor, _RecordOutputIterator>::
computeWindowQuery(RecordOutputIterator iter, const BBox& window) const {
   // Copy the window into keys for faster access.
   const Key windowXMin = window.lower[0];
   const Key windowYMin = window.lower[1];
   const Key windowZMin = window.lower[2];
   const Key windowXMax = window.upper[0];
   const Key windowYMax = window.upper[1];
   const Key windowZMax = window.upper[2];
   ConstIterator recordsEnd = _records.end();
   std::size_t count = 0;
   for (ConstIterator i = _records.begin(); i != recordsEnd; ++i) {
      const MultiKey& mk = _multiKeyAccessor(*i);
      if (mk[0] >= windowXMin &&
            mk[0] <= windowXMax &&
            mk[1] >= windowYMin &&
            mk[1] <= windowYMax &&
            mk[2] >= windowZMin &&
            mk[2] <= windowZMax) {
         *(iter++) = (*i);
         ++count;
      }
   }
   return count;
}


template < typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor,
         typename _RecordOutputIterator >
inline
std::size_t
OctreeLeaf<_Record, _MultiKey, _Key, _MultiKeyAccessor, _RecordOutputIterator>::
computeWindowQueryCheckDomain(RecordOutputIterator iter,
                              const BBox& window) const {
   std::size_t count = 0;
   if (window.isIn(Base::_domain)) {
      for (ConstIterator i = _records.begin(); i != _records.end(); ++i) {
         *(iter++) = (*i);
      }
      count += _records.size();
   }
   else {
      // copy the window into keys for faster access.
      const Key windowXMin = window.lower[0];
      const Key windowYMin = window.lower[1];
      const Key windowZMin = window.lower[2];
      const Key windowXMax = window.upper[0];
      const Key windowYMax = window.upper[1];
      const Key windowZMax = window.upper[2];
      ConstIterator recordsEnd = _records.end();
      for (ConstIterator i = _records.begin(); i != recordsEnd; ++i) {
         const MultiKey& mk = _multiKeyAccessor(*i);
         if (mk[0] >= windowXMin &&
               mk[0] <= windowXMax &&
               mk[1] >= windowYMin &&
               mk[1] <= windowYMax &&
               mk[2] >= windowZMin &&
               mk[2] <= windowZMax) {
            *(iter++) = (*i);
            ++count;
         }
      }
   }
   return count;
}


//
// File IO
//


template < typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor,
         typename _RecordOutputIterator >
inline
void
OctreeLeaf<_Record, _MultiKey, _Key, _MultiKeyAccessor, _RecordOutputIterator>::
put(std::ostream& out) const {
   for (ConstIterator i = _records.begin(); i != _records.end(); ++i) {
      out << **i << '\n';
   }
}


template < typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor,
         typename _RecordOutputIterator >
inline
void
OctreeLeaf<_Record, _MultiKey, _Key, _MultiKeyAccessor, _RecordOutputIterator>::
print(std::ostream& out, std::string tabbing) const {
   out << tabbing << Base::_domain << '\n';
   for (ConstIterator i = _records.begin(); i != _records.end(); ++i) {
      out << tabbing << **i << '\n';
   }
}


//-----------------------------Octree-----------------------------------


//
// File I/O
//


template < typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor,
         typename _RecordOutputIterator >
inline
void
Octree<_Record, _MultiKey, _Key, _MultiKeyAccessor, _RecordOutputIterator>::
put(std::ostream& out) const {
   out << getSize() << " records"
       << '\n'
       << "domain = " << getDomain()
       << '\n';
   _root->put(out);
}


template < typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor,
         typename _RecordOutputIterator >
inline
void
Octree<_Record, _MultiKey, _Key, _MultiKeyAccessor, _RecordOutputIterator>::
print(std::ostream& out) const {
   out << getSize() << " records" << '\n';
   out << "domain = " << getDomain() << '\n';

   std::string empty;
   _root->print(out, empty);
}

} // namespace geom
