// -*- C++ -*-

#if !defined(__geom_CellSearchAllNeighbors_ipp__)
#error This file is an implementation detail of the class CellSearchAllNeighbors.
#endif

#include <string.h>

namespace geom{

template<unsigned int _Dim>
inline
void
CellSearchAllNeighbors<_Dim>::
allNeighbors(_InputIterator begin, _InputIterator end) {
   // Clear results from previous queries.
   records.clear();
   packedNeighbors.clear();
   neighborDelimiters.clear();
   _recordData.clear();
   // Dispense with the trivial case.
   if (begin == end) {
      return;
   }
   // Copy the records and compute the locations.
   RecLoc ril;
   while (begin != end) {
      ril.record = begin->first;

      for ( std::size_t i = 0; i != _Dim; ++i )
	ril.location[i] = begin->second[i];

      _recordData.push_back(ril);

      begin++;
   }
   // Bound the locations.
   BoundedBox box = {_recordData[0].location, _recordData[0].location};
   for (std::size_t i = 1; i != _recordData.size(); ++i) {
      box.add(_recordData[i].location);
   }

   // Expand to avoid errors in converting locations to cell multi-indices.
   box.offset(std::sqrt(std::numeric_limits<double>::epsilon()) *
              (1. + max(box.upper - box.lower)));
   // Set the lower corner and the cell array extents.
   _lowerCorner = box.lower;
   _upperCorner = box.upper;

   _cellArray.rebuild(computeExtentsAndSizes(_recordData.size(), box));
   // Sort by the cell indices and define the cells.
   cellSort();
   // Record the order of the records.
   records.resize(_recordData.size());
   for (std::size_t i = 0; i != records.size(); ++i) {
      records[i] = _recordData[i].record;
   }
   // Resize the delimiters.
   neighborDelimiters.resize(_recordData.size() + 1);
   // Perform neighbor queries for each record.
   for (std::size_t i = 0; i != _recordData.size(); ++i) {
      neighborDelimiters[i] = packedNeighbors.size();
      // Record the neighbors in the vector of packed neighbors.
      recordNeighbors(i);
   }
   // Set the last delimiter.
   neighborDelimiters.back() = packedNeighbors.size();
}


template<unsigned int _Dim>
inline
void
CellSearchAllNeighbors<_Dim>::
recordNeighbors(const std::size_t index) {
   typedef array::SimpleMultiIndexRange<_Dim-1> Range;
   typedef array::SimpleMultiIndexRangeIterator<_Dim-1> RangeIterator;

   // The record and the location for the query.
   const Record record = _recordData[index].record;
   const Point x = _recordData[index].location;
   // The window for the ORQ has corners at x - searchRadius and
   // x + searchRadius. Convert the corners to cell array indices.
   IndexList lo = locationToIndices(x - searchRadius);
   IndexList hi = locationToIndices(x + searchRadius);
   // Iterate over cells in all except the first dimension.
   typename Range::IndexList extents, bases;
   for (std::size_t i = 0; i != _Dim - 1; ++i) {
      extents[i] = hi[i+1] - lo[i+1] + 1;
      bases[i] = lo[i+1];
   }
   IndexList start, stop;
   start[0] = lo[0];
   stop[0] = hi[0] + 1;
   const double squaredRadius = searchRadius * searchRadius;
   const Range range = {extents, bases};
   const RangeIterator end = RangeIterator::end(range);
   for (RangeIterator i = RangeIterator::begin(range); i != end; ++i) {
      for (std::size_t d = 0; d != _Dim - 1; ++d) {
         start[d+1] = stop[d+1] = (*i)[d];
      }
      // Iterate over the records in the row.
      const ConstIterator jEnd = _cellArray(stop);
      for (ConstIterator j = _cellArray(start); j != jEnd; ++j) {
         if (j->record != record && 
             squaredDistance(x, j->location) < squaredRadius) {
            packedNeighbors.push_back(j->record);
         }
      }
   }
}


template<unsigned int _Dim>
inline
typename CellSearchAllNeighbors<_Dim>::IndexList
CellSearchAllNeighbors<_Dim>::
locationToIndices(const Point& x) {
   IndexList index;
   for (std::size_t i = 0; i != _Dim; ++i) {
      index[i] = std::min(_cellArray.extents()[i] - 1,
                          Index(std::max(0., x[i] - _lowerCorner[i]) *
                                _inverseCellLengths[i]));
   }
   // Adjust to map to a non-empty cell.
   if (index[0] == _cellArray.extents()[0] - 1) {
      --index[0];
   }
   return index;
}

template<unsigned int _Dim>
inline
typename CellSearchAllNeighbors<_Dim>::Index
CellSearchAllNeighbors<_Dim>::
containerIndex(const Point& x) {
   IndexList index;
   for (std::size_t i = 0; i != _Dim; ++i) {
      index[i] = Index((x[i] - _lowerCorner[i]) * _inverseCellLengths[i]);
      //cout << x[i] << "\t" << _lowerCorner[i] << "\t" << _inverseCellLengths[i] << "\t" << index[i] << endl;
   }
   return _cellArray.arrayIndex(index);
}


template<unsigned int _Dim>
inline
void
CellSearchAllNeighbors<_Dim>::
cellSort() {
   // Calculate the cell container index for each record.
   _cellIndices.resize(_recordData.size());
   // Count the number of records in each cell.
   _cellCounts.resize(_cellArray.size());
   // Initialize to zero.
   memset(&_cellCounts[0], 0, _cellCounts.size() * sizeof(std::size_t));
   for (std::size_t i = 0; i != _recordData.size(); ++i) {
      _cellIndices[i] = containerIndex(_recordData[i].location);
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

#if 0
// No longer used.
template<unsigned int _Dim>
inline
void
CellSearchAllNeighbors<_Dim>::
plainSort() {
   // Calculate the cell container index for each record.
   for (std::size_t i = 0; i != _recordData.size(); ++i) {
      _recordData[i].cellIndex = containerIndex(_recordData[i].location);
   }
   // Sort by the cell indices.
   CompareIndex compare;
   std::sort(_recordData.begin(), _recordData.end(), compare);
   // Define the cells.
   {
      ConstIterator iter = _recordData.begin();
      for (std::size_t i = 0; i != _cellArray.size(); ++i) {
         _cellArray[i] = iter;
         while (iter != _recordData.end() && iter->cellIndex <= i) {
            ++iter;
         }
      }
   }
}
#endif


// Compute the array extents and the sizes for the cells.
template<unsigned int _Dim>
inline
typename CellSearchAllNeighbors<_Dim>::IndexList
CellSearchAllNeighbors<_Dim>::
computeExtentsAndSizes(const std::size_t numberOfCells, const BoundedBox& domain) {
   assert(numberOfCells > 0);

   // Work from the the least to greatest Cartesian extent to compute the
   // grid extents.
   IndexList extents;
   IndexList order;
   Point ext = domain.upper - domain.lower;
   ads::computeOrder(ext.begin(), ext.end(), order.begin());
   for (std::size_t i = 0; i != _Dim; ++i) {
      // Normalize the domain to numberOfCells content.
      double content = product(ext);
      assert(content != 0);
      const double factor = std::pow(numberOfCells / content,
                                     double(1.0 / (_Dim - i)));
      for (std::size_t j = i; j != _Dim; ++j) {
         ext[order[j]] *= factor;
      }
      // The current index;
      const std::size_t n = order[i];
      // Add 0.5 and truncate to round to the nearest integer.
      ext[n] = extents[n] = std::max(std::size_t(1),
                                     std::size_t(ext[n] + 0.5));
   }

   // From the domain and the cell array extents, compute the cell size.
   for (std::size_t n = 0; n != _Dim; ++n) {
      const double d = domain.upper[n] -
         domain.lower[n];
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

//    for ( std::size_t n = 0; n != _Dim; ++n ) {
//      cout << extents[n] << endl;
//    }

   return extents;
}


} // namespace geom
