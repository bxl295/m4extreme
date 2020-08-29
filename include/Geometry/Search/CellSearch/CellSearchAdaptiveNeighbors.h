// -*- C++ -*-

/*!
  \file CellSearchAdaptiveNeighbors.h
  \brief A class for computing neighbors within a specified radius.
*/

#if !defined(__geom_CellSearchAdaptiveNeighbors_h__)
#define __geom_CellSearchAdaptiveNeighbors_h__

#include "packages/geom/kernel/BBox.h"
#include "packages/array/SimpleMultiArray.h"
#include "packages/array/SimpleMultiIndexRangeIterator.h"
#include "packages/ads/functor/Dereference.h"
#include "packages/ads/algorithm/sort.h"

#include "Set/SetLib.h"
#include "Geometry/Search/Search.h"

#include <cstring>

namespace geom {

//! A class for computing neighbors within a specified radius.
/*!
  \param _D The space dimension.
  \param _Record The record type, which is most likely a pointer to a class
  or an iterator into a container.
  \param _Location A functor that takes the record type as an argument
  and returns the location for the record.

  \par The Location functor.
  Suppose that one has the following Node structure for representing a node
  in a 3-D simulation. The \c coords member points to some externally allocated
  memory.
  \code
  struct Node {
     double* coords;
  };
  \endcode
  Next suppose that the record type is a pointer to Node. Below is a functor
  that converts a record to a Cartesian point (\c std::tr1::array<double,3>).
  \code
  struct Location :
     public std::unary_function<Node*, std::tr1::array<double,3> > {
     result_type
     operator()(argument_type r) {
        result_type location = {{r->coords[0], r->coords[1], r->coords[2]}};
        return location;
     }
  };
  \endcode

  \par Constructor.
  In most circumstances, one constructs the neighbors search data structure
  without any arguments. One must pass a location functor
  only if it does not have a default constructor.
  \code
  const double searchRadius = 1;
  geom::CellSearchAdaptiveNeighbors<3, Node*, Location> neighborSearch;
  \endcode

  \par Neighbor queries.
  First use the \c initialize() member function to register the records. 
  These will be sorted and stored in the cell array. One must re-initialize
  whenever the locations of the records are modified.
  Next use the \c neighborsQuery() member function to find the neighbors
  of specified points. Specifically, it finds the records within a ball
  with given center and radius.
  \code
  std::vector<Node> nodes;
  std::vector<geom::Ball<double, 3> > balls;
  ...
  neighborSearch.initialize(&nodes[0], &nodes[0] + nodes.size());
  std::vector<Node*> neighbors;
  // For each query ball.
  for (std::size_t i = 0; i != balls.size(); ++i) {
     neighborSearch.neighborQuery(balls[i].center, balls[i].radius, &neighbors);
     // For each neighbor of this point.
     for (std::size_t j = 0; j != neighbors.size(); ++j) {
        // Do something with the neighbor.
        foo(neighbors[j]);
     }
  }
  \endcode

  \par Neighbor queries and %ORQ data structures.
  One performs neighbor queries by first bounding the search ball with an
  axis-oriented bounding box. One then uses an orthogonal range query 
  (%ORQ) data structure to find all of the records in the bounding box. 
  Octrees, K-d trees, and cell arrays are common %ORQ data structures.
  The tree data structures are heirarchical; the leaves store small sets
  of records with a fixed (small) maximum size. A dense cell array covers
  the bounding box of the records with a uniform array of cells. For most
  record distributions, cell arrays are more efficient than tree data
  structures. One can directly access and iterate over the cells that
  intersect the query window.

  \par Efficient cell representation.
  The simplest way to store a cell array would be use an array of 
  variabled-sized containers. For example, a multi-dimensional array 
  of \c std::vector's.
  \code
   array::SimpleMultiArray<std::vector<Record>, D> cellArray;
  \endcode
  This would be very inefficient both in terms of storage and cache
  utilization. We store pairs of the records and their locations in
  a packed vector. (We store the location to avoid the cost of 
  recomputing it from the record.) The record/location pairs are 
  sorted according to their container index in the multi-array.
  (Note that the sorting may be done in linear time.) The cell array
  is simply a multi-array of iterators into the vector of 
  record/location pairs. The iterator for a given cell points to the
  first record in that cell. One uses the following cell to obtain
  the end iterator.

  \par Cell size.
  Choosing a good cell size is important for the performance of %ORQ's 
  with cell arrays. If the cells are much larger than the query windows
  then one will have to examine many more candidate records than the 
  number of neighbors for any particular query. If the cells are much 
  smaller than the query window, then the cost of iterating over cells
  may dominate. If the search radius is fixed, a common approach is
  to set the cell lengths equal to the search radius. However, this 
  is problematic because, depending on the distribution of the record
  locations, there may be many more cells than records. The solution
  is to set the number of cells to be approximately equal to the number
  of records and then compute appropriate cell dimensions based on 
  that constraint.

  \par 3-D Performance.
  We perform neighbor queries for a set of records uniformly distributed
  in a cube whose volume is equal to the number of records. We perform a 
  query around the center of each record and choose a 
  search radius of 1.4 to yield roughly 10 neighbors per query. Below
  we show performance results as we vary the number of records. The test
  is conducted on a MacBook Pro with a 2.8 GHz Intel Core 2 Duo processor
  with 4 GB of 1067 MHz DDR3 RAM. Note that for the largest test case, 
  cache misses become the dominant cost.
  <table>
  <tr>
  <th> Records
  <th> Average # of neighbors
  <th> Initialization time
  <th> Query time per neighbor
  <tr>
  <td> 1,000
  <td> 10.7
  <td> 0.2 milliseconds
  <td> 59 nanoseconds
  <tr>
  <td> 10,000
  <td> 11.7
  <td> 6 milliseconds
  <td> 53 nanoseconds
  <tr>
  <td> 100,000
  <td> 12.1
  <td> 72 milliseconds
  <td> 59 nanoseconds
  <tr>
  <td> 1,000,000
  <td> 12.3
  <td> 322 milliseconds
  <td> 251 nanoseconds
  </table>
*/
template<std::size_t _D, typename T>
class CellSearchAdaptiveNeighbors : public Geometry::Search<T> {
   //
   // Constants.
   //
public:

   //! The space dimension.
   static const std::size_t D = _D;

   //
   // Types.
   //
public:
   typedef typename Geometry::Search<T>  BASE;
   typedef typename BASE::point_type  point_type;
   typedef typename BASE::Record      Record;
   typedef std::map<Record, point_type*> pointset_type;

   //! A Cartesian point.
   typedef std::tr1::array<double, D> Point;

protected:

   //! Bounding box.
   typedef geom::BBox<double, D> BBox;
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
      Point location;
   };

   typedef typename std::vector<RecLoc>::const_iterator ConstIterator;

   //
   // Data
   //

public:

   //! The Cartesian location of the lower corner of the cell array.
   Point _lowerCorner;
   Point _upperCorner;

private:

   //! The records along with the cell indices and locations.
   std::vector<RecLoc> _recordData;
   //! The array of cells.
   /*! The array is padded by one empty slice in the first dimension. 
     This makes it easier to iterate over a range of cells. */
   array::SimpleMultiArray<ConstIterator, D> _cellArray;

   //! The inverse cell lengths.
   /*! This is used for converting locations to cell multi-indices. */
   Point _inverseCellLengths;

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
   CellSearchAdaptiveNeighbors(const CellSearchAdaptiveNeighbors&);

   // Assignment operator not implemented.
   CellSearchAdaptiveNeighbors&
   operator=(const CellSearchAdaptiveNeighbors&);

   //--------------------------------------------------------------------------
   //! \name Constructors.
   // @{
public:

   //! Construct from the search radius.
   CellSearchAdaptiveNeighbors(const pointset_type & points_) :
      _recordData(),
      _cellArray(),
      // Fill with invalid values.
      _lowerCorner(ext::filled_array<Point>
                          (std::numeric_limits<double>::quiet_NaN())),
      _inverseCellLengths(ext::filled_array<Point>
                          (std::numeric_limits<double>::quiet_NaN())) {
     initialize(points_.begin(), points_.end());
   }

   //! Find the records that are in the specified ball.
   /*! Store the neighbors in the supplied container, which will first be 
     cleared. */
   void
   operator () (const point_type & p, double radius,
		std::vector<Record> & neighbors) {
      neighbors.clear();
      neighborQuery(p, radius, std::back_inserter(neighbors));
   }

   const set<Record> *
   operator () (Record cy) const {
     assert(false); // Not implemented;
   }

   // @}
   //--------------------------------------------------------------------------
   //! \name Window Queries.
   // @{
private:

   typedef typename pointset_type::const_iterator _InputIterator;
  
   //! Initialize with the given sequence of records.   
   void
   initialize(_InputIterator begin, _InputIterator end);

   //! Find the records that are in the specified ball.
   template<typename _OutputIterator>
   void
   neighborQuery(const point_type & p, double radius, _OutputIterator neighbors);

protected:

   //! Convert a location to a valid cell array multi-index.
   IndexList
   locationToIndices(const Point& x);

   //! Convert a location to a container index.
   Index
   containerIndex(const Point& x);

   //! Sort _recordData by cell container indices.
   void
   cellSort();

   //! Compute the array extents and the sizes for the cells.
   IndexList
   computeExtentsAndSizes(std::size_t numberOfCells, const BBox& domain);

   // @}
};

} // namespace geom

#define __geom_CellSearchAdaptiveNeighbors_ipp__
#include "CellSearchAdaptiveNeighbors.ipp"
#undef __geom_CellSearchAdaptiveNeighbors_ipp__

// // 3-D specialization.
// #define __geom_CellSearchAdaptiveNeighbors3_ipp__
// #include "CellSearchAdaptiveNeighbors3.ipp"
// #undef __geom_CellSearchAdaptiveNeighbors3_ipp__

#endif
