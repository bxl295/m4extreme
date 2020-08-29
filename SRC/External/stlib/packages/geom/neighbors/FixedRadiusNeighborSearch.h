// -*- C++ -*-

/*!
  \file geom/neighbors/FixedRadiusNeighborSearch.h
  \brief A class for fixed-radius neighbor search on points.
*/

#if !defined(__geom_FixedRadiusNeighborSearch_h__)
#define __geom_FixedRadiusNeighborSearch_h__

#include "../orq/CellArrayStatic.h"

#include "../../array/StaticArrayOfArrays.h"

namespace geom {

//! Fixed-radius neighbor search on points.
/*!
  Inherit from geom::CellArrayStatic so that we have access to the
  multi-key accessor.

  We perform tests on a MacBook Pro with a 2.8 GHz Intel Core 2 Duo processor
  and 4 GB of 1067 MHz DDR3 memory. Below are some performance timings for
  various mathematical operations.

  <table>
  <tr>
  <th> Function
  <th> Time (ns)
  <tr>
  <td> multiplication
  <td> 3.29094
  <tr>
  <td> division
  <td> 8.492
  <tr>
  <td> sqrt()
  <td> 6.99952
  <tr>
  <td> exp()
  <td> 13.3119
  <tr>
  <td> log()
  <td> 23.4017
  </table>

  We perform neighbors queries on a set of uniformly randomly distributed
  points in a 3-D equilateral cube. The volume of the cube matches the number
  of records (points). First we perform neighbor queries for a search radius
  of 5. We vary the number of records.

  <table>
  <tr>
  <th> # Records
  <th> # Reported
  <th> Construction (s)
  <th> Reporting (s)
  <th> Time per query (\f$\mu\f$s)
  <th> Time per reported record (ns)
  <tr>
  <td> 100
  <td> 9598
  <td> 2.6084e-05
  <td> 0.000190526
  <td> 1.90526
  <td> 19.8506
  <tr>
  <td> 1,000
  <td> 269394
  <td> 8.0808e-05
  <td> 0.00858815
  <td> 8.58815
  <td> 31.8795
  <tr>
  <td> 10,000
  <td> 3976238
  <td> 0.000626642
  <td> 0.156945
  <td> 15.6945
  <td> 39.4706
  <tr>
  <td> 100,000
  <td> 46255506
  <td> 0.00829673
  <td> 2.47815
  <td> 24.7815
  <td> 53.5752
  <tr>
  <td> 1,000,000
  <td> 494567104
  <td> 0.228114
  <td> 93.455
  <td> 93.455
  <td> 188.963
  </table>

  New we perform neighbor queries for 10,000 records, varying the search radius.
  <table>
  <tr>
  <th> %Search Radius
  <th> # Reported
  <th> Construction (s)
  <th> Reporting (s)
  <th> Time per query (\f$\mu\f$s)
  <th> Time per reported record (ns)
  <tr>
  <td> 1
  <td> 39376
  <td> 0.000915934
  <td> 0.00774989
  <td> 0.774989
  <td> 196.818
  <tr>
  <td> 2
  <td> 300734
  <td> 0.000807395
  <td> 0.024397
  <td> 2.4397
  <td> 81.1247
  <tr>
  <td> 3
  <td> 961394
  <td> 0.00089317
  <td> 0.055341
  <td> 5.5341
  <td> 57.5633
  <tr>
  <td> 4
  <td> 2155226
  <td> 0.000787288
  <td> 0.100605
  <td> 10.0605
  <td> 46.6796
  <tr>
  <td> 5
  <td> 3976238
  <td> 0.000884905
  <td> 0.157609
  <td> 15.7609
  <td> 39.6377
  <tr>
  <td> 10
  <td> 23288132
  <td> 0.000677505
  <td> 0.624287
  <td> 62.4287
  <td> 26.8071
  <tr>
  <td> 20
  <td> 84828450
  <td> 0.000828621
  <td> 1.43156
  <td> 143.156
  <td> 16.8759
  </table>
*/
template<std::size_t N,
         typename _Record,
         typename _MultiKey =
         typename std::iterator_traits<_Record>::value_type,
         typename _Key = typename _MultiKey::value_type,
         typename _MultiKeyAccessor = ads::Dereference<_Record> >
class FixedRadiusNeighborSearch :
   public CellArrayStatic<N, _Record, _MultiKey, _Key, _MultiKeyAccessor> {
   //
   // Private types.
   //
private:

   typedef CellArrayStatic<N, _Record, _MultiKey, _Key, _MultiKeyAccessor> Base;

   //
   // Public types.
   //
public:

   //! A pointer to the record type.
   typedef _Record Record;
   //! The multi-key type.
   typedef _MultiKey MultiKey;
   //! The key type.
   typedef _Key Key;

   //! A Cartesian point.
   typedef typename Base::Point Point;

   //
   // Protected types.
   //
protected:

   //! Bounding box.
   typedef typename Base::BBox BBox;
   //! A single index.
   typedef typename Base::Index Index;
   //! A multi-index.
   typedef typename Base::IndexList IndexList;

   //
   // Data
   //
private:

   //! The first record. This is used to compute record indices.
   const Record _recordsBegin;
   //! The search radius.
   const double _radius;
   //! The squared search radius.
   const double _squaredRadius;
   //! The bounding box for the ORQ.
   BBox _boundingBox;
   //! The records in the box.
   std::vector<Record> _recordsInBox;

protected:

   using Base::_multiKeyAccessor;

private:

   //
   // Not implemented
   //

   // Default constructor not implemented.
   FixedRadiusNeighborSearch();

   // Copy constructor not implemented
   FixedRadiusNeighborSearch(const FixedRadiusNeighborSearch&);

   // Assignment operator not implemented
   FixedRadiusNeighborSearch&
   operator=(const FixedRadiusNeighborSearch&);

   //--------------------------------------------------------------------------
   //! \name Constructors.
   // @{
public:

   //! Construct from a range of records and the search radius.
   /*!
     \pre There must be a non-zero number of records.

     \param first The first record.
     \param last The last record.
     \param radius The search radius.

     \c first and \c last must be random access iterators because they are
     used to compute record indices.
   */
   template<class _RecordRandomAccessIterator>
   FixedRadiusNeighborSearch(_RecordRandomAccessIterator first,
                             _RecordRandomAccessIterator last,
                             const double radius) :
      Base(first, last),
      _recordsBegin(first),
      _radius(radius),
      _squaredRadius(radius* radius),
      _boundingBox(),
      _recordsInBox() {
      assert(radius >= 0);
   }

   // @}
   //--------------------------------------------------------------------------
   /*! \name Accesors.
     Functionality inherited from CellArrayStatic.
   */
   // @{
public:

   //! Return the number of records.
   using Base::getSize;

   // @}
   //--------------------------------------------------------------------------
   /*! \name Manipulators.
     Functionality inherited from CellArrayStatic.
   */
   // @{
public:

   //! Rebuild for new or modified records.
   using Base::rebuild;

   // @}
   //--------------------------------------------------------------------------
   //! \name Searching.
   // @{
public:

   //! Find the neighbors of the specified record.
   /*!
     Report the indices of the records that are within the search radius.
   */
   template<typename _IndexOutputIterator>
   void
   findNeighbors(_IndexOutputIterator iter, std::size_t recordIndex);

   // @}
};

} // namespace geom

#define __geom_FixedRadiusNeighborSearch_ipp__
#include "FixedRadiusNeighborSearch.ipp"
#undef __geom_FixedRadiusNeighborSearch_ipp__

#endif
