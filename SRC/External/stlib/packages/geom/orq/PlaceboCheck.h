// -*- C++ -*-

/*!
  \file PlaceboCheck.h
  \brief A placebo class.
*/

#if !defined(__geom_PlaceboCheck_h__)
#define __geom_PlaceboCheck_h__

#include "Placebo.h"

namespace geom {

//! A placebo for ORQ's in N-D
/*!
  Stores a vector of records.
*/
template < std::size_t N,
         typename _Record,
         typename _MultiKey =
         typename std::iterator_traits<_Record>::value_type,
         typename _Key = typename _MultiKey::value_type,
         typename _MultiKeyAccessor = ads::Dereference<_Record> >
class PlaceboCheck :
   public Placebo<N, _Record, _MultiKey, _Key, _MultiKeyAccessor> {
private:

   //
   // Private types.
   //

   typedef Placebo<N, _Record, _MultiKey, _Key, _MultiKeyAccessor> Base;

public:

   //
   // Public types.
   //

   //! A pointer to the record type.
   typedef _Record Record;
   //! The multi-key type.
   typedef _MultiKey MultiKey;
   //! The key type.
   typedef _Key Key;

   //! A Cartesian point.
   typedef typename Base::Point Point;
   //! Bounding box.
   typedef typename Base::BBox BBox;

protected:

   using Base::_multiKeyAccessor;

private:

   //
   // Not implemented.
   //

   //! Copy constructor not implemented
   PlaceboCheck(const PlaceboCheck&);

   //! Assignment operator not implemented
   PlaceboCheck&
   operator=(const PlaceboCheck&);

   //
   // Member variables.
   //

   using Base::_records;

public:

   //-------------------------------------------------------------------------
   //! \name Constructors.
   //@{

   //! Default constructor.
   PlaceboCheck() :
      Base() {}

   //! Reserve storage for \c size records.
   explicit
   PlaceboCheck(const std::size_t size) :
      Base(size) {}

   //! Construct from a range of records.
   /*!
     \param first the beginning of a range of records.
     \param last the end of a range of records.
   */
   template <typename InputIterator>
   PlaceboCheck(InputIterator first, InputIterator last) :
      Base(first, last) {}

   //@}
   //-------------------------------------------------------------------------
   //! \name Window queries.
   //@{

   //! Get the records in the window.  Return the # of records inside.
   template<typename OutputIterator>
   std::size_t
   computeWindowQuery(OutputIterator iter, const BBox& window) const;

   //@}
   //-------------------------------------------------------------------------
   //! \name Insert records.
   //@{

   //! Insert a single record or a range of records.
   using Base::insert;

   //@}
   //-------------------------------------------------------------------------
   //! \name Memory usage.
   //@{

   //! Return the total memory usage.
   using Base::getMemoryUsage;

   //@}
   //-------------------------------------------------------------------------
   //! \name Accessors.
   //@{

   //
   // Inherited from ORQ.
   //

   //! Return the number of records.
   using Base::getSize;

   //! Return true if the grid is empty.
   using Base::isEmpty;

   //
   // Inherited from Placebo.
   //

   //! Return the number of records to return with each query.
   using Base::getQuerySize;

   //@}
   //-------------------------------------------------------------------------
   //! \name Manipulators.
   //@{

   //! Set the number of records to return with each query.
   using Base::setQuerySize;

   //! Shuffle the record pointers.
   using Base::shuffle;

   //@}
   //-------------------------------------------------------------------------
   //! \name File I/O.
   //@{

   //! Print the records.
   using Base::put;

   //@}

protected:

   //! Increment the number of records.
   using Base::incrementSize;

   //! Return a starting point for the window query.
   using Base::getStartingPoint;
};

//! Write to a file stream.
/*! \relates PlaceboCheck */
template < std::size_t N,
         typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor >
inline
std::ostream&
operator<<(std::ostream& out,
           const PlaceboCheck<N, _Record, _MultiKey, _Key, _MultiKeyAccessor>& x) {
   x.put(out);
   return out;
}

} // namespace geom

#define __geom_PlaceboCheck_ipp__
#include "PlaceboCheck.ipp"
#undef __geom_PlaceboCheck_ipp__

#endif
