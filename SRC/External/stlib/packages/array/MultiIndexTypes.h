// -*- C++ -*-

/*!
  \file array/MultiIndexTypes.h
  \brief Types for multi-array indexing.
*/

#if !defined(__array_MultiIndexTypes_h__)
#define __array_MultiIndexTypes_h__

#include "IndexTypes.h"
#include "MultiArrayStorage.h"

namespace array {

//! Types for %array indexing.
template<std::size_t _Dimension>
class MultiIndexTypes {
   //
   // Constants.
   //
public:

   //! The number of dimensions.
   static const std::size_t Dimension = _Dimension;

   //
   // Public types.
   //
public:

   // Types for STL compliance.

   //! The size type.
   typedef IndexTypes::size_type size_type;
   //! Pointer difference type.
   typedef IndexTypes::difference_type difference_type;

   // Other types.

   //! An %array index is a signed integer.
   typedef IndexTypes::Index Index;
   //! A list of sizes.
   typedef std::tr1::array<size_type, _Dimension> SizeList;
   //! A list of indices.
   typedef std::tr1::array<Index, _Dimension> IndexList;
   //! The storage order.
   typedef MultiArrayStorage<_Dimension> Storage;
};

} // namespace array

#endif
