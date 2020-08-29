// -*- C++ -*-

/*!
  \file array/IndexTypes.h
  \brief Types for array indexing.
*/

#if !defined(__array_IndexTypes_h__)
#define __array_IndexTypes_h__

#include <cstddef>

namespace array {

//! Types for array indexing.
class IndexTypes {
   //
   // Public types.
   //
public:

   // Types for STL compliance.

   //! The size type.
   typedef std::size_t size_type;
   //! Pointer difference type.
   typedef std::ptrdiff_t difference_type;

   // Other types.

   //! An array index is a signed integer.
   typedef difference_type Index;
};

} // namespace array

#endif
