// -*- C++ -*-

/*!
  \file array/MultiArrayStorage.h
  \brief
*/

#if !defined(__array_MultiArrayStorage_h__)
#define __array_MultiArrayStorage_h__

#include "../ext/array.h"

#ifdef DEBUG_stlib
#include <algorithm>
#endif

namespace array {

//! Class to indicate row-major storage in a multidimensional %array such as MultiArray.
struct RowMajor {
};

//! Class to indicate column-major storage in a multidimensional %array such as MultiArray.
struct ColumnMajor {
};

//! MultiArrayStorage orders for the %array.
/*!
  The storage order is a permutation of 0 .. _Dimension - 1. In 3-D, {2, 1, 0}
  follows the C convention and {0, 1, 2} follows the fortran convention.
*/
template<std::size_t _Dimension>
class
      MultiArrayStorage : public std::tr1::array<std::size_t, _Dimension> {
   //
   // Constants.
   //
public:

   //! The number of dimensions.
   static const std::size_t Dimension = _Dimension;

   //
   // Types.
   //
private:

   typedef std::tr1::array<std::size_t, _Dimension> Base;

   //
   // Using.
   //
public:

   using Base::begin;
   using Base::end;

   //--------------------------------------------------------------------------
   //! \name Constructors etc.
   //@{
public:

   // The default copy constructor, assignment operator, and destructor are fine.

   //! Default constructor. Uninitialized memory.
   MultiArrayStorage() :
      Base() {
   }

   //! Row-major storage. The last coordinate varies fastest.
   MultiArrayStorage(RowMajor /*dummy*/) {
      for (std::size_t i = 0; i != Dimension; ++i) {
         (*this)[i] = Dimension - 1 - i;
      }
   }

   //! Column-major storage. The first coordinate varies fastest.
   MultiArrayStorage(ColumnMajor /*dummy*/) {
      for (std::size_t i = 0; i != Dimension; ++i) {
         (*this)[i] = i;
      }
   }

   //! Construct from the storage order %array.
   MultiArrayStorage(const std::tr1::array<std::size_t, _Dimension>& storage) :
      Base(storage) {
#ifdef DEBUG_stlib
      // Check that the storage order is a permutation of 0 .. Dimension - 1.
      Base tmp = *this;
      std::sort(tmp.begin(), tmp.end());
      for (std::size_t i = 0; i != Dimension; ++i) {
         assert(tmp[i] == i);
      }
#endif
   }

   //@}
};

} // namespace array

#endif
