// -*- C++ -*-

#if !defined(__levelSet_Patch_h__)
#define __levelSet_Patch_h__

#include "count.h"
#include "../array/EquilateralArrayRef.h"

namespace levelSet {

//! A patch is a multi-array that has equal extents in each dimension.
/*!
  \param _T The value type.
  \param _D The dimension
  \param _N The extent in each dimension.
*/
template<typename _T, std::size_t _D, std::size_t _N>
class Patch :
      public array::EquilateralArrayRef<_T, _D, _N> {
   //
   // Types.
   //
private:

   //! The base class.
   typedef array::EquilateralArrayRef<_T, _D, _N> Base;

   //
   // Member data.
   //
public:

   typename Base::value_type fillValue;

   //--------------------------------------------------------------------------
   /*! \name Constructors etc.
     We use the default copy constructor, assignment operator, and
     destructor. */
   // @{
public:

   //! Default constructor. Unrefined with NaN for the fill value.
   Patch() :
      Base(),
      fillValue(std::numeric_limits<typename Base::value_type>::quiet_NaN()) {
   }

   //! Construct from a pointer to the array data.
   Patch(typename Base::iterator data) :
      Base(data),
      fillValue(0) {
   }

   // @}
   //--------------------------------------------------------------------------
   //! \name Accessors.
   // @{
public:

   bool
   isValid() const {
      // Either the data pointer must be valid, or there must be a valid
      // (non-zero) fill value.
      return Base::data() == 0 && fillValue != 0 ||
         Base::data() != 0 && fillValue == 0;
   }

   //! Return true if the patch is refined.
   bool
   isRefined() const {
      return Base::data();
   }

   //! Return true if the patch should be coarsened.
   bool
   shouldBeCoarsened() const {
      return isRefined() && allSame(Base::begin(), Base::end());
   }

   // @}
   //--------------------------------------------------------------------------
   //! \name Manipulators.
   // @{
public:

   //! Reference an externally allocated array.
   void
   refine(const typename Base::iterator data) {
      Base::setData(data);
      fillValue = 0;
   }

   //! Coarsen the patch. Set the fill value to the first grid value.
   /*! \pre The patch must be refined. */
   void
   coarsen() {
      assert(Base::data());
      fillValue = *Base::begin();
      Base::setData(0);
   }

   //! Clear the patch. Set the data to 0 and the fill value to NaN.
   void
   clear() {
      Base::setData(0);
      fillValue = std::numeric_limits<typename Base::value_type>::quiet_NaN();
   }

private:

   // Hide setData().
   void
   setData();

   // @}
};


} // namespace levelSet

#endif
