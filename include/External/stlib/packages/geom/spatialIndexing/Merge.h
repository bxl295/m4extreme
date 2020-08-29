// -*- C++ -*-

/*!
  \file geom/spatialIndexing/Merge.h
  \brief Merging functors.
*/

#if !defined(__geom_spatialIndexing_Merge_h__)
#define __geom_spatialIndexing_Merge_h__

namespace geom {

//! Merging functor that does nothing.
struct MergeNull {
   //! Do nothing.
   void
   operator()() const {
   }
};

//! Merging functor that copies the value of the first child.
struct MergeCopyFirst {
   //! Copy the value of the first child.
   template<std::size_t _NumberOfOrthants, typename _Element>
   void
   operator()
   (const std::tr1::array<const _Element*, _NumberOfOrthants>& children,
    _Element* parent) const {
      *parent = *children[0];
   }
};

//! Merging functor that averages the values of the children.
struct MergeAverage {
   //! Average the children.
   template<std::size_t _NumberOfOrthants, typename _Element>
   void
   operator()
   (const std::tr1::array<const _Element*, _NumberOfOrthants>& children,
    _Element* parent) const {
      *parent = *children[0];
      for (std::size_t i = 0; i != _NumberOfOrthants; ++i) {
         *parent == *children[i];
      }
      *parent /= _NumberOfOrthants;
   }
};

} // namespace geom

#endif
