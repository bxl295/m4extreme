// -*- C++ -*-

#if !defined(__levelSet_MolecularSurfaces_h__)
#define __levelSet_MolecularSurfaces_h__

#include "GridVirtual.h"
#include "marchingSimplices.h"
#include "powerDistance.h"
#include "negativeDistance.h"
#include "solventExcluded.h"

namespace levelSet {


//! Compute the surfaces for a molecule.
template<typename _T, std::size_t _D, std::size_t _N>
class MolecularSurfaces {
   //
   // Types.
   //
public:
   //! A Cartesian point.
   typedef std::tr1::array<_T, _D> Point;

protected:
   //! A bounding box.
   typedef geom::BBox<_T, _D> BBox;
   //! The virtual grid.
   typedef GridVirtual<_T, _D, _N> Grid;
   //! A single index.
   typedef typename Grid::Index Index;
   //! A list of indices.
   typedef typename Grid::IndexList IndexList;
   //! A patch in the virtual grid.
   typedef typename Grid::Patch Patch;

   //
   // Member data.
   //
public:

   //! The list of atoms.
   const std::vector<geom::Ball<_T, _D> > atoms;
   //! The probe radius.
   const _T probeRadius;
   //! The target grid spacing.
   const _T targetGridSpacing;

   //--------------------------------------------------------------------------
   /*! \name Constructors etc.
     The copy constructor and assignment operator are disabled as there is
     constant member data. We use the synthesized destructor.
   */
   // @{
public:

   //! Construct from the list of atoms, the probe radius, and the target grid spacing.
   /*! 
    \pre The list of atoms must be non-empty. The probe radius and the target
    grid spacing must be positive. */
   MolecularSurfaces(const std::vector<geom::Ball<_T, _D> >& atoms_,
                     const _T probeRadius_, const _T targetGridSpacing_) :
      atoms(atoms_),
      probeRadius(probeRadius_),
      targetGridSpacing(targetGridSpacing_) {
      // Check the input.
      assert(! atoms.empty());
      assert(probeRadius > 0);
      assert(targetGridSpacing > 0);
   }

private:

   //! The copy constructor is not implemented.
   MolecularSurfaces(const MolecularSurfaces&);

   //! The assignment operator is not implemented.
   MolecularSurfaces& 
   operator=(MolecularSurfaces&);

   // @}
   //--------------------------------------------------------------------------
   //! \name Calculate the surfaces.
   // @{
public:

   //! Calculate the the van der Waals surface.
   /*! Report the vertices of the boundary tesselation. */
   template<typename _OutputIterator>
   std::pair<_T, _T>
   vanDerWaals(_OutputIterator vertices) const {
      return unionOfBalls(atoms, vertices);
   }

   //! Calculate the the solvent accessible surface.
   /*! Report the vertices of the boundary tesselation. */
   template<typename _OutputIterator>
   std::pair<_T, _T>
   solventAccessible(_OutputIterator vertices) const;

   //! Calculate the the solvent excluded surface.
   /*! Report the vertices of the boundary tesselation. */
   template<typename _OutputIterator>
   std::pair<_T, _T>
   solventExcluded(_OutputIterator vertices) const;

   //! Calculate the the solventd-excluded cavities.
   /*! Report the vertices of the boundary tesselation. */
   template<typename _OutputIterator>
   std::pair<_T, _T>
   solventExcludedCavities(_OutputIterator vertices) const;

protected:

   //! Calculate the surface of the union of the balls.
   /*! Report the vertices of the boundary tesselation. */
   template<typename _OutputIterator>
   std::pair<_T, _T>
   unionOfBalls(const std::vector<geom::Ball<_T, _D> >& balls,
                _OutputIterator vertices) const;

   //! Calculate the patch dependencies for computing the negative distance.
   void
   negativeDistanceDependencies
   (const Grid& grid, const std::vector<geom::Ball<_T, _D> >& balls,
    array::StaticArrayOfArrays<std::size_t>* dependencies) const;

   //! Calculate an appropriate domain for the balls.
   /*! Bound the balls and then expand by the target grid spacing. */
   BBox
   calculateDomain(const std::vector<geom::Ball<_T, _D> >& balls) const;

   // @}
};


} // namespace levelSet

#define __levelSet_MolecularSurfaces_ipp__
#include "MolecularSurfaces.ipp"
#undef __levelSet_MolecularSurfaces_ipp__

#endif
