// -*- C++ -*-

#if !defined(__levelSet_h__)
#define __levelSet_h__

#include "solventExcluded.h"

//! Classes and functions for level set calculations.
namespace levelSet {
}

/*!
\mainpage Level Sets

\section levelSetIntroduction Introduction

\par
Consult \ref levelSetSethian1999 "Level Set Methods and Fast Marching Methods" 
and \ref levelSetOsher2003 "Level Set Methods and Dynamic Implicit Surfaces."
for an overview of level set methods.

\par Solvent-Accessible Manifold.
Each atom in a molecule may be modeled as a ball. For each, the radius is 
determined by the
<a href="http://en.wikipedia.org/wiki/Van_der_Waals_radius">van der Waals
volume</a>. The <em>van der Waals manifold</em> is the union of the balls.
The solvent water molecule is modeled with a single ball. The radius is
typically 1.4 Angstroms. Consider a solvent ball rolling over the surface
of the van der Waals manifold. The 
<a href="http://en.wikipedia.org/wiki/Accessible_surface_area">
solvent-accessible manifold</a> is bounded by the locus points of the center 
of the solvent. Equivalently, it is the union of the molecular balls where
each has been expanded by the solvent radius.

\par Solvent-Excluded Manifold.
The solvent ball is in a valid position if it does not intersect any of the
molecular balls. The solvent-excluded manifold is the negation of the union
of all valid solvent positions. Equivalently, one may construct it by
offsetting the surface of the solvent-accessible manifold by the negative
of the solvent radius. (Of course, when doing this one must account for 
self-intersections and take the entropy-satisfying solution.)

\par Cavities.
The <em>solvent-excluded cavities</em> are the difference between the 
solvent-excluded manifold and the van der Waals manifold. These are gaps
between the atoms where the solvent cannot reach. Consider a solvent that 
is restricted to valid (non-intersecting) positions. It is <em>trapped</em> if 
there is no continuous path to a position outside the convex hull of the 
van der Waals manifold. The union of all trapped solvents forms the 
<em>solvent-accessible cavities</em>. Roughly speaking, these are the cavities
where solvent molecules may be trapped.

\par Classes.
The levelSet::GridUniform class is a grid with uniform spacing.
The levelSet::Grid class uses a single level of adaptive mesh refinement.
The patches are elements of the levelSet::Patch class. The 
levelSet::GridVirtual class has information about patches, but does not
store them.

\par Functions.
There are functions for performing the following calculations.
- \ref levelSetBoolean "Boolean" operations.
- \ref levelSetCount "Count" the number of known or unknown values in a grid.
- \ref levelSetFlood "Flood fill" operations to set unknown values to a
signed constant.
- Level sets related to the \ref levelSetSolventExcluded "solvent excluded surface".
- Level sets for the \ref levelSetPositiveDistance "positive"
and \ref levelSetNegativeDistance "negative" distance to the union of a
set of balls.

\par Special values.
Use NaN (Not a Number, with value \c std::numeric_limits<T>::quiet_NaN())
to denote an unknown value. Use <tt>x != x</tt> to test whether the
variable \c x is NaN. This curious predicate arises from the axiom that NaN is 
not equal to any other number, including itself. Use \f$\infty\f$,
with value \c std::numeric_limits<T>::infinity(), to 
denote large unknown positive distances. Use \f$-\infty\f$ for large
unknown negative distances.

\section levelSetRelated Related Work

\par
\ref levelSetCan2006 "Efficient molecular surface generation using level-set methods."
The authors use level-set methods to calculate molecular surfaces and cavities.
To calculate the SES, they start by seeding the distance at the atom centers.
The distance is propogated outward using the Fast Marching Method. Instead of
using finite differencing, they propagate the index of the closest atom center
and perform Euclidean distance calculations. This yields an approximation of
the true distance. The distance is calculated back inwards the distance of
the probe radius to obtain the SES. The level-set methods are carried out
on a dense, uniform grid. This, along with a rough characterization of
the surface, limits the accuracy of the method.

\par
\ref levelSetDias2010 "CUDA-based Triangulations of Convolution Molecular Surfaces."
The authors use the GPU to compule Blinn surfaces for molecules. The level set
for the surface is computed on the GPU with one thead per atom. Then they
use a Marching Cubes algorithm with linear interpolation on the GPU to 
triangulate the surface.
*/

/*!
\page levelSetBibliography Bibliography

<!--------------------------------------------------------------------------->
-# \anchor levelSetSethian1999
J. A. Sethian. "Level Set Methods and Fast Marching Methods" Cambridge.
<!--------------------------------------------------------------------------->
-# \anchor levelSetOsher2003
Stanley Osher, Ronald Fedkiw. "Level Set Methods and Dynamic Implicit
Surfaces." Springer.
<!--------------------------------------------------------------------------->
-# \anchor levelSetCan2006
Tolga Can, Chao-I Chen, and Yuan-Fang Wang. 
"Efficient molecular surface generation using level-set methods."
Journal of Molecular Graphics & Modelling, Vol. 25, Issue 4, Pages 442-454,
Dec. 2006.
<!--------------------------------------------------------------------------->
-# \anchor levelSetDias2010
Sergio Dias, Kuldeep Bora, and Abel Gomes.
"CUDA-based Triangulations of Convolution Molecular Surfaces."
HPDC'10, June 20–25, 2010, Chicago, Illinois, USA.
*/
#endif
