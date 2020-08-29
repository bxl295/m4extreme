// -*- C++ -*-

/*!
  \file SimplexJac.h
  \brief Implements operations for the Jacobian matrix of a simplex.
*/

#if !defined(__geom_SimplexJac_h__)
#define __geom_SimplexJac_h__

#include "geometry.h"

#include "../../../ads/tensor/SquareMatrix.h"

#if _MSC_VER >= 1600
#include <tuple>
#else
#include <tr1/tuple>
#endif

namespace geom {

namespace {

const std::tr1::tuple < ads::SquareMatrix<1>, ads::SquareMatrix<2>,
      ads::SquareMatrix<3> >
      simplexJacIdentityToReference
      (ads::SquareMatrix<1>(1.0),
       ads::SquareMatrix<2>(1.0, - std::sqrt(3.0) / 3.0,
                            0.0, 2.0 * std::sqrt(3.0) / 3.0),
       ads::SquareMatrix<3>(1.0, - std::sqrt(3.0) / 3.0,   - std::sqrt(6.0) / 6.0,
                            0.0, 2.0 * std::sqrt(3.0) / 3.0, - std::sqrt(6.0) / 6.0,
                            0.0, 0.0,                      std::sqrt(6.0) / 2.0));

std::tr1::array<ads::SquareMatrix<1>, 1>
simplexJacGradientMatrix1 = {{ads::SquareMatrix<1>(1.0)}};

std::tr1::array<ads::SquareMatrix<2>, 2>
simplexJacGradientMatrix2 = {{
      ads::SquareMatrix<2>
      (-1.0, - std::sqrt(3.0) / 3.0,
      0.0, 0.0),
      ads::SquareMatrix<2>
      (0.0, 0.0,
      -1.0, - std::sqrt(3.0) / 3.0)
   }
};

std::tr1::array<ads::SquareMatrix<3>, 3>
simplexJacGradientMatrix3 = {{
      ads::SquareMatrix<3>
      (-1.0, - std::sqrt(3.0) / 3.0, - std::sqrt(6.0) / 6.0,
      0.0, 0.0, 0.0,
      0.0, 0.0, 0.0),
      ads::SquareMatrix<3>
      (0.0, 0.0, 0.0,
      -1.0, - std::sqrt(3.0) / 3.0, - std::sqrt(6.0) / 6.0,
      0.0, 0.0, 0.0),
      ads::SquareMatrix<3>
      (0.0, 0.0, 0.0,
      0.0, 0.0, 0.0,
      -1.0, - std::sqrt(3.0) / 3.0, - std::sqrt(6.0) / 6.0)
   }
};

const std::tr1::tuple < std::tr1::array<ads::SquareMatrix<1>, 1>,
      std::tr1::array<ads::SquareMatrix<2>, 2>,
      std::tr1::array<ads::SquareMatrix<3>, 3> >
      simplexJacGradientMatrix(simplexJacGradientMatrix1,
                               simplexJacGradientMatrix2,
                               simplexJacGradientMatrix3);

const std::tr1::array<double, 3>
simplexJacDeterminantIdentityToReference = {{1.0, 2.0 * std::sqrt(3.0) / 3.0, std::sqrt(2.0)}};

const std::tr1::array<double, 3>
simplexJacDimensionFactorial = {{1, 2, 6}};

}

//! Base class with const data.
template<std::size_t N>
class SimplexJacStaticData {
protected:

   //! An NxN matrix.
   typedef ads::SquareMatrix<N> Matrix;

   //! The Jacobian matrix that maps the identity simplex to the reference simplex.
   /*! In 2-D, it maps the identity triangle:
     (0,0), (1,0), (1/2,sqrt(3)/2)
     to the reference triangle:
     (0,0), (1,0), (0,1).
     In 3-D it maps the identity tetrahedron:
     (0,0,0), (1,0,0), (1/2,sqrt(3)/2,0), (1/2,sqrt(3)/6,sqrt(2)/sqrt(3))
     to the reference tetrahedron:
     (0,0,0), (1,0,0), (0,1,0), (0,0,1).
   */
   const Matrix&
   identityToReference() const {
      return std::tr1::get < N - 1 > (simplexJacIdentityToReference);
   }

   //! The gradient of the Jacobian matrix.
   const std::tr1::array<Matrix, N>&
   gradientMatrix() const {
      return std::tr1::get < N - 1 > (simplexJacGradientMatrix);
   }

   //! The determinant of identityToReference;
   double
   determinantIdentityToReference() const {
      return simplexJacDeterminantIdentityToReference[N-1];
   }

   //! The factorial of the dimension.
   double
   dimensionFactorial() const {
      return simplexJacDimensionFactorial[N-1];
   }
};

// CONTINUE: Get rid of T.
//! Implements operations for the Jacobian matrix of a simplex.
/*!
  \param N is the dimension.
  \param T is the number type.  By default it is double.

  <b>The Jacobian Matrix</b>

  Consider a simplex \f$T\f$ with vertices
  \f$ \{ \mathbf{x}_0, \ldots \mathbf{x}_{N-1} \} \f$.  We call this
  the physical simplex.
  The identity simplex \f$T_I\f$ may be mapped to the physical simplex
  \f$T\f$ by an affine transformation and a translation.
  \f[ T = S T_I + \mathbf{x}_0 \f]
  \f$S\f$ is the Jacobian matrix of the transformation.

  In 2-D, the identity triangle has vertices:
  \f$(0,0)\f$, \f$(1,0)\f$ and \f$(1/2,\sqrt{3}/2)\f$.
  In 3-D, the identity tetrahedron has vertices:
  \f$(0,0,0)\f$, \f$(1,0,0)\f$, \f$(1/2,\sqrt{3}/2,0)\f$ and
  \f$(1/2,\sqrt{3}/6,\sqrt{2} / \sqrt{3})\f$.

  The logical simplex, \f$T_L\f$ is the simplex whose vertices are the origin
  and unit displacements in each coordinate direction.
  In 2-D, the logical triangle has vertices:
  \f$(0,0)\f$, \f$(1,0)\f$ and \f$(0,1)\f$.
  In 3-D, the logical tetrahedron has vertices:
  \f$(0,0,0)\f$, \f$(1,0,0)\f$, \f$(0,1,0)\f$, \f$(0,0,1)\f$.
  It is easy to map the logical simplex to the physical simplex.
  \f[ T = A T_L + \mathbf{x}_0 \f]
  The columns of \f$A\f$ are the displacements of the vertices from
  the first vertex \f$ \mathbf{x}_0 \f$.  In 3-D, this is
  \f[
  A =
  \left(
  \begin{array}{ccc}
  x_1 - x_0 & x_2 - x_0 & x_3 - x_0 \\
  y_1 - y_0 & y_2 - y_0 & y_3 - y_0 \\
  z_1 - z_0 & z_2 - z_0 & z_3 - z_0
  \end{array}
  \right).
  \f]
  It is also easy to map the logical simplex to the identity simplex.
  \f[
  T_I = W T_L
  \f]
  The columns of \f$W\f$ are the vertices of \f$T_I\f$.  In 3-D, this is
  \f[
  W =
  \left(
  \begin{array}{ccc}
  1 & 1/2 & 1/2 \\
  0 & \sqrt{3} / 2 & \sqrt{3} / 6 \\
  0 & 0 & \sqrt{2} / \sqrt{3}
  \end{array}
  \right).
  \f]
  By combining the former transformation with the inverse of the latter,
  we can map the identity simplex to the physical simplex.
  \f[
  T = A W^{-1} T_I + \mathbf{x}_0
  \f]
  The Jacobian matrix of the transformation is \f$ S = A W^{-1} \f$.

  <b>Usage</b>

  Construct a \c SimplexJac with the default constructor or from a
  \c SimplexJac::Simplex.
  \code
  typedef geom::SimplexJac<3> TetJac;
  typedef TetJac::Vertex Vertex;
  typedef TetJac::Simplex Tetrahedron;
  // Default constructor.
  TetJac tet;
  \endcode
  \code
  // The identity tetrahedron.
  Tetrahedron t(Vertex(0, 0, 0),
                Vertex(1, 0, 0),
	        Vertex(1./2, std::sqrt(3.)/2, 0),
	        Vertex(1./2, std::sqrt(3.)/6, std::sqrt(2./3.)));
  TetJac tet(t);
  \endcode
  The \c Simplex constructor calls \c set() to enable evaluation of the
  determinant, the content and their gradients.

  To evaluate the determinant or the content of the simplex, first call
  \c setFunction() to set the Jacobian matrix and then use the getDeterminant()
  and computeContent() member functions.
  \code
  tet.setFunction(t);
  std::cout << "Identity tetrahedron:\n"
            << "determinant = " << tet.getDeterminant()
            << "\nvolume = " << tet.computeContent()
            << '\n';
  \endcode
  To evaluate the determinant and content and/or their gradients, first call
  \c set() to set the Jacobian matrix and its gradient.  Then use the
  member functions to access the appropriate qantities.
  \code
  tet.set(t);
  std::cout << "Identity tetrahedron:\n"
            << "determinant = " << tet.getDeterminant()
            << "\ngrad determinant = " << tet.getGradientDeterminant()
            << "\nvolume = " << tet.computeContent()
            << "\ngrad volume = " << tet.computeGradientContent()
            << '\n';
  \endcode
*/
template < std::size_t N, typename T = double >
class SimplexJac : public SimplexJacStaticData<N> {
private:

   typedef SimplexJacStaticData<N> Base;

public:

   //
   // Public types.
   //

   //! The number type.
   typedef T Number;

   //! The class for a vertex.
   typedef std::tr1::array<Number, N> Vertex;

   //! The simplex type.
   typedef std::tr1::array < Vertex, N + 1 > Simplex;

   //! An NxN matrix.
   typedef ads::SquareMatrix<N, Number> Matrix;

private:

   //
   // Member data.
   //

   // The Jacobian matrix that maps the identity simplex to this
   // simplex (after the first vertex has been translated to the origin.)
   Matrix _matrix;

   // The determinant of _matrix.
   Number _determinant;

   // The gradient of the determinant of the Jacobian matrix.
   Vertex _gradientDeterminant;

   using Base::identityToReference;
   using Base::gradientMatrix;
   using Base::determinantIdentityToReference;
   using Base::dimensionFactorial;

public:

   //--------------------------------------------------------------------------
   //! \name Constructors etc.
   //! @{

   //! Default constructor.  Un-initialized memory.
   SimplexJac() :
      _matrix(),
      _determinant(),
      _gradientDeterminant() {}

   //! Copy constructor.
   SimplexJac(const SimplexJac& other) :
      _matrix(other._matrix),
      _determinant(other._determinant),
      _gradientDeterminant(other._gradientDeterminant) {}


   //! Construct from a simplex.
   SimplexJac(const Simplex& s) {
      set(s);
   }

   //! Assignment operator.
   SimplexJac&
   operator=(const SimplexJac& other) {
      if (&other != this) {
         _matrix = other._matrix;
         _determinant = other._determinant;
         _gradientDeterminant = other._gradientDeterminant;
      }
      return *this;
   }

   //! Trivial destructor.
   ~SimplexJac() {}

   //! @}
   //--------------------------------------------------------------------------
   //! \name Accessors
   //! @{

   //! Return a const reference to the Jacobian matrix.
   const Matrix&
   getMatrix() const {
      return _matrix;
   }

   //! Return a const reference to the gradient of the Jacobian matrix.
   const std::tr1::array<Matrix, N>&
   getGradientMatrix() const {
      return gradientMatrix();
   }

   //! Return the determinant of the Jacobian matrix.
   Number
   getDeterminant() const {
      return _determinant;
   }

   //! Return a const reference to the gradient of the determinant of the Jacobian matrix.
   const Vertex&
   getGradientDeterminant() const {
      return _gradientDeterminant;
   }

   //! Return the content (hypervolume) of the simplex.
   Number
   computeContent() const {
      return _determinant / dimensionFactorial() /
      determinantIdentityToReference();
   }

   //! Calculate the gradient of the content (hypervolume) of the simplex.
   void
   computeGradientContent(Vertex* grad) const {
      *grad = _gradientDeterminant;
      *grad /= (dimensionFactorial() * determinantIdentityToReference());
   }

   //! Return the gradient of the content (hypervolume) of the simplex.
   Vertex
   computeGradientContent() const {
      Vertex grad;
      computeGradientContent(&grad);
      return grad;
   }

   //! @}
   //--------------------------------------------------------------------------
   //! \name Manipulators
   //! @{

   //! Set the vertices.  Calculate the Jacobian matrix and determinant.
   void
   setFunction(const Simplex& s);

   //! Set the vertices.  Calculate the Jacobian matrix and the determinant and its gradient.
   void
   set(const Simplex& s);

   // CONTINUE: These kind of functions do not know the orientation of the
   // simplex.  A 3-2 mesh could become tangled and still have high quality.
   //! Set the vertices.  Calculate the Jacobian matrix and determinant.
   /*!
     This first projects the simplex to N-D and then call the above
     setFunction().
   */
   void
   setFunction(const std::tr1::array < std::tr1::array < Number, N + 1 > , N + 1 > & s);

   //! Set the vertices.  Calculate the Jacobian matrix and the determinant and its gradient.
   /*!
     This first projects the simplex to N-D and then call the above set().
   */
   void
   set(const std::tr1::array < std::tr1::array < Number, N + 1 > , N + 1 > & s);

   //! @}

};

} // namespace geom

#define __geom_SimplexJac_ipp__
#include "SimplexJac.ipp"
#undef __geom_SimplexJac_ipp__

#endif
