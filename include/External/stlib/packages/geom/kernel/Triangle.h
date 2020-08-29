// -*- C++ -*-

/*!
  \file Triangle.h
  \brief Implements a class for a triangle.
*/

#if !defined(__geom_Triangle_h__)
#define __geom_Triangle_h__

#include "Plane.h"

namespace geom {

//! A class for a triangle in N dimensions.
/*!
  \param N is the dimension.
  \param T is the number type.  By default it is double.
*/
template < std::size_t N, typename T = double >
class Triangle {
   //
   // Public types.
   //

public:

   //! The representation for a point.
   typedef std::tr1::array<T, N> Point;

   //
   // Data
   //

private:

   std::tr1::array<Point, 3> _vertices;

public:

   //--------------------------------------------------------------------------
   //! \name Constructors etc.
   // @{

   //! Default constructor.  Uninitialized memory.
   Triangle() {}

   //! Construct from three points.
   Triangle(const Point& a, const Point& b, const Point& c) {
      _vertices[0] = a;
      _vertices[1] = b;
      _vertices[2] = c;
   }

   //! Copy constructor.
   Triangle(const Triangle& other) {
      _vertices[0] = other._vertices[0];
      _vertices[1] = other._vertices[1];
      _vertices[2] = other._vertices[2];
   }

   //! Assignment operator.
   const Triangle&
   operator=(const Triangle& other);

   //! Trivial destructor.
   ~Triangle() {}

   // @}
   //--------------------------------------------------------------------------
   //! \name Accesors.
   // @{

   //! Get the specified vertex.
   const Point&
   getVertex(const std::size_t n) const {
      return _vertices[n];
   }

   // @}
   //--------------------------------------------------------------------------
   //! \name Manipulators.
   // @{

   //! Set the specified vertex.
   void
   setVertex(const std::size_t n, const Point& x) const {
      _vertices[n] = x;
   }

   // @}
   //--------------------------------------------------------------------------
   //! \name Translations.
   // @{

   //! Translate the triangle by +p.
   Triangle&
   operator+=(const Point& p);

   //! Translate the triangle by -p.
   Triangle&
   operator-=(const Point& p);

   // @}
};

//
// Mathematical Functions.
//


//! Return the supporting plane of a triangle in 3 dimensions.
/*! \relates Triangle */
template<typename T>
Plane<T>
buildSupportingPlane(const Triangle<3, T>& triangle);


//
// File I/O operators.
//


//! Read the three vertices.
/*! \relates Triangle */
template<std::size_t N, typename T>
std::istream&
operator>>(std::istream& in, Triangle<N, T>& t);


//! Write the three vertices.
/*! \relates Triangle */
template<std::size_t N, typename T>
std::ostream&
operator<<(std::ostream& out, const Triangle<N, T>& t);


} // namespace geom

#define __geom_Triangle_ipp__
#include "Triangle.ipp"
#undef __geom_Triangle_ipp__

#endif
