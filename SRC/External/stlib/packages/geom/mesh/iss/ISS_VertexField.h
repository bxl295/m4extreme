// -*- C++ -*-

/*!
  \file ISS_VertexField.h
  \brief Interpolation on an indexed simplex set.
*/

#if !defined(__geom_ISS_VertexField_h__)
#define __geom_ISS_VertexField_h__

// Interpolation functions.
#include "../../../numerical/interpolation/simplex.h"

#if _MSC_VER >= 1600
#include <array>
#else
#include <tr1/array>
#endif

namespace geom {

//! Indexed simplex set with fields at the vertices.
/*!
  \param ISS is the indexed simplex set.
  \param F is the field type.  By default it is the number type of the mesh.

  This class stores a reference to an indexed simplex set and a pointer
  to the field data.
*/
template < class ISS,
         typename F = typename ISS::Number >
class ISS_VertexField {
   //
   // Private types.
   //

private:

   //! The indexed simplex set.
   typedef ISS IssType;
   //! The size type.
   typedef typename IssType::SizeType SizeType;

   //
   // Public types.
   //

public:

   //! The number type.
   typedef typename IssType::Number Number;
   //! A vertex.
   typedef typename IssType::Vertex Vertex;

   //
   // Field types.
   //

   //! The field type.
   typedef F Field;

   //
   // Data.
   //

private:

   //! The indexed simplex set.
   const IssType& _iss;
   //! The field values at the vertices.
   const Field* _fields;
   //! The simplex of positions.
   mutable std::tr1::array < Vertex, ISS::M + 1 > _pos;
   //! The simplex of vertex field values.
   mutable std::tr1::array < Field, ISS::M + 1 > _val;

   //
   // Not implemented.
   //

private:

   //! Default constructor not implemented.
   ISS_VertexField();

   //! Assignment operator not implemented.
   ISS_VertexField&
   operator=(const ISS_VertexField& x);

public:

   //--------------------------------------------------------------------------
   //! \name Constructors etc.
   //! @{

   //! Construct from the indexed simplex set and the fields.
   /*!
     \param iss is the indexed simplex set.
     \param fields is the array of fields.
   */
   ISS_VertexField(const IssType& iss, const Field* fields) :
      _iss(iss),
      _fields(fields),
      _pos(),
      _val() {
   }

   //! Copy constructor.
   ISS_VertexField(const ISS_VertexField& other) :
      _iss(other._iss),
      _fields(other._fields),
      _pos(),
      _val() {}

   //! Destructor has no effect on the indexed simplex set.
   ~ISS_VertexField() {}

   //! @}
   //--------------------------------------------------------------------------
   //! \name Mathematical Functions.
   //! @{

   //! Return the interpolated field for the n_th simplex and the point \c x.
   const Field&
   interpolate(std::size_t n, const Vertex& x) const;

   //! @}
};

} // namespace geom

#define __geom_ISS_VertexField_ipp__
#include "ISS_VertexField.ipp"
#undef __geom_ISS_VertexField_ipp__

#endif
