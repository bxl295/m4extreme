// -*- C++ -*-

#if !defined(__geom_Line_2_ipp__)
#error This file is an implementation detail of the class Line_2.
#endif

namespace geom {

//
// Assignment operator
//


template<typename T>
inline
Line_2<T>&
Line_2<T>::operator=(const Line_2& other) {
   // Avoid assignment to self
   if (&other != this) {
      _segment = other._segment;
      _normal = other._normal;
   }
   // Return *this so assignments can chain
   return *this;
}


//
// Mathematical functions
//


template<typename T>
inline
void
Line_2<T>::
computeIntersection(Point q1, Point q2, Point* intersectionPoint) const {
#ifdef DEBUG_stlib
   assert(computeSignedDistance(q1) * computeSignedDistance(q2) <= 0);
#endif
   q1 -= getPointOn();
   q2 -= getPointOn();

   const Number p1 = dot(q1, getTangent());
   const Number p2 = dot(q2, getTangent());
   const Number h1 = dot(q1, getNormal());
   const Number h2 = dot(q2, getNormal());

   *intersectionPoint = getPointOn() +
                        ((p1 * h2 - p2 * h1) / (h2 - h1)) * getTangent();
}

} // namespace geom
