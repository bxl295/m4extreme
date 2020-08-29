// -*- C++ -*-

#if !defined(__geom_Point_ipp__)
#error This file is an implementation detail of Point.
#endif

#include <cassert>
#include <cmath>

namespace geom {

//
// Math operators.
//


// Compute an orthogonal vector.
template<typename T>
inline
void
computeAnOrthogonalVector(std::tr1::array<T, 3> vector,
                          std::tr1::array<T, 3>* orthogonal) {
   normalize(&vector);

   // One of (1,0,0) and (0,1,0) is independent to the vector.
   std::tr1::array<T, 3> x = {{1, 0, 0}};
   const T d0 = dot(vector, x);
   x[0] = 0;
   x[1] = 1;
   const T d1 = dot(vector, x);
   // If (1,0,0) is a better independent direction than (0,1,0).
   if (std::abs(d0) < std::abs(d1)) {
      (*orthogonal)[0] = 1;
      (*orthogonal)[1] = 0;
      (*orthogonal)[2] = 0;
      vector *= d0;
   }
   // Otherwise (0,1,0) is a better independent direction than (1,0,0).
   else {
      (*orthogonal)[0] = 0;
      (*orthogonal)[1] = 1;
      (*orthogonal)[2] = 0;
      vector *= d1;
   }
   *orthogonal -= vector;
}


//
// Angles between points
//

template<typename T>
inline
int
computeSignOfTurn(const std::tr1::array<T, 2>& p,
const std::tr1::array<T, 2>& q,
const std::tr1::array<T, 2>& r) {
   const T disc = discriminant(q - p, r - q);
   if (disc > 0) {
      return -1;
   }
   else if (disc < 0) {
      return 1;
   }
   return 0;
}


template<typename T>
inline
int
computeApproximateSignOfTurn(const std::tr1::array<T, 2>& p,
const std::tr1::array<T, 2>& q,
const std::tr1::array<T, 2>& r) {
   std::tr1::array<T, 2> u = q - p;
   std::tr1::array<T, 2> v = r - q;
   normalize(&u);
   normalize(&v);
   const T disc = discriminant(u, v);
   if (disc > 1e-5) {
      return -1;
   }
   else if (disc < -1e-5) {
      return 1;
   }
   return 0;
}


template<typename T>
inline
T
computePseudoAngle(const std::tr1::array<T, 2>& vec) {
   const T dxpdy = std::abs(vec[0]) + std::abs(vec[1]);
   T theta = (dxpdy == 0) ? 0 : vec[1] / dxpdy;
   if (vec[0] < 0) {
      theta = 2 - theta;
   }
   else if (vec[1] < 0) {
      theta = 4 + theta;
   }
   return theta;
}


template<std::size_t N, typename T>
inline
T
computeAngle(const std::tr1::array<T, N>& a, const std::tr1::array<T, N>& b) {
   const T denominator = magnitude(a) * magnitude(b);
#ifdef DEBUG_stlib
   assert(denominator != 0);
#endif
   return std::acos(dot(a, b) / denominator);
}


//
// Rotations
//

template<typename T>
inline
void
rotatePiOver2(std::tr1::array<T, 2>* p) {
   const T temp = (*p)[0];
   (*p)[0] = -(*p)[1];
   (*p)[1] = temp;
}


template<typename T>
inline
void
rotateMinusPiOver2(std::tr1::array<T, 2>* p) {
   const T temp = (*p)[0];
   (*p)[0] = (*p)[1];
   (*p)[1] = -temp;
}

} // namespace geom
