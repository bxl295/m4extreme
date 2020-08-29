// -*- C++ -*-

#if !defined(__geom_Segment_ipp__)
#error This file is an implementation detail of the class Segment.
#endif

namespace geom {

//
// Constructors
//


template<std::size_t N, typename T>
inline
Segment<N, T>&
Segment<N, T>::
operator=(const Segment& other) {
   // Avoid assignment to self
   if (&other != this) {
      _source = other._source;
      _target = other._target;
   }
   // Return *this so assignments can chain
   return *this;
}


//
// Arithmetic member operators
//


template<std::size_t N, typename T>
inline
Segment<N, T>&
Segment<N, T>::
operator+=(const Point& p) {
   _source += p;
   _target += p;
   return (*this);
}


template<std::size_t N, typename T>
inline
Segment<N, T>&
Segment<N, T>::
operator-=(const Point& p) {
   _source -= p;
   _target -= p;
   return (*this);
}


} // namespace geom
