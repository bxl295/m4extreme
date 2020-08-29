// -*- C++ -*-

#if !defined(__numerical_constants_Exponentiation_ipp__)
#error This file is an implementation detail of Exponentiation.
#endif

namespace numerical {


//! 0^0 is indeterminate, so we do not provide a result.
/*!
  Trying to compute 0^0 will cause a compilation error.
*/
template<>
struct Exponentiation<int, 0, 0> {
};

// We must specialize for each integer type because a specialization may not
// depend on a template parameter.

//! 0^E = 0 for nonzero E.
template<char Exponent>
struct Exponentiation<char, 0, Exponent> {
   static const char Result = 0;
};
template<signed char Exponent>
struct Exponentiation<signed char, 0, Exponent> {
   static const signed char Result = 0;
};
template<unsigned char Exponent>
struct Exponentiation<unsigned char, 0, Exponent> {
   static const unsigned char Result = 0;
};
template<short Exponent>
struct Exponentiation<short, 0, Exponent> {
   static const short Result = 0;
};
template<unsigned short Exponent>
struct Exponentiation<unsigned short, 0, Exponent> {
   static const unsigned short Result = 0;
};
template<int Exponent>
struct Exponentiation<int, 0, Exponent> {
   static const int Result = 0;
};
template<unsigned Exponent>
struct Exponentiation<unsigned, 0, Exponent> {
   static const unsigned Result = 0;
};
template<long Exponent>
struct Exponentiation<long, 0, Exponent> {
   static const long Result = 0;
};
template<unsigned long Exponent>
struct Exponentiation<unsigned long, 0, Exponent> {
   static const unsigned long Result = 0;
};


//! B^0 = 1 for nonzero B.
template<char Base>
struct Exponentiation<char, Base, 0> {
   static const char Result = 1;
};
template<signed char Base>
struct Exponentiation<signed char, Base, 0> {
   static const signed char Result = 1;
};
template<unsigned char Base>
struct Exponentiation<unsigned char, Base, 0> {
   static const unsigned char Result = 1;
};
template<short Base>
struct Exponentiation<short, Base, 0> {
   static const short Result = 1;
};
template<unsigned short Base>
struct Exponentiation<unsigned short, Base, 0> {
   static const unsigned short Result = 1;
};
template<int Base>
struct Exponentiation<int, Base, 0> {
   static const int Result = 1;
};
template<unsigned Base>
struct Exponentiation<unsigned, Base, 0> {
   static const unsigned Result = 1;
};
template<long Base>
struct Exponentiation<long, Base, 0> {
   static const long Result = 1;
};
template<unsigned long Base>
struct Exponentiation<unsigned long, Base, 0> {
   static const unsigned long Result = 1;
};


} // namespace numerical
