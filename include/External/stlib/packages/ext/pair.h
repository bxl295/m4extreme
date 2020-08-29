// -*- C++ -*-

/*!
  \file ext/pair.h
  \brief Functions for pair.
*/

#if !defined(__ext_pair_h__)
#define __ext_pair_h__

#include <utility>
#include <iostream>

/*!
\page extPair Extensions to std::pair

Here we functions for \ref extPairFile "file I/O" to extend the functionality
of the std::pair struct [\ref extAustern1999 "Austern, 1999"].
*/

namespace std {

//----------------------------------------------------------------------------
//! \defgroup extPairFile Pair File I/O
//@{

//! Write the first and second element, separated by a space.
/*!
  Format:
  x.first x.second
*/
template<typename _T1, typename _T2>
inline
ostream&
operator<<(ostream& out, const pair<_T1, _T2>& x) {
   out << x.first << ' ' << x.second;
   return out;
}

//! Read the first and second element.
template<typename _T1, typename _T2>
inline
istream&
operator>>(istream& in, pair<_T1, _T2>& x) {
   in >> x.first >> x.second;
   return in;
}

//@}

}

#endif
