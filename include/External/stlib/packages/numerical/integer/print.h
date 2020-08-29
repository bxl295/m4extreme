// -*- C++ -*-

/*!
  \file numerical/integer/print.h
  \brief Select an integer type with at least the specified number of print.
*/

#if !defined(__numerical_integer_print_h__)
#define __numerical_integer_print_h__

#include "../../loki/TypeManip.h"

#include <limits>
#include <iostream>

#include <cassert>

namespace numerical {

//-----------------------------------------------------------------------------
//! \defgroup numerical_integer_print Print the bits of an integer.
//@{

//! Print the bits of the integer.
template<typename _Integer>
void
printBits(std::ostream& out, _Integer x);


//! Print the specified bits of the integer.
template<typename _Integer>
void
printBits(std::ostream& out, _Integer x, std::size_t indexBeginning,
          std::size_t indexEnd);

//\@}

} // namespace numerical

#define __numerical_integer_print_ipp__
#include "print.ipp"
#undef __numerical_integer_print_ipp__

#endif
