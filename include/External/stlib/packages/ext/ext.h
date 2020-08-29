// -*- C++ -*-

#if !defined(__ext_h__)
#define __ext_h__

#include "ext/array.h"
#include "ext/pair.h"
#include "ext/vector.h"

//! The namespace for the C++ Standard Library.
namespace std {
//! The namespace for the C++ Standard Library extensions in Technical Report 1.
namespace tr1 {
}
}

//! Some functions that extend the capabilities of the C++ standard library are defined in this namespace.
namespace ext {
}

/*!
\mainpage Standard Library Extensions

\section extIntroduction Introduction

This package provides extensions to the following standard library
classes:
- \ref extArray
- \ref extPair
- \ref extVector

For the most part functions are injected into the std or std::tr1 namespace.
There is no need to specify the namespace when using the functions.
The function will be resolved with
<a href="http://en.wikipedia.org/wiki/Argument_dependent_name_lookup">
argument dependent name lookup</a>. In the example below we find the minimum
element in a std::vector and calculate the dot product of two
std::vector's.
\verbatim
std::vector<double> x, y;
...
const double minValue = min(x);
const double d = dot(x, y); \endverbatim
This is more concise than the following code.
\verbatim
const double minValue = *std::min_element(x.begin(), x.end())
double d = 0;
for (std::size_t i = 0; i != x.size(); ++i) {
   d += x[i] * y[i];
} \endverbatim
*/

/*!
\page extBibliography Bibliography

-# \anchor extAustern1999
Matthew H. Austern. "Generic Programming and the STL: Using and Extending
the C++ Standard Template Library." Addison-Wesley.
-# \anchor extBecker2007
Pete Becker. "The C++ Standard Library Extensions." Addison-Wesley.
*/

#endif
