// -*- C++ -*-

#if !defined(__hj_debug_h__)
#define __hj_debug_h__

#include "../array/IndexTypes.h"

#include <tr1/array>

namespace hj {

namespace debug {

typedef array::IndexTypes::Index Index;

//
// 2-D.
//

// CONTINUE int to ptrdiff_t. Use std::tr1::array
std::tr1::array<std::tr1::array<Index, 2>, 4> _adj2 = {{
   {{ 1,  0}},
   {{-1,  0}},
   {{ 0,  1}},
   {{ 0, -1}}
   }};

std::tr1::array<std::tr1::array<Index, 2>, 4> _diag2 = {{
   {{ 1,  1}},
   {{ 1, -1}},
   {{-1,  1}},
   {{-1, -1}}
   }};

//! Return true if the difference is in an adjacent direction.
inline
bool
is_adjacent(const std::tr1::array<Index, 2>& di) {
   for (std::size_t i = 0; i != _adj2.size(); ++i) {
      if (di == _adj2[i]) {
         return true;
      }
   }
   return false;
}

//! Return true if the difference is in a diagonal direction.
inline
bool
is_diagonal(const std::tr1::array<Index, 2>& di) {
   for (std::size_t i = 0; i != _diag2.size(); ++i) {
      if (di == _diag2[i]) {
         return true;
      }
   }
   return false;
}

//
// 3-D.
//

std::tr1::array<std::tr1::array<Index, 3>, 6> _adj3 = {{
   {{ 1,  0,  0}},
   {{-1,  0,  0}},
   {{ 0,  1,  0}},
   {{ 0, -1,  0}},
   {{ 0,  0,  1}},
   {{ 0,  0, -1}}
   }};

std::tr1::array<std::tr1::array<Index, 3>, 12> _diag3 = {{
   {{ 1,  1,  0}},
   {{ 1, -1,  0}},
   {{-1,  1,  0}},
   {{-1, -1,  0}},
   {{ 0,  1,  1}},
   {{ 0,  1, -1}},
   {{ 0, -1,  1}},
   {{ 0, -1, -1}},
   {{ 1,  0,  1}},
   {{-1,  0,  1}},
   {{ 1,  0, -1}},
   {{-1,  0, -1}}
   }};

//! Return true if the difference is in an adjacent direction.
inline
bool
is_adjacent(const std::tr1::array<Index, 3>& di) {
   for (std::size_t i = 0; i != _adj3.size(); ++i) {
      if (di == _adj3[i]) {
         return true;
      }
   }
   return false;
}

//! Return true if the difference is in a diagonal direction.
inline
bool
is_diagonal(const std::tr1::array<Index, 3>& di) {
   for (std::size_t i = 0; i != _diag3.size(); ++i) {
      if (di == _diag3[i]) {
         return true;
      }
   }
   return false;
}

} // namespace debug

} // namespace hj

#endif
