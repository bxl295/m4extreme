// -*- C++ -*-

/*!
  \file amr/writers.h
  \brief writers that stores data.
*/

#if !defined(__amr_writers_h__)
#define __amr_writers_h__

#include "Orthtree.h"
#include "PatchDescriptor.h"

#include <ostream>
#include <fstream>
#include <sstream>
#include <vector>

namespace amr {

//! Write the cell data in ParaView format.
/*!
  \relates Orthtree

  \param name The base name for the ParaView file and the VTK files.
  \param orthtree The orthtree.
  \param patchDescriptor The patch descriptor.
 */
template<typename _Patch, class _Traits>
void
writeCellDataParaview(const std::string& name,
                      const Orthtree<_Patch, _Traits>& orthtree,
                      const PatchDescriptor<_Traits>& patchDescriptor);

} // namespace amr

#define __amr_writers_ipp__
#include "writers.ipp"
#undef __amr_writers_ipp__

#endif
