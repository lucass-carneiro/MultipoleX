#ifndef MULTIPOLEX_TYPE_ALIASES_HPP
#define MULTIPOLEX_TYPE_ALIASES_HPP

#include "variable_desc.hpp"

#include <cctk.h>

#include <vector>

namespace MultipoleX {

using real_vec = std::vector<CCTK_REAL>;

// TODO: WTH is chimera of a type?
using ylm_t = std::vector<std::vector<std::vector<std::vector<CCTK_REAL> > > >;

using variable_desc_vec = std::vector<variable_desc>;

} // namespace MultipoleX

#endif // MULTIPOLEX_TYPE_ALIASES_HPP