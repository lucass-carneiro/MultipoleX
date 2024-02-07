#ifndef MULTIPOLEX_VARIABLE_DESC_HPP
#define MULTIPOLEX_VARIABLE_DESC_HPP

#include <cctk.h>

#include <string>

namespace MultipoleX {

struct variable_desc {
  int index;
  int imag_index;
  CCTK_INT spin_weight;
  std::string name;
};

} // namespace MultipoleX

#endif // MULTIPOLEX_VARIABLE_DESC_HPP