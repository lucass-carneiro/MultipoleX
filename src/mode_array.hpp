#ifndef MULTIPOLEX_MODE_ARRAY_HPP
#define MULTIPOLEX_MODE_ARRAY_HPP

#include "type_aliases.hpp"

#include <cctk.h>

namespace MultipoleX {

/*
 * A simple array class to hold complex modes abs(m) <= l, l <= lmax, for nradii
 * radii for nvars variables
 */
class mode_array {
public:
  mode_array(int nvars, int nradii, int lmax);

  CCTK_REAL &operator()(int v, int ri, int l, int m, bool is_im);
  const CCTK_REAL &operator()(int v, int ri, int l, int m, bool is_im) const;

  int get_nvars() const;
  int get_nradii() const;
  int get_lmax() const;

private:
  std::size_t mode_idx(int v, int ri, int l, int m, int is_im) const;

  const int nvars;
  const int nradii;
  const int lmax;
  real_vec modes;
};

} // namespace MultipoleX

#endif // MULTIPOLEX_MODE_ARRAY_HPP