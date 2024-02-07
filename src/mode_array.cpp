#include "mode_array.hpp"

#include <cctk.h>

#include <cassert>

MultipoleX::mode_array::mode_array(int nvars, int nradii, int lmax)
    : nvars{nvars}, nradii{nradii}, lmax{lmax},
      modes{size_t(nvars * nradii * (lmax + 1) * (lmax + 1) * 2)} {}

CCTK_REAL &MultipoleX::mode_array::operator()(int v, int ri, int l, int m,
                                              bool is_im) {
  return modes.at(mode_idx(v, ri, l, m, is_im));
}

const CCTK_REAL &MultipoleX::mode_array::operator()(int v, int ri, int l, int m,
                                                    bool is_im) const {
  return modes.at(mode_idx(v, ri, l, m, is_im));
}

int MultipoleX::mode_array::get_nvars() const { return nvars; }

int MultipoleX::mode_array::get_nradii() const { return nradii; }

int MultipoleX::mode_array::get_lmax() const { return lmax; }

std::size_t MultipoleX::mode_array::mode_idx(int v, int ri, int l, int m,
                                             int is_im) const {
  assert(v >= 0 && v < nvars);
  assert(ri >= 0 && ri < nradii);
  assert(l >= 0 && l <= lmax);
  assert(m <= l && -m <= l);
  return size_t(v * nradii * (lmax + 1) * (lmax + 1) * 2 +
                ri * (lmax + 1) * (lmax + 1) * 2 + (l * l + (m + l)) * 2 +
                is_im);
}
