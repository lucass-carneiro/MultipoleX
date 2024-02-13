#ifndef MULTIPOLEX_MULTIPOLE_HPP
#define MULTIPOLEX_MULTIPOLE_HPP

#include <cctk.h>
#include <cctk_Arguments.h>

#include <string>
#include <vector>

// MultipoleX_Calc
//      This is the main scheduling file.  Because we are completely local here
//      and do not use cactus arrays etc, we schedule only one function and then
//      like program like one would in C, C++ with this function taking the
//      place of int main(void).
//
//      This function calls functions to accomplish 3 things:
//        1) Interpolate psi4 onto a sphere
//        2) Integrate psi4 with the ylm's over that sphere
//        2) Output the mode decomposed psi4
extern "C" void MultipoleX_Calc(CCTK_ARGUMENTS);

namespace MultipoleX {

// information about variable which we decompose
struct variable_desc {
  int index;
  int imag_index;
  CCTK_INT spin_weight;
  std::string name;
};

// a simple array class to hold complex modes abs(m) <= l, l <= lmax, for
// nradii radii for nvars variables
class mode_array {
public:
  mode_array(int nvars, int nradii, int lmax);

  CCTK_REAL &operator()(int v, int ri, int l, int m, bool is_im);

  const CCTK_REAL &operator()(int v, int ri, int l, int m, bool is_im) const;

  int get_nvars() const;
  int get_nradii() const;
  int get_lmax() const;

private:
  size_t mode_idx(int v, int ri, int l, int m, int is_im) const;

  const int nvars;
  const int nradii;
  const int lmax;
  std::vector<CCTK_REAL> modes;
};

} // namespace MultipoleX

#endif // MULTIPOLEX_MULTIPOLE_HPP