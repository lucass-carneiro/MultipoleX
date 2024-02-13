#ifndef MULTIPOLEX_INTERPOLATE_HPP
#define MULTIPOLEX_INTERPOLATE_HPP

#include <cctk.h>
#include <cctk_Arguments.h>

#include <vector>

namespace MultipoleX {

// Multipole_Interp:
//      This function interpolates psi4 onto the sphere in cartesian
//      coordinates as created by MultipoleX::CoordSetup.
void Multipole_Interp(CCTK_ARGUMENTS, std::vector<CCTK_REAL> const &x,
                      std::vector<CCTK_REAL> const &y,
                      std::vector<CCTK_REAL> const &z, int real_idx,
                      int imag_idx, std::vector<CCTK_REAL> &psi4r,
                      std::vector<CCTK_REAL> &psi4i);

} // namespace MultipoleX

#endif // MULTIPOLEX_INTERPOLATE_HPP