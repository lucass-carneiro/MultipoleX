#ifndef MULTIPOLEX_INTERPOLATE_HPP
#define MULTIPOLEX_INTERPOLATE_HPP

#include "type_aliases.hpp"

#include <cctk_Arguments.h>

namespace MultipoleX {

/*
 * This function interpolates psi4 onto the sphere in cartesian coordinates as
 * created by MultipoleX_CoordSetup.
 */
void Interp(CCTK_ARGUMENTS, const real_vec &x, const real_vec &y,
            const real_vec &z, int real_idx, int imag_idx, real_vec &psi4r,
            real_vec &psi4i);

} // namespace MultipoleX

#endif // MULTIPOLEX_INTERPOLATE_HPP