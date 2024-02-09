#ifndef MULTIPOLEX_SPHERICALHARMONIC_HPP
#define MULTIPOLEX_SPHERICALHARMONIC_HPP

#include <cctk.h>

#include "type_aliases.hpp"

namespace MultipoleX {

void HarmonicSetup(int s, int l, int m, const real_vec &th, const real_vec &ph,
                   real_vec &reY, real_vec &imY);

void SphericalHarmonic(int s, int l, int m, CCTK_REAL th, CCTK_REAL ph,
                       CCTK_REAL *reY, CCTK_REAL *imY);

} // namespace MultipoleX

#endif // MULTIPOLEX_SPHERICALHARMONIC_HPP