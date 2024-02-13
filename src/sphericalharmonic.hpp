#ifndef MULTIPOLEX_SPHERICALHARMONIC_HPP
#define MULTIPOLEX_SPHERICALHARMONIC_HPP

#include <cctk.h>

#include <vector>

namespace MultipoleX {

void HarmonicSetup(int s, int l, int m, std::vector<CCTK_REAL> const &th,
                   std::vector<CCTK_REAL> const &ph,
                   std::vector<CCTK_REAL> &reY, std::vector<CCTK_REAL> &imY);

void SphericalHarmonic(int s, int l, int m, CCTK_REAL th, CCTK_REAL ph,
                       CCTK_REAL *reY, CCTK_REAL *imY);

} // namespace MultipoleX

#endif // MULTIPOLEX_SPHERICALHARMONIC_HPP
