#ifndef MULTIPOLEX_INTEGRATE_HPP
#define MULTIPOLEX_INTEGRATE_HPP

#include "type_aliases.hpp"

#include <cctk.h>

namespace MultipoleX {

CCTK_REAL Midpoint2DIntegral(CCTK_REAL const *f, int nx, int ny, CCTK_REAL hx,
                             CCTK_REAL hy);

CCTK_REAL Trapezoidal2DIntegral(CCTK_REAL const *f, int nx, int ny,
                                CCTK_REAL hx, CCTK_REAL hy);

CCTK_REAL Simpson2DIntegral(CCTK_REAL const *f, int nx, int ny, CCTK_REAL hx,
                            CCTK_REAL hy);

CCTK_REAL DriscollHealy2DIntegral(CCTK_REAL const *f, int nx, int ny,
                                  CCTK_REAL hx, CCTK_REAL hy);

void Integrate(const real_vec &array1r, const real_vec &array1i,
               const real_vec &array2r, const real_vec &array2i,
               const real_vec &th, const real_vec &pph, CCTK_REAL *out_valr,
               CCTK_REAL *out_vali);

} // namespace MultipoleX

#endif