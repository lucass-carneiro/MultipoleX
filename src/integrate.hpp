
#ifndef MULTIPOLEX_INTEGRATE_HPP
#define MULTIPOLEX_INTEGRATE_HPP

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

} // namespace MultipoleX

#endif // MULTIPOLEX_INTEGRATE_HPP