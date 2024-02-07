#include "coordinates.hpp"

#include <cmath>

void MultipoleX::CoordSetup(real_vec &xhat, real_vec &yhat, real_vec &zhat,
                            real_vec &th, real_vec &ph) {
  using std::cos;
  using std::sin;

  DECLARE_CCTK_PARAMETERS;

  // Add an offset for midpoint integration.
  CCTK_REAL is_midpoint = 0.0;
  if (CCTK_Equals(integration_method, "midpoint")) {
    is_midpoint = 1.0;
  }
  const CCTK_REAL dth = M_PI / (ntheta + is_midpoint);
  const CCTK_REAL dph = 2 * M_PI / (nphi + is_midpoint);
  for (int it = 0; it <= ntheta; it++) {
    for (int ip = 0; ip <= nphi; ip++) {
      const int i = Multipole_Index(it, ip, ntheta);

      /*
       * Check for when midpoint enabled:
       *   dth = M_PI/(ntheta+1) -> ntheta = M_PI/dth - 1
       *   th[i] = i * dth + 0.5*dth
       *  Therefore:
       *   th[0]      = 0.5*dth      <- GOOD.
       *   th[ntheta] = ntheta*dth + 0.5*dth
       *              = ((M_PI/dth)-1)*dth + 0.5*dth
       *              = M_PI - dth + 0.5*dth
       *              = M_PI - 0.5*dth <- GOOD.
       *  Similarly for ph.
       *
       * Check for when midpoint disabled:
       *   dth = M_PI/ntheta -> ntheta = M_PI/dth
       *   th[i] = i * dth
       *  Therefore:
       *   th[0]      = 0      <- GOOD.
       *   th[ntheta] = ntheta*dth
       *              = M_PI/dth*dth
       *              = M_PI     <- GOOD.
       * Similarly for ph.
       */

      th[i] = it * dth + 0.5 * dth * is_midpoint;
      ph[i] = ip * dph + 0.5 * dph * is_midpoint;
      xhat[i] = cos(ph[i]) * sin(th[i]);
      yhat[i] = sin(ph[i]) * sin(th[i]);
      zhat[i] = cos(th[i]);
    }
  }
}

void MultipoleX::ScaleCartesian(const CCTK_REAL r, const real_vec &xhat,
                                const real_vec &yhat, const real_vec &zhat,
                                real_vec &x, real_vec &y, real_vec &z) {
  for (size_t i = 0; i < xhat.size(); i++) {
    x[i] = r * xhat[i];
    y[i] = r * yhat[i];
    z[i] = r * zhat[i];
  }
}