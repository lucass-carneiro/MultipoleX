#ifndef MULTIPOLEX_COORDINATES_HPP
#define MULTIPOLEX_COORDINATES_HPP

#include "type_aliases.hpp"

#include <cctk.h>

namespace MultipoleX {

void CoordSetup(real_vec &xhat, real_vec &yhat, real_vec &zhat, real_vec &th,
                real_vec &ph);

void ScaleCartesian(CCTK_REAL r, const real_vec &xhat, const real_vec &yhat,
                    const real_vec &zhat, real_vec &x, real_vec &y,
                    real_vec &z);

} // namespace MultipoleX

#endif // MULTIPOLEX_COORDINATES_HPP