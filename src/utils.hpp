#ifndef MULTIPOLEX_UTILS_HPP
#define MULTIPOLEX_UTILS_HPP

#include "multipole.hpp"

#include <cctk.h>

#include <vector>
#include <string>

namespace MultipoleX {

enum mp_coord { mp_theta, mp_phi };

void OutputArrayToFile(CCTK_ARGUMENTS, const std::string &name,
                       std::vector<CCTK_REAL> const &th,
                       std::vector<CCTK_REAL> const &ph,
                       std::vector<CCTK_REAL> const &x,
                       std::vector<CCTK_REAL> const &y,
                       std::vector<CCTK_REAL> const &z,
                       std::vector<CCTK_REAL> const &data);

void Output1D(CCTK_ARGUMENTS, const std::string &name,
              std::vector<CCTK_REAL> const &th,
              std::vector<CCTK_REAL> const &ph, mp_coord coord,
              std::vector<CCTK_REAL> const &data);

void OutputComplexToFile(CCTK_ARGUMENTS, const std::string &name,
                         CCTK_REAL redata, CCTK_REAL imdata);

void OutputComplexToH5File(CCTK_ARGUMENTS,
                           const std::vector<variable_desc> &vars,
                           const CCTK_REAL radii[],
                           const MultipoleX::mode_array &modes);

void CoordSetup(std::vector<CCTK_REAL> &xhat, std::vector<CCTK_REAL> &yhat,
                std::vector<CCTK_REAL> &zhat, std::vector<CCTK_REAL> &th,
                std::vector<CCTK_REAL> &ph);

void ScaleCartesian(CCTK_REAL r, std::vector<CCTK_REAL> const &xhat,
                    std::vector<CCTK_REAL> const &yhat,
                    std::vector<CCTK_REAL> const &zhat,
                    std::vector<CCTK_REAL> &x, std::vector<CCTK_REAL> &y,
                    std::vector<CCTK_REAL> &z);

static inline int Index(int it, int ip, int ntheta) {
  return it + (ntheta + 1) * ip;
}

void Integrate(std::vector<CCTK_REAL> const &array1r,
               std::vector<CCTK_REAL> const &array1i,
               std::vector<CCTK_REAL> const &array2r,
               std::vector<CCTK_REAL> const &array2i,
               std::vector<CCTK_REAL> const &th,
               std::vector<CCTK_REAL> const &pph, CCTK_REAL *out_valr,
               CCTK_REAL *out_vali);

} // namespace MultipoleX

#endif // MULTIPOLEX_UTILS_HPP
