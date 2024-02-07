#ifndef MULTIPOLEX_OUTPUT_HPP
#define MULTIPOLEX_OUTPUT_HPP

#include "type_aliases.hpp"
#include "mode_array.hpp"

#include <cctk.h>
#include <cctk_Arguments.h>

#include <string>

namespace MultipoleX {

enum mp_coord { mp_theta, mp_phi };

void OutputArrayToFile(CCTK_ARGUMENTS, const std::string &name,
                       const real_vec &th, const real_vec &ph,
                       const real_vec &x, const real_vec &y, const real_vec &z,
                       const real_vec &data);

void Output1D(CCTK_ARGUMENTS, const std::string &name, const real_vec &th,
              const real_vec &ph, mp_coord coord, const real_vec &data);

void OutputComplexToFile(CCTK_ARGUMENTS, const string &name, CCTK_REAL redata,
                         CCTK_REAL imdata);

// TODO: Maybe the binary output routines should be updated
void OutputComplexToH5File(CCTK_ARGUMENTS, const variable_desc_vec &vars,
                           const CCTK_REAL radii[], const mode_array &modes);

output_modes(CCTK_ARGUMENTS, const vector<Multipole::variable_desc> &vars,
             const CCTK_REAL radii[], const Multipole::mode_array &modes) {

} // namespace MultipoleX

#endif // MULTIPOLEX_OUTPUT_HPP