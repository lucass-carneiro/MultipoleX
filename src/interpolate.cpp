#include "interpolate.hpp"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <util_Table.h>

static void report_interp_error(int ierr) {
  if (ierr < 0) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "CCTK_InterpGridArrays returned error code %d", ierr);
  }
}

/* TODO: This function uses the cactus interpolation infrastructure, so it
 * should work out of the box, but maybe it won't
 */
void MultipoleX::Interp(CCTK_ARGUMENTS, const real_vec &xs, const real_vec &ys,
                        const real_vec &zs, int real_idx, int imag_idx,
                        real_vec &sphere_real, real_vec &sphere_imag) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /*
   * Need parameters for the following:
   * ntheta (dtheta = pi/(ntheta)
   * nphi (dphi = 2pi/(nphi)
   * r (radius of sphere)
   * NOTE: depending on the interval of integration, denominator above may
   * need to be modified to avoid double counting
   */

  CCTK_INT num_input_arrays = imag_idx == -1 ? 1 : 2;
  CCTK_INT num_output_arrays = imag_idx == -1 ? 1 : 2;
  const CCTK_INT num_dims = 3;
  int ierr = -1;

  const void *interp_coords[num_dims] = {(const void *)xs.data(),
                                         (const void *)ys.data(),
                                         (const void *)zs.data()};

  const CCTK_INT input_array_indices[2] = {real_idx, imag_idx};

  const CCTK_INT output_array_types[2] = {CCTK_VARIABLE_REAL,
                                          CCTK_VARIABLE_REAL};

  void *output_arrays[2] = {(void *)sphere_real.data(),
                            (void *)sphere_imag.data()};

  const int operator_handle = CCTK_InterpHandle(interpolator_name);

  int param_table_handle = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);
  ierr = Util_TableSetFromString(param_table_handle, interpolator_pars);

  const int coord_system_handle = CCTK_CoordSystemHandle(coord_system);

  ierr = CCTK_InterpGridArrays(
      cctkGH, num_dims, operator_handle, param_table_handle,
      coord_system_handle,
      CCTK_MyProc(cctkGH) == 0 ? (ntheta + 1) * (nphi + 1)
                               : 0, // Only the 0 processor needs the points
      CCTK_VARIABLE_REAL, interp_coords, num_input_arrays, input_array_indices,
      num_output_arrays, output_array_types, output_arrays);

  report_interp_error(ierr);

  if (imag_idx == -1) {
    for (int i = 0; i < (ntheta + 1) * (nphi + 1); i++) {
      sphere_imag[i] = 0;
    }
  }

  Util_TableDestroy(param_table_handle);
}
