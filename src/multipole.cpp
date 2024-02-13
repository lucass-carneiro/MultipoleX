#include "multipole.hpp"
#include "interpolate.hpp"
#include "utils.hpp"
#include "sphericalharmonic.hpp"

#include <cctk_Parameters.h>
#include <cctk_Functions.h>
#include <util_Table.h>

#include <cassert>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>

#include <cassert>

MultipoleX::mode_array::mode_array(int nvars, int nradii, int lmax)
    : nvars{nvars}, nradii{nradii}, lmax{lmax},
      modes(nvars * nradii * (lmax + 1) * (lmax + 1) * 2) {}

CCTK_REAL &MultipoleX::mode_array::operator()(int v, int ri, int l, int m,
                                              bool is_im) {
  return modes.at(mode_idx(v, ri, l, m, is_im));
}

const CCTK_REAL &MultipoleX::mode_array::operator()(int v, int ri, int l, int m,
                                                    bool is_im) const {
  return modes.at(mode_idx(v, ri, l, m, is_im));
}

int MultipoleX::mode_array::get_nvars() const { return nvars; }
int MultipoleX::mode_array::get_nradii() const { return nradii; }
int MultipoleX::mode_array::get_lmax() const { return lmax; }

size_t MultipoleX::mode_array::mode_idx(int v, int ri, int l, int m,
                                        int is_im) const {
  assert(v >= 0 && v < nvars);
  assert(ri >= 0 && ri < nradii);
  assert(l >= 0 && l <= lmax);
  assert(m <= l && -m <= l);
  return size_t(v * nradii * (lmax + 1) * (lmax + 1) * 2 +
                ri * (lmax + 1) * (lmax + 1) * 2 + (l * l + (m + l)) * 2 +
                is_im);
}

static void fill_variable(int idx, const char *optstring, void *callback_arg) {
  assert(idx >= 0);
  assert(callback_arg != 0);

  std::vector<MultipoleX::variable_desc> &vs =
      *(std::vector<MultipoleX::variable_desc> *)callback_arg;

  MultipoleX::variable_desc v;

  v.index = idx;

  // Default values if there is no option string or if the options are
  // not present
  v.imag_index = -1;
  v.spin_weight = 0;
  v.name = std::string(CCTK_VarName(v.index));

  if (optstring != 0) {
    int table = Util_TableCreateFromString(optstring);

    if (table >= 0) {
      const int buffer_length = 256;
      char buffer[buffer_length];

      Util_TableGetInt(table, &v.spin_weight, "sw");
      if (Util_TableGetString(table, buffer_length, buffer, "cmplx") >= 0) {
        v.imag_index = CCTK_VarIndex(buffer);
      }
      if (Util_TableGetString(table, buffer_length, buffer, "name") >= 0) {
        v.name = std::string(buffer);
      }

      const int ierr = Util_TableDestroy(table);
      if (ierr) {
        CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                    "Could not destroy table: %d", ierr);
      }
    }
  }
  vs.push_back(v);
}

static void
parse_variables_string(const std::string &var_string,
                       std::vector<MultipoleX::variable_desc> &vars) {
  int ierr = CCTK_TraverseString(var_string.c_str(), fill_variable, &vars,
                                 CCTK_GROUP_OR_VAR);
  assert(ierr >= 0);
}

static void output_modes(CCTK_ARGUMENTS,
                         const std::vector<MultipoleX::variable_desc> &vars,
                         const CCTK_REAL radii[],
                         const MultipoleX::mode_array &modes) {
  using std::ios;
  using std::ostringstream;
  using std::setiosflags;
  using std::setprecision;

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (output_ascii) {
    if (CCTK_MyProc(cctkGH) == 0) {
      for (int v = 0; v < modes.get_nvars(); v++) {
        for (int i = 0; i < modes.get_nradii(); i++) {
          const CCTK_REAL rad = radii[i];
          for (int l = 0; l <= modes.get_lmax(); l++) {
            for (int m = -l; m <= l; m++) {
              ostringstream name;
              name << "mp_" << vars[v].name << "_l" << l << "_m" << m << "_r"
                   << setiosflags(ios::fixed) << setprecision(2) << rad
                   << ".asc";
              MultipoleX::OutputComplexToFile(CCTK_PASS_CTOC, name.str(),
                                              modes(v, i, l, m, 0),
                                              modes(v, i, l, m, 1));
            }
          }
        }
      }
    }
  }
  if (output_hdf5) {
    if (CCTK_MyProc(cctkGH) == 0) {
      MultipoleX::OutputComplexToH5File(CCTK_PASS_CTOC, vars, radii, modes);
    }
  }
}

void output_1D(CCTK_ARGUMENTS, const MultipoleX::variable_desc &v,
               CCTK_REAL rad, const std::vector<CCTK_REAL> &th,
               const std::vector<CCTK_REAL> &ph,
               const std::vector<CCTK_REAL> &real,
               const std::vector<CCTK_REAL> &imag) {
  using std::ios;
  using std::ostringstream;
  using std::setiosflags;
  using std::setprecision;

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if (CCTK_MyProc(cctkGH) == 0 && output_ascii) {
    if (out_1d_every != 0 && cctk_iteration % out_1d_every == 0) {
      ostringstream real_base;
      real_base << "mp_" << std::string(CCTK_VarName(v.index)) << "_r"
                << setiosflags(ios::fixed) << setprecision(2) << rad;
      MultipoleX::Output1D(CCTK_PASS_CTOC,
                           real_base.str() + std::string(".th.asc"), th, ph,
                           MultipoleX::mp_theta, real);
      MultipoleX::Output1D(CCTK_PASS_CTOC,
                           real_base.str() + std::string(".ph.asc"), th, ph,
                           MultipoleX::mp_phi, real);

      if (v.imag_index != -1) {
        ostringstream imag_base;
        imag_base << "mp_" << std::string(CCTK_VarName(v.imag_index)) << "_r"
                  << setiosflags(ios::fixed) << setprecision(2) << rad;
        MultipoleX::Output1D(CCTK_PASS_CTOC,
                             imag_base.str() + std::string(".th.asc"), th, ph,
                             MultipoleX::mp_theta, imag);
        MultipoleX::Output1D(CCTK_PASS_CTOC,
                             imag_base.str() + std::string(".ph.asc"), th, ph,
                             MultipoleX::mp_phi, imag);
      }
    }
  }
}

static bool int_in_array(int a, const std::vector<int> &array) {
  for (size_t i = 0; i < array.size(); i++) {
    if (array[i] == a)
      return true;
  }
  return false;
}

static int find_int_in_array(int a, const std::vector<int> &array) {
  for (size_t i = 0; i < array.size(); i++) {
    if (array[i] == a)
      return i;
  }
  return -1;
}

static void get_spin_weights(const std::vector<MultipoleX::variable_desc> &vars,
                             std::vector<int> &spin_weights) {
  for (size_t i = 0; i < vars.size(); i++) {
    if (!int_in_array(vars[i].spin_weight, spin_weights)) {
      spin_weights.push_back(vars[i].spin_weight);
    }
  }
}

// For backward compatibility we allow the user to set l_mode instead
// of l_max, but if it is left at the default of -1, l_max is used.
static int get_l_max() {
  DECLARE_CCTK_PARAMETERS;
  return l_mode == -1 ? l_max : l_mode;
}

static void setup_harmonics(
    const std::vector<int> &spin_weights, int lmax, std::vector<CCTK_REAL> &th,
    std::vector<CCTK_REAL> &ph, int array_size,
    std::vector<std::vector<std::vector<std::vector<CCTK_REAL> > > > &reY,
    std::vector<std::vector<std::vector<std::vector<CCTK_REAL> > > > &imY) {
  th.resize(array_size);
  ph.resize(array_size);

  reY.resize(spin_weights.size());
  imY.resize(spin_weights.size());
  for (size_t si = 0; si < spin_weights.size(); si++) {
    int sw = spin_weights[si];

    reY[si].resize(lmax + 1);
    imY[si].resize(lmax + 1);
    for (int l = 0; l <= lmax; l++) {
      reY[si][l].resize(2 * lmax + 1);
      imY[si][l].resize(2 * lmax + 1);
      for (int m = -l; m <= lmax; m++) {
        reY[si][l][m + l].resize(array_size);
        imY[si][l][m + l].resize(array_size);
        MultipoleX::HarmonicSetup(sw, l, m, th, ph, reY[si][l][m + l],
                                  imY[si][l][m + l]);
      }
    }
  }
}

extern "C" void MultipoleX_ParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MultipoleX_ParamCheck;
  DECLARE_CCTK_PARAMETERS;

  if (l_mode != -1) {
    CCTK_WARN(CCTK_WARN_ALERT,
              "The parameter l_mode is deprecated. Use l_max instead.  For "
              "compatibility, l_max = l_mode is being used.");
  }

  if (!CCTK_Equals(mode_type, "deprecated")) {
    CCTK_WARN(CCTK_WARN_ALERT, "The parameter mode_type is deprecated and is "
                               "no longer used.  All modes will be computed.");
  }

  if (l_min != -1) {
    CCTK_WARN(CCTK_WARN_ALERT,
              "The parameter l_min is deprecated and is no longer used.  Modes "
              "from l = 0 will be computed.");
  }

  if (m_mode != -100) {
    CCTK_WARN(
        CCTK_WARN_ALERT,
        "The parameter m_mode is deprecated. All m modes will be computed.");
  }
}

extern "C" void MultipoleX_Calc(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MultipoleX_Calc;
  DECLARE_CCTK_PARAMETERS;

  static std::vector<CCTK_REAL> xs, ys, zs;
  static std::vector<CCTK_REAL> xhat, yhat, zhat;
  static std::vector<CCTK_REAL> th, ph;
  static std::vector<CCTK_REAL> real, imag;
  static std::vector<std::vector<std::vector<std::vector<CCTK_REAL> > > > reY;
  static std::vector<std::vector<std::vector<std::vector<CCTK_REAL> > > > imY;
  static std::vector<MultipoleX::variable_desc> vars;
  static std::vector<int> spin_weights;

  static bool initialized = false;

  const int array_size = (ntheta + 1) * (nphi + 1);

  if (out_every == 0 || cctk_iteration % out_every != 0)
    return;

  int lmax = get_l_max();

  if (!initialized) {
    real.resize(array_size);
    imag.resize(array_size);
    th.resize(array_size);
    ph.resize(array_size);
    xs.resize(array_size);
    ys.resize(array_size);
    zs.resize(array_size);
    xhat.resize(array_size);
    yhat.resize(array_size);
    zhat.resize(array_size);

    parse_variables_string(std::string(variables), vars);
    get_spin_weights(vars, spin_weights);
    MultipoleX::CoordSetup(xhat, yhat, zhat, th, ph);
    setup_harmonics(spin_weights, lmax, th, ph, array_size, reY, imY);
    initialized = true;
  }

  MultipoleX::mode_array modes(vars.size(), nradii, lmax);
  for (size_t v = 0; v < vars.size(); v++) {
    // assert(vars[v].spin_weight == -2);

    int si = find_int_in_array(vars[v].spin_weight, spin_weights);
    assert(si != -1);

    for (int i = 0; i < nradii; i++) {
      // Compute x^i = r * \hat x^i
      MultipoleX::ScaleCartesian(radius[i], xhat, yhat, zhat, xs, ys, zs);

      // Interpolate Psi4r and Psi4i
      MultipoleX::Multipole_Interp(CCTK_PASS_CTOC, xs, ys, zs, vars[v].index,
                                   vars[v].imag_index, real, imag);
      for (int l = 0; l <= lmax; l++) {
        for (int m = -l; m <= l; m++) {
          // Integrate sYlm (real + i imag) over the sphere at radius r
          MultipoleX::Integrate(reY[si][l][m + l], imY[si][l][m + l], real,
                                imag, th, ph, &modes(v, i, l, m, 0),
                                &modes(v, i, l, m, 1));

        } // loop over m
      }   // loop over l
      output_1D(CCTK_PASS_CTOC, vars[v], radius[i], th, ph, real, imag);
    } // loop over radii
  }   // loop over variables
  output_modes(CCTK_PASS_CTOC, vars, radius, modes);
}
