#include "type_aliases.hpp"
#include "mode_array.hpp"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <vector>

// TODO: Checked arguments
extern "C" void Multipole_Calc(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  static real_vec xs, ys, zs;
  static real_vec xhat, yhat, zhat;
  static real_vec th, ph;
  static real_vec real, imag;
  static ylm_t reY;
  static ylm_t imY;
  static variable_desc_vec vars;
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

    parse_variables_string(string(variables), vars);
    get_spin_weights(vars, spin_weights);
    Multipole_CoordSetup(xhat, yhat, zhat, th, ph);
    setup_harmonics(spin_weights, lmax, th, ph, array_size, reY, imY);
    initialized = true;
  }

  mode_array modes(vars.size(), nradii, lmax);
  for (size_t v = 0; v < vars.size(); v++) {
    // assert(vars[v].spin_weight == -2);

    int si = find_int_in_array(vars[v].spin_weight, spin_weights);
    assert(si != -1);

    for (int i = 0; i < nradii; i++) {
      // Compute x^i = r * \hat x^i
      Multipole_ScaleCartesian(radius[i], xhat, yhat, zhat, xs, ys, zs);

      // Interpolate Psi4r and Psi4i
      Multipole_Interp(CCTK_PASS_CTOC, xs, ys, zs, vars[v].index,
                       vars[v].imag_index, real, imag);
      for (int l = 0; l <= lmax; l++) {
        for (int m = -l; m <= l; m++) {
          // Integrate sYlm (real + i imag) over the sphere at radius r
          Multipole_Integrate(reY[si][l][m + l], imY[si][l][m + l], real, imag,
                              th, ph, &modes(v, i, l, m, 0),
                              &modes(v, i, l, m, 1));

        } // loop over m
      }   // loop over l
      output_1D(CCTK_PASS_CTOC, vars[v], radius[i], th, ph, real, imag);
    } // loop over radii
  }   // loop over variables
  output_modes(CCTK_PASS_CTOC, vars, radius, modes);
}

bool int_in_array(int a, const vector<int> &array) {
  for (size_t i = 0; i < array.size(); i++) {
    if (array[i] == a)
      return true;
  }
  return false;
}

int find_int_in_array(int a, const vector<int> &array) {
  for (size_t i = 0; i < array.size(); i++) {
    if (array[i] == a)
      return i;
  }
  return -1;
}

static void get_spin_weights(const vector<Multipole::variable_desc> &vars,
                             vector<int> &spin_weights) {
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

static void setup_harmonics(const vector<int> &spin_weights, int lmax,
                            vector<CCTK_REAL> &th, vector<CCTK_REAL> &ph,
                            int array_size,
                            vector<vector<vector<vector<CCTK_REAL> > > > &reY,
                            vector<vector<vector<vector<CCTK_REAL> > > > &imY) {
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
        Multipole_HarmonicSetup(sw, l, m, th, ph, reY[si][l][m + l],
                                imY[si][l][m + l]);
      }
    }
  }
}