#include "utils.hpp"
#include "integrate.hpp"
#include "interpolate.hpp"
#include "sphericalharmonic.hpp"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cstdio>
#include <cstring>
#include <cmath>
#include <cassert>

static CCTK_REAL test_integral(int n,
                               CCTK_REAL (*integration_fn)(const CCTK_REAL *,
                                                           int, int, CCTK_REAL,
                                                           CCTK_REAL),
                               const int is_midpoint) {
  const int nx = n;
  const int ny = n;
  const int array_size = (nx + 1) * (ny + 1);

  CCTK_REAL *f = new CCTK_REAL[array_size];

  const CCTK_REAL dx = 1. / (nx + is_midpoint);
  const CCTK_REAL dy = 1. / (ny + is_midpoint);

  for (int ix = 0; ix <= nx; ix++) {
    for (int iy = 0; iy <= ny; iy++) {
      const int i = MultipoleX::Index(ix, iy, nx);

      const CCTK_REAL x = ix * dx + 0.5 * dx * is_midpoint;
      const CCTK_REAL y = iy * dy + 0.5 * dy * is_midpoint;
      const CCTK_REAL PI = acos(-1.0);

      f[i] = x * pow(y, 2) * pow(cos(2 * PI * y), 2) * pow(sin(2 * PI * x), 2);
    }
  }

  const CCTK_REAL result = integration_fn(f, nx, ny, dx, dy);
  delete[] f;
  return result;
}

static CCTK_REAL test_pi_symmetric_sphere_integral(
    CCTK_REAL (*integration_fn)(const CCTK_REAL *, int, int, CCTK_REAL,
                                CCTK_REAL),
    const int is_midpoint) {
  const int n = 100;
  const int nth = n;
  const int nph = n;
  const int array_size = (nth + 1) * (nph + 1);
  const CCTK_REAL PI = acos(-1.0);

  CCTK_REAL *f = new CCTK_REAL[array_size];

  const CCTK_REAL dth = PI / (nth + is_midpoint);
  const CCTK_REAL dph = 2 * PI / (nph + is_midpoint);

  for (int ith = 0; ith <= nth; ith++) {
    for (int iph = 0; iph <= nph; iph++) {
      const int i = MultipoleX::Index(ith, iph, nth);

      const CCTK_REAL th = ith * dth + 0.5 * dth * is_midpoint;
      const CCTK_REAL ph = iph * dph + 0.5 * dph * is_midpoint;

      f[i] = -(cos(ph) * sqrt(5 / PI) * pow(cos(th / 2.), 3) * sin(th / 2.));
    }
  }

  const CCTK_REAL result = integration_fn(f, nth, nph, dth, dph);
  delete[] f;
  return result;
}

static CCTK_REAL integration_convergence_order(
    CCTK_REAL (*integration_fn)(const CCTK_REAL *, int, int, CCTK_REAL,
                                CCTK_REAL),
    CCTK_REAL *store_low, CCTK_REAL *store_high, const int is_midpoint) {
  const int n1 = 100;
  const int n2 = 200;
  const CCTK_REAL PI = acos(-1.0);
  const CCTK_REAL result1 = test_integral(100, integration_fn, is_midpoint);
  *store_low = result1;
  const CCTK_REAL result2 = test_integral(200, integration_fn, is_midpoint);
  *store_high = result2;
  const CCTK_REAL exact = 1. / 24 + 1. / (64 * pow(PI, 2));
  const CCTK_REAL error1 = fabs(result1 - exact);
  const CCTK_REAL error2 = fabs(result2 - exact);
  return log10(error1 / error2) / log10((CCTK_REAL)n2 / n1);
}

extern "C" void MultipoleX_TestIntegrationConvergence(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MultipoleX_TestIntegrationConvergence;

  *test_simpson_convergence_order = integration_convergence_order(
      &MultipoleX::Trapezoidal2DIntegral, test_simpson_result_low,
      test_simpson_result_high, 0);
  *test_trapezoidal_convergence_order = integration_convergence_order(
      &MultipoleX::Trapezoidal2DIntegral, test_trapezoidal_result_low,
      test_trapezoidal_result_high, 0);
  *test_midpoint_convergence_order = integration_convergence_order(
      &MultipoleX::Midpoint2DIntegral, test_midpoint_result_low,
      test_midpoint_result_high, 1);
}

extern "C" void MultipoleX_TestIntegrationSymmetry(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MultipoleX_TestIntegrationSymmetry;

  *test_simpson_pi_symmetry =
      test_pi_symmetric_sphere_integral(&MultipoleX::Trapezoidal2DIntegral, 0);
  *test_midpoint_pi_symmetry =
      test_pi_symmetric_sphere_integral(&MultipoleX::Midpoint2DIntegral, 1);
  *test_trapezoidal_pi_symmetry =
      test_pi_symmetric_sphere_integral(&MultipoleX::Trapezoidal2DIntegral, 0);
  *test_driscollhealy_pi_symmetry = test_pi_symmetric_sphere_integral(
      &MultipoleX::DriscollHealy2DIntegral, 0);
  printf("Pi symmetry Simpson integral: %.19g\n", *test_simpson_pi_symmetry);
  printf("Pi symmetry midpoint integral: %.19g\n", *test_midpoint_pi_symmetry);
  printf("Pi symmetry trapezoidal integral: %.19g\n",
         *test_trapezoidal_pi_symmetry);
  printf("Pi symmetry Driscoll and Healy integral: %.19g\n",
         *test_driscollhealy_pi_symmetry);
}

void MultipoleX_TestOrthonormality(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MultipoleX_TestOrthonormality;
  DECLARE_CCTK_PARAMETERS;

  /* Campute Cartesian coordinates of points on the sphere */
  int array_size = (ntheta + 1) * (nphi + 1);

  std::vector<CCTK_REAL> th(array_size);
  std::vector<CCTK_REAL> ph(array_size);
  std::vector<CCTK_REAL> xhat(array_size);
  std::vector<CCTK_REAL> yhat(array_size);
  std::vector<CCTK_REAL> zhat(array_size);

  MultipoleX::CoordSetup(xhat, yhat, zhat, th, ph);

  const int max_l_modes = 10;

  /* Populate spherical-harmonic array */
  std::vector<std::vector<std::vector<std::vector<CCTK_REAL> > > > reY(1);
  std::vector<std::vector<std::vector<std::vector<CCTK_REAL> > > > imY(1);

  for (int sw = 0; sw <= 0; sw++) {
    reY[sw].resize(max_l_modes);
    imY[sw].resize(max_l_modes);
    for (int l = 0; l < max_l_modes; l++) {
      reY[sw][l].resize(2 * max_l_modes + 1);
      imY[sw][l].resize(2 * max_l_modes + 1);
      for (int m = -l; m <= l; m++) {
        reY[sw][l][m + l].resize(array_size, 0.);
        imY[sw][l][m + l].resize(array_size, 0.);

        MultipoleX::HarmonicSetup(sw, l, m, th, ph, reY[sw][l][m + l],
                                  imY[sw][l][m + l]);
      }
    }
  }

  /* Loop over l and m, assign Ylm to (rel,imag), and compute the scalar
     product with all spherical harmonics (loop over li, mi) */
  // [0..max_l_modes) has N=max_l_modes^2
  // comparing each mode with each other but skipping the duplicates
  // gives N*(N+1)/2
  // only 1 spin-weight mode is tested
  const int N = max_l_modes * max_l_modes;
  int idx = 0;
  for (int sw = 0; sw <= 0; sw++) {
    for (int l = 0; l < max_l_modes; l++) {
      for (int m = -l; m <= l; m++) {

        /* Compute scalar product of (real,imag) and all the Ylimi */
        for (int li = 0; li < max_l_modes; li++) {
          for (int mi = -li; mi <= li; mi++) {
            // only handle lower triangle in ((l,m),(li,mi)) space
            if (l * l + l + m < li * li + li + mi)
              continue;

            CCTK_REAL real_lm = 0.0, imag_lm = 0.0;
            MultipoleX::Integrate(reY[sw][li][mi + li], imY[sw][li][mi + li],
                                  reY[sw][l][m + l], imY[sw][l][m + l], th, ph,
                                  &real_lm, &imag_lm);

            assert(idx < 1 * N * (N + 1) / 2);
            test_orthonormality[idx++] =
                sqrt(real_lm * real_lm + imag_lm * imag_lm);
          }
        }
      }
    }
  }
  assert(idx == 1 * N * (N + 1) / 2);

  return;
}
