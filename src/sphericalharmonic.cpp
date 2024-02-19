// TODO: Maybe use SSHT instead of computing our own harmonics

#include "sphericalharmonic.hpp"

#include <loop.hxx>

#include <cctk.h>
#include <cctk_Parameters.h>
#include <cctk_Arguments.h>

#include <cassert>
#include <cmath>
#include <iostream>

// TODO: Use M_Pi?
static const CCTK_REAL PI = acos(-1.0);

static double factorial(int n) {
  double returnval = 1;
  for (int i = n; i >= 1; i--) {
    returnval *= i;
  }
  return returnval;
}

static inline double combination(int n, int m) {
  // Binomial coefficient is undefined if these conditions do not hold
  assert(n >= 0);
  assert(m >= 0);
  assert(m <= n);

  return factorial(n) / (factorial(m) * factorial(n - m));
}

static inline int imin(int a, int b) { return a < b ? a : b; }

static inline int imax(int a, int b) { return a > b ? a : b; }

void MultipoleX::SphericalHarmonic(int s, int l, int m, CCTK_REAL th,
                                   CCTK_REAL ph, CCTK_REAL *reY,
                                   CCTK_REAL *imY) {
  //  assert(s == -2 && l == 2 && m == 2);
  //  *reY = 1.0/2.0 * sqrt(5/PI) * pow(cos(th/2), 4) * cos(2*ph);
  //  *imY = 1.0/2.0 * sqrt(5/PI) * pow(cos(th/2), 4) * sin(2*ph);
  double all_coeff = 0, sum = 0;
  all_coeff = pow(-1.0, m);
  all_coeff *= sqrt(factorial(l + m) * factorial(l - m) * (2 * l + 1) /
                    (4. * PI * factorial(l + s) * factorial(l - s)));
  sum = 0.;
  for (int i = imax(m - s, 0); i <= imin(l + m, l - s); i++) {
    double sum_coeff = combination(l - s, i) * combination(l + s, i + s - m);
    sum += sum_coeff * pow(-1.0, l - i - s) * pow(cos(th / 2.), 2 * i + s - m) *
           pow(sin(th / 2.), 2 * (l - i) + m - s);
  }
  *reY = all_coeff * sum * cos(m * ph);
  *imY = all_coeff * sum * sin(m * ph);
}

void MultipoleX::HarmonicSetup(int s, int l, int m,
                               std::vector<CCTK_REAL> const &th,
                               std::vector<CCTK_REAL> const &ph,
                               std::vector<CCTK_REAL> &reY,
                               std::vector<CCTK_REAL> &imY) {
  for (size_t i = 0; i < th.size(); i++) {
    MultipoleX::SphericalHarmonic(s, l, m, th[i], ph[i], &reY[i], &imY[i]);
  }
}

// Fill a grid function with a given spherical harmonic
extern "C" void MultipoleX_SetHarmonic(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_MultipoleX_SetHarmonic;
  DECLARE_CCTK_PARAMETERS;

  using std::acos;
  using std::atan2;
  using std::sqrt;

  grid.loop_all<0, 0, 0>(
      grid.nghostzones,
      [=](const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const auto x = vcoordx(p.I);
        const auto y = vcoordy(p.I);
        const auto z = vcoordz(p.I);
        const auto r = sqrt(x * x + y * y + z * z);
        const auto theta = acos(z / r);
        const auto phi = atan2(y, x);

        CCTK_REAL re = 0;
        CCTK_REAL im = 0;

        MultipoleX::SphericalHarmonic(test_sw, test_l, test_m, theta, phi, &re,
                                      &im);

        const auto fac = test_mode_proportional_to_r ? r : 1.0;

        harmonic_re(p.I) = re * fac;
        harmonic_im(p.I) = im * fac;
      });
}
