#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

// TODO: Change to checkd arguments
extern "C" void Multipole_ParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

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