#ifndef ALEPH_DEDX_VALID_H
#define ALEPH_DEDX_VALID_H

#include <cmath>

namespace FCCAnalyses {
namespace AlephDedx {

// Single source of the dE/dx validity rule: a failed leg stores the track's
// omega as its value (verbatim copy), so valid = value != omega, with finite
// positive value and error. dQdx.type (pad-leg status only) is not used.
inline bool dEdxValid(float value, float error, float trackOmega) {
  return std::isfinite(value) && value > 0.f && value != trackOmega &&
         std::isfinite(error) && error > 0.f;
}

} // namespace AlephDedx
} // namespace FCCAnalyses

#endif
