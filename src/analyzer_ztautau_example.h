#ifndef ALEPH_ZTAUTAU_EXAMPLE_H
#define ALEPH_ZTAUTAU_EXAMPLE_H

/*
  Helpers of Ztautau_example.py, a Z -> tau tau selection with one charged
  track per thrust hemisphere: hemisphere split of the tracks and of the
  energy-flow neutrals, track kinematics, the energy-flow type of a track, and
  the tau decay mode from the MCParticles decay chain. Lengths cm, momenta and
  energies GeV, angles rad, field T.

  Every per-hemisphere result is an array of size 2: entry 0 = the hemisphere
  the thrust axis points into (p . T > 0), entry 1 = the other one.
*/

#include <cmath>
#include <stdexcept>

#include <ROOT/RVec.hxx>

#include "edm4hep/MCParticleData.h"
#include "edm4hep/ParticleIDData.h"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/TrackData.h"
#include "edm4hep/TrackState.h"

#include "aleph_units.h"
#include "analyzer_alephvars.h"

namespace FCCAnalyses {
namespace AlephZtautau {

using ROOT::VecOps::RVec;

constexpr int kNHem = 2;
using AlephVars::kNotAvailable;     // -999: value of a quantity that is not available
using AlephVars::kNotAvailableInt;

// Hemisphere of a direction: 0 if its projection on the thrust axis is positive, else 1.
inline int hemisphereOf(double px, double py, double pz, double tx, double ty, double tz) {
  return px * tx + py * ty + pz * tz > 0. ? 0 : 1;
}

// ---------------------------------------------------------------------------
// Hemisphere split of the selected tracks
// ---------------------------------------------------------------------------
struct HemTracks {
  RVec<int> n;        // tracks per hemisphere
  RVec<int> origIdx;  // Tracks index of the track of a hemisphere with exactly one track, else -1
};

// states / origIdx: a selected track collection and its Tracks indices (same
// order); direction = perigee direction (phi, tanLambda).
inline HemTracks splitTracks(const RVec<edm4hep::TrackState>& states, const RVec<int>& origIdx,
                             float tx, float ty, float tz) {
  if (states.size() != origIdx.size())
    throw std::runtime_error("AlephZtautau::splitTracks: track states and origIdx differ in size");
  HemTracks out{RVec<int>(kNHem, 0), RVec<int>(kNHem, -1)};
  for (size_t k = 0; k < states.size(); ++k) {
    const auto& s = states[k];
    const int h = hemisphereOf(std::cos(s.phi), std::sin(s.phi), s.tanLambda, tx, ty, tz);
    ++out.n[h];
    out.origIdx[h] = origIdx[k];
  }
  for (int h = 0; h < kNHem; ++h)
    if (out.n[h] != 1) out.origIdx[h] = -1;
  return out;
}

// ---------------------------------------------------------------------------
// Energy-flow neutrals per hemisphere
// ---------------------------------------------------------------------------
struct HemNeutrals {
  RVec<float> energy;  // summed energy per hemisphere [GeV]
  RVec<int> count;     // number per hemisphere
};

// Particles of energy-flow type `type` with energy >= eMin; pid index-parallel to particles.
inline HemNeutrals neutralsPerHem(const RVec<edm4hep::ReconstructedParticleData>& particles,
                                  const RVec<edm4hep::ParticleIDData>& pid, int type, double eMin,
                                  float tx, float ty, float tz) {
  if (pid.size() != particles.size())
    throw std::runtime_error("AlephZtautau::neutralsPerHem: particles and ParticleID differ in size");
  HemNeutrals out{RVec<float>(kNHem, 0.f), RVec<int>(kNHem, 0)};
  for (size_t i = 0; i < particles.size(); ++i) {
    if (pid[i].type != type || particles[i].energy < eMin) continue;
    const auto& m = particles[i].momentum;
    const int h = hemisphereOf(m.x, m.y, m.z, tx, ty, tz);
    out.energy[h] += particles[i].energy;
    ++out.count[h];
  }
  return out;
}

// ---------------------------------------------------------------------------
// Track kinematics and energy-flow type
// ---------------------------------------------------------------------------
struct TrackKin {
  RVec<float> p, pt, costheta, phi;  // at the perigee [GeV, rad]
  RVec<int> q;                       // charge [e]
  RVec<float> d0, z0;                // impact parameters w.r.t. the beam spot [cm]
  RVec<float> d0Err, z0Err;          // sqrt of the d0 / z0 covariance diagonal [cm]
};

// For each Tracks index in origIdx (-1 = none: kNotAvailable / 0), from the
// stored (ALEPH-convention) track state: momentum from the curvature in the
// field bz, charge -sign(omega), and d0, z0 re-referenced to the beam spot (bsx, bsy,
// bsz) to first order, d0 in the LCIO sign convention (as AlephTrkAux::pushTrk).
inline TrackKin trackKinematics(const RVec<int>& origIdx, const RVec<edm4hep::TrackData>& tracks,
                                const RVec<edm4hep::TrackState>& states, double bz,
                                double bsx, double bsy, double bsz) {
  TrackKin out;
  for (int o : origIdx) {
    if (o < 0 || (size_t)o >= tracks.size() || tracks[o].trackStates_begin >= states.size()) {
      for (auto* a : {&out.p, &out.pt, &out.costheta, &out.phi, &out.d0, &out.z0, &out.d0Err, &out.z0Err})
        a->push_back(kNotAvailable);
      out.q.push_back(0);
      continue;
    }
    const auto& s = states[tracks[o].trackStates_begin];
    const double om = std::abs(s.omega);
    const double pt = om > 0. ? AlephUnits::kPtPerTeslaCm * bz / om : 0.;
    const double pz = pt * s.tanLambda, p = std::hypot(pt, pz);
    out.p.push_back(p);
    out.pt.push_back(pt);
    out.costheta.push_back(p > 0. ? pz / p : kNotAvailable);
    out.phi.push_back(s.phi);
    out.q.push_back(s.omega < 0. ? 1 : (s.omega > 0. ? -1 : 0));
    // flipped convention (d0 = -D0, omega -> -omega): perigee at +d0 * (-sin phi, cos phi)
    const double d0f = -s.D0, omf = -s.omega;
    const double cphi = std::cos(s.phi), sphi = std::sin(s.phi);
    const double sbs = bsx * cphi + bsy * sphi;
    out.d0.push_back(d0f + bsx * sphi - bsy * cphi - 0.5 * omf * sbs * sbs);
    out.z0.push_back(s.Z0 - bsz + s.tanLambda * sbs);
    // packed lower triangle of (d0, phi, omega, z0, tanLambda): d0 = 0, z0 = 9
    out.d0Err.push_back(s.covMatrix[0] > 0. ? std::sqrt(s.covMatrix[0]) : kNotAvailable);
    out.z0Err.push_back(s.covMatrix[9] > 0. ? std::sqrt(s.covMatrix[9]) : kNotAvailable);
  }
  return out;
}

// ParticleID.type of each track's energy-flow particle, -1 = none. rpOfTrack:
// Tracks index -> particle index (AlephTrkAux::rpIndexByTrack); pid
// index-parallel to those particles.
inline RVec<int> trackEFType(const RVec<int>& origIdx, const RVec<int>& rpOfTrack,
                             const RVec<edm4hep::ParticleIDData>& pid) {
  RVec<int> out;
  for (int o : origIdx) {
    const int r = (o >= 0 && (size_t)o < rpOfTrack.size()) ? rpOfTrack[o] : -1;
    out.push_back(r >= 0 && (size_t)r < pid.size() ? pid[r].type : -1);
  }
  return out;
}

// ---------------------------------------------------------------------------
// Tau decay mode from the MCParticles decay chain (simulation)
// ---------------------------------------------------------------------------
enum TauMode : int {
  kModeNoTau = -1,     // MCParticles[0] and [1] are not a tau pair (not a tau sample)
  kModeOther = 0,      // none of the classes below (e.g. five charged prongs)
  kModeElectron = 1,   // one charged prong, an electron
  kModeMuon = 2,       // one charged prong, a muon
  kModeHadron = 3,     // one charged hadron, no pi0
  kModeHadronPi0 = 4,  // one charged hadron, one pi0
  kModeHadron2Pi0 = 5, // one charged hadron, two or more pi0
  kMode3Prong = 6      // three charged prongs, any number of pi0
};

// generatorStatus = FKIN history code: 10000 * KS + mother row for generator
// particles, 0 for GEANT secondaries (ALEPH kxl7al.F, fyfkin.F)
constexpr int kGeneratorStatusMin = 10000;
constexpr int kMaxDecayDepth = 100;  // guard against a cyclic daughter relation

struct TauDecay {
  int nCharged = 0, nPi0 = 0;
  int prongPdg = 0;              // PDG code of the (last) charged prong
  double vis[3] = {0., 0., 0.};  // visible momentum: final-state particles except neutrinos [GeV]
  int mode = kModeOther;
};

// Follows the generator daughters of particle i. Final = no generator daughter,
// or a pi0, K0S, K0L, Lambda or photon (their decay products are not counted).
inline void walkDecay(int i, const RVec<edm4hep::MCParticleData>& mc, const RVec<int>& daughters,
                      TauDecay& d, int depth) {
  if (depth > kMaxDecayDepth) throw std::runtime_error("AlephZtautau::walkDecay: decay chain deeper than kMaxDecayDepth");
  const auto& p = mc[i];
  RVec<int> gd;
  for (unsigned k = p.daughters_begin; k < p.daughters_end && k < daughters.size(); ++k) {
    const int j = daughters[k];
    if (j < 0 || (size_t)j >= mc.size())
      throw std::runtime_error("AlephZtautau::walkDecay: MCParticles daughter out of range");
    if (mc[j].generatorStatus >= kGeneratorStatusMin) gd.push_back(j);
  }
  const int apdg = std::abs(p.PDG);
  const bool stop = apdg == 111 || apdg == 310 || apdg == 130 || apdg == 3122 || apdg == 22;
  if (depth > 0 && (gd.empty() || stop)) {
    if (apdg == 12 || apdg == 14 || apdg == 16) return;
    d.vis[0] += p.momentum.x;
    d.vis[1] += p.momentum.y;
    d.vis[2] += p.momentum.z;
    if (apdg == 111) {
      ++d.nPi0;
    } else if (std::abs(p.charge) > 0.5f) {
      ++d.nCharged;
      d.prongPdg = p.PDG;
    }
    return;
  }
  for (int j : gd) walkDecay(j, mc, daughters, d, depth + 1);
}

inline TauDecay tauDecay(int tau, const RVec<edm4hep::MCParticleData>& mc, const RVec<int>& daughters) {
  TauDecay d;
  walkDecay(tau, mc, daughters, d, 0);
  if (d.nCharged == 1) {
    const int a = std::abs(d.prongPdg);
    if (a == 11) d.mode = kModeElectron;
    else if (a == 13) d.mode = kModeMuon;
    else d.mode = d.nPi0 == 0 ? kModeHadron : (d.nPi0 == 1 ? kModeHadronPi0 : kModeHadron2Pi0);
  } else if (d.nCharged == 3) {
    d.mode = kMode3Prong;
  }
  return d;
}

struct TauTruth {
  RVec<int> mode, nPi0;  // per hemisphere
  RVec<float> visP;      // magnitude of the visible momentum [GeV]
};

// The taus are MCParticles[0] and [1]; the one whose visible momentum has the
// larger projection on the thrust axis (tx, ty, tz) is given to hemisphere 0.
inline TauTruth tauTruthPerHem(const RVec<edm4hep::MCParticleData>& mc, const RVec<int>& daughters,
                               float tx, float ty, float tz) {
  TauTruth out{RVec<int>(kNHem, kModeNoTau), RVec<int>(kNHem, kNotAvailableInt), RVec<float>(kNHem, kNotAvailable)};
  if (mc.size() < 2 || std::abs(mc[0].PDG) != 15 || mc[0].PDG != -mc[1].PDG) return out;
  const TauDecay dec[2] = {tauDecay(0, mc, daughters), tauDecay(1, mc, daughters)};
  const double p0 = dec[0].vis[0] * tx + dec[0].vis[1] * ty + dec[0].vis[2] * tz;
  const double p1 = dec[1].vis[0] * tx + dec[1].vis[1] * ty + dec[1].vis[2] * tz;
  const int tauOfHem0 = p0 >= p1 ? 0 : 1;
  for (int h = 0; h < kNHem; ++h) {
    const TauDecay& d = dec[h == 0 ? tauOfHem0 : 1 - tauOfHem0];
    out.mode[h] = d.mode;
    out.nPi0[h] = d.nPi0;
    out.visP[h] = std::sqrt(d.vis[0] * d.vis[0] + d.vis[1] * d.vis[1] + d.vis[2] * d.vis[2]);
  }
  return out;
}

}  // namespace AlephZtautau
}  // namespace FCCAnalyses

#endif  // ALEPH_ZTAUTAU_EXAMPLE_H
