#ifndef ALEPH_SVNEW_ANALYZERS_H
#define ALEPH_SVNEW_ANALYZERS_H

// Standalone secondary-vertex finder on the secondary-track collection, run
// after the V0 module has claimed its tracks (V0-first pipeline).

#include "ROOT/RVec.hxx"
#include "TVector3.h"
#include "edm4hep/TrackState.h"
#include "FCCAnalyses/VertexingUtils.h"
#include "FCCAnalyses/VertexFitterSimple.h"

#include <algorithm>
#include <numeric>
#include <vector>

#include "aleph_units.h"
#include "analyzer_trkaux.h"

namespace FCCAnalyses {
namespace AlephSVNew {

using ROOT::VecOps::RVec;
using AlephTrkAux::fitTracksCm;

constexpr double SVN_MPI = AlephMasses::kPiCh;

// Adopted secondary-vertex selection: the single source for these values.
constexpr double SVN_CHI2 = 10.;        // normalised vertex chi2 (seed and growth)
constexpr double SVN_DIS_LO = 0.03;     // PV displacement window, low edge [cm]
constexpr double SVN_DIS_HI = 3.;       // PV displacement window, high edge [cm]
constexpr double SVN_SIGL_MAX = 0.10;   // longitudinal vertex sigma guard [cm]
constexpr int SVN_MAX_TRK = 8;          // maximum tracks per candidate (growth cap)
constexpr double SVN_TRK_CHI2 = 5.;     // per-track chi2 contribution cap (<=0 off)
constexpr double SVN_COS_POINT = 0.7;   // minimum cosPointing

// V0-track masking modes for findSVs
constexpr int SVN_MASK_NONE = 0;        // mask nothing (unmasked control twin)
constexpr int SVN_MASK_MODE = 1;        // adopted: mask the tight-claimed V0 tracks

// sigma along direction u from the edm4hep lower-triangular covMatrix
// (xx, yx, yy, zx, zy, zz); fit runs in the cm-as-mm homothety, positions cm.
template <typename Cov>
inline double sigmaAlong(const Cov& c, const TVector3& u) {
  double x = u.x(), y = u.y(), z = u.z();
  double var = c[0] * x * x + c[2] * y * y + c[5] * z * z +
               2. * (c[1] * x * y + c[3] * x * z + c[4] * y * z);
  return (var > 0.) ? std::sqrt(var) : 0.;
}

struct SVCand {
  VertexingUtils::FCCAnalysesVertex vtx;
  std::vector<int> trk;
  double chi2;   // normalised
  double mass;
};

// One multi-track fit + the group-level quality cuts; `group` is a scratch
// buffer reused across calls.
inline bool svFitGroup(const RVec<edm4hep::TrackState>& np_tracks,
                       const std::vector<int>& idx, double solenoidBz,
                       RVec<edm4hep::TrackState>& group,
                       VertexingUtils::FCCAnalysesVertex& out) {
  group.clear();
  for (int k : idx) group.push_back(np_tracks[k]);
  auto v = fitTracksCm(group, np_tracks, solenoidBz);
  if ((int)v.updated_track_momentum_at_vertex.size() != (int)idx.size()) return false;
  double chi2 = v.vertex.chi2;  // normalised
  if (!(chi2 == chi2) || chi2 >= SVN_CHI2) return false;
  // per-track compatibility: every track must individually fit the vertex
  if (v.reco_chi2.size() == idx.size())
    for (float rc : v.reco_chi2)
      if (rc > SVN_TRK_CHI2) return false;
  out = v;
  return true;
}

inline bool svPassWindows(const VertexingUtils::FCCAnalysesVertex& v,
                          const TVector3& pv, double& mass_out) {
  TVector3 x(v.vertex.position[0], v.vertex.position[1], v.vertex.position[2]);
  TVector3 d = x - pv;
  double dis = d.Mag();
  if (dis < SVN_DIS_LO || dis > SVN_DIS_HI) return false;
  TVector3 psum(0., 0., 0.);
  double esum = 0.;
  for (const auto& tp : v.updated_track_momentum_at_vertex) {
    psum += tp;
    esum += std::sqrt(tp.Mag2() + SVN_MPI * SVN_MPI);
  }
  const double pmag = psum.Mag();
  if (pmag <= 0.) return false;
  // pointing: the SV momentum must not be anti-aligned with the flight line
  if (d.Dot(psum) / (dis * pmag) < SVN_COS_POINT) return false;
  // collinear-degeneracy guard: fit constrained along the track bundle?
  if (sigmaAlong(v.vertex.covMatrix, psum.Unit()) > SVN_SIGL_MAX) return false;
  double m2 = esum * esum - psum.Mag2();
  mass_out = (m2 > 0.) ? std::sqrt(m2) : 0.;
  return true;
}

// Two-track seed pass over the WHOLE collection, computed once per event and
// read by every masking mode: masking only ever removes tracks, and both the
// seed list and the pair-compatibility table are then filtered, never refitted.
struct SVSeeds {
  std::vector<SVCand> seeds;  // window-passing pairs, (i, j) ascending
  std::vector<char> pairok;   // nTr x nTr: the pair fits a common vertex
  int nTr = 0;
};

inline SVSeeds svSeedPass(const RVec<edm4hep::TrackState>& np_tracks,
                          const VertexingUtils::FCCAnalysesVertex& PV,
                          double solenoidBz) {
  SVSeeds out;
  const int nTr = np_tracks.size();
  out.nTr = nTr;
  if (nTr < 2) return out;
  TVector3 pv(PV.vertex.position[0], PV.vertex.position[1], PV.vertex.position[2]);

  RVec<edm4hep::TrackState> group;
  group.reserve(nTr);
  out.pairok.assign((size_t)nTr * nTr, 0);
  for (int i = 0; i < nTr - 1; ++i) {
    for (int j = i + 1; j < nTr; ++j) {
      VertexingUtils::FCCAnalysesVertex v;
      if (!svFitGroup(np_tracks, {i, j}, solenoidBz, group, v)) continue;
      out.pairok[(size_t)i * nTr + j] = out.pairok[(size_t)j * nTr + i] = 1;
      double m;
      if (!svPassWindows(v, pv, m)) continue;
      out.seeds.push_back({v, {i, j}, v.vertex.chi2, m});
    }
  }
  return out;
}

// findSVs: v0s/v0_tight mask V0-claimed tracks; mask_mode 0 = none, 1 = tight
// only. Returns FCCAnalysesV0 (pdgAbs = 0, invM = N-pion mass); reco_ind
// indexes the SECONDARY collection, NOT the original Tracks index space of
// v0n_trk1/trk2 (map through sec2origIdx to compare). Every selection value is
// a constant defined above.
inline VertexingUtils::FCCAnalysesV0 findSVs(
    const RVec<edm4hep::TrackState>& np_tracks,
    const VertexingUtils::FCCAnalysesVertex& PV,
    const VertexingUtils::FCCAnalysesV0& v0s,
    const RVec<int>& v0_tight,
    int mask_mode,
    double solenoidBz,
    const SVSeeds& pass) {

  VertexingUtils::FCCAnalysesV0 result;
  const int nTr = np_tracks.size();
  if (nTr < 2 || pass.nTr != nTr) return result;

  TVector3 pv(PV.vertex.position[0], PV.vertex.position[1], PV.vertex.position[2]);

  std::vector<bool> masked(nTr, false);
  if (mask_mode > 0) {
    for (size_t iv = 0; iv < v0s.vtx.size(); ++iv) {
      if (mask_mode == 1 && (iv >= v0_tight.size() || v0_tight[iv] != 1)) continue;
      for (int ti : v0s.vtx[iv].reco_ind)
        if (ti >= 0 && ti < nTr) masked[ti] = true;
    }
  }

  // reused across every fit: cleared, not reallocated, on each call
  RVec<edm4hep::TrackState> group;
  group.reserve(nTr);
  std::vector<int> trial;
  trial.reserve(nTr);
  const std::vector<char>& pairok = pass.pairok;

  // ---- growth of ONE candidate ------------------------------------------
  // repeatedly attach the available (unblocked, pair-linked) track giving the
  // best refit chi2, while windows/guards still pass.
  auto growCand = [&](const SVCand& seed, const std::vector<bool>& blocked) {
    SVCand c = seed;
    bool grew = true;
    while (grew && (int)c.trk.size() < SVN_MAX_TRK) {
      grew = false;
      SVCand best = c;
      for (int k = 0; k < nTr; ++k) {
        if (blocked[k]) continue;
        if (std::find(c.trk.begin(), c.trk.end(), k) != c.trk.end()) continue;
        bool linked = false;
        for (int m0 : c.trk)
          if (pairok[(size_t)m0 * nTr + k]) { linked = true; break; }
        if (!linked) continue;
        trial = c.trk;
        trial.push_back(k);
        VertexingUtils::FCCAnalysesVertex v;
        if (!svFitGroup(np_tracks, trial, solenoidBz, group, v)) continue;
        double m;
        if (!svPassWindows(v, pv, m)) continue;
        if (!grew || v.vertex.chi2 < best.chi2) {
          best = {v, trial, v.vertex.chi2, m};
          grew = true;
        }
      }
      if (grew) c = best;
    }
    return c;
  };

  // ---- seed ordering: best normalised chi2 claims first -------------------
  // seeds on a masked track are dropped first, so the surviving order is the
  // one a seed pass restricted to the unmasked tracks would have produced.
  std::vector<size_t> order;
  order.reserve(pass.seeds.size());
  for (size_t s = 0; s < pass.seeds.size(); ++s)
    if (!masked[pass.seeds[s].trk[0]] && !masked[pass.seeds[s].trk[1]])
      order.push_back(s);
  std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
    return pass.seeds[a].chi2 < pass.seeds[b].chi2;
  });

  std::vector<bool> used(nTr, false);
  std::vector<bool> blocked = masked;
  for (size_t s : order) {
    if (used[pass.seeds[s].trk[0]] || used[pass.seeds[s].trk[1]]) continue;
    const SVCand c = growCand(pass.seeds[s], blocked);
    for (int t : c.trk) { used[t] = true; blocked[t] = true; }
    result.vtx.push_back(c.vtx);
    result.pdgAbs.push_back(0);
    result.invM.push_back(c.mass);
  }
  return result;
}

// svn-specific getters (generic candChi2/candDxyz/candP/candCosPointing/
// candPointSig come from AlephTruth/AlephV0New on the shared struct).
inline RVec<int> candNtracks(const VertexingUtils::FCCAnalysesV0& svs) {
  RVec<int> out;
  for (const auto& v : svs.vtx) out.push_back((int)v.reco_ind.size());
  return out;
}

// SV position components relative to the PV (comp 0/1/2 = x/y/z), cm.
inline RVec<float> candDcomp(const VertexingUtils::FCCAnalysesV0& svs,
                             const VertexingUtils::FCCAnalysesVertex& PV,
                             int comp) {
  RVec<float> out;
  for (const auto& v : svs.vtx)
    out.push_back((float)(v.vertex.position[comp] - PV.vertex.position[comp]));
  return out;
}

inline RVec<float> candSigL(const VertexingUtils::FCCAnalysesV0& svs) {
  RVec<float> out;
  for (const auto& v : svs.vtx) {
    TVector3 psum(0., 0., 0.);
    for (const auto& tp : v.updated_track_momentum_at_vertex) psum += tp;
    out.push_back(psum.Mag() > 0.
                      ? (float)sigmaAlong(v.vertex.covMatrix, psum.Unit())
                      : -1.f);
  }
  return out;
}

// Flat candidate<->track association: candTrkSV[k] = candidate index,
// candTrkIdx[k] = its track index in the SECONDARY collection (NOT the
// v0n_trk1/2 original-Tracks space; map through sec2origIdx to compare).
inline RVec<int> candTrkSV(const VertexingUtils::FCCAnalysesV0& svs) {
  RVec<int> out;
  for (size_t i = 0; i < svs.vtx.size(); ++i)
    for (size_t k = 0; k < svs.vtx[i].reco_ind.size(); ++k) out.push_back((int)i);
  return out;
}

inline RVec<int> candTrkIdx(const VertexingUtils::FCCAnalysesV0& svs) {
  RVec<int> out;
  for (const auto& v : svs.vtx)
    for (int t : v.reco_ind) out.push_back(t);
  return out;
}

}  // namespace AlephSVNew
}  // namespace FCCAnalyses

#endif
