#ifndef ANALYZER_JOINTVTX_H
#define ANALYZER_JOINTVTX_H

/*
  Stage1 adapter of the joint PV + SV vertexer with V0s: builds the Event of
  jointvtx_event.h from stage1 columns and exposes the Result as branches.

  Per-track branches run over the whole Tracks collection; tracks outside the
  baseline list read JV_FILL_INT / JV_FILL_FLOAT. Covariances are packed
  xx, yx, yy, zx, zy, zz. The per-SV derived quantities come from the
  getters of the existing SV branches, on the converted vertices
  (svsAsFCC / pvAsFCC).
*/

#include <algorithm>
#include <cmath>
#include <vector>

#include <Eigen/Dense>
#include <ROOT/RVec.hxx>
#include "TVector3.h"

#include "edm4hep/TrackData.h"
#include "edm4hep/TrackState.h"
#include "FCCAnalyses/VertexingUtils.h"

#include "jointvtx_common.h"
#include "jointvtx_event.h"
#include "analyzer_pvnew.h"
#include "analyzer_truth.h"
#include "analyzer_v0new.h"
#include "analyzer_svnew.h"

namespace FCCAnalyses {
namespace AlephJointVtx {

using ROOT::VecOps::RVec;

// beam-spot setting of the nominal settings
constexpr BeamSpotSetting JV_BS_SETTING = BeamSpotSetting::S1;

constexpr int JV_FILL_INT = -1;
constexpr float JV_FILL_FLOAT = -999.f;

// smallest |omega| of the momentum, 1/cm
constexpr double JV_OMEGA_MIN = 1e-9;

namespace detail {

constexpr int kPack[6][2] = {{0, 0}, {1, 0}, {1, 1}, {2, 0}, {2, 1}, {2, 2}};

inline int hitCount(const edm4hep::TrackData& t, const RVec<int>& hits, unsigned k) {
  const unsigned i = t.subdetectorHitNumbers_begin + k;
  return (i < t.subdetectorHitNumbers_end && i < hits.size()) ? hits[i] : 0;
}

inline std::vector<int> localIndex(const Event& ev) {
  std::vector<int> loc(std::max(ev.nTracksAll, 0), -1);
  for (int i = 0; i < ev.size(); ++i)
    if (ev.origIdx[i] >= 0 && ev.origIdx[i] < ev.nTracksAll) loc[ev.origIdx[i]] = i;
  return loc;
}

inline int origOf(const Event& ev, int local) {
  return (local >= 0 && local < ev.size()) ? ev.origIdx[local] : JV_FILL_INT;
}

template <typename T, typename F>
inline RVec<F> overTracks(const Event& ev, const std::vector<T>& perLocal, F fill) {
  RVec<F> out(std::max(ev.nTracksAll, 0), fill);
  const int n = std::min<int>(ev.size(), perLocal.size());
  for (int i = 0; i < n; ++i)
    if (ev.origIdx[i] >= 0 && ev.origIdx[i] < ev.nTracksAll) out[ev.origIdx[i]] = static_cast<F>(perLocal[i]);
  return out;
}

}  // namespace detail

// states / statesFlipped / tracks / origIdx: the baseline selection (raw and
// sign-flipped states, Tracks index); jetTrackIdx: Tracks indices of the jet
// constituents; v0Trk1 / v0Trk2 / v0Tight: the existing V0 finder's candidates.
inline Event buildEvent(const RVec<edm4hep::TrackState>& states,
                        const RVec<edm4hep::TrackState>& statesFlipped,
                        const RVec<edm4hep::TrackData>& tracks,
                        const RVec<int>& subdetHits,
                        const RVec<int>& origIdx,
                        int nTracksAll,
                        const RVec<RVec<int>>& jetTrackIdx,
                        const RVec<double>& jetTheta,
                        const RVec<double>& jetPhi,
                        const RVec<double>& jetP,
                        double bsx, double bsy, double bsz,
                        double thrustX, double thrustY, double thrustZ,
                        const RVec<int>& v0Trk1,
                        const RVec<int>& v0Trk2,
                        const RVec<int>& v0Tight,
                        double ptConst = JV_PT_CONST) {
  Event ev;
  const int n = static_cast<int>(states.size());
  ev.nTracksAll = nTracksAll;
  ev.ptConst = ptConst;
  ev.bs = Vec3(bsx, bsy, bsz);
  ev.thrust = Vec3(thrustX, thrustY, thrustZ);
  for (int k = 0; k < 2; ++k) {
    if (k < static_cast<int>(jetTheta.size())) ev.jetTheta[k] = jetTheta[k];
    if (k < static_cast<int>(jetPhi.size())) ev.jetPhi[k] = jetPhi[k];
    if (k < static_cast<int>(jetP.size())) ev.jetP[k] = jetP[k];
  }
  std::vector<int> jetOf(std::max(nTracksAll, 0), -1);
  for (int j = 0; j < std::min<int>(2, jetTrackIdx.size()); ++j)
    for (int ti : jetTrackIdx[j])
      if (ti >= 0 && ti < nTracksAll) jetOf[ti] = j;

  const AlephPVNew::detail::TrackSet ts = AlephPVNew::detail::convert(statesFlipped);
  ev.fpar = ts.par;
  ev.fcov = ts.cov;
  ev.par.reserve(n);
  ev.cov.reserve(n);
  for (int i = 0; i < n; ++i) {
    const auto& s = states[i];
    const double om = s.omega, tl = s.tanLambda;
    Vec5 pr;
    pr << static_cast<double>(s.D0), static_cast<double>(s.phi), om, static_cast<double>(s.Z0), tl;
    ev.par.push_back(pr);
    std::array<double, 15> c15;
    for (int k = 0; k < 15; ++k) c15[k] = static_cast<double>(s.covMatrix[k]);
    ev.cov.push_back(c15);
    ev.q.push_back(om > 0. ? -1 : 1);
    ev.p.push_back(ptConst / std::max(std::abs(om), JV_OMEGA_MIN) * std::sqrt(1. + tl * tl));
    const int o = origIdx[i];
    ev.origIdx.push_back(o);
    ev.jet.push_back((o >= 0 && o < nTracksAll) ? jetOf[o] : -1);
    ev.nvdet.push_back(detail::hitCount(tracks[i], subdetHits, 0));
    ev.nitc.push_back(detail::hitCount(tracks[i], subdetHits, 1));
    ev.ntpc.push_back(detail::hitCount(tracks[i], subdetHits, 2));
  }
  const std::vector<int> loc = detail::localIndex(ev);
  const auto toLocal = [&](int o) { return (o >= 0 && o < nTracksAll) ? loc[o] : -1; };
  const size_t nv = std::min({v0Trk1.size(), v0Trk2.size(), v0Tight.size()});
  for (size_t k = 0; k < nv; ++k) {
    ev.v0T1.push_back(toLocal(v0Trk1[k]));
    ev.v0T2.push_back(toLocal(v0Trk2[k]));
    ev.v0Tight.push_back(v0Tight[k]);
  }
  return ev;
}

// ---- per event
// 0 = ok; otherwise the failed step (ResultStatus, jointvtx_event.h): no vertices, per-track fill values
inline int eventStatus(const Result& r) { return r.status; }
inline float pvComp(const Result& r, int i) { return static_cast<float>(r.pv[i]); }
inline float pvCov(const Result& r, int k) {
  return static_cast<float>(r.pvC(detail::kPack[k][0], detail::kPack[k][1]));
}
inline int nSV(const Result& r) { return static_cast<int>(r.svs.size()); }
inline int nV0(const Result& r) { return static_cast<int>(r.v0s.size()); }

// ---- per SV
inline RVec<float> svPos(const Result& r, int i) {
  RVec<float> out;
  for (const auto& s : r.svs) out.push_back(static_cast<float>(s.x[i]));
  return out;
}
inline RVec<float> svCov(const Result& r, int k) {
  RVec<float> out;
  for (const auto& s : r.svs) out.push_back(static_cast<float>(s.C(detail::kPack[k][0], detail::kPack[k][1])));
  return out;
}
inline RVec<float> svMass(const Result& r) {
  RVec<float> out;
  for (const auto& s : r.svs) out.push_back(static_cast<float>(s.mass));
  return out;
}
inline RVec<int> svNtracks(const Result& r) {
  RVec<int> out;
  for (const auto& s : r.svs) out.push_back(static_cast<int>(s.trk.size()));
  return out;
}
// per constituent track, flat over the SVs
inline RVec<int> svTrkSV(const Result& r) {
  RVec<int> out;
  for (size_t k = 0; k < r.svs.size(); ++k)
    for (size_t j = 0; j < r.svs[k].trk.size(); ++j) out.push_back(static_cast<int>(k));
  return out;
}
inline RVec<int> svTrkOrigIdx(const Result& r, const Event& ev) {
  RVec<int> out;
  for (const auto& s : r.svs)
    for (int t : s.trk) out.push_back(detail::origOf(ev, t));
  return out;
}

// The vertices as the existing SV getters read them (float positions and
// covariances, momentum = the summed track momenta).
inline VertexingUtils::FCCAnalysesVertex pvAsFCC(const Result& r) {
  VertexingUtils::FCCAnalysesVertex out;
  edm4hep::VertexData vd;
  vd.position = edm4hep::Vector3f(float(r.pv[0]), float(r.pv[1]), float(r.pv[2]));
  std::array<float, 6> cm{};
  for (int k = 0; k < 6; ++k) cm[k] = float(r.pvC(detail::kPack[k][0], detail::kPack[k][1]));
  vd.covMatrix = cm;
  out.vertex = vd;
  out.ntracks = static_cast<int>(std::count(r.assign.begin(), r.assign.end(), 0));
  out.mc_ind = -1;
  return out;
}
inline VertexingUtils::FCCAnalysesV0 svsAsFCC(const Result& r) {
  VertexingUtils::FCCAnalysesV0 out;
  for (const auto& s : r.svs) {
    VertexingUtils::FCCAnalysesVertex v;
    edm4hep::VertexData vd;
    vd.position = edm4hep::Vector3f(float(s.x[0]), float(s.x[1]), float(s.x[2]));
    std::array<float, 6> cm{};
    for (int k = 0; k < 6; ++k) cm[k] = float(s.C(detail::kPack[k][0], detail::kPack[k][1]));
    vd.covMatrix = cm;
    v.vertex = vd;
    v.ntracks = static_cast<int>(s.trk.size());
    v.mc_ind = -1;
    for (int t : s.trk) v.reco_ind.push_back(t);
    v.updated_track_momentum_at_vertex.push_back(TVector3(s.P[0], s.P[1], s.P[2]));
    out.vtx.push_back(v);
    out.pdgAbs.push_back(0);
    out.invM.push_back(s.mass);
  }
  return out;
}

// ---- per V0
inline RVec<int> v0TrkOrigIdx(const Result& r, const Event& ev, int leg) {
  RVec<int> out;
  for (const auto& v : r.v0s) out.push_back(detail::origOf(ev, leg == 0 ? v.t1 : v.t2));
  return out;
}
inline RVec<int> v0Pdg(const Result& r) {
  RVec<int> out;
  for (const auto& v : r.v0s) out.push_back(v.pdg);
  return out;
}
inline RVec<int> v0Hyp(const Result& r) {
  RVec<int> out;
  for (const auto& v : r.v0s) out.push_back(v.hyp);
  return out;
}
inline RVec<float> v0Pos(const Result& r, int i) {
  RVec<float> out;
  for (const auto& v : r.v0s) out.push_back(static_cast<float>(v.x[i]));
  return out;
}
inline RVec<float> v0Cov(const Result& r, int k) {
  RVec<float> out;
  for (const auto& v : r.v0s) out.push_back(static_cast<float>(v.C(detail::kPack[k][0], detail::kPack[k][1])));
  return out;
}
inline RVec<float> v0Mass(const Result& r) {
  RVec<float> out;
  for (const auto& v : r.v0s) out.push_back(static_cast<float>(v.m));
  return out;
}
inline RVec<float> v0Weight(const Result& r) {
  RVec<float> out;
  for (const auto& v : r.v0s) out.push_back(static_cast<float>(v.w));
  return out;
}
inline RVec<float> v0P(const Result& r) {
  RVec<float> out;
  for (const auto& v : r.v0s) out.push_back(static_cast<float>(v.p));
  return out;
}
inline RVec<float> v0Pcomp(const Result& r, int i) {
  RVec<float> out;
  for (const auto& v : r.v0s) out.push_back(static_cast<float>(v.P[i]));
  return out;
}

// ---- per track, over Tracks
inline RVec<int> trkAssign(const Result& r, const Event& ev) { return detail::overTracks(ev, r.assign, JV_FILL_INT); }
inline RVec<int> trkV0(const Result& r, const Event& ev) { return detail::overTracks(ev, r.tv0, JV_FILL_INT); }
inline RVec<float> trkWeight(const Result& r, const Event& ev) { return detail::overTracks(ev, r.wass, JV_FILL_FLOAT); }
inline RVec<float> trkWeightPV(const Result& r, const Event& ev) { return detail::overTracks(ev, r.wpv, JV_FILL_FLOAT); }
inline RVec<float> trkIPSig(const Result& r, const Event& ev) { return detail::overTracks(ev, r.s3, JV_FILL_FLOAT); }

}  // namespace AlephJointVtx
}  // namespace FCCAnalyses

#endif
