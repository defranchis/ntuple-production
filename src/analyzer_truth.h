#ifndef ALEPHTRUTH_H
#define ALEPHTRUTH_H

/*
  Truth-matching utilities for V0 and secondary-vertex studies.
  Decay graphs are recovered geometrically (the podio parent/daughter links are
  empty): daughters are trackable MC particles produced at the mother endpoint
  within tol. MC positions are in CM, not mm.
*/

#include <cmath>
#include <map>
#include <set>
#include <utility>

#include <ROOT/RVec.hxx>
#include "TVector3.h"

#include "edm4hep/MCParticleData.h"
#include "edm4hep/TrackData.h"
#include "edm4hep/TrackState.h"
#include "podio/ObjectID.h"

#include "FCCAnalyses/VertexingUtils.h"
#include "FCCAnalyses/VertexFinderLCFIPlus.h"

namespace FCCAnalyses {
namespace AlephTruth {

using ROOT::VecOps::RVec;

// Many-to-many track<->MC maps from the trackMCLink collection. Link
// convention: _trackMCLink_from indexes Tracks, _trackMCLink_to indexes
// MCParticles.

inline RVec<RVec<int>> buildMCToTracks(size_t n_mc,
                                       const RVec<podio::ObjectID>& link_from,
                                       const RVec<podio::ObjectID>& link_to) {
  if (link_from.size() != link_to.size())
    throw std::runtime_error("AlephTruth::buildMCToTracks: link size mismatch");
  RVec<RVec<int>> out(n_mc);
  for (size_t i = 0; i < link_to.size(); ++i) {
    int mc_idx = link_to[i].index;
    int trk_idx = link_from[i].index;
    if (mc_idx < 0 || mc_idx >= (int)n_mc) continue;
    out[mc_idx].push_back(trk_idx);
  }
  return out;
}

inline RVec<RVec<int>> buildTrackToMCs(size_t n_tracks,
                                       const RVec<podio::ObjectID>& link_from,
                                       const RVec<podio::ObjectID>& link_to) {
  if (link_from.size() != link_to.size())
    throw std::runtime_error("AlephTruth::buildTrackToMCs: link size mismatch");
  RVec<RVec<int>> out(n_tracks);
  for (size_t i = 0; i < link_from.size(); ++i) {
    int trk_idx = link_from[i].index;
    int mc_idx = link_to[i].index;
    if (trk_idx < 0 || trk_idx >= (int)n_tracks) continue;
    out[trk_idx].push_back(mc_idx);
  }
  return out;
}

// Mother-anchored true V0 finding (Ks -> pi+pi-, Lambda -> p pi).

struct TrueV0s {
  RVec<int>   pdg;        // signed: 310, +-3122
  RVec<int>   mother_idx; // index in MCParticles
  RVec<int>   dau1;       // MC index; for Lambda dau1 = baryon (p/pbar)
  RVec<int>   dau2;
  RVec<float> p;          // mother |p| [GeV]
  RVec<float> costheta;   // mother pz/|p|
  RVec<float> px;         // mother momentum components [GeV]
  RVec<float> py;
  RVec<float> pz;
  RVec<float> fd;         // flight distance |endpoint - vertex| of mother [cm]
  RVec<float> dpv;        // decay point distance from gen PV (MCParticles[0].vertex) [cm]
  RVec<float> vx;         // decay point components [cm]
  RVec<float> vy;
  RVec<float> vz;
  RVec<int>   nmatched;   // # daughters with >=1 linked track (0/1/2)
};

// A particle Geant actually tracked: Geant secondary (status 0) or generator
// final-state handed to Geant (10000+N); excludes generator-decayed
// intermediates (1x0000).
inline bool isTrackable(int genStatus) {
  return genStatus == 0 || (genStatus >= 10000 && genStatus < 100000);
}

inline TrueV0s findTrueV0s(const RVec<edm4hep::MCParticleData>& mc,
                           const RVec<RVec<int>>& mcToTracks,
                           double tol = 1e-4) {
  TrueV0s out;
  if (mc.empty()) return out;
  TVector3 genPV(mc[0].vertex.x, mc[0].vertex.y, mc[0].vertex.z);

  // charged particles Geant tracked, indexed once
  std::vector<int> charged;
  charged.reserve(64);
  for (size_t k = 0; k < mc.size(); ++k)
    if (mc[k].charge != 0 && isTrackable(mc[k].generatorStatus))
      charged.push_back(k);

  for (size_t i = 0; i < mc.size(); ++i) {
    int apdg = std::abs(mc[i].PDG);
    if (apdg != 310 && apdg != 3122) continue;

    TVector3 ep(mc[i].endpoint.x, mc[i].endpoint.y, mc[i].endpoint.z);

    // collect charged daughters at the endpoint
    std::vector<int> daus;
    for (int k : charged) {
      if ((int)i == k) continue;
      double dx = mc[k].vertex.x - ep.X();
      double dy = mc[k].vertex.y - ep.Y();
      double dz = mc[k].vertex.z - ep.Z();
      if (dx * dx + dy * dy + dz * dz < tol * tol) daus.push_back(k);
    }
    if (daus.size() != 2) continue;

    const auto& d1 = mc[daus[0]];
    const auto& d2 = mc[daus[1]];
    if (d1.charge * d2.charge >= 0) continue;

    int a1 = std::abs(d1.PDG), a2 = std::abs(d2.PDG);
    int dau_b = -1, dau_m = -1; // baryon (or first pion), meson
    int pdg_signed = 0;
    if (apdg == 310) {
      if (!(a1 == 211 && a2 == 211)) continue;
      pdg_signed = 310;
      dau_b = daus[0];
      dau_m = daus[1];
    } else {
      // Lambda -> p pi- (pdg 3122), Lambdabar -> pbar pi+ (pdg -3122)
      int ip = (a1 == 2212) ? 0 : (a2 == 2212 ? 1 : -1);
      if (ip < 0) continue;
      int im = 1 - ip;
      int am = std::abs(mc[daus[im]].PDG);
      if (am != 211) continue;
      dau_b = daus[ip];
      dau_m = daus[im];
      pdg_signed = (mc[dau_b].PDG > 0) ? 3122 : -3122;
    }

    // dedupe duplicated mother records decaying at the same point
    bool dup = false;
    for (size_t j = 0; j < out.pdg.size(); ++j) {
      if (std::abs(out.pdg[j]) != apdg) continue;
      int m = out.mother_idx[j];
      double dx = mc[m].endpoint.x - ep.X();
      double dy = mc[m].endpoint.y - ep.Y();
      double dz = mc[m].endpoint.z - ep.Z();
      if (dx * dx + dy * dy + dz * dz < tol * tol) { dup = true; break; }
    }
    if (dup) continue;

    TVector3 pv3(mc[i].momentum.x, mc[i].momentum.y, mc[i].momentum.z);
    TVector3 prod(mc[i].vertex.x, mc[i].vertex.y, mc[i].vertex.z);

    int nm = 0;
    if (!mcToTracks[dau_b].empty()) ++nm;
    if (!mcToTracks[dau_m].empty()) ++nm;

    out.pdg.push_back(pdg_signed);
    out.mother_idx.push_back(i);
    out.dau1.push_back(dau_b);
    out.dau2.push_back(dau_m);
    out.p.push_back(pv3.Mag());
    out.costheta.push_back(pv3.Mag() > 0 ? pv3.Z() / pv3.Mag() : 0.);
    out.px.push_back(pv3.X());
    out.py.push_back(pv3.Y());
    out.pz.push_back(pv3.Z());
    out.fd.push_back((ep - prod).Mag());
    out.dpv.push_back((ep - genPV).Mag());
    out.vx.push_back(ep.X());
    out.vy.push_back(ep.Y());
    out.vz.push_back(ep.Z());
    out.nmatched.push_back(nm);
  }
  return out;
}

// Index recovery: V0 candidates hold reco_ind into SecondaryTracks_looseBS;
// map those entries back to original Tracks indices.

// Mirror of AlephSelection::select_tracks_baseline returning, for each
// selected entry (same order), the ORIGINAL index in the Tracks collection.
inline RVec<int> selectedBaselineOriginalIndices(
    const RVec<edm4hep::TrackData>& tracks_in,
    const RVec<edm4hep::TrackState>& trackstates_in,
    const RVec<edm4hep::TrackState>& selected_states_check) {
  RVec<int> out;
  for (size_t it = 0; it < tracks_in.size(); ++it) {
    const auto& track = tracks_in[it];
    if (track.ndf == 0) continue;
    if (track.chi2 / track.ndf > 10.) continue;
    auto n_trackstates = track.trackStates_end - track.trackStates_begin;
    if (n_trackstates != 1)
      throw std::runtime_error("AlephTruth: expected exactly one TrackState per Track");
    const auto& ts = trackstates_in[track.trackStates_begin];
    const auto& cov = ts.covMatrix;
    if (cov[0] <= 1e-12 || cov[2] <= 1e-12 || cov[9] <= 1e-12) continue;
    if (!std::isfinite(cov[0]) || !std::isfinite(cov[2]) || !std::isfinite(cov[9])) continue;
    out.push_back(it);
  }
  if (out.size() != selected_states_check.size())
    throw std::runtime_error("AlephTruth: baseline-selection mirror out of sync with analyzer.h");
  return out;
}

// For each entry of SecondaryTracks (flipped param space), find its position in
// the flipped selected-baseline list and return the original Tracks index.
inline RVec<int> secondaryToOriginalTrack(
    const RVec<edm4hep::TrackState>& secondaries,
    const RVec<edm4hep::TrackState>& selected_flipped,
    const RVec<int>& selected_orig_idx) {
  RVec<int> out;
  out.reserve(secondaries.size());
  for (const auto& s : secondaries) {
    int found = -1;
    for (size_t k = 0; k < selected_flipped.size(); ++k) {
      if (VertexingUtils::compare_Tracks(s, selected_flipped[k])) { found = selected_orig_idx[k]; break; }
    }
    if (found < 0)
      throw std::runtime_error("AlephTruth: secondary track not found in selected baseline");
    out.push_back(found);
  }
  return out;
}

// Pair-index recovery for V0 candidates: the compiled get_V0s leaves
// FCCAnalysesVertex.reco_ind empty, so replicate its booking loop (same
// windows and exclusivity). classifyV0s throws on any pdg/invM mismatch.

struct V0Pairs {
  RVec<int>    i1, i2;   // indices into the secondaries collection, booking order
  RVec<int>    pdgAbs;
  RVec<double> invM;
};

inline V0Pairs rerunV0Pairing(const RVec<edm4hep::TrackState>& np_tracks,
                              const VertexingUtils::FCCAnalysesVertex& PV,
                              double solenoidBz,
                              double chi2_cut = 10.) {
  V0Pairs out;
  const int nTr = np_tracks.size();
  if (nTr < 2) return out;
  // loose windows from get_V0s_ALEPH (Gamma invM_high=-1 -> never booked)
  const double Ks_lo = 0.1, Ks_hi = 1.4, Ks_dis = 0.1, Ks_cos = 0.999;
  const double L_lo = 0.1, L_hi = 1.4, L_dis = 0.1, L_cos = 0.999;
  const double G_hi = -1., G_dis = 0.9, G_cos = 0.999;

  RVec<bool> isInV0(nTr, false);
  RVec<edm4hep::TrackState> tr_pair(2);
  VertexingUtils::FCCAnalysesVertex V0_vtx;

  auto book = [&](int i, int j, int pdg, double m) {
    isInV0[i] = true;
    isInV0[j] = true;
    out.i1.push_back(i);
    out.i2.push_back(j);
    out.pdgAbs.push_back(pdg);
    out.invM.push_back(m);
  };

  for (int i = 0; i < nTr - 1; ++i) {
    if (isInV0[i]) continue; // exclusive_tracks
    tr_pair[0] = np_tracks[i];
    for (int j = i + 1; j < nTr; ++j) {
      if (isInV0[j]) continue;
      if (tr_pair[0].omega * np_tracks[j].omega > 0) continue; // same charge
      tr_pair[1] = np_tracks[j];
      RVec<double> cand = VertexFinderLCFIPlus::get_V0candidate(
          V0_vtx, tr_pair, PV, true, chi2_cut, solenoidBz);
      if (cand[0] == -1) continue;
      if (cand[0] > Ks_lo && cand[0] < Ks_hi && cand[4] > Ks_dis && cand[5] > Ks_cos)
        book(i, j, 310, cand[0]);
      if (cand[1] > L_lo && cand[1] < L_hi && cand[4] > L_dis && cand[5] > L_cos)
        book(i, j, 3122, cand[1]);
      if (cand[2] > L_lo && cand[2] < L_hi && cand[4] > L_dis && cand[5] > L_cos)
        book(i, j, 3122, cand[2]);
      if (cand[3] < G_hi && cand[4] > G_dis && cand[5] > G_cos)
        book(i, j, 22, cand[3]);
    }
  }
  return out;
}

// Per true V0: how many of its daughters (0-2) have >=1 linked track that
// survives into the secondary-track set.
inline RVec<int> daughtersInSecondaries(const TrueV0s& tv,
                                        const RVec<RVec<int>>& mcToTracks,
                                        const RVec<int>& sec2orig) {
  // membership flag indexed by original-track index (bounds-checked on read)
  int nmax = -1;
  for (int o : sec2orig) nmax = std::max(nmax, o);
  std::vector<char> inSec(nmax + 1, 0);
  for (int o : sec2orig)
    if (o >= 0) inSec[o] = 1;
  RVec<int> out;
  for (size_t k = 0; k < tv.pdg.size(); ++k) {
    int n = 0;
    for (int dau : {tv.dau1[k], tv.dau2[k]}) {
      for (int t : mcToTracks[dau])
        if (t >= 0 && t < (int)inSec.size() && inSec[t]) { ++n; break; }
    }
    out.push_back(n);
  }
  return out;
}

// For V0 collections whose vertices carry a filled reco_ind (e.g. the
// AlephV0New module): build the pair association directly, no replica needed.
inline V0Pairs pairsFromRecoInd(const VertexingUtils::FCCAnalysesV0& v0s) {
  V0Pairs out;
  for (size_t c = 0; c < v0s.vtx.size(); ++c) {
    const auto& ri = v0s.vtx[c].reco_ind;
    if (ri.size() != 2)
      throw std::runtime_error("AlephTruth::pairsFromRecoInd: reco_ind not size 2");
    out.i1.push_back(ri[0]);
    out.i2.push_back(ri[1]);
    out.pdgAbs.push_back(v0s.pdgAbs[c]);
    out.invM.push_back(v0s.invM[c]);
  }
  return out;
}

// Truth classification of reco V0 candidates (event-level V0s_event order).

struct V0TruthInfo {
  RVec<int>   cls;          // 0 combinatorial, 1 true Ks, 2 true Lambda,
                            // 3 gamma-conversion, 4 half-match (one true daughter + other track)
  RVec<int>   true_idx;     // index into TrueV0s arrays, -1 if none
  RVec<int>   pair_mult;    // # candidates built from this exact track pair (multi-hypothesis booking)
  RVec<int>   track_shared; // 1 if either track also appears in ANOTHER candidate pair
  RVec<float> alpha;        // Armenteros-Podolanski longitudinal asymmetry
  RVec<float> qt;           // Armenteros-Podolanski transverse momentum [GeV]
  RVec<int>   trk1;         // original Tracks indices of the pair
  RVec<int>   trk2;
};

inline bool contains(const RVec<int>& v, int x) {
  for (int e : v) if (e == x) return true;
  return false;
}

inline V0TruthInfo classifyV0s(const VertexingUtils::FCCAnalysesV0& v0s,
                               const V0Pairs& vp,
                               const RVec<edm4hep::TrackState>& secondaries,
                               const RVec<int>& sec2orig,
                               const RVec<RVec<int>>& trackToMCs,
                               const RVec<edm4hep::MCParticleData>& mc,
                               const TrueV0s& tv,
                               double tol = 1e-4) {
  V0TruthInfo out;
  size_t n = v0s.vtx.size();

  // cross-check: the replicated booking must reproduce the compiled candidates
  if (vp.pdgAbs.size() != n)
    throw std::runtime_error("AlephTruth: rerunV0Pairing candidate count mismatch");
  for (size_t c = 0; c < n; ++c)
    if (vp.pdgAbs[c] != v0s.pdgAbs[c] || vp.invM[c] != v0s.invM[c])
      throw std::runtime_error("AlephTruth: rerunV0Pairing pdg/invM mismatch");

  // original track indices per candidate
  std::vector<std::pair<int, int>> pairs(n, {-1, -1});
  std::map<std::pair<int, int>, int> pair_count;
  std::map<int, int> track_count; // distinct pairs a track appears in
  for (size_t c = 0; c < n; ++c) {
    int t1 = sec2orig.at(vp.i1[c]);
    int t2 = sec2orig.at(vp.i2[c]);
    auto key = std::minmax(t1, t2);
    pairs[c] = {t1, t2};
    pair_count[key]++;
  }
  // track sharing across DIFFERENT pairs
  for (auto& pc : pair_count) {
    track_count[pc.first.first]++;
    track_count[pc.first.second]++;
  }

  for (size_t c = 0; c < n; ++c) {
    int t1 = pairs[c].first, t2 = pairs[c].second;
    const RVec<int>& M1 = trackToMCs.at(t1);
    const RVec<int>& M2 = trackToMCs.at(t2);

    int cls = 0, true_idx = -1;

    // true V0: both daughters matched (either track order)
    for (size_t k = 0; k < tv.pdg.size(); ++k) {
      bool direct = contains(M1, tv.dau1[k]) && contains(M2, tv.dau2[k]);
      bool swapped = contains(M1, tv.dau2[k]) && contains(M2, tv.dau1[k]);
      if (direct || swapped) {
        cls = (std::abs(tv.pdg[k]) == 310) ? 1 : 2;
        true_idx = k;
        break;
      }
    }

    // gamma conversion: e+e- pair from a common displaced point
    if (cls == 0) {
      for (int m1 : M1) {
        if (cls) break;
        if (std::abs(mc[m1].PDG) != 11) continue;
        for (int m2 : M2) {
          if (std::abs(mc[m2].PDG) != 11) continue;
          if (mc[m1].charge * mc[m2].charge >= 0) continue;
          double dx = mc[m1].vertex.x - mc[m2].vertex.x;
          double dy = mc[m1].vertex.y - mc[m2].vertex.y;
          double dz = mc[m1].vertex.z - mc[m2].vertex.z;
          if (dx * dx + dy * dy + dz * dz < tol * tol) { cls = 3; break; }
        }
      }
    }

    // half-match: exactly one track is a true-V0 daughter
    if (cls == 0) {
      for (size_t k = 0; k < tv.pdg.size(); ++k) {
        bool one = contains(M1, tv.dau1[k]) || contains(M1, tv.dau2[k]);
        bool two = contains(M2, tv.dau1[k]) || contains(M2, tv.dau2[k]);
        if (one != two) { cls = 4; true_idx = k; break; }
      }
    }

    // Armenteros-Podolanski from updated momenta at the fitted vertex, in
    // fitted pair order (i1,i2). SecondaryTracks are flipD0_copy'ed (raw ALEPH
    // omega carries -charge), so physical charge = +sign(omega).
    float alpha = -99., qt = -99.;
    const auto& upd = v0s.vtx[c].updated_track_momentum_at_vertex;
    if (upd.size() == 2) {
      TVector3 pa = upd[0], pb = upd[1];
      TVector3 ptot = pa + pb;
      if (ptot.Mag() > 0) {
        double la = pa.Dot(ptot) / ptot.Mag();
        double lb = pb.Dot(ptot) / ptot.Mag();
        double q_a = (secondaries.at(vp.i1[c]).omega > 0) ? 1. : -1.;
        double lplus = (q_a > 0) ? la : lb;
        double lminus = (q_a > 0) ? lb : la;
        alpha = (lplus + lminus != 0.) ? (lplus - lminus) / (lplus + lminus) : -99.;
        qt = pa.Cross(ptot.Unit()).Mag();
      }
    }

    out.cls.push_back(cls);
    out.true_idx.push_back(true_idx);
    auto key = std::minmax(t1, t2);
    out.pair_mult.push_back(pair_count[key]);
    out.track_shared.push_back((track_count[t1] > 1 || track_count[t2] > 1) ? 1 : 0);
    out.alpha.push_back(alpha);
    out.qt.push_back(qt);
    out.trk1.push_back(t1);
    out.trk2.push_back(t2);
  }
  return out;
}

// Per true V0: found by get_V0s at all / found under the correct hypothesis.
inline RVec<int> trueV0FoundAny(const TrueV0s& tv, const V0TruthInfo& info) {
  RVec<int> out(tv.pdg.size(), 0);
  for (size_t c = 0; c < info.cls.size(); ++c)
    if ((info.cls[c] == 1 || info.cls[c] == 2) && info.true_idx[c] >= 0)
      out[info.true_idx[c]] = 1;
  return out;
}

inline RVec<int> trueV0FoundCorrect(const TrueV0s& tv, const V0TruthInfo& info,
                                    const VertexingUtils::FCCAnalysesV0& v0s) {
  RVec<int> out(tv.pdg.size(), 0);
  for (size_t c = 0; c < info.cls.size(); ++c) {
    if (info.true_idx[c] < 0) continue;
    if (info.cls[c] != 1 && info.cls[c] != 2) continue;
    if (v0s.pdgAbs[c] == std::abs(tv.pdg[info.true_idx[c]]))
      out[info.true_idx[c]] = 1;
  }
  return out;
}

// Simple event-order candidate kinematics (independent of jet assignment).
inline RVec<float> candDxyz(const VertexingUtils::FCCAnalysesV0& v0s,
                            const VertexingUtils::FCCAnalysesVertex& PV) {
  RVec<float> out;
  TVector3 pv(PV.vertex.position[0], PV.vertex.position[1], PV.vertex.position[2]);
  for (const auto& v : v0s.vtx) {
    TVector3 x(v.vertex.position[0], v.vertex.position[1], v.vertex.position[2]);
    out.push_back((x - pv).Mag());
  }
  return out;
}

// Fitted-vertex position component, axis 0/1/2 = x/y/z, in cm.
inline RVec<float> candVtxPos(const VertexingUtils::FCCAnalysesV0& v0s, int axis) {
  RVec<float> out;
  for (const auto& v : v0s.vtx) out.push_back(v.vertex.position[axis]);
  return out;
}

inline RVec<float> candChi2(const VertexingUtils::FCCAnalysesV0& v0s) {
  RVec<float> out;
  for (const auto& v : v0s.vtx) out.push_back(v.vertex.chi2);
  return out;
}

inline RVec<float> candP(const VertexingUtils::FCCAnalysesV0& v0s) {
  RVec<float> out;
  for (const auto& v : v0s.vtx) {
    TVector3 p(0., 0., 0.);
    for (const auto& tp : v.updated_track_momentum_at_vertex) p += tp;
    out.push_back(p.Mag());
  }
  return out;
}

// Component comp (0/1/2 = x/y/z) of the SAME summed vertex momentum whose
// magnitude candP returns, so sqrt(px^2+py^2+pz^2) reproduces candP exactly.
inline RVec<float> candPcomp(const VertexingUtils::FCCAnalysesV0& v0s, int comp) {
  RVec<float> out;
  for (const auto& v : v0s.vtx) {
    TVector3 p(0., 0., 0.);
    for (const auto& tp : v.updated_track_momentum_at_vertex) p += tp;
    out.push_back(p[comp]);
  }
  return out;
}

inline RVec<float> candCosPointing(const VertexingUtils::FCCAnalysesV0& v0s,
                                   const VertexingUtils::FCCAnalysesVertex& PV) {
  RVec<float> out;
  TVector3 pv(PV.vertex.position[0], PV.vertex.position[1], PV.vertex.position[2]);
  for (const auto& v : v0s.vtx) {
    TVector3 x(v.vertex.position[0], v.vertex.position[1], v.vertex.position[2]);
    TVector3 p(0., 0., 0.);
    for (const auto& tp : v.updated_track_momentum_at_vertex) p += tp;
    TVector3 d = x - pv;
    out.push_back((d.Mag() > 0 && p.Mag() > 0) ? d.Dot(p) / (d.Mag() * p.Mag()) : -2.);
  }
  return out;
}

} // namespace AlephTruth
} // namespace FCCAnalyses

#endif
