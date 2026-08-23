#ifndef ALEPH_TRKAUX_H
#define ALEPH_TRKAUX_H

/*
  Track-level auxiliaries, vertex-fit glue and LUND-truth helpers shared by the
  V0, SV, phi->KK and D* finders. Lengths cm, momenta GeV; the per-track
  quantities are aligned with a selected trackstate collection through its
  original-Tracks index map.
*/

#include <algorithm>
#include <cmath>
#include <vector>

#include <ROOT/RVec.hxx>
#include "TVector3.h"

#include "edm4hep/TrackState.h"
#include "edm4hep/MCParticleData.h"
#include "edm4hep/ParticleIDData.h"
#include "FCCAnalyses/VertexingUtils.h"
#include "FCCAnalyses/VertexFitterSimple.h"

#include "aleph_units.h"

namespace FCCAnalyses {
namespace AlephTrkAux {

using ROOT::VecOps::RVec;

// Momentum at the perigee (pT = kPtPerTeslaCm*Bz/|omega|); used for pre-fit
// mass windows and for daughters that enter no fit.
inline TVector3 perigeeMomentum(const edm4hep::TrackState& t, double Bz) {
  double om = std::abs(t.omega);
  if (om <= 0.) return TVector3(0., 0., 0.);
  double pt = AlephUnits::kPtPerTeslaCm * Bz / om;
  return TVector3(pt * std::cos(t.phi), pt * std::sin(t.phi), pt * t.tanLambda);
}

// The one vertex-fit entry point of every finder: fit `group` (alltracks passed
// so reco_ind is filled), then undo the cm-as-mm homothety of the momenta,
// which come out 10x too small, ONCE at the source.
inline VertexingUtils::FCCAnalysesVertex fitTracksCm(
    const RVec<edm4hep::TrackState>& group,
    const RVec<edm4hep::TrackState>& alltracks, double Bz) {
  auto v = VertexFitterSimple::VertexFitter_Tk(
      0, group, alltracks, false, 0., 0., 0., 0., 0., 0., Bz, false);
  for (auto& tp : v.updated_track_momentum_at_vertex) tp *= 10.;
  return v;
}

// Armenteros-Podolanski variables of a track pair: qt and
// alpha = (pL+ - pL-)/(pL+ + pL-). q1sign = physical charge of p1; for a
// same-charge pair the labels are conventional, so pass +1 to order by p1.
// alpha_null is returned when the longitudinal momenta cancel.
inline void apVars(const TVector3& p1, const TVector3& p2, double q1sign,
                   double& alpha, double& qt, double alpha_null = 0.) {
  TVector3 p = p1 + p2;
  double pmag = p.Mag();
  qt = p1.Cross(p.Unit()).Mag();
  double la = p1.Dot(p) / pmag, lb = p2.Dot(p) / pmag;
  double lplus = (q1sign > 0) ? la : lb, lminus = (q1sign > 0) ? lb : la;
  alpha = (lplus + lminus != 0.) ? (lplus - lminus) / (lplus + lminus)
                                 : alpha_null;
}

// Sum of two packed lower-triangular position covariances (xx, yx, yy, zx, zy,
// zz) as a full 3x3 matrix.
template <typename CovA, typename CovB>
inline void sumCovPacked(const CovA& ca, const CovB& cb, double C[3][3]) {
  C[0][0] = double(ca[0]) + cb[0];
  C[0][1] = C[1][0] = double(ca[1]) + cb[1];
  C[1][1] = double(ca[2]) + cb[2];
  C[0][2] = C[2][0] = double(ca[3]) + cb[3];
  C[1][2] = C[2][1] = double(ca[4]) + cb[4];
  C[2][2] = double(ca[5]) + cb[5];
}

// Membership table indexed by original-Tracks index: the per-event
// replacement for repeated std::find scans over an index list.
inline std::vector<char> memberMask(const RVec<int>& set_orig) {
  int mx = -1;
  for (int o : set_orig) mx = std::max(mx, o);
  std::vector<char> m((size_t)(mx + 1), 0);
  for (int o : set_orig)
    if (o >= 0) m[o] = 1;
  return m;
}

inline bool inMask(const std::vector<char>& m, int o) {
  return o >= 0 && (size_t)o < m.size() && m[o] != 0;
}

// 3D compatibility significance of two vertex positions: sqrt of the chi2 of
// d = x1 - x2 under the summed position covariances (packed lower triangle
// xx, yx, yy, zx, zy, zz). -1 when the summed covariance is not invertible.
template <typename CovA, typename CovB>
inline float vertexDistSig(const TVector3& d, const CovA& ca, const CovB& cb) {
  double C[3][3];
  sumCovPacked(ca, cb, C);
  double det = C[0][0] * (C[1][1] * C[2][2] - C[1][2] * C[2][1])
             - C[0][1] * (C[1][0] * C[2][2] - C[1][2] * C[2][0])
             + C[0][2] * (C[1][0] * C[2][1] - C[1][1] * C[2][0]);
  if (!(std::abs(det) > 0.) || !std::isfinite(det)) return -1.;
  double inv[3][3];
  inv[0][0] = (C[1][1] * C[2][2] - C[1][2] * C[2][1]) / det;
  inv[0][1] = (C[0][2] * C[2][1] - C[0][1] * C[2][2]) / det;
  inv[0][2] = (C[0][1] * C[1][2] - C[0][2] * C[1][1]) / det;
  inv[1][0] = (C[1][2] * C[2][0] - C[1][0] * C[2][2]) / det;
  inv[1][1] = (C[0][0] * C[2][2] - C[0][2] * C[2][0]) / det;
  inv[1][2] = (C[0][2] * C[1][0] - C[0][0] * C[1][2]) / det;
  inv[2][0] = (C[1][0] * C[2][1] - C[1][1] * C[2][0]) / det;
  inv[2][1] = (C[0][1] * C[2][0] - C[0][0] * C[2][1]) / det;
  inv[2][2] = (C[0][0] * C[1][1] - C[0][1] * C[1][0]) / det;
  double dv[3] = {d.X(), d.Y(), d.Z()};
  double s2 = 0.;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) s2 += dv[i] * inv[i][j] * dv[j];
  return (s2 > 0. && std::isfinite(s2)) ? std::sqrt(s2) : -1.;
}

// Component `which` (0 = VDET, 1 = ITC, 2 = TPC) of the per-track
// subdetectorHitNumbers block. -1 when the block is missing or too short.
inline RVec<int> subdetHits(const RVec<int>& orig_idx,
                            const RVec<unsigned int>& begin,
                            const RVec<unsigned int>& end,
                            const RVec<int>& values, int which) {
  RVec<int> out;
  for (int o : orig_idx) {
    int v = -1;
    if (o >= 0 && o < (int)begin.size() && o < (int)end.size()) {
      unsigned int b = begin[o], e = end[o];
      if (b + which < e && b + which < values.size()) v = values[b + which];
    }
    out.push_back(v);
  }
  return out;
}

// 1 for each entry whose original-Tracks index appears in `set_orig`; backs
// the "daughter was in the fitted primary set" flag.
inline RVec<int> flagInSet(const RVec<int>& orig_idx, const RVec<int>& set_orig) {
  RVec<int> out;
  const std::vector<char> in = memberMask(set_orig);
  for (int o : orig_idx) out.push_back(inMask(in, o) ? 1 : 0);
  return out;
}

// Original-Tracks indices of the daughters of the selected V0 candidates
// (select with the tight flag); consumed by the pairing veto of the later finders.
inline RVec<int> claimedOrigIdx(const RVec<int>& d1, const RVec<int>& d2,
                                const RVec<int>& keep) {
  RVec<int> out;
  for (size_t i = 0; i < d1.size() && i < d2.size(); ++i) {
    if (i < keep.size() && !keep[i]) continue;
    if (d1[i] >= 0) out.push_back(d1[i]);
    if (d2[i] >= 0) out.push_back(d2[i]);
  }
  return out;
}

inline RVec<float> trackChi2Ndf(const RVec<int>& orig_idx,
                                const RVec<float>& chi2,
                                const RVec<int>& ndf) {
  RVec<float> out;
  for (int o : orig_idx) {
    float v = -1.f;
    if (o >= 0 && o < (int)chi2.size() && o < (int)ndf.size() && ndf[o] != 0)
      v = chi2[o] / float(ndf[o]);
    out.push_back(v);
  }
  return out;
}

// ---------------------------------------------------------------------------
// Particle-flow join: original track index -> ReconstructedParticle -> PF type.
// ---------------------------------------------------------------------------

// PF type code of the ParticleID collection (index-parallel to RecoParticles);
// the same source pfcand_isChargedHad reads.
constexpr int kPFChargedHad = 0;

// ReconstructedParticle index of every original-Tracks index, -1 when no RP
// points at the track. tracks_begin/tracks_end are offsets into the RP->Track
// relation, and begin != end is the "has a track" test; a neutral's begin is
// in range but meaningless.
template <typename Idx>
inline RVec<int> rpIndexByTrack(const RVec<Idx>& tracks_begin,
                                const RVec<Idx>& tracks_end,
                                const RVec<int>& rpTrackIndex, size_t nTracks) {
  RVec<int> out(nTracks, -1);
  const size_t n = std::min(tracks_begin.size(), tracks_end.size());
  for (size_t i = 0; i < n; ++i)
    for (size_t k = tracks_begin[i]; k < (size_t)tracks_end[i]; ++k) {
      if (k >= rpTrackIndex.size()) break;
      const int t = rpTrackIndex[k];
      if (t >= 0 && (size_t)t < nTracks && out[t] < 0) out[t] = (int)i;
    }
  return out;
}

// Tri-state PF charged-hadron flag of a candidate leg: 1 = its linked
// ReconstructedParticle is a PF charged hadron, 0 = the RP is e/mu/other,
// -1 = the track has no linked RP.
inline RVec<int> legIsChargedHad(const RVec<int>& orig_idx,
                                 const RVec<int>& rp_of_track,
                                 const RVec<edm4hep::ParticleIDData>& pid) {
  RVec<int> out;
  for (int o : orig_idx) {
    int v = -1;
    if (o >= 0 && o < (int)rp_of_track.size()) {
      const int r = rp_of_track[o];
      if (r >= 0 && r < (int)pid.size())
        v = (pid[r].type == kPFChargedHad) ? 1 : 0;
    }
    out.push_back(v);
  }
  return out;
}

// ---------------------------------------------------------------------------
// Per-track membership over the ORIGINAL Tracks frame: which sets and stored
// candidates each track belongs to, so the offline 1/n de-duplication weight
// and any set algebra need no re-running of the finders.
// ---------------------------------------------------------------------------
enum TrkMemberBit : int {
  kTrkPV        = 1 << 0,  // in the fitted primary-vertex track set
  kTrkV0        = 1 << 1,  // daughter of any stored V0 candidate (loose tier)
  kTrkV0Tight   = 1 << 2,  // daughter of a tight V0 candidate
  kTrkPhi       = 1 << 3,  // leg of any stored phi->KK candidate
  kTrkPhiWp     = 1 << 4,  // leg of a phi candidate passing the wp flag
  kTrkD0        = 1 << 5,  // leg of any stored D0 candidate
  kTrkDstar     = 1 << 6,  // leg of any stored D* candidate, slow pion included
  kTrkDstarTight= 1 << 7,  // leg of a D* candidate passing the tight flag
  kTrkSV        = 1 << 8,  // track of a masked-mode (svn) secondary vertex
  kTrkBaseline  = 1 << 9   // baseline-selected track
};

struct TrackTags {
  RVec<int> member;  // OR of TrkMemberBit, indexed by original track index
  RVec<int> nCand;   // stored V0 + phi + D0 + D* candidates using the track
};

namespace detail {

inline void tagBit(RVec<int>& m, const RVec<int>& idx, int bit) {
  for (int o : idx)
    if (o >= 0 && o < (int)m.size()) m[o] |= bit;
}

inline void tagBitIf(RVec<int>& m, const RVec<int>& idx, const RVec<int>& flag,
                     int bit) {
  for (size_t k = 0; k < idx.size() && k < flag.size(); ++k)
    if (flag[k] && idx[k] >= 0 && idx[k] < (int)m.size()) m[idx[k]] |= bit;
}

inline void countLegs(RVec<int>& n, const RVec<int>& idx) {
  for (int o : idx)
    if (o >= 0 && o < (int)n.size()) ++n[o];
}

}  // namespace detail

// Every index list is in the ORIGINAL Tracks frame except sv_trk_idx, which is
// in the secondary frame and is walked through sec2orig. A module that did not
// run passes empty lists.
inline TrackTags trackTags(size_t nTracks,
                           const RVec<int>& baseline_orig,
                           const RVec<int>& prim_orig,
                           const RVec<int>& sv_trk_idx,
                           const RVec<int>& sec2orig,
                           const RVec<int>& v0_d1, const RVec<int>& v0_d2,
                           const RVec<int>& v0_tight,
                           const RVec<int>& phi_t1, const RVec<int>& phi_t2,
                           const RVec<int>& phi_wp,
                           const RVec<int>& d0_k, const RVec<int>& d0_pi,
                           const RVec<int>& ds_k, const RVec<int>& ds_pi,
                           const RVec<int>& ds_pis, const RVec<int>& ds_tight) {
  TrackTags out;
  out.member = RVec<int>(nTracks, 0);
  out.nCand = RVec<int>(nTracks, 0);

  detail::tagBit(out.member, baseline_orig, kTrkBaseline);
  detail::tagBit(out.member, prim_orig, kTrkPV);
  for (int s : sv_trk_idx)
    if (s >= 0 && s < (int)sec2orig.size()) {
      const int o = sec2orig[s];
      if (o >= 0 && o < (int)nTracks) out.member[o] |= kTrkSV;
    }

  for (const RVec<int>* leg : {&v0_d1, &v0_d2}) {
    detail::tagBit(out.member, *leg, kTrkV0);
    detail::tagBitIf(out.member, *leg, v0_tight, kTrkV0Tight);
    detail::countLegs(out.nCand, *leg);
  }
  for (const RVec<int>* leg : {&phi_t1, &phi_t2}) {
    detail::tagBit(out.member, *leg, kTrkPhi);
    detail::tagBitIf(out.member, *leg, phi_wp, kTrkPhiWp);
    detail::countLegs(out.nCand, *leg);
  }
  for (const RVec<int>* leg : {&d0_k, &d0_pi}) {
    detail::tagBit(out.member, *leg, kTrkD0);
    detail::countLegs(out.nCand, *leg);
  }
  for (const RVec<int>* leg : {&ds_k, &ds_pi, &ds_pis}) {
    detail::tagBit(out.member, *leg, kTrkDstar);
    detail::tagBitIf(out.member, *leg, ds_tight, kTrkDstarTight);
    detail::countLegs(out.nCand, *leg);
  }
  return out;
}

// ---------------------------------------------------------------------------
// Per-daughter storage block, one instance per candidate role. `p`/`costheta`
// are the momentum actually used for that daughter (at the fitted vertex, or
// at the perigee for a leg that entered no fit). `pool` stays -1 for finders
// that do not stage their tracks.
// ---------------------------------------------------------------------------
struct TrkBlock {
  RVec<int>   origIdx;   // index into the original Tracks collection
  RVec<int>   q;         // physical charge, +sign(omega)
  RVec<float> p, costheta;
  RVec<float> d0, z0;    // perigee parameters [cm]
  RVec<float> sigd0;     // sqrt(cov[0]) [cm]
  RVec<int>   nvdet, nitc;
  RVec<float> chi2ndf;   // track fit chi2/ndf
  RVec<int>   isprim;    // 1 = in the fitted primary set
  RVec<int>   pool;      // 0 = primary set, 1 = secondary set, 2 = neither
};

inline void reserveTrk(TrkBlock& b, size_t n) {
  b.origIdx.reserve(n); b.q.reserve(n); b.p.reserve(n); b.costheta.reserve(n);
  b.d0.reserve(n); b.z0.reserve(n); b.sigd0.reserve(n); b.nvdet.reserve(n);
  b.nitc.reserve(n); b.chi2ndf.reserve(n); b.isprim.reserve(n); b.pool.reserve(n);
}

// Per-entry auxiliaries of one selected trackstate collection; a short or
// absent vector yields the -1 "unknown" sentinel.
struct TrkAux {
  const RVec<int>*   orig_idx;
  const RVec<int>*   nvdet;
  const RVec<int>*   nitc;
  const RVec<float>* chi2ndf;
  const RVec<int>*   isprim;
  const RVec<int>*   pool;

  int ival(const RVec<int>* src, int k) const {
    return (src && k >= 0 && k < (int)src->size()) ? (*src)[k] : -1;
  }
  float fval(const RVec<float>* src, int k) const {
    return (src && k >= 0 && k < (int)src->size()) ? (*src)[k] : -1.f;
  }
};

inline void pushTrk(TrkBlock& b, const RVec<edm4hep::TrackState>& tracks, int k,
                    const TVector3& p, const TrkAux& aux) {
  b.origIdx.push_back(aux.ival(aux.orig_idx, k));
  b.q.push_back((tracks[k].omega > 0) ? 1 : -1);
  b.p.push_back(p.Mag());
  b.costheta.push_back(p.Mag() > 0. ? p.Z() / p.Mag() : -99.);
  b.d0.push_back(tracks[k].D0);
  b.z0.push_back(tracks[k].Z0);
  b.sigd0.push_back(tracks[k].covMatrix[0] > 0. ? std::sqrt(tracks[k].covMatrix[0]) : -1.);
  b.nvdet.push_back(aux.ival(aux.nvdet, k));
  b.nitc.push_back(aux.ival(aux.nitc, k));
  b.chi2ndf.push_back(aux.fval(aux.chi2ndf, k));
  b.isprim.push_back(aux.ival(aux.isprim, k));
  b.pool.push_back(aux.ival(aux.pool, k));
}

// ---------------------------------------------------------------------------
// MC truth. The MCParticles daughter relations are EMPTY here; provenance
// comes from generatorStatus = 10000*KS + line, with `line` the 1-BASED LUND
// line of the MOTHER (MCParticles entry i is LUND line i+1). Geant
// secondaries carry generatorStatus 0 and have no LUND mother.
// ---------------------------------------------------------------------------
inline int lundMotherIdx(int genStatus) {
  return (genStatus > 0) ? (genStatus % 10000) - 1 : -1;
}

// heaviest quark in a hadron PDG code: 5 (b), 4 (c), else 0
inline int heavyFlavour(int pdg) {
  int a = std::abs(pdg);
  if (a >= 1000000000) return 0;  // nucleus
  int nq3 = (a / 10) % 10, nq2 = (a / 100) % 10, nq1 = (a / 1000) % 10;
  int mx = std::max({nq1, nq2, nq3});
  if (mx >= 5) return (mx == 5) ? 5 : 0;  // top-flavoured codes not expected
  return (mx == 4) ? 4 : 0;
}

// LUND mother index of every MCParticles entry, -1 when there is none.
inline std::vector<int> lundMothers(const RVec<edm4hep::MCParticleData>& mc) {
  const int n = mc.size();
  std::vector<int> moth(n, -1);
  for (int i = 0; i < n; ++i) {
    int m = lundMotherIdx(mc[i].generatorStatus);
    moth[i] = (m >= 0 && m < n) ? m : -1;
  }
  return moth;
}

// LUND mother relation inverted once per event; children stay in ascending
// MCParticles order, as a direct scan over the mother array would visit them.
inline std::vector<std::vector<int>> lundChildren(const std::vector<int>& moth) {
  std::vector<std::vector<int>> ch(moth.size());
  for (size_t k = 0; k < moth.size(); ++k)
    if (moth[k] >= 0) ch[moth[k]].push_back((int)k);
  return ch;
}

// 5 if a b hadron is found anywhere in the ancestry, else 4 if a c hadron, else 0.
inline int ancestorFlavour(const std::vector<int>& moth,
                           const RVec<edm4hep::MCParticleData>& mc, int start) {
  int orig = 0;
  for (int a = start, guard = 0; a >= 0 && guard < 200; a = moth[a], ++guard) {
    int hf = heavyFlavour(mc[a].PDG);
    if (hf == 5) return 5;
    if (hf == 4 && orig == 0) orig = 4;  // keep looking: a b ancestor wins
  }
  return orig;
}

// MC particles linked to an original-Tracks index, all links: the relation
// carries several unweighted partners for some tracks, in no particular order.
inline RVec<int> mcLinksOfTrack(const RVec<RVec<int>>& trackToMCs, int trk, int nmc) {
  RVec<int> out;
  if (trk < 0 || trk >= (int)trackToMCs.size()) return out;
  for (int m : trackToMCs[trk])
    if (m >= 0 && m < nmc) out.push_back(m);
  return out;
}

// Truth of one candidate leg: appends the PDG of MC particle `m` (first link
// by default; callers overwrite with the link that matched a true decay) and
// of its LUND mother (0 when absent).
inline void pushLegTruth(const RVec<edm4hep::MCParticleData>& mc, int m,
                         RVec<int>& mcpdg, RVec<int>& mothpdg) {
  const int n = mc.size();
  const int mo = (m >= 0 && m < n) ? lundMotherIdx(mc[m].generatorStatus) : -1;
  mcpdg.push_back((m >= 0 && m < n) ? mc[m].PDG : 0);
  mothpdg.push_back((mo >= 0 && mo < n) ? mc[mo].PDG : 0);
}
// Replaces the entry just pushed for the current leg with MC particle `m`;
// call at most once per leg, after pushLegTruth.
inline void setLegTruth(const RVec<edm4hep::MCParticleData>& mc, int m,
                        RVec<int>& mcpdg, RVec<int>& mothpdg) {
  mcpdg.pop_back(); mothpdg.pop_back();
  pushLegTruth(mc, m, mcpdg, mothpdg);
}

// first link of a leg, -1 when unlinked
inline int firstOr(const RVec<int>& links) { return links.empty() ? -1 : links[0]; }

// 1 for each of the `n_true` true particles reconstructed by at least one
// class-1 candidate carrying `flag` (a loose/tight label; an empty flag counts
// every class-1 candidate). true_idx = per-candidate index into the true list.
inline RVec<int> foundFlags(size_t n_true, const RVec<int>& true_idx,
                            const RVec<int>& cls, const RVec<int>& flag) {
  RVec<int> out(n_true, 0);
  for (size_t k = 0; k < cls.size(); ++k) {
    if (cls[k] != 1) continue;
    if (k < flag.size() && !flag[k]) continue;
    int t = true_idx[k];
    if (t >= 0 && t < (int)out.size()) out[t] = 1;
  }
  return out;
}

}  // namespace AlephTrkAux
}  // namespace FCCAnalyses

#endif  // ALEPH_TRKAUX_H
