#ifndef ALEPHPHIKK_H
#define ALEPHPHIKK_H

/*
  phi(1020) -> K+ K- reconstruction from the FULL baseline-selected track list;
  no dE/dx quantity enters any selection. Input = flipD0_copy'ed trackstates
  (physical charge = +sign(omega)); lengths cm, momenta GeV.
*/

#include <algorithm>
#include <cmath>

#include <ROOT/RVec.hxx>
#include "TVector3.h"

#include "edm4hep/TrackState.h"
#include "edm4hep/MCParticleData.h"
#include "FCCAnalyses/VertexingUtils.h"
#include "FCCAnalyses/VertexFitterSimple.h"

#include "aleph_units.h"
#include "analyzer_trkaux.h"
#include "analyzer_v0new.h"

namespace FCCAnalyses {
namespace AlephPhiKK {

using ROOT::VecOps::RVec;

// shared track auxiliaries and LUND-truth helpers (see analyzer_trkaux.h)
using AlephTrkAux::perigeeMomentum;
using AlephTrkAux::vertexDistSig;
using AlephTrkAux::apVars;
using AlephTrkAux::fitTracksCm;
using AlephTrkAux::TrkBlock;
using AlephTrkAux::TrkAux;
using AlephTrkAux::pushTrk;
using AlephTrkAux::lundMothers;
using AlephTrkAux::lundChildren;
using AlephTrkAux::ancestorFlavour;
using AlephTrkAux::pushLegTruth;
using AlephTrkAux::foundFlags;

constexpr double M_K   = AlephMasses::kK;    // charged kaon
constexpr double M_PHI = AlephMasses::kPhi;
constexpr int    PDG_PHI = 333;
constexpr int    PDG_K   = 321;

// daughter momentum and energy in the phi rest frame
inline double pstarKK() {
  return std::sqrt(0.25 * M_PHI * M_PHI - M_K * M_K);  // 0.1269181 GeV
}
constexpr double estarKK() { return 0.5 * M_PHI; }

// Storage defaults: deliberately loose, the sample is meant to be re-cut
// offline. Near-collinear KK pairs leave the vertex unconstrained along the
// flight direction, hence the wide DPV_FID.
constexpr double M_LO = 0.98, M_HI = 1.10;   // stored KK mass window [GeV]
constexpr double PRE_MARGIN = 0.02;          // pre-fit window margin [GeV]
constexpr double CHI2_CUT = 25.;             // loose sanity only (ndf = 1)
constexpr double DPV_FID = 50.;              // |vtx-PV| storage fiducial [cm]
// applied to the PERIGEE momentum, while the stored trk p is the at-vertex one
constexpr double P_MIN_DEF = 0.3;            // per-track |p| floor [GeV]

// Working points stored as flags, evaluated on the stored quantities (post-fit
// mass, at-vertex daughter momenta). Charge-blind by design: signal =
// wp && !same_sign, control = wp && same_sign under identical cuts.
constexpr double WP_DM = 0.012, TIGHT_DM = 0.005;  // |m - m_phi| [GeV]
constexpr double WP_PDAU = 1.0;                    // daughter |p| [GeV]
constexpr double WP_DPV = 1.0;                     // |vtx-PV| [cm]
constexpr double TIGHT_SIGD0 = 0.01;               // daughter sigma(d0) [cm]

// Armenteros-Podolanski band variable for the equal-mass locus
// (alpha/alpha_max)^2 + (qT/p*)^2, centred on alpha = 0; 1 on the exact
// phi -> K+K- ellipse. For equal masses it reparametrises the mass window.
inline double phiBandEll(double alpha, double qt, double pmag) {
  const double ps = pstarKK();
  double beta = pmag / std::sqrt(pmag * pmag + M_PHI * M_PHI);
  if (beta <= 0.) return -1.;
  double amax = ps / (beta * estarKK());
  return std::sqrt(alpha * alpha / (amax * amax) + qt * qt / (ps * ps));
}

// ---------------------------------------------------------------------------
// Candidate container. One entry per accepted track pair; nothing is claimed
// exclusively, so a track may appear in several candidates.
// ---------------------------------------------------------------------------
struct PhiKKCands {
  RVec<float> invM;       // KK invariant mass at the fitted vertex [GeV]
  RVec<float> p, px, py, pz;  // pair momentum at the fitted vertex [GeV]
  RVec<float> alpha, qt;  // Armenteros-Podolanski (qt in GeV)
  RVec<float> bandEll;    // equal-mass AP band variable (1 = exact locus)
  RVec<float> chi2;       // vertex fit chi2, normalised (ndf = 1)
  RVec<float> vx, vy, vz; // fitted vertex position [cm]
  RVec<float> dpv;        // |vertex - PV| [cm]
  RVec<float> dpvSig;     // 3D significance of the same, -1 if undefined
  RVec<int>   same_sign;  // 1 = same-charge pair (combinatorial control)
  // per-daughter blocks; trk1 is the higher-|p| track of the pair
  TrkBlock    trk1, trk2;
  RVec<int>   wp, tight;  // working-point flags, charge-blind (see the
                          // constants above)
};

// ---------------------------------------------------------------------------
// The finder.
//   tracks     : flipD0_copy'ed baseline-selected trackstates
//   orig_idx   : original-Tracks index of each entry of `tracks`
//   nvdet/nitc : VDET/ITC hit counts of each entry (-1 = unknown)
//   chi2ndf    : track fit chi2/ndf of each entry (-1 = unknown)
//   isprim     : 1 if the entry is in the fitted primary set (stored, never cut on)
//   PV         : fitted primary vertex (positions numerically cm)
//   veto_orig  : original-Tracks indices excluded from the pairing pool (the
//                tight-claimed V0 daughters; empty when the V0 module is off)
// Every selection value is a constant defined above.
// ---------------------------------------------------------------------------
inline PhiKKCands findPhiKK(
    const RVec<edm4hep::TrackState>& tracks,
    const RVec<int>& orig_idx,
    const RVec<int>& nvdet,
    const RVec<int>& nitc,
    const RVec<float>& chi2ndf,
    const RVec<int>& isprim,
    const VertexingUtils::FCCAnalysesVertex& PV,
    double solenoidBz,
    const RVec<int>& veto_orig) {

  PhiKKCands out;
  const int nTr = tracks.size();
  if (nTr < 2) return out;

  const TrkAux aux{&orig_idx, &nvdet, &nitc, &chi2ndf, &isprim, nullptr};
  auto oidx = [&](int k) {
    return (k >= 0 && k < (int)orig_idx.size()) ? orig_idx[k] : -1;
  };

  // momentum prefilter: one pass, then pairs are formed among survivors.
  // eK = kaon-hypothesis energy of the perigee momentum, for the pre-fit mass.
  const std::vector<char> vetoed = AlephTrkAux::memberMask(veto_orig);
  std::vector<int> good;
  good.reserve(nTr);
  std::vector<TVector3> pper(nTr);
  std::vector<double> eK(nTr, 0.);
  for (int i = 0; i < nTr; ++i) {
    const auto& t = tracks[i];
    if (AlephTrkAux::inMask(vetoed, oidx(i))) continue;
    pper[i] = perigeeMomentum(t, solenoidBz);
    if (pper[i].Mag() < P_MIN_DEF) continue;
    eK[i] = std::sqrt(pper[i].Mag2() + M_K * M_K);
    good.push_back(i);
  }
  if (good.size() < 2) return out;

  TVector3 pv(PV.vertex.position[0], PV.vertex.position[1], PV.vertex.position[2]);
  const double pre_lo = M_LO - PRE_MARGIN, pre_hi = M_HI + PRE_MARGIN;

  RVec<edm4hep::TrackState> tr_pair(2);
  for (size_t a = 0; a + 1 < good.size(); ++a) {
    const int i = good[a];
    for (size_t b = a + 1; b < good.size(); ++b) {
      const int j = good[b];
      const bool ss = (tracks[i].omega * tracks[j].omega > 0);

      // pre-fit KK mass from the perigee momenta: removes the bulk of the pair
      // combinatorics before any fit, at the cost of a small tail at the low
      // mass edge (the fit can still move pT by more than the margin).
      TVector3 p_pre = pper[i] + pper[j];
      double e_pre = eK[i] + eK[j];
      double m_pre = std::sqrt(std::max(0., e_pre * e_pre - p_pre.Mag2()));
      if (m_pre < pre_lo || m_pre > pre_hi) continue;

      tr_pair[0] = tracks[i];
      tr_pair[1] = tracks[j];
      auto v = fitTracksCm(tr_pair, tracks, solenoidBz);
      if (v.updated_track_momentum_at_vertex.size() != 2) continue;

      double chi2 = v.vertex.chi2;
      if (!(chi2 == chi2)) continue;
      if (chi2 >= CHI2_CUT) continue;

      TVector3 pa = v.updated_track_momentum_at_vertex[0];
      TVector3 pb = v.updated_track_momentum_at_vertex[1];
      double m = AlephV0New::invMass(pa, M_K, pb, M_K);
      if (m < M_LO || m > M_HI) continue;

      TVector3 p = pa + pb;
      double pmag = p.Mag();
      if (!(pmag > 0.)) continue;
      // physical charge = +sign(omega) for the flipD0_copy'ed collection; a
      // same-sign pair has conventional labels, so it is ordered by track order
      int qi = (tracks[i].omega > 0) ? 1 : -1;
      double alpha, qt;
      apVars(pa, pb, ss ? 1 : qi, alpha, qt);
      double ell = phiBandEll(alpha, qt, pmag);

      TVector3 x(v.vertex.position[0], v.vertex.position[1], v.vertex.position[2]);
      TVector3 d = x - pv;
      double dis = d.Mag();
      if (dis > DPV_FID) continue;
      float dsig = vertexDistSig(d, v.vertex.covMatrix, PV.vertex.covMatrix);

      // daughter 1 = higher-|p| track, so the per-daughter branches have a
      // definite meaning without needing the charge to disambiguate
      int i1 = i, i2 = j;
      TVector3 q1 = pa, q2 = pb;
      if (pb.Mag() > pa.Mag()) { i1 = j; i2 = i; q1 = pb; q2 = pa; }

      out.invM.push_back(m);
      out.p.push_back(pmag);
      out.px.push_back(p.X()); out.py.push_back(p.Y()); out.pz.push_back(p.Z());
      out.alpha.push_back(alpha);
      out.qt.push_back(qt);
      out.bandEll.push_back(ell);
      out.chi2.push_back(chi2);
      out.vx.push_back(x.X()); out.vy.push_back(x.Y()); out.vz.push_back(x.Z());
      out.dpv.push_back(dis);
      out.dpvSig.push_back(dsig);
      out.same_sign.push_back(ss ? 1 : 0);
      pushTrk(out.trk1, tracks, i1, q1, aux);
      pushTrk(out.trk2, tracks, i2, q2, aux);

      const double dm = std::abs(m - M_PHI);
      const bool pdau = (q1.Mag() > WP_PDAU && q2.Mag() > WP_PDAU);
      const bool prompt = (dis < WP_DPV);
      const double sd1 = out.trk1.sigd0.back(), sd2 = out.trk2.sigd0.back();
      out.wp.push_back((dm < WP_DM && pdau && prompt) ? 1 : 0);
      out.tight.push_back((dm < TIGHT_DM && pdau && prompt &&
                           sd1 < TIGHT_SIGD0 && sd2 < TIGHT_SIGD0) ? 1 : 0);
    }
  }
  return out;
}

// ---------------------------------------------------------------------------
// MC truth. The MCParticles daughter relations are EMPTY here; provenance
// comes from generatorStatus = 10000*KS + line, with `line` the 1-BASED LUND
// line of the MOTHER (MCParticles entry i is LUND line i+1). Geant
// secondaries carry generatorStatus 0 and have no LUND mother.
// ---------------------------------------------------------------------------
struct TruePhis {
  RVec<int>   idx;        // MCParticles index of the phi
  RVec<int>   dauPlus;    // MC index of the K+
  RVec<int>   dauMinus;   // MC index of the K-
  RVec<int>   mothPdg;    // PDG of the immediate LUND mother (0 = none)
  RVec<int>   origin;     // 5 = from b hadron, 4 = from c hadron, 0 = other
  RVec<float> p, pt, costheta, px, py, pz;
  RVec<float> vx, vy, vz; // production point [cm]
  RVec<float> dauPlus_p, dauMinus_p;
  RVec<int>   nmatched;   // daughters with >= 1 linked track (0/1/2)
};

inline TruePhis findTruePhis(const RVec<edm4hep::MCParticleData>& mc,
                             const RVec<RVec<int>>& mcToTracks) {
  TruePhis out;
  const int n = mc.size();
  if (n == 0) return out;
  std::vector<int> moth = lundMothers(mc);
  const std::vector<std::vector<int>> child = lundChildren(moth);
  for (int i = 0; i < n; ++i) {
    if (mc[i].PDG != PDG_PHI) continue;
    int kp = -1, km = -1;
    bool extra = false;
    for (int k : child[i]) {
      if (mc[k].PDG == PDG_K) { if (kp < 0) kp = k; else extra = true; }
      else if (mc[k].PDG == -PDG_K) { if (km < 0) km = k; else extra = true; }
      else extra = true;
    }
    if (kp < 0 || km < 0 || extra) continue;  // not phi -> K+K-

    int mp = (moth[i] >= 0) ? mc[moth[i]].PDG : 0;
    int orig = ancestorFlavour(moth, mc, moth[i]);
    TVector3 pv3(mc[i].momentum.x, mc[i].momentum.y, mc[i].momentum.z);
    int nm = 0;
    if (kp < (int)mcToTracks.size() && !mcToTracks[kp].empty()) ++nm;
    if (km < (int)mcToTracks.size() && !mcToTracks[km].empty()) ++nm;

    out.idx.push_back(i);
    out.dauPlus.push_back(kp);
    out.dauMinus.push_back(km);
    out.mothPdg.push_back(mp);
    out.origin.push_back(orig);
    out.p.push_back(pv3.Mag());
    out.pt.push_back(pv3.Perp());
    out.costheta.push_back(pv3.Mag() > 0. ? pv3.Z() / pv3.Mag() : 0.);
    out.px.push_back(pv3.X()); out.py.push_back(pv3.Y()); out.pz.push_back(pv3.Z());
    out.vx.push_back(mc[i].vertex.x);
    out.vy.push_back(mc[i].vertex.y);
    out.vz.push_back(mc[i].vertex.z);
    out.dauPlus_p.push_back(TVector3(mc[kp].momentum.x, mc[kp].momentum.y, mc[kp].momentum.z).Mag());
    out.dauMinus_p.push_back(TVector3(mc[km].momentum.x, mc[km].momentum.y, mc[km].momentum.z).Mag());
    out.nmatched.push_back(nm);
  }
  return out;
}

// Per-candidate truth. `cls` is 1 when BOTH daughter tracks link to kaons
// sharing the same true phi -> K+K- mother, else 0; the daughter MC PDG and
// its LUND mother PDG make the kaon-sample composition measurable offline.
struct PhiKKTruth {
  RVec<int> cls;          // 1 = true phi -> K+K-
  RVec<int> truephi_idx;  // index into the TruePhis lists, -1 if none
  RVec<int> trk1_mcpdg, trk2_mcpdg;      // PDG of the linked MC particle (0 = unlinked)
  RVec<int> trk1_mothpdg, trk2_mothpdg;  // PDG of its LUND mother (0 = none)
};

inline PhiKKTruth classifyPhiKK(const PhiKKCands& c,
                                const RVec<RVec<int>>& trackToMCs,
                                const RVec<edm4hep::MCParticleData>& mc,
                                const TruePhis& tp) {
  PhiKKTruth out;
  for (size_t k = 0; k < c.invM.size(); ++k) {
    int m1, m2, mo1, mo2;
    pushLegTruth(trackToMCs, mc, c.trk1.origIdx[k],
                 out.trk1_mcpdg, out.trk1_mothpdg, m1, mo1);
    pushLegTruth(trackToMCs, mc, c.trk2.origIdx[k],
                 out.trk2_mcpdg, out.trk2_mothpdg, m2, mo2);
    const int pdg1 = out.trk1_mcpdg.back(), pdg2 = out.trk2_mcpdg.back();
    int ti = -1;
    if (std::abs(pdg1) == PDG_K && std::abs(pdg2) == PDG_K && pdg1 == -pdg2 &&
        mo1 >= 0 && mo1 == mo2)
      for (size_t t = 0; t < tp.idx.size(); ++t)
        if (tp.idx[t] == mo1) { ti = (int)t; break; }
    out.truephi_idx.push_back(ti);
    out.cls.push_back(ti >= 0 ? 1 : 0);
  }
  return out;
}

// 1 if the true phi at position t was reconstructed by at least one candidate.
inline RVec<int> truePhiFound(const TruePhis& tp, const PhiKKTruth& info) {
  return foundFlags(tp.idx.size(), info.truephi_idx, info.cls, {});
}

} // namespace AlephPhiKK
} // namespace FCCAnalyses

#endif
