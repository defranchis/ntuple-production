#ifndef ALEPHV0NEW_H
#define ALEPHV0NEW_H

/*
  Standalone improved V0 reconstruction; output type
  VertexingUtils::FCCAnalysesV0. Input = flipD0_copy'ed trackstates (params AND
  covariance ALEPH->LCIO transformed); all fit-chain positions are in cm.
*/

#include <algorithm>
#include <cmath>
#include <numeric>

#include <ROOT/RVec.hxx>
#include "TVector3.h"

#include "aleph_units.h"
#include "analyzer_trkaux.h"
#include "dedx_valid.h"
#include "edm4hep/TrackState.h"
#include "FCCAnalyses/VertexingUtils.h"
#include "FCCAnalyses/VertexFitterSimple.h"

namespace FCCAnalyses {
namespace AlephV0New {

using ROOT::VecOps::RVec;
using AlephTrkAux::apVars;
using AlephTrkAux::fitTracksCm;

constexpr double m_pi_ = AlephMasses::kPiCh;
constexpr double m_p_  = AlephMasses::kProton;
constexpr double MKS   = AlephMasses::kKs;
constexpr double MLAM  = AlephMasses::kLambda;

// ---------------------------------------------------------------------------
// Cut package: single named source. TIGHT = adopted package (findV0s defaults;
// candTight re-evaluates it offline), LOOSE = ML-training superset tier. Mass
// windows, chi2 and displacement window are COMMON to both tiers.
// ---------------------------------------------------------------------------
constexpr double KS_M_LO = 0.40, KS_M_HI = 0.60;
constexpr double LAM_M_LO = 1.08, LAM_M_HI = 1.20;
constexpr double DIS_LO = 0.1, DIS_HI = 150.;
constexpr double CHI2_CUT = 10.;
constexpr double TIGHT_COS_KS_LOWP = 0.999, TIGHT_COS_KS_MIDP = 0.9995,
                 TIGHT_COS_KS_HIGHP = 0.9999;
// Lambda tight pointing, p-tiered in cos; deliberately NOT a mirror of the
// Ks ladder.
constexpr double TIGHT_COS_LAM_LOWP = 0.99995, TIGHT_COS_LAM_MIDP = 0.9999,
                 TIGHT_COS_LAM_HIGHP = 0.9999;
constexpr double TIGHT_QT_MIN_LAM = 0.04;
constexpr double AP_BAND_KS = 0.05, AP_LAM_LO = 0.10, AP_LAM_HI = 0.20;
constexpr double TIGHT_NSIG_KS_LOWP = 3., TIGHT_NSIG_KS_HIGHP = 4.;
constexpr double LOOSE_COS_POINT = 0.999;
constexpr double LOOSE_QT_MIN_LAM = 0.02;
constexpr double LOOSE_NSIG_KS = 6.;
constexpr double LOOSE_LAM_BAND_LO = 0.20, LOOSE_LAM_BAND_HI = 0.40;
// LOOSE Lambda band-distance acceptance: the stored loose tier keeps
// |ell-1| / thr(p) below this fraction of the ramp half-width (see
// lamBandThrLoose).
constexpr double LOOSE_LAM_BAND_FRAC = 0.8;
constexpr double LAM_P_LO = 8., LAM_P_HI = 20.;
// Lambda AP-band ellipse resolution sigma_ell(p), quadrature model;
// "nsig" is in units of this 68% width.
constexpr double SIG_ELL_LAM_A = 0.01622, SIG_ELL_LAM_B = 0.0033748,
                 SIG_ELL_LAM_C = 0.00015544;
// Ks AP-band ellipse resolution, linear model — single source for the
// tight/loose band cuts AND the bandSig pull.
constexpr double SIG_ELL_KS_A = 0.007, SIG_ELL_KS_B = 0.0015;
constexpr double TIGHT_LAM_NSIG = 3.;
// Mass resolution sigma_m(p) [GeV], quadrature model — used by candMassSig.
constexpr double SIG_M_KS_A = 2.658e-3, SIG_M_KS_B = 0.5214e-3,
                 SIG_M_KS_C = 0.01418e-3;
constexpr double SIG_M_LAM_A = 1.045e-3, SIG_M_LAM_B = 0.2357e-3,
                 SIG_M_LAM_C = 0.005511e-3;

// Shared per-hypothesis acceptance helpers, used by findV0s (both tiers) and
// candTight (offline re-evaluation of the booked hypothesis).
inline double ksPointThr(double pmag, double lowp, double midp, double highp) {
  return (pmag < 2.) ? lowp : (pmag < 4.) ? midp : highp;
}
// Lambda tight pointing tiers — same tier boundaries as Ks, different low-p
// value; single source for findV0s AND candTight.
inline double lamPointThr(double pmag, double lowp = TIGHT_COS_LAM_LOWP,
                          double midp = TIGHT_COS_LAM_MIDP,
                          double highp = TIGHT_COS_LAM_HIGHP) {
  return (pmag < 2.) ? lowp : (pmag < 4.) ? midp : highp;
}
inline double ksBandEll(double alpha, double qt, double pmag) {
  // frozen tuned values of the Ks AP band; do not re-derive
  const double PSTAR_K = 0.20582, ESTAR_K = 0.248806;
  double beta = pmag / std::sqrt(pmag * pmag + MKS * MKS);
  double amax = PSTAR_K / (beta * ESTAR_K);
  return std::sqrt(std::pow(alpha / amax, 2) + std::pow(qt / PSTAR_K, 2));
}
inline double sigmaEllKs(double pmag) {
  return SIG_ELL_KS_A + SIG_ELL_KS_B * pmag;
}
inline double ksBandThr(double pmag, double floor_, double nsig_lo, double nsig_hi) {
  // resolution-scaled width; floor_ acts as the low-p floor
  double nsig = (pmag < 15.) ? nsig_lo : nsig_hi;
  return std::max(floor_, nsig * sigmaEllKs(pmag));
}
inline double lamBandEll(double alpha, double qt, double pmag) {
  // frozen tuned values of the Lambda AP band; do not re-derive
  const double PSTAR_L = 0.1005, ALPHA0_L = 0.69157;
  double beta = pmag / std::sqrt(pmag * pmag + MLAM * MLAM);
  double amp = 2. * PSTAR_L / (beta * MLAM);
  return std::sqrt(std::pow((std::abs(alpha) - ALPHA0_L) / amp, 2) +
                   std::pow(qt / PSTAR_L, 2));
}
inline double lamBandThr(double pmag, double lo, double hi) {
  return (pmag < LAM_P_LO) ? lo : (pmag < LAM_P_HI) ? lo + (hi - lo) * (pmag - LAM_P_LO) / (LAM_P_HI - LAM_P_LO) : hi;
}
inline double sigmaEllLam(double pmag) {
  return std::sqrt(SIG_ELL_LAM_A * SIG_ELL_LAM_A +
                   std::pow(SIG_ELL_LAM_B * pmag, 2) +
                   std::pow(SIG_ELL_LAM_C * pmag * pmag, 2));
}
// TIGHT Lambda AP band: resolution-scaled, floored at floor_ and capped at the
// NOMINAL loose ramp edge (so the tight package is config-independent). The
// tight-inside-loose invariant is enforced on the loose side.
inline double lamBandThrTight(double pmag, double floor_ = AP_LAM_LO,
                              double nsig = TIGHT_LAM_NSIG) {
  return std::min(std::max(floor_, nsig * sigmaEllLam(pmag)),
                  lamBandThr(pmag, LOOSE_LAM_BAND_LO, LOOSE_LAM_BAND_HI));
}

// LOOSE Lambda AP band: BAND-DISTANCE convention — acceptance is a fixed
// fraction of the ramp half-width (|ell-1| < LOOSE_LAM_BAND_FRAC * thr(p)),
// floored at the tight threshold (tight_floor mirrors the caller's tight-clause
// floor) so the loose tier stays a superset of the tight one at every p.
inline double lamBandThrLoose(double pmag, double lo, double hi,
                              double tight_floor = AP_LAM_LO) {
  return std::max(LOOSE_LAM_BAND_FRAC * lamBandThr(pmag, lo, hi),
                  lamBandThrTight(pmag, tight_floor));
}

// TIGHT (adopted) package for ONE hypothesis: mass window, p-tiered pointing
// and AP band, plus the qT veto for Lambda. Single source for the finder tier
// and for the offline candTight flag.
inline bool ksTight(double m, double cp, double pmag, double alpha, double qt) {
  bool ok = (m > KS_M_LO && m < KS_M_HI) &&
            cp > ksPointThr(pmag, TIGHT_COS_KS_LOWP, TIGHT_COS_KS_MIDP,
                            TIGHT_COS_KS_HIGHP);
  if (ok)
    ok = std::abs(ksBandEll(alpha, qt, pmag) - 1.) <
         ksBandThr(pmag, AP_BAND_KS, TIGHT_NSIG_KS_LOWP, TIGHT_NSIG_KS_HIGHP);
  return ok;
}

inline bool lamTight(double m, double cp, double pmag, double alpha, double qt) {
  bool ok = (m > LAM_M_LO && m < LAM_M_HI) && cp > lamPointThr(pmag) &&
            qt > TIGHT_QT_MIN_LAM;
  if (ok)
    ok = std::abs(lamBandEll(alpha, qt, pmag) - 1.) < lamBandThrTight(pmag);
  return ok;
}

// momenta of the two tracks at the fitted vertex, already rescaled to the true
// GeV scale inside findV0s
inline void pairMomenta(const VertexingUtils::FCCAnalysesVertex& v,
                        TVector3& p1, TVector3& p2) {
  p1 = v.updated_track_momentum_at_vertex[0];
  p2 = v.updated_track_momentum_at_vertex[1];
}

inline double invMass(const TVector3& p1, double m1, const TVector3& p2, double m2) {
  double e1 = std::sqrt(p1.Mag2() + m1 * m1);
  double e2 = std::sqrt(p2.Mag2() + m2 * m2);
  TVector3 p = p1 + p2;
  double e = e1 + e2;
  return std::sqrt(std::max(0., e * e - p.Mag2()));
}

// ---------------------------------------------------------------------------
// Standard selection of ONE fitted pair, from quantities derived from the fit.
// tier: 0 rejected, 1 loose, 2 tight; pdg/m = booked hypothesis when tier > 0.
// ---------------------------------------------------------------------------
struct V0Sel { int tier; int pdg; double m; };

inline V0Sel evalV0Selection(double chi2, double dis, double cp, double pmag,
                             double alpha, double qt, double mks, double mlam) {
  V0Sel s{0, 310, mks};
  if (chi2 >= CHI2_CUT || !(chi2 == chi2)) return s;
  if (dis < DIS_LO || dis > DIS_HI) return s;
  if (pmag <= 0) return s;

  // TIGHT (adopted) package first; arbitration among the tight-passing
  // hypotheses only, so the tight subset is EXACTLY what the module would
  // output with the loose tier switched off.
  bool okKs = ksTight(mks, cp, pmag, alpha, qt);
  bool okLam = lamTight(mlam, cp, pmag, alpha, qt);
  bool tight = okKs || okLam;
  if (!tight) {
    // LOOSE training tier: flat pointing, widened AP bands, relaxed
    // Lambda qT veto; windows/chi2/displacement common.
    bool inWinKs = (mks > KS_M_LO && mks < KS_M_HI);
    bool inWinLam = (mlam > LAM_M_LO && mlam < LAM_M_HI);
    okKs = inWinKs && cp > LOOSE_COS_POINT;
    if (okKs)
      okKs = std::abs(ksBandEll(alpha, qt, pmag) - 1.) <
             ksBandThr(pmag, AP_BAND_KS, LOOSE_NSIG_KS, LOOSE_NSIG_KS);
    okLam = inWinLam && cp > LOOSE_COS_POINT && qt > LOOSE_QT_MIN_LAM;
    if (okLam)
      okLam = std::abs(lamBandEll(alpha, qt, pmag) - 1.) <
              lamBandThrLoose(pmag, LOOSE_LAM_BAND_LO, LOOSE_LAM_BAND_HI, AP_LAM_LO);
    if (!okKs && !okLam) return s;
  }
  double dks = std::abs(mks - MKS) / (0.5 * (KS_M_HI - KS_M_LO));
  double dlam = std::abs(mlam - MLAM) / (0.5 * (LAM_M_HI - LAM_M_LO));
  if (okKs && (!okLam || dks <= dlam)) { s.pdg = 310;  s.m = mks; }
  else                                 { s.pdg = 3122; s.m = mlam; }
  s.tier = tight ? 2 : 1;
  return s;
}

// ---------------------------------------------------------------------------
// The finder. np_tracks = flipD0_copy'ed non-primary trackstates, PV = fitted
// primary vertex (positions in cm). TWO-TIER: only tight-failing pairs enter the
// LOOSE tier, and tight candidates claim tracks first, so candTight==1 selects
// exactly the tight-only output. Returns candidates in claim order (tight block
// first, chi2 ascending within a tier); pdgAbs = best hypothesis (310 or 3122),
// invM its mass.
// ---------------------------------------------------------------------------
inline VertexingUtils::FCCAnalysesV0 findV0s(
    const RVec<edm4hep::TrackState>& np_tracks,
    const VertexingUtils::FCCAnalysesVertex& PV,
    double solenoidBz) {

  VertexingUtils::FCCAnalysesV0 result;
  const int nTr = np_tracks.size();
  if (nTr < 2) return result;

  TVector3 pv(PV.vertex.position[0], PV.vertex.position[1], PV.vertex.position[2]);

  struct Cand {
    VertexingUtils::FCCAnalysesVertex vtx;
    int i, j;
    int pdg;       // best hypothesis
    double m;      // mass under best hypothesis
    double chi2;
    bool tight;    // passed the tight (adopted) package
  };
  std::vector<Cand> cands;

  RVec<edm4hep::TrackState> tr_pair(2);
  for (int i = 0; i < nTr - 1; ++i) {
    tr_pair[0] = np_tracks[i];
    for (int j = i + 1; j < nTr; ++j) {
      if (np_tracks[i].omega * np_tracks[j].omega > 0) continue; // same charge
      tr_pair[1] = np_tracks[j];

      auto v = fitTracksCm(tr_pair, np_tracks, solenoidBz);
      if (v.updated_track_momentum_at_vertex.size() != 2) continue;
      double chi2 = v.vertex.chi2; // normalised, ndf=1

      // displacement (cm) + pointing
      TVector3 x(v.vertex.position[0], v.vertex.position[1], v.vertex.position[2]);
      TVector3 d = x - pv;
      double dis = d.Mag();

      TVector3 p1, p2;
      pairMomenta(v, p1, p2);
      TVector3 p = p1 + p2;
      double pmag = p.Mag();
      double cp = (dis > 0. && pmag > 0.) ? d.Dot(p) / (dis * pmag) : -2.;
      // physical charge = +sign(omega) for the flipD0_copy'ed collection (raw
      // ALEPH omega carries -charge, the flip restores +charge)
      double q1 = (np_tracks[i].omega > 0) ? 1. : -1.;
      double alpha = -99., qt = -99.;
      if (pmag > 0.) apVars(p1, p2, q1, alpha, qt);

      // hypothesis masses: Ks(pipi), Lambda(p pi) with proton = higher-|p| track
      // (in a Lambda decay the baryon carries most of the momentum)
      double mks = invMass(p1, m_pi_, p2, m_pi_);
      double mlam = (p1.Mag() > p2.Mag()) ? invMass(p1, m_p_, p2, m_pi_)
                                          : invMass(p1, m_pi_, p2, m_p_);

      V0Sel sel = evalV0Selection(chi2, dis, cp, pmag, alpha, qt, mks, mlam);
      if (sel.tier == 0) continue;

      cands.push_back({v, i, j, sel.pdg, sel.m, chi2, sel.tier == 2});
    }
  }

  // quality-ranked global claiming: tight candidates claim first (preserving
  // the tight-only output), then loose; best chi2 first within a tier
  std::vector<size_t> order(cands.size());
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
    if (cands[a].tight != cands[b].tight) return cands[a].tight;
    return cands[a].chi2 < cands[b].chi2;
  });

  std::vector<bool> used(nTr, false);
  for (size_t k : order) {
    const Cand& c = cands[k];
    if (used[c.i] || used[c.j]) continue;
    used[c.i] = true;
    used[c.j] = true;
    result.vtx.push_back(c.vtx);
    result.pdgAbs.push_back(c.pdg);
    result.invM.push_back(c.m);
  }
  return result;
}

// ---------------------------------------------------------------------------
// Truth-free per-candidate diagnostics (work on data; reco_ind is filled by
// this module, momenta are already at the true GeV scale).
// ---------------------------------------------------------------------------

inline RVec<float> candAlpha(const VertexingUtils::FCCAnalysesV0& v0s,
                             const RVec<edm4hep::TrackState>& secondaries) {
  RVec<float> out;
  for (const auto& v : v0s.vtx) {
    if (v.reco_ind.size() < 2 || v.updated_track_momentum_at_vertex.size() < 2 ||
        v.reco_ind[0] < 0 || v.reco_ind[0] >= (int)secondaries.size()) {
      out.push_back(-99.);
      continue;
    }
    double q1 = (secondaries[v.reco_ind[0]].omega > 0) ? 1. : -1.;
    double alpha, qt;
    apVars(v.updated_track_momentum_at_vertex[0],
           v.updated_track_momentum_at_vertex[1], q1, alpha, qt, -99.);
    out.push_back(alpha);
  }
  return out;
}

// Offline tight-package flag: 1 if the candidate's BOOKED hypothesis passes the
// adopted tight package, 0 if it entered via the loose training tier. Uses the
// shared helpers/constants of findV0s; assumes no variant override. On data.
inline RVec<int> candTight(const VertexingUtils::FCCAnalysesV0& v0s,
                           const VertexingUtils::FCCAnalysesVertex& PV,
                           const RVec<edm4hep::TrackState>& secondaries) {
  RVec<int> out;
  TVector3 pv(PV.vertex.position[0], PV.vertex.position[1], PV.vertex.position[2]);
  for (size_t c = 0; c < v0s.vtx.size(); ++c) {
    const auto& v = v0s.vtx[c];
    if (v.reco_ind.size() < 2 || v.updated_track_momentum_at_vertex.size() < 2 ||
        v.reco_ind[0] < 0 || v.reco_ind[0] >= (int)secondaries.size()) {
      out.push_back(0);
      continue;
    }
    TVector3 p1 = v.updated_track_momentum_at_vertex[0];
    TVector3 p2 = v.updated_track_momentum_at_vertex[1];
    TVector3 p = p1 + p2;
    double pmag = p.Mag();
    TVector3 x(v.vertex.position[0], v.vertex.position[1], v.vertex.position[2]);
    TVector3 d = x - pv;
    double dis = d.Mag();
    if (pmag <= 0 || dis <= 0) { out.push_back(0); continue; }
    double cp = d.Dot(p) / (dis * pmag);
    double q1 = (secondaries[v.reco_ind[0]].omega > 0) ? 1. : -1.;
    double alpha, qt;
    apVars(p1, p2, q1, alpha, qt);
    double m = v0s.invM[c];
    bool ok = (v0s.pdgAbs[c] == 310) ? ksTight(m, cp, pmag, alpha, qt)
                                     : lamTight(m, cp, pmag, alpha, qt);
    out.push_back(ok ? 1 : 0);
  }
  return out;
}

// ML-input pulls of the BOOKED hypothesis, in resolution units. Both SIGNED;
// -999 = undefined candidate. bandSig = (bandEll - 1) / sigma_ell(p),
// massSig = (invM - m_hyp) / sigma_m(p). Other cut variables are stored raw.
inline RVec<float> candBandSig(const VertexingUtils::FCCAnalysesV0& v0s,
                               const RVec<edm4hep::TrackState>& secondaries) {
  RVec<float> out;
  for (size_t c = 0; c < v0s.vtx.size(); ++c) {
    const auto& v = v0s.vtx[c];
    if (v.reco_ind.size() < 2 || v.updated_track_momentum_at_vertex.size() < 2 ||
        v.reco_ind[0] < 0 || v.reco_ind[0] >= (int)secondaries.size()) {
      out.push_back(-999.);
      continue;
    }
    TVector3 p1 = v.updated_track_momentum_at_vertex[0];
    TVector3 p2 = v.updated_track_momentum_at_vertex[1];
    TVector3 p = p1 + p2;
    double pmag = p.Mag();
    if (pmag <= 0) { out.push_back(-999.); continue; }
    double q1 = (secondaries[v.reco_ind[0]].omega > 0) ? 1. : -1.;
    double alpha, qt;
    apVars(p1, p2, q1, alpha, qt);
    if (v0s.pdgAbs[c] == 310)
      out.push_back((ksBandEll(alpha, qt, pmag) - 1.) / sigmaEllKs(pmag));
    else
      out.push_back((lamBandEll(alpha, qt, pmag) - 1.) / sigmaEllLam(pmag));
  }
  return out;
}

inline RVec<float> candMassSig(const VertexingUtils::FCCAnalysesV0& v0s) {
  RVec<float> out;
  for (size_t c = 0; c < v0s.vtx.size(); ++c) {
    const auto& v = v0s.vtx[c];
    if (v.updated_track_momentum_at_vertex.size() < 2) {
      out.push_back(-999.);
      continue;
    }
    TVector3 p = v.updated_track_momentum_at_vertex[0] +
                 v.updated_track_momentum_at_vertex[1];
    double pmag = p.Mag(), p2 = pmag * pmag;
    bool isKs = (v0s.pdgAbs[c] == 310);
    double a = isKs ? SIG_M_KS_A : SIG_M_LAM_A;
    double b = isKs ? SIG_M_KS_B : SIG_M_LAM_B;
    double cc = isKs ? SIG_M_KS_C : SIG_M_LAM_C;
    double sig = std::sqrt(a * a + b * b * p2 + cc * cc * p2 * p2);
    out.push_back((v0s.invM[c] - (isKs ? MKS : MLAM)) / sig);
  }
  return out;
}

// Pointing significance: chi2-like significance of the displacement component
// PERPENDICULAR to the candidate momentum (all in cm). d = candidate vertex -
// reference vertex, p = candidate momentum, cV/cR = packed lower-triangular
// position covariances (xx,yx,yy,zx,zy,zz); reference = PV for candPointSig, an
// SV for candSVPointing. Returns -1 for degenerate or singular geometry.
template <typename CovV, typename CovR>
inline float pointSigTransverse(const TVector3& d, const TVector3& p,
                                const CovV& cV, const CovR& cR) {
  if (p.Mag() <= 0 || d.Mag() <= 0) return -1.;
  TVector3 ph = p.Unit();
  TVector3 u1 = ph.Orthogonal().Unit();
  TVector3 u2 = ph.Cross(u1);
  double C[3][3];
  AlephTrkAux::sumCovPacked(cV, cR, C);
  auto quad = [&](const TVector3& a, const TVector3& b) {
    double s = 0.;
    double av[3] = {a.X(), a.Y(), a.Z()}, bv[3] = {b.X(), b.Y(), b.Z()};
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) s += av[i] * C[i][j] * bv[j];
    return s;
  };
  double c11 = quad(u1, u1), c22 = quad(u2, u2), c12 = quad(u1, u2);
  double det = c11 * c22 - c12 * c12;
  if (det <= 0. || c11 <= 0. || c22 <= 0.) return -1.;
  double d1 = d.Dot(u1), d2 = d.Dot(u2);
  double sig2 = (d1 * (c22 * d1 - c12 * d2) + d2 * (c11 * d2 - c12 * d1)) / det;
  return sig2 > 0. ? std::sqrt(sig2) : 0.;
}

inline RVec<float> candPointSig(const VertexingUtils::FCCAnalysesV0& v0s,
                                const VertexingUtils::FCCAnalysesVertex& PV) {
  RVec<float> out;
  TVector3 pv(PV.vertex.position[0], PV.vertex.position[1], PV.vertex.position[2]);
  for (const auto& v : v0s.vtx) {
    TVector3 x(v.vertex.position[0], v.vertex.position[1], v.vertex.position[2]);
    TVector3 p(0., 0., 0.);
    for (const auto& tp : v.updated_track_momentum_at_vertex) p += tp;
    out.push_back(pointSigTransverse(x - pv, p, v.vertex.covMatrix,
                                     PV.vertex.covMatrix));
  }
  return out;
}

inline RVec<float> candQt(const VertexingUtils::FCCAnalysesV0& v0s) {
  RVec<float> out;
  for (const auto& v : v0s.vtx) {
    TVector3 pa = v.updated_track_momentum_at_vertex[0];
    TVector3 p = pa + v.updated_track_momentum_at_vertex[1];
    out.push_back(pa.Cross(p.Unit()).Mag());
  }
  return out;
}

// ---------------------------------------------------------------------------
// vertex-fit covariance exposure + per-daughter joins
// ---------------------------------------------------------------------------

// Vertex-fit covariance component ic of every candidate (packed lower triangle:
// 0=xx 1=yx 2=yy 3=zx 4=zy 5=zz, cm^2 — same packing as Vertex_refit_cov_*).
// Works for the V0 and the new-SV module output (both are FCCAnalysesV0).
inline RVec<float> candCovComp(const VertexingUtils::FCCAnalysesV0& v0s,
                               int ic) {
  RVec<float> out;
  for (const auto& v : v0s.vtx) out.push_back(v.vertex.covMatrix[ic]);
  return out;
}

// Daughter k (0/1) of every candidate as an ORIGINAL Tracks index: reco_ind
// (secondary space) walked through sec2orig. -1 when unavailable; truth-free.
inline RVec<int> candDaughterOrigIdx(const VertexingUtils::FCCAnalysesV0& v0s,
                                     const RVec<int>& sec2orig, int k) {
  RVec<int> out;
  for (const auto& v : v0s.vtx) {
    int idx = -1;
    if (k < (int)v.reco_ind.size()) {
      int s = v.reco_ind[k];
      if (s >= 0 && s < (int)sec2orig.size()) idx = sec2orig[s];
    }
    out.push_back(idx);
  }
  return out;
}

// ---------------------------------------------------------------------------
// Pointing of every V0 candidate at the nearest SV (largest cos between the
// candidate momentum and the SV->candidate line). FEATURE ONLY: no feedback into
// selection. SVs sharing a daughter track are excluded, both legs walked through
// sec2orig into the ORIGINAL Tracks space (unmapped -1 never matches, veto fails
// open). Sentinels cos=-2, sig=-1, idx=-1 when no usable SV remains; pointSig is
// ALSO -1 on a singular covariance, so test "no SV" on idx, never on pointSig.
// ---------------------------------------------------------------------------
struct V0SVPointing {
  RVec<float> cosPoint;  // cos(candidate momentum, SV->candidate vector)
  RVec<float> pointSig;  // transverse pointing significance wrt that SV
  RVec<int>   svIdx;     // index of that SV in the SV collection
};

inline V0SVPointing candSVPointing(const VertexingUtils::FCCAnalysesV0& v0s,
                                   const VertexingUtils::FCCAnalysesV0& svs,
                                   const RVec<int>& sec2orig) {
  V0SVPointing out;
  const int nSec = (int)sec2orig.size();
  auto toOrig = [&](int s) { return (s >= 0 && s < nSec) ? sec2orig[s] : -1; };

  std::vector<std::vector<int>> sv_orig(svs.vtx.size());
  for (size_t s = 0; s < svs.vtx.size(); ++s)
    for (int t : svs.vtx[s].reco_ind) {
      int o = toOrig(t);
      if (o >= 0) sv_orig[s].push_back(o);
    }

  for (const auto& v : v0s.vtx) {
    int o1 = (v.reco_ind.size() > 0) ? toOrig(v.reco_ind[0]) : -1;
    int o2 = (v.reco_ind.size() > 1) ? toOrig(v.reco_ind[1]) : -1;
    TVector3 x(v.vertex.position[0], v.vertex.position[1], v.vertex.position[2]);
    TVector3 p(0., 0., 0.);
    for (const auto& tp : v.updated_track_momentum_at_vertex) p += tp;

    int best = -1;
    double best_cos = -2.;
    if (p.Mag() > 0.) {
      for (size_t s = 0; s < svs.vtx.size(); ++s) {
        bool shared = false;
        for (int o : sv_orig[s])
          if (o == o1 || o == o2) { shared = true; break; }
        if (shared) continue;
        const auto& sv = svs.vtx[s].vertex;
        TVector3 d = x - TVector3(sv.position[0], sv.position[1], sv.position[2]);
        double dm = d.Mag();
        if (dm <= 0.) continue;
        double cp = d.Dot(p) / (dm * p.Mag());
        if (best < 0 || cp > best_cos) { best_cos = cp; best = (int)s; }
      }
    }
    if (best < 0) {
      out.cosPoint.push_back(-2.);
      out.pointSig.push_back(-1.);
      out.svIdx.push_back(-1);
      continue;
    }
    const auto& sv = svs.vtx[best].vertex;
    TVector3 d = x - TVector3(sv.position[0], sv.position[1], sv.position[2]);
    out.cosPoint.push_back(best_cos);
    out.pointSig.push_back(pointSigTransverse(d, p, v.vertex.covMatrix,
                                              sv.covMatrix));
    out.svIdx.push_back(best);
  }
  return out;
}

// Measurement index of every original track index, or -1: the dE/dx join built
// ONCE per collection per event, in place of a scan per requested track. A
// track measured twice keeps the first measurement, as the scan did.
// An entry is -1 when the measurement fails the shared validity gate (value ==
// omega of the track = the failed-leg sentinel, or non-finite/non-positive
// value or error), so value and error branches share one lookup.
inline RVec<int> dedxIndexByTrack(const RVec<float>& value,
                                  const RVec<float>& error,
                                  const RVec<int>& meas_track_idx,
                                  const RVec<edm4hep::TrackState>& trackStates) {
  RVec<int> out(trackStates.size(), -1);
  std::vector<char> seen(out.size(), 0);
  const size_t nm = std::min({value.size(), error.size(),
                              meas_track_idx.size()});
  for (size_t j = 0; j < nm; ++j) {
    const int t = meas_track_idx[j];
    if (t < 0 || t >= (int)out.size() || seen[t]) continue;
    seen[t] = 1;
    if (AlephDedx::dEdxValid(value[j], error[j], trackStates[t].omega))
      out[t] = (int)j;
  }
  return out;
}

// Per-track quantity by ORIGINAL track index through that join; -1 when the
// track has no valid measurement.
inline RVec<float> trackQuantityByIndex(const RVec<int>& want,
                                        const RVec<float>& values,
                                        const RVec<int>& meas_of_track) {
  RVec<float> out;
  for (int w : want) {
    float val = -1.f;
    if (w >= 0 && w < (int)meas_of_track.size()) {
      const int j = meas_of_track[w];
      if (j >= 0 && j < (int)values.size()) val = values[j];
    }
    out.push_back(val);
  }
  return out;
}

} // namespace AlephV0New
} // namespace FCCAnalyses

#endif
