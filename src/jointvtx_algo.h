#ifndef JOINTVTX_ALGO_H
#define JOINTVTX_ALGO_H

/*
  Decision logic of the joint PV + SV vertexer: runEvent(Event, Settings) -> Result.

  Order of one event (that of the reference implementation):
    momentum scale, tracks, lifetime-sign axes
    stage 0   track selection: finite states, TPC gate, removal of the existing finder's tight V0
              daughters, production-radius bound from the hits
    stage 1   reference PV (claims no tracks), lifetime-signed IP significances s3 at it
    stage 2-3 candidate supply (vertex-density maxima -> candidates -> scores -> merging)
    V0 class  V0 candidates, built once at the reference PV
    fit list  S0 tracks in the PV window + candidate tracks + V0-candidate tracks; local indices
    stage 4   annealed joint fit, drop loop, final selection
    stage 5   seeded second pass (pairs of left-over tracks -> growth -> stage 4 again)
    split     B -> D cascade split of the final SVs
    outputs   hard assignment, SV masses, s3 at the final PV, V0 list
  Every per-event quantity is a local of runEvent (EventState); nothing is static.
  runEventT takes a tracer whose hooks see the state after each step (NoTrace: none).
  V0 class (v0Class): the V0 candidates are built once and reused by every joint fit; each weight-matrix call
  recomputes their marginal weights against the PV-form denominators and scales the vertex weights of their tracks,
  so the V0 state is that of the last call; the joint fit ends with one more call at the last temperature.
  A fit or solve that fails where the reference raises (singular or non-finite linear algebra; its run ends) marks
  the event failed: runEvent stops after that step and returns a Result without vertices whose status names the
  step. A failure the reference catches keeps its branch: the jet-projected PV falls back to the teardown fit.
*/

#include <algorithm>
#include <cmath>
#include <vector>

#include <Eigen/Dense>

#include "jointvtx_common.h"
#include "jointvtx_event.h"
#include "jointvtx_settings.h"
#include "jointvtx_track.h"
#include "jointvtx_lookup.h"
#include "jointvtx_seedless.h"
#include "jointvtx_fit.h"
#include "jointvtx_v0.h"

namespace FCCAnalyses {
namespace AlephJointVtx {

constexpr double JV_SCORE_OFF = -1e8;      // a score threshold below this disables the score cut

namespace detail {
constexpr bool sameText(const char* a, const char* b) {
  while (*a != '\0' && *a == *b) {
    ++a;
    ++b;
  }
  return *a == *b;
}
}  // namespace detail

// configuration switches of the algorithm: only the nominal value of each is implemented
static_assert(JV_JOINT, "JV_JOINT: only true is implemented");
static_assert(detail::sameText(JV_SV_MODE, "finder"), "JV_SV_MODE: only \"finder\" is implemented");
static_assert(detail::sameText(JV_PV_REF, "jetproj_hard"), "JV_PV_REF: only \"jetproj_hard\" is implemented");
static_assert(detail::sameText(JV_SUPPLY, "vd"), "JV_SUPPLY: only \"vd\" is implemented");
static_assert(detail::sameText(JV_SEED_MODE, "signed"), "JV_SEED_MODE: only \"signed\" is implemented");
static_assert(!JV_SEED_KIN, "JV_SEED_KIN: only false is implemented");
static_assert(!JV_VETO_STAGE1_V0, "JV_VETO_STAGE1_V0: only false is implemented");
static_assert(detail::sameText(JV_GROW_RULE, "symmetric"), "JV_GROW_RULE: only \"symmetric\" is implemented");
static_assert(!JV_GROW_SYM_REFIT, "JV_GROW_SYM_REFIT: only false is implemented");
static_assert(JV_GROW_MIN_VDET == 0, "JV_GROW_MIN_VDET: only 0 is implemented");
static_assert(detail::sameText(JV_CHI2_MODE_PV, "point"), "JV_CHI2_MODE_PV: only \"point\" is implemented");
static_assert(detail::sameText(JV_CHI2_MODE_SV, "vcov"), "JV_CHI2_MODE_SV: only \"vcov\" is implemented");
static_assert(detail::sameText(JV_WEIGHT_FORM, "student"), "JV_WEIGHT_FORM: only \"student\" is implemented");
static_assert(detail::sameText(JV_PRIOR_FORM, "lorentz"), "JV_PRIOR_FORM: only \"lorentz\" is implemented");
static_assert(!JV_PRIOR_UPDATE, "JV_PRIOR_UPDATE: only false is implemented");
static_assert(detail::sameText(JV_ASSOC_RULE, "any"), "JV_ASSOC_RULE: only \"any\" is implemented");
static_assert(!JV_ANNEAL_CAND, "JV_ANNEAL_CAND: only false is implemented");
static_assert(detail::sameText(JV_MERGE_WHEN, "each"), "JV_MERGE_WHEN: only \"each\" is implemented");
static_assert(detail::sameText(JV_PV_FINAL, "soft"), "JV_PV_FINAL: only \"soft\" is implemented");
static_assert(detail::sameText(JV_SV_FINAL, "frozen"), "JV_SV_FINAL: only \"frozen\" is implemented");
static_assert(detail::sameText(JV_S5_MODE, "seeded"), "JV_S5_MODE: only \"seeded\" is implemented");
static_assert(!JV_ITER_UPDATE_S3, "JV_ITER_UPDATE_S3: only false is implemented");
static_assert(detail::sameText(JV_SPLIT_PICK, "valid"), "JV_SPLIT_PICK: only \"valid\" is implemented");
static_assert(JV_V0_POINT, "JV_V0_POINT: only true is implemented");
static_assert(!JV_V0_PT_POINT, "JV_V0_PT_POINT: only false is implemented");

// candidate or secondary vertex: position, covariance, tracks (event indices before the fit list is
// built, local fit indices after), chi2 of its fit; cF: C is column-major in the reference (a fit result),
// which fixes the summation order of u C u
struct Vtx {
  Vec3 x = Vec3::Zero();
  Mat3 C = Mat3::Zero();
  std::vector<int> trk;
  double chi2 = 0.;
  bool cF = true;
};

// call sites of the weight matrix inside the joint fit
enum WmSite : int { kWmAnneal = 0, kWmPvTmin = 1, kWmDrop = 2, kWmFinal = 3, kWmOther = 4 };

// working state of one event
struct EventState {
  int n = 0;
  std::vector<Trk> t;            // fitter-convention states, covariances and hit counts
  std::vector<double> p;         // momenta at the momentum constant of the settings
  double pscale = 1.;            // momentum constant of the settings over that of the event's momenta
  double ptHalfEv = JV_PT_HALF;  // half of the event's constant: momentumAt(.., ptHalfEv) * pscale
  double ptHalfCfg = JV_PT_HALF; // half of the settings' constant (V0 and conversion kernels)
  std::vector<Vec3> axes;        // lifetime-sign axes
  std::vector<int> jetOf;
  // stage 0
  std::vector<char> ok;          // track enters the vertexer
  std::vector<char> tightv0;     // daughter of a tight V0 of the existing finder
  std::vector<double> rmax;      // bound on the transverse production radius [cm]
  // reference PV and significances
  Vec3 pvx = Vec3::Constant(JV_NAN);
  Mat3 pvC = Mat3::Constant(JV_NAN);
  std::vector<int> pvtrk;
  std::vector<double> s3;
  // candidates and fit list
  std::vector<Vtx> cands;
  std::vector<int> fitIdx;       // event index of each local fit index
  std::vector<int> loc;          // local fit index of each event track, -1 outside
  std::vector<Trk> tf;           // tracks of the fit list
  std::vector<double> s3f;       // s3 of the fit tracks (reference PV)
  // joint fit
  Vec3 px = Vec3::Constant(JV_NAN);
  Mat3 pC = Mat3::Constant(JV_NAN);
  std::vector<Vtx> V;
  MatX W;                        // (1 + number of SVs) x fit tracks
  std::vector<double> pvw;       // weights of the last PV refit (empty: row 0 of W)
  std::vector<int> hard;         // per fit track: argmax of W, -1 at or below hardMin
  V0Set v0;                      // V0 class: candidates and the state of the last weight-matrix call
  bool failed = false;           // a fit or solve failed where the reference raises
};

namespace detail {

// hooks of a traced run; the arguments are those of the trace records
struct NoTrace {
  void stage0(const EventState&) {}
  void refPv(const EventState&) {}
  void s3(const EventState&) {}
  void vdMaxima(const VdMaxima&) {}
  void candidate(const std::vector<int>&, const Vtx*) {}
  void score(int, const std::vector<int>&, double, double) {}
  void merge(int, const std::vector<Vtx>&) {}
  void supply(const EventState&) {}
  void v0cand(const EventState&) {}
  void fitIdx(const EventState&) {}
  void wm(int, double, const std::vector<Vtx>&, const MatX&, const std::vector<double>*) {}
  void finalOk(const std::vector<int>&, bool) {}
  void stage4(const Vec3&, const Mat3&, const std::vector<Vtx>&, const MatX&) {}
  void iter(int, const std::vector<int>&, int, const std::vector<Vtx>&) {}
  void split(const std::vector<Vtx>&, const std::vector<int>&) {}
  void out(const Result&) {}
  void ties(int, int) {}   // site, number of exactly tied sort keys
};

// np.max over a column: NaN propagates
inline double colMaxNaN(const MatX& W, int j) {
  double m = W(0, j);
  for (int r = 1; r < W.rows(); ++r) {
    if (std::isnan(m)) break;
    const double v = W(r, j);
    if (std::isnan(v) || v > m) m = v;
  }
  return m;
}

// momenta of the settings' constant; tracks of the kernels; lifetime-sign axes
inline void initState(const Event& ev, const Settings& s, EventState& st) {
  const int n = ev.size();
  st.n = n;
  st.pscale = s.ptConst / ev.ptConst;
  st.ptHalfEv = 0.5 * ev.ptConst;
  st.ptHalfCfg = 0.5 * s.ptConst;
  st.p = ev.p;
  if (st.pscale != 1.0)
    for (int i = 0; i < n; ++i) st.p[i] = ev.p[i] * st.pscale;
  st.t.resize(n);
  for (int i = 0; i < n; ++i) {
    st.t[i].p = ev.fpar[i];
    st.t[i].C = ev.fcov[i];
    st.t[i].nvdet = ev.nvdet[i];
    st.t[i].nitc = ev.nitc[i];
  }
  st.jetOf = ev.jet;
  const std::vector<double> th(ev.jetTheta.begin(), ev.jetTheta.end()), ph(ev.jetPhi.begin(), ev.jetPhi.end());
  st.axes = jetAxes(th, ph, ev.thrust);
}

// np.min of v over idx (NaN propagates)
inline double minAt(const std::vector<double>& v, const std::vector<int>& idx) {
  double m = v[idx[0]];
  for (size_t k = 1; k < idx.size(); ++k) m = npMinimum(m, v[idx[k]]);
  return m;
}

// momentum of track i at the arc of its xy point nearest to x, with the momentum constant of the settings
inline Vec3 trackMomentum(const EventState& st, int i, const Vec3& x) {
  Vec3 p = momentumAt(st.t[i].p, arcNear(st.t[i].p, x), st.ptHalfEv);
  if (st.pscale != 1.0)
    for (int a = 0; a < 3; ++a) p[a] = p[a] * st.pscale;
  return p;
}

inline InvMass massAt(const EventState& st, const std::vector<int>& trk, const Vec3& x,
                      const double* masses = nullptr) {
  return invMass(st.t, trk, x, masses, st.ptHalfEv, st.pscale);
}

// acceptance window of a vertex (tracks trk at x, covariance C) around the reference PV
inline bool accept(const EventState& st, const Settings& s, const std::vector<int>& trk, const Vec3& x,
                   const Mat3& C, bool cF, const Vec3& pvx, const Mat3& pvC) {
  if (!s.accept) return true;
  AcceptCfg c;
  c.disLo = s.accDisLo;
  c.disHi = s.accDisHi;
  c.cosMin = s.accCos;
  c.siglMax = s.accSigl;
  c.fsig2 = s.accFsig2;
  c.maxTrk = s.accMaxTrk;
  const Vec3 P = massAt(st, trk, x).P;
  return acceptWindow(static_cast<int>(trk.size()), x, C, pvx, pvC, P, c, cF);
}

// features of a vertex for the logistic scores; s3: lifetime-signed significances at the reference PV
inline CandFeats vertexFeats(const EventState& st, const std::vector<int>& trk, const Vec3& x, const Mat3& C,
                             bool cF, double chi2, const Vec3& pvx, const Mat3& pvC, const std::vector<double>& s3) {
  CandFeatsIn in;
  in.x = x;
  in.pvx = pvx;
  in.C = C;
  in.pvC = pvC;
  in.covFortran = cF;
  in.chi2 = chi2;
  const InvMass mp = massAt(st, trk, x);
  in.mass = mp.m;
  in.P = mp.P;
  std::vector<Vec5> states;
  for (int k : trk) {
    in.s3.push_back(s3[k]);
    in.p.push_back(st.p[k]);
    in.nvdet.push_back(st.t[k].nvdet);
    in.jet.push_back(st.jetOf[k]);
    states.push_back(st.t[k].p);
  }
  in.drmax = maxEtaPhiDr(states);
  in.jetAxes = st.axes;
  return candFeats(in);
}

// logistic score: which 0 = candidate coefficients, 1 = final-vertex coefficients
inline double vertexScore(const EventState& st, const std::vector<int>& trk, const Vec3& x, const Mat3& C, bool cF,
                          double chi2, const Vec3& pvx, const Mat3& pvC, const std::vector<double>& s3, int which) {
  const ScoreInputs in = scoreInputs(vertexFeats(st, trk, x, C, cF, chi2, pvx, pvC, s3));
  return logisticScore(which == 0 ? JV_SCORE_CAND : JV_SCORE_FINAL, in);
}

// final selection of a secondary vertex v (tracks trk, event indices) against the PV pv: two tracks and the final
// score of the tracks' chi2 sum at v.x reaching slFinalMin
template <class Tr>
inline bool finalOkVertex(EventState& st, const Settings& s, const std::vector<int>& trk, const Vtx& v,
                          const Vtx& pv, const std::vector<double>& s3, Tr& tr) {
  if (s.slFinalMin < JV_SCORE_OFF) return true;
  std::vector<double> ch, L;
  if (chi2ToPoint(st.t, trk, v.x, ch, L) != kLinOk || JOINTVTX_FAIL_AT(kFsFinalChi2)) {
    st.failed = true;
    return false;
  }
  if (trk.size() < 2) return false;
  const double sc = vertexScore(st, trk, v.x, v.C, v.cF, npSum(ch), pv.x, pv.C, s3, 1);
  tr.score(1, trk, sc, s.slFinalMin);
  return sc >= s.slFinalMin;
}

// lifetime prior of the tracks of significances s3f at temperature T (Lorentz form)
inline std::vector<double> lifetimePriors(const Settings& s, const std::vector<double>& s3f, double T) {
  std::vector<double> pr(s3f.size(), 1.);
  if (!s.lifetimePrior) return pr;
  for (size_t k = 0; k < s3f.size(); ++k) {
    const double sp = npMaximum(s3f[k], 0.) / s.priorS;
    pr[k] = 1. / (1. + sp * sp);
    if (s.priorT) pr[k] = std::pow(pr[k], 1. / T);
  }
  return pr;
}

// stage 0: usable tracks (finite state and covariance, TPC gate), tight-V0 daughters, hit-radius bound
inline void stage0(const Event& ev, const Settings& s, EventState& st) {
  const int n = st.n;
  st.ok.assign(n, 1);
  for (int i = 0; i < n; ++i) st.ok[i] = allFinite(ev.fpar[i]) && allFinite(ev.fcov[i]);
  if (s.gateMinTpc > 0)
    for (int i = 0; i < n; ++i) st.ok[i] = st.ok[i] && ev.ntpc[i] >= s.gateMinTpc;
  st.tightv0.assign(n, 0);
  const size_t nv = std::min(ev.v0T1.size(), std::min(ev.v0T2.size(), ev.v0Tight.size()));
  for (size_t k = 0; k < nv; ++k) {
    const int a = ev.v0T1[k], b = ev.v0T2[k];
    if (ev.v0Tight[k] && a >= 0 && b >= 0 && a < n && b < n) st.tightv0[a] = st.tightv0[b] = 1;
  }
  if (s.removeTightV0)
    for (int i = 0; i < n; ++i) st.ok[i] = st.ok[i] && !st.tightv0[i];
  st.rmax.assign(n, JV_INF);
  if (s.hitRadius)
    for (int i = 0; i < n; ++i) st.rmax[i] = maxRadiusFromHits(ev.nvdet[i], ev.nitc[i]);
}

// stage 1: reference PV of the stage-0 tracks inside the pre-selection window: jet-projected PV (hard form),
// with the beam-spot-constrained teardown fit as the fallback when it fails or has fewer than two tracks; the PV
// tracks are those within pvChi2max of it. The reference fits the teardown in every event before the jet-projected
// PV, whose hard form does not use it; here it runs only as the fallback, which changes no value
inline void referencePv(const Event& ev, const Settings& s, EventState& st) {
  std::vector<int> pre;
  for (int i = 0; i < st.n; ++i)
    if (st.ok[i] && preWindow(ev.par[i], ev.bs, s.pvWindowD0, s.pvWindowZ0)) pre.push_back(i);
  BeamSpot bsp;
  bsp.c = ev.bs;
  bsp.sig = s.bsSig;
  if (pre.size() >= 2) {
    std::vector<Vec5> par;
    std::vector<Vec15> cov;
    std::vector<int> jet;
    for (int i : pre) {
      par.push_back(ev.par[i]);
      Vec15 c;
      for (int k = 0; k < 15; ++k) c[k] = ev.cov[i][k];
      cov.push_back(c);
      jet.push_back(ev.jet[i]);
    }
    const std::vector<double> th(ev.jetTheta.begin(), ev.jetTheta.end());
    const std::vector<double> ph(ev.jetPhi.begin(), ev.jetPhi.end());
    const std::vector<double> jp(ev.jetP.begin(), ev.jetP.end());
    const JpPV r = jetProjPV(par, cov, jet, th, ph, jp, bsp);
    if (r.ok()) {
      std::vector<double> ch, L;
      if (chi2ToPoint(st.t, pre, r.x, ch, L) != kLinOk || JOINTVTX_FAIL_AT(kFsRefPvChi2)) st.failed = true;
      st.pvx = r.x;
      st.pvC = r.C;
      st.pvtrk.clear();
      for (size_t k = 0; k < pre.size(); ++k)
        if (ch[k] < s.pvChi2max) st.pvtrk.push_back(pre[k]);
      return;
    }
  }
  const TeardownFit f = teardown(st.t, pre, ev.bs, &bsp, s.pvChi2max);
  if (!f.ok() || JOINTVTX_FAIL_AT(kFsTeardown)) st.failed = true;
  st.pvx = f.x;
  st.pvC = f.C;
  st.pvtrk = f.kept;
}

// lifetime-signed 3D IP significances of all tracks at a PV (all NaN and the failure status if a track cannot be
// linearised)
inline LinStatus significances(const EventState& st, const Vec3& pv, const Mat3& pvC, std::vector<double>& s3) {
  IpInfo r;
  const LinStatus ls = ipInfo(st.t, pv, pvC, st.axes, st.jetOf, r);
  if (ls != kLinOk) {
    s3.assign(st.n, JV_NAN);
    return ls;
  }
  s3 = r.s3;
  return kLinOk;
}

// V0-like pair (opposite charges, pion-pair mass near the K0S, p-pi near the Lambda, or e-e below the
// conversion window); fixed windows
inline bool v0LikePair(const Event& ev, const EventState& st, int i, int j, const Vec3& x) {
  return v0Like(ev.q[i], ev.q[j], trackMomentum(st, i, x), trackMomentum(st, j, x));
}

// candidate from a vertex-density maximum: fit the claimed tracks from x0, drop the worst track while it is above
// slTrimChi2, the vertex lies beyond a track's production-radius bound or there are more than accMaxTrk tracks;
// then the fit quality, mass, flight-significance, direction, V0-like, acceptance and score cuts
template <class Tr>
inline bool candidate(const Event& ev, const Settings& s, EventState& st, std::vector<int> trk, const Vec3& x0,
                      const Vec3& pvx, const Mat3& pvC, Vtx& out, Tr& tr) {
  Vec3 x = x0;
  VtxFit f;
  while (trk.size() >= 2) {
    f = fitVertex(st.t, trk, x, nullptr, nullptr, JV_CAND_FIT_NITER, JV_CAND_FIT_TOL);
    if (!f.ok() || JOINTVTX_FAIL_AT(kFsCandFit)) {
      st.failed = true;
      return false;
    }
    if (!allFinite(f.x)) return false;
    x = f.x;
    const int j = argmaxNaN(f.chi2);
    const double rad = std::hypot(x[0], x[1]);
    if (trk.size() > 2 && (f.chi2[j] > s.slTrimChi2 || rad >= minAt(st.rmax, trk) ||
                           static_cast<int>(trk.size()) > s.accMaxTrk)) {
      trk.erase(trk.begin() + j);
      continue;
    }
    break;
  }
  if (trk.size() < 2) return false;
  if (npSum(f.chi2) / static_cast<double>(2 * static_cast<int>(trk.size()) - 3) >= s.vtxChi2ndf ||
      std::hypot(x[0], x[1]) >= minAt(st.rmax, trk))
    return false;
  while (trk.size() > 2 && massAt(st, trk, x).m >= s.massMax) {
    trk.erase(trk.begin() + argmaxNaN(f.chi2));
    f = fitVertex(st.t, trk, x, nullptr, nullptr, JV_CAND_FIT_NITER, JV_CAND_FIT_TOL);
    if (!f.ok() || JOINTVTX_FAIL_AT(kFsCandMassFit)) {
      st.failed = true;
      return false;
    }
    x = f.x;
  }
  if (dsig2(x, f.C, pvx, pvC) < detail::powSquare(s.slFsig)) return false;
  const InvMass mp = massAt(st, trk, x);
  const Vec3 d = x - pvx;
  if (dotSeq(d, mp.P) <= 0 || mp.m >= s.massMax) return false;
  if (trk.size() == 2 && s.v0Veto && v0LikePair(ev, st, trk[0], trk[1], x)) return false;
  if (!accept(st, s, trk, x, f.C, true, pvx, pvC)) return false;
  const double chi2 = npSum(f.chi2);
  if (s.slScoreMin > JV_SCORE_OFF) {
    const double sc = vertexScore(st, trk, x, f.C, true, chi2, pvx, pvC, st.s3, 0);
    tr.score(0, trk, sc, s.slScoreMin);
    if (sc < s.slScoreMin) return false;
  }
  out.x = x;
  out.C = f.C;
  out.trk = trk;
  out.chi2 = chi2;
  out.cF = true;
  return true;
}

// merging of candidates, restarted after every change: a pair sharing at least mergeShare of the smaller track
// list or closer than mergeSig is refitted with the union of the tracks; the union replaces both if its chi2 / ndf
// is below vtxChi2ndf, otherwise the shared tracks stay with the larger candidate and the smaller one is refitted
// (it keeps its chi2) or dropped below two tracks
inline void mergeCandidates(EventState& st, const Settings& s, std::vector<Vtx>& cands) {
  bool changed = true;
  while (changed && cands.size() > 1) {
    changed = false;
    for (size_t a = 0; a < cands.size() && !changed; ++a) {
      for (size_t b = a + 1; b < cands.size(); ++b) {
        const Vtx& A = cands[a];
        const Vtx& B = cands[b];
        int sh = 0;
        for (int k : A.trk)
          if (std::find(B.trk.begin(), B.trk.end(), k) != B.trk.end()) ++sh;
        const bool close = dsig2(A.x, A.C, B.x, B.C) < detail::powSquare(s.mergeSig);
        const int nA = static_cast<int>(A.trk.size()), nB = static_cast<int>(B.trk.size());
        const bool sharedBig = sh >= s.mergeShare * std::min(nA, nB);
        if (!(s.mergeRequireShare ? (sharedBig && close) : (sharedBig || close))) continue;
        std::vector<int> un = A.trk;
        un.insert(un.end(), B.trk.begin(), B.trk.end());
        std::sort(un.begin(), un.end());
        un.erase(std::unique(un.begin(), un.end()), un.end());
        Vec3 x0;
        for (int i = 0; i < 3; ++i) x0[i] = 0.5 * (A.x[i] + B.x[i]);
        const VtxFit f = fitVertex(st.t, un, x0, nullptr, nullptr, JV_CAND_FIT_NITER, JV_CAND_FIT_TOL);
        if (!f.ok() || JOINTVTX_FAIL_AT(kFsMergeFit)) {
          st.failed = true;
          return;
        }
        const double tot = npSum(f.chi2);
        const double ndf = static_cast<double>(std::max(2 * static_cast<int>(un.size()) - 3, 1));
        if (tot / ndf < s.vtxChi2ndf && allFinite(f.x)) {
          Vtx m;
          m.x = f.x;
          m.C = f.C;
          m.trk = un;
          m.chi2 = tot;
          cands[a] = m;
          cands.erase(cands.begin() + static_cast<std::ptrdiff_t>(b));
          changed = true;
          break;
        } else if (sh > 0) {
          const size_t is = nA < nB ? a : b, ib = nA < nB ? b : a;
          Vtx& sm = cands[is];
          const std::vector<int>& big = cands[ib].trk;
          std::vector<int> keep;
          for (int k : sm.trk)
            if (std::find(big.begin(), big.end(), k) == big.end()) keep.push_back(k);
          sm.trk = keep;
          if (sm.trk.size() < 2) {
            cands.erase(cands.begin() + static_cast<std::ptrdiff_t>(is));
          } else {
            const VtxFit g = fitVertex(st.t, sm.trk, sm.x, nullptr, nullptr, JV_CAND_FIT_NITER, JV_CAND_FIT_TOL);
            if (!g.ok() || JOINTVTX_FAIL_AT(kFsMergeRefit)) {
              st.failed = true;
              return;
            }
            sm.x = g.x;
            sm.C = g.C;
            sm.cF = true;
          }
          changed = true;
          break;
        }
      }
    }
  }
}

// candidates from the vertex-density maxima of the stage-0 tracks, in order of decreasing density: a maximum
// resolved from the PV with at least two claiming tracks gives a candidate
template <class Tr>
inline void supplyVd(const Event& ev, const Settings& s, EventState& st, const Vec3& pvx, const Mat3& pvC,
                     std::vector<Vtx>& out, Tr& tr) {
  out.clear();
  std::vector<int> idx;
  for (int i = 0; i < st.n; ++i)
    if (st.ok[i]) idx.push_back(i);
  VdParams P;
  P.k0 = s.vdK0;
  P.tube = s.vdTube;
  P.gridChi2 = s.vdGridChi2;
  P.vmin = s.vdVmin;
  P.nClimb = s.vdNClimb;
  P.maxStarts = s.vdMaxStarts;
  P.fwdMargin = s.vdFwdMargin;
  P.memberChi2 = s.vdMemberChi2;
  P.gridPvChi2 = s.vdGridPvChi2;
  P.gridCos = s.vdGridCos;
  P.R0 = s.vdR0;
  P.slChi2 = s.slChi2;
  P.slSymmetric = s.slSymmetric;
  P.slForward = s.slForward;
  std::vector<VdStart> starts;
  VdSupplyInfo info;
  VdMaxima M;
  const LinStatus vs = vdSupply(st.t, idx, st.rmax, pvx, pvC, P, starts, info, &M);
  if (vs != kLinOk || (info.maxima && JOINTVTX_FAIL_AT(kFsVdSupply))) {
    st.failed = true;
    return;
  }
  if (!info.maxima) return;
  tr.vdMaxima(M);
  for (const VdStart& q : starts) {
    Vtx v;
    const bool got = candidate(ev, s, st, q.trk, q.x, pvx, pvC, v, tr);
    if (st.failed) return;
    tr.candidate(q.trk, got ? &v : nullptr);
    if (got) out.push_back(v);
  }
}

// stages 2-3: candidate supply (vertex-density maxima), merging, then the acceptance window
template <class Tr>
inline void supply(const Event& ev, const Settings& s, EventState& st, Tr& tr) {
  const Vec3 pvx = st.pvx;
  const Mat3 pvC = st.pvC;
  supplyVd(ev, s, st, pvx, pvC, st.cands, tr);
  if (st.failed) return;
  const int nIn = static_cast<int>(st.cands.size());
  mergeCandidates(st, s, st.cands);
  if (st.failed) return;
  tr.merge(nIn, st.cands);
  std::vector<Vtx> kept;
  for (const Vtx& v : st.cands)
    if (accept(st, s, v.trk, v.x, v.C, v.cF, pvx, pvC)) kept.push_back(v);
  st.cands = kept;
}

// V0 class: candidates at the reference PV with its denominators and significances; their tracks enter the fit list
inline void buildV0Set(const Event& ev, const Settings& s, EventState& st, std::vector<char>& inc) {
  std::vector<double> Aref;
  if (v0RefDenominators(s, st.t, st.ok, st.pvx, st.s3, Aref) != kLinOk || JOINTVTX_FAIL_AT(kFsV0RefDen)) {
    st.failed = true;
    return;
  }
  V0CandInfo inf;
  const std::vector<V0Cand> c = AlephJointVtx::v0Candidates(s, st.t, ev.q, st.jetOf, st.ok, st.rmax, st.pvx, st.pvC,
                                                            Aref, st.ptHalfCfg, st.s3, &inf);
  if (inf.chi2Failed || inf.fitFailed) {
    st.failed = true;
    return;
  }
  for (const V0Cand& k : c) inc[k.a] = inc[k.b] = 1;
  st.v0.init(s, c);
}

// fit list: S0 tracks inside the PV pre-selection window, candidate tracks, V0-candidate tracks
inline void preWindowMask(const Event& ev, const Settings& s, const EventState& st, std::vector<char>& inc) {
  for (int i = 0; i < st.n; ++i)
    if (st.ok[i] && preWindow(ev.par[i], ev.bs, s.pvWindowD0, s.pvWindowZ0)) inc[i] = 1;
}

inline void buildFitList(EventState& st, const std::vector<char>& inc) {
  st.fitIdx.clear();
  st.loc.assign(st.n, -1);
  for (int i = 0; i < st.n; ++i)
    if (inc[i]) {
      st.loc[i] = static_cast<int>(st.fitIdx.size());
      st.fitIdx.push_back(i);
    }
  for (Vtx& v : st.cands)
    for (int& k : v.trk) k = st.loc[k];
  st.tf.clear();
  st.s3f.clear();
  for (int i : st.fitIdx) {
    st.tf.push_back(st.t[i]);
    st.s3f.push_back(st.s3[i]);
  }
}

// hard assignment of the columns of W: argmax over the vertices (first NaN), -1 when the column maximum (NaN
// propagating) is at most hardMin
inline std::vector<int> hardOf(const MatX& W, double hardMin) {
  const int nf = static_cast<int>(W.cols());
  std::vector<int> h(nf, -1);
  for (int j = 0; j < nf; ++j) {
    const int a = argmaxNaN(W.col(j));
    h[j] = (colMaxNaN(W, j) <= hardMin) ? -1 : a;
  }
  return h;
}

// V0 marginals of the last weight-matrix call for the trace (none without V0 candidates)
inline const std::vector<double>* v0Marginals(const Settings& s, const EventState& st) {
  return s.v0Class && st.v0.n() > 0 ? &st.v0.P : nullptr;
}

inline std::vector<double> rowOf(const MatX& W, int r) {
  std::vector<double> v(static_cast<size_t>(W.cols()));
  for (int j = 0; j < W.cols(); ++j) v[j] = W(r, j);
  return v;
}

// soft weights of the fit tracks to the vertices V (0 = PV) at temperature T (Student-t form): PV numerators from
// the track-point chi2 times the lifetime prior; SV numerators from the chi2 with the vertex covariance, zero unless
// chi2 < assocChi2, the vertex inside the track's radius bound and the track direction forward (svForward); PV
// denominator = outlier constant (chi2CutNeg for tracks of negative lifetime sign) + numerator sum, SV denominator =
// SV outlier constant (chi2CutSv) + numerator sum. Wv: weights of the last fit per vertex
inline MatX weightMatrix(const Settings& s, EventState& st, const std::vector<Vtx>& V,
                         const std::vector<std::vector<double>>& Wv, double T) {
  const int nt = static_cast<int>(st.tf.size());
  const int nv = static_cast<int>(V.size());
  MatX num = MatX::Zero(nv, nt);
  const Vec3 pv = V[0].x;
  std::vector<double> ch, dens, cosf, Lo;
  for (int vi = 0; vi < nv; ++vi) {
    const Vec3 d = V[vi].x - pv;
    const double dn = detail::norm3(d);
    Vec3 fdir(1., 0., 0.);
    if (dn > 0)
      for (int k = 0; k < 3; ++k) fdir[k] = d[k] / dn;
    const int mode = vi == 0 ? 0 : 1;
    Mat3 H = Mat3::Zero();
    if (mode && (npInv<3>(V[vi].C, H) == kLinSingular || JOINTVTX_FAIL_AT(kFsWeightInv))) st.failed = true;
    if (chi2Vertex(st.tf, V[vi].x, H, Wv[vi], mode, s.looEps, fdir, ch, dens, cosf, Lo) != kLinOk ||
        JOINTVTX_FAIL_AT(kFsWeightChi2)) {
      st.failed = true;
      ch.assign(nt, JV_NAN);
    }
    if (vi == 0) {
      for (int j = 0; j < nt; ++j) num(0, j) = studentF(ch[j], T, s.studentNu);
      continue;
    }
    const double rv = std::hypot(V[vi].x[0], V[vi].x[1]);
    for (int j = 0; j < nt; ++j) {
      const double eff = ch[j] + s.densK * dens[j];
      double nn = studentF(eff + s.svPen, T, s.studentNu);
      bool allow = ch[j] < s.assocChi2 && rv < st.rmax[st.fitIdx[j]];
      if (!std::isnan(s.svForward)) allow = allow && cosf[j] >= s.svForward;
      if (s.svMinVdet > 0) allow = allow && st.tf[j].nvdet >= s.svMinVdet;
      if (s.svMinP > 0) allow = allow && st.p[st.fitIdx[j]] >= s.svMinP;
      if (!allow) nn = 0.;
      if (s.svPriorS > 0) {
        const double sp = npMaximum(st.s3f[j], 0.);
        nn = nn * (sp * sp / (sp * sp + detail::powSquare(s.svPriorS)));
      }
      num(vi, j) = nn;
    }
  }
  const std::vector<double> pr = lifetimePriors(s, st.s3f, T);
  for (int j = 0; j < nt; ++j) num(0, j) = num(0, j) * pr[j];
  const double cut = studentF(s.chi2Cut, T, s.studentNu);
  const double cutNeg = std::isnan(s.chi2CutNeg) ? cut : studentF(s.chi2CutNeg, T, s.studentNu);
  const bool svCut = !std::isnan(s.chi2CutSv);
  const double cutSv = svCut ? studentF(s.chi2CutSv, T, s.studentNu) : 0.;
  MatX W(nv, nt);
  const bool v0 = s.v0Class && st.v0.n() > 0;
  std::vector<double> A(v0 ? nt : 0);
  for (int j = 0; j < nt; ++j) {
    double sumn = num(0, j);
    for (int vi = 1; vi < nv; ++vi) sumn += num(vi, j);
    const double c0 = st.s3f[j] < 0 ? cutNeg : cut;
    W(0, j) = num(0, j) / (c0 + sumn);
    for (int vi = 1; vi < nv; ++vi) W(vi, j) = num(vi, j) / ((svCut ? cutSv : c0) + sumn);
    if (v0) A[j] = c0 + sumn;
  }
  if (v0) st.v0.apply(W, A, T, num, V);
  return W;
}

// stage 4: annealed joint fit of the PV (start pvx, pvC) and the candidates (local track indices). Per temperature
// of temps + [last]: weights at the current positions; PV refit with the weights at max(T, pvTmin); SV refits
// (weights above svRefitWmin, weight sum >= 1.5, not frozen); merging of SVs into the PV or into each other
// (weight-averaged position, covariance and weights of the first). Then at most three rounds dropping SVs with fewer
// than minSvTracks hard tracks or failing the acceptance window or the final score (each round refits the PV).
// Outputs: PV, SVs with their hard tracks, W, and the weights of the last PV refit
template <class Tr>
inline void jointFit(const Event& ev, const Settings& s, EventState& st, const std::vector<Vtx>& cands,
                     const Vec3& pvx, const Mat3& pvC, Tr& tr) {
  const int nt = static_cast<int>(st.tf.size());
  std::vector<Vtx> V;
  Vtx p0;
  p0.x = pvx;
  p0.C = pvC;
  p0.cF = false;
  V.push_back(p0);
  std::vector<std::vector<double>> Wv(1, std::vector<double>(nt, 0.));
  for (const Vtx& c : cands) {
    Vtx v;
    v.x = c.x;
    v.C = c.C;
    v.trk = c.trk;
    v.cF = false;
    V.push_back(v);
    std::vector<double> w0(nt, 0.);
    for (int k : c.trk) w0[k] = 1.;
    Wv.push_back(w0);
  }
  const Mat3 Binv = beamInfo(s.bsSig);
  const double msig2 = detail::powSquare(s.mergeSig);
  const int nT = static_cast<int>(s.temps.size());
  const double Tlast = s.temps[nT - 1];
  MatX W;
  std::vector<double> pvw;
  auto refitPv = [&](double T) {
    MatX Wt;
    const MatX* Wp = &W;
    if (!std::isnan(s.pvTmin) && T < s.pvTmin) {
      Wt = weightMatrix(s, st, V, Wv, s.pvTmin);
      tr.wm(kWmPvTmin, s.pvTmin, V, Wt, v0Marginals(s, st));
      Wp = &Wt;
    }
    std::vector<int> sel;
    std::vector<double> w;
    for (int j = 0; j < nt; ++j)
      if ((*Wp)(0, j) > 1e-4) {
        sel.push_back(j);
        w.push_back((*Wp)(0, j));
      }
    std::vector<double> Ls(nt, 0.);
    const VtxFit f = fitVertexNb(st.tf, sel, w, V[0].x, Binv, ev.bs, Ls, JV_CAND_FIT_NITER, JV_CAND_FIT_TOL);
    if (!f.ok() || JOINTVTX_FAIL_AT(kFsPvRefit)) st.failed = true;
    if (allFinite(f.x)) {
      V[0].x = f.x;
      V[0].C = f.C;
      V[0].cF = true;
      for (int j = 0; j < nt; ++j) {
        const double v = (*Wp)(0, j);
        Wv[0][j] = v <= 1e-4 ? 0. : v;
      }
    }
    pvw = Wv[0];
  };
  for (int it = 0; it <= nT; ++it) {
    const double T = it < nT ? s.temps[it] : Tlast;
    const bool frozen = it >= nT - (s.freezeLast ? s.freezeN : 0);
    W = weightMatrix(s, st, V, Wv, T);
    tr.wm(kWmAnneal, T, V, W, v0Marginals(s, st));
    if (it == nT) break;
    refitPv(T);
    if (!frozen) {
      for (size_t vi = 1; vi < V.size(); ++vi) {
        std::vector<int> sel;
        std::vector<double> w;
        for (int j = 0; j < nt; ++j)
          if (W(vi, j) > s.svRefitWmin) {
            sel.push_back(j);
            w.push_back(W(vi, j));
          }
        if (npSum(w) < 1.5 || sel.size() < 2) continue;
        std::vector<double> Ls(nt, 0.);
        const VtxFit f = fitVertexNb(st.tf, sel, w, V[vi].x, Mat3::Zero(), Vec3::Zero(), Ls, JV_CAND_FIT_NITER,
                                     JV_CAND_FIT_TOL);
        if (!f.ok() || JOINTVTX_FAIL_AT(kFsSvRefit)) st.failed = true;
        if (allFinite(f.x)) {
          V[vi].x = f.x;
          V[vi].C = f.C;
          V[vi].cF = true;
          std::vector<double> wn(nt, 0.);
          for (size_t k = 0; k < sel.size(); ++k) wn[sel[k]] = w[k];
          Wv[vi] = wn;
        }
      }
    }
    std::vector<char> keep(V.size(), 1);
    for (size_t a = 1; a < V.size(); ++a) {
      if (!keep[a]) continue;
      if (dsig2(V[a].x, V[a].C, V[0].x, V[0].C) < msig2) {
        keep[a] = 0;
        continue;
      }
      for (size_t b = a + 1; b < V.size(); ++b) {
        if (keep[b] && dsig2(V[a].x, V[a].C, V[b].x, V[b].C) < msig2) {
          const double wa = npSum(rowOf(W, static_cast<int>(a))), wb = npSum(rowOf(W, static_cast<int>(b)));
          const double den = std::max(wa + wb, 1e-9);
          for (int k = 0; k < 3; ++k) V[a].x[k] = (wa * V[a].x[k] + wb * V[b].x[k]) / den;
          keep[b] = 0;
        }
      }
    }
    std::vector<Vtx> Vk;
    std::vector<std::vector<double>> Wk;
    for (size_t a = 0; a < V.size(); ++a)
      if (keep[a]) {
        Vk.push_back(V[a]);
        Wk.push_back(Wv[a]);
      }
    V.swap(Vk);
    Wv.swap(Wk);
  }
  for (int round = 0; round < 3; ++round) {
    const std::vector<int> hard = hardOf(W, s.hardMin);
    std::vector<size_t> good = {0};
    for (size_t vi = 1; vi < V.size(); ++vi) {
      std::vector<int> te;
      for (int j = 0; j < nt; ++j)
        if (hard[j] == static_cast<int>(vi)) te.push_back(st.fitIdx[j]);
      if (static_cast<int>(te.size()) < s.minSvTracks) continue;
      if (s.acceptFinal && !accept(st, s, te, V[vi].x, V[vi].C, V[vi].cF, V[0].x, V[0].C)) continue;
      const bool ok = finalOkVertex(st, s, te, V[vi], V[0], st.s3, tr);
      tr.finalOk(te, ok);
      if (ok) good.push_back(vi);
    }
    if (good.size() == V.size()) break;
    std::vector<Vtx> Vg;
    std::vector<std::vector<double>> Wg;
    for (size_t g : good) {
      Vg.push_back(V[g]);
      Wg.push_back(Wv[g]);
    }
    V.swap(Vg);
    Wv.swap(Wg);
    W = weightMatrix(s, st, V, Wv, Tlast);
    tr.wm(kWmDrop, Tlast, V, W, v0Marginals(s, st));
    refitPv(Tlast);
    W = weightMatrix(s, st, V, Wv, Tlast);
    tr.wm(kWmDrop, Tlast, V, W, v0Marginals(s, st));
  }
  if (s.v0Class) {
    // the last call may have been the PV refit at pvTmin: V0 state of the last temperature
    W = weightMatrix(s, st, V, Wv, Tlast);
    tr.wm(kWmFinal, Tlast, V, W, v0Marginals(s, st));
  }
  const std::vector<int> hard = hardOf(W, s.hardMin);
  for (size_t vi = 1; vi < V.size(); ++vi) {
    V[vi].trk.clear();
    for (int j = 0; j < nt; ++j)
      if (hard[j] == static_cast<int>(vi)) V[vi].trk.push_back(j);
  }
  st.px = V[0].x;
  st.pC = V[0].C;
  st.V.assign(V.begin() + 1, V.end());
  st.W = W;
  st.pvw = pvw;
  tr.stage4(st.px, st.pC, st.V, st.W);
}

// seed pairs of the second pass: two-track fits of the tracks cand (event indices) within seedDr in eta-phi, with
// chi2 < seedChi2, inside both radius bounds, resolved from the PV, flying forward, not V0-like and inside the
// acceptance window; ordered by chi2 (stable)
inline std::vector<Vtx> seedPairs(const Event& ev, const Settings& s, EventState& st, const std::vector<int>& cand,
                                  const Vec3& pvx, const Mat3& pvC) {
  std::vector<Vtx> out;
  for (size_t a = 0; a < cand.size(); ++a)
    for (size_t b = a + 1; b < cand.size(); ++b) {
      const int i = cand[a], j = cand[b];
      if (s.accept && etaPhiDr(st.t[i].p, st.t[j].p) > s.seedDr) continue;
      if (s.seedSameJet && st.jetOf[i] != st.jetOf[j]) continue;
      const PairFit r = pairFit(st.t, i, j, st.rmax[i], st.rmax[j], nullptr, s.seedCrossings);
      if (r.status != kLinOk || JOINTVTX_FAIL_AT(kFsSeedPair)) {
        st.failed = true;
        return {};
      }
      if (!r.found) continue;
      if (!(r.chi2 < s.seedChi2)) continue;
      if (std::hypot(r.x[0], r.x[1]) >= std::min(st.rmax[i], st.rmax[j])) continue;
      if (dsig2(r.x, r.C, pvx, pvC) < detail::powSquare(s.seedFsig)) continue;
      const std::vector<int> ij = {i, j};
      const Vec3 d = r.x - pvx;
      if (dotSeq(d, massAt(st, ij, r.x).P) <= 0) continue;
      if (s.v0Veto && v0Like(ev.q[i], ev.q[j], trackMomentum(st, i, r.x), trackMomentum(st, j, r.x), s.ksWin, s.lamWin,
                             s.eeWin))
        continue;
      if (s.seedWindow && !accept(st, s, ij, r.x, r.C, true, pvx, pvC)) continue;
      Vtx v;
      v.x = r.x;
      v.C = r.C;
      v.trk = ij;
      v.chi2 = r.chi2;
      out.push_back(v);
    }
  std::stable_sort(out.begin(), out.end(), [](const Vtx& A, const Vtx& B) { return A.chi2 < B.chi2; });
  return out;
}

// growth of the seeds with the pool tracks (event indices) nearer to the seed than to the PV (chi2 < growChi2),
// nearest first; a track is kept if its chi2 < growChi2, chi2 / ndf < vtxChi2ndf, the vertex lies inside every
// radius bound and the mass stays below massMax; at most accMaxTrk tracks. Then merging and the acceptance window
template <class Tr>
inline std::vector<Vtx> growSeeds(const Settings& s, EventState& st, const std::vector<Vtx>& seeds,
                                  const std::vector<int>& pool, const Vec3& pvx, const Mat3& pvC, Tr& tr) {
  std::vector<Vtx> cands;
  std::vector<double> chi2pv, L;
  if (chi2ToPoint(st.t, pool, pvx, chi2pv, L) != kLinOk || JOINTVTX_FAIL_AT(kFsGrowChi2Pv)) {
    st.failed = true;
    return {};
  }
  std::vector<char> used(st.n, 0);
  for (const Vtx& seed : seeds) {
    Vtx v = seed;
    if (used[v.trk[0]] && used[v.trk[1]]) continue;
    std::vector<double> ch, Lp;
    if (chi2ToPoint(st.t, pool, v.x, ch, Lp) != kLinOk || JOINTVTX_FAIL_AT(kFsGrowChi2Seed)) {
      st.failed = true;
      return {};
    }
    const double rad = std::hypot(v.x[0], v.x[1]);
    std::vector<int> pk;
    std::vector<double> pc;
    for (size_t k = 0; k < pool.size(); ++k) {
      const bool in = std::find(v.trk.begin(), v.trk.end(), pool[k]) != v.trk.end();
      if (ch[k] < s.growChi2 && rad < st.rmax[pool[k]] && !in && ch[k] < chi2pv[k]) {
        pk.push_back(pool[k]);
        pc.push_back(ch[k]);
      }
    }
    const std::vector<int> order = stableArgsort(pc);
    int nties = 0;
    for (size_t k = 1; k < order.size(); ++k)
      if (pc[order[k]] == pc[order[k - 1]]) ++nties;
    if (nties) tr.ties(0, nties);
    for (int o : order) {
      const int k = pk[o];
      if (s.accept && static_cast<int>(v.trk.size()) >= s.accMaxTrk) break;
      std::vector<int> trial = v.trk;
      trial.push_back(k);
      const VtxFit f = fitVertex(st.t, trial, v.x, nullptr, nullptr, JV_CAND_FIT_NITER, JV_CAND_FIT_TOL);
      if (!f.ok() || JOINTVTX_FAIL_AT(kFsGrowFit)) {
        st.failed = true;
        return {};
      }
      const double ndf = static_cast<double>(2 * static_cast<int>(trial.size()) - 3);
      const double tot = npSum(f.chi2);
      if (f.chi2.back() < s.growChi2 && tot / ndf < s.vtxChi2ndf && std::hypot(f.x[0], f.x[1]) < minAt(st.rmax, trial) &&
          massAt(st, trial, f.x).m < s.massMax) {
        v.x = f.x;
        v.C = f.C;
        v.trk = trial;
        v.chi2 = tot;
        v.cF = true;
      }
    }
    for (int k : v.trk) used[k] = 1;
    cands.push_back(v);
  }
  const int nIn = static_cast<int>(cands.size());
  mergeCandidates(st, s, cands);
  if (st.failed) return {};
  tr.merge(nIn, cands);
  std::vector<Vtx> kept;
  for (const Vtx& v : cands)
    if (accept(st, s, v.trk, v.x, v.C, v.cF, pvx, pvC)) kept.push_back(v);
  return kept;
}

// stage 5: up to nIter passes over the tracks left unassigned or significantly displaced but on the PV (s3 at the
// reference PV above iterSig), with s3 above seedSig: seed pairs -> growth -> candidate score -> joint fit of the
// current SVs and the new candidates from the current PV
template <class Tr>
inline void secondPass(const Event& ev, const Settings& s, EventState& st, Tr& tr) {
  const int nf = static_cast<int>(st.fitIdx.size());
  for (int itr = 0; itr < s.nIter; ++itr) {
    const std::vector<int> hard = hardOf(st.W, s.hardMin);
    std::vector<int> left;
    for (int j = 0; j < nf; ++j) {
      const double sj = st.s3[st.fitIdx[j]];
      bool lm = hard[j] < 0 || (hard[j] == 0 && sj > s.iterSig);
      if (s.v0Class && st.v0.n() > 0) lm = lm && st.v0.share[j] <= s.hardMin;   // tracks held by a V0 seed no SV
      if (lm && sj > s.seedSig) left.push_back(st.fitIdx[j]);
    }
    std::vector<Vtx> none;
    if (left.size() < 2) {
      tr.iter(itr, left, 0, none);
      break;
    }
    const std::vector<Vtx> seeds = seedPairs(ev, s, st, left, st.px, st.pC);
    if (st.failed) return;
    if (seeds.empty()) {
      tr.iter(itr, left, 0, none);
      break;
    }
    std::vector<Vtx> grown = growSeeds(s, st, seeds, left, st.px, st.pC, tr);
    if (st.failed) return;
    const double th5 = s.slScoreSeeded ? s.slScoreMin : -1e9;
    std::vector<Vtx> newc;
    if (th5 > JV_SCORE_OFF) {
      for (const Vtx& v : grown) {
        const double sc = vertexScore(st, v.trk, v.x, v.C, v.cF, v.chi2, st.px, st.pC, st.s3, 0);
        tr.score(0, v.trk, sc, s.slScoreMin);
        if (sc >= th5) newc.push_back(v);
      }
    } else {
      newc = grown;
    }
    tr.iter(itr, left, static_cast<int>(seeds.size()), newc);
    if (newc.empty()) break;
    std::vector<Vtx> allc = st.V;
    for (const Vtx& v : newc) {
      Vtx u;
      u.x = v.x;
      u.C = v.C;
      for (int k : v.trk)
        if (st.loc[k] >= 0) u.trk.push_back(st.loc[k]);
      allc.push_back(u);
    }
    const Vec3 px = st.px;
    const Mat3 pC = st.pC;
    jointFit(ev, s, st, allc, px, pC, tr);
    if (st.failed) return;
  }
}

// hard assignment of the final weights
inline void hardAssign(const Settings& s, EventState& st) { st.hard = hardOf(st.W, s.hardMin); }

// B -> D cascade split: each final SV with splitNmin..splitNmax hard tracks is split into the bipartition (both
// halves >= 2 tracks, track 0 always in the first half) of lowest summed chi2 among those with chi2 / ndf below
// vtxChi2ndf per half, halves resolved from each other (splitSig) and from the PV (splitFsig) and the far half
// downstream of the near one (splitOrder); accepted if it lowers the one-vertex chi2 by more than splitDchi2, each
// half lies inside its tracks' radius bounds, flies forward and (splitFinal) passes the final score. The second
// half is appended as a new SV
template <class Tr>
inline void splitCascades(const Event&, const Settings& s, EventState& st, Tr& tr) {
  std::vector<Vtx>& V = st.V;
  std::vector<int>& hard = st.hard;
  const Vec3 px = st.px;
  const Mat3 pC = st.pC;
  Vtx pvV;
  pvV.x = px;
  pvV.C = pC;
  const size_t nv0 = V.size();
  const int nf = static_cast<int>(hard.size());
  auto ndfOf = [](size_t n) { return static_cast<double>(std::max(2 * static_cast<int>(n) - 3, 1)); };
  auto geomOk = [&](const Vec3& x1, const Mat3& C1, const Vec3& x2, const Mat3& C2) {
    if (dsig2(x1, C1, x2, C2) <= detail::powSquare(s.splitSig) ||
        dsig2(x1, C1, px, pC) <= detail::powSquare(s.splitFsig) || dsig2(x2, C2, px, pC) <= detail::powSquare(s.splitFsig))
      return false;
    const double e1 = detail::norm3(Vec3(x1 - px)), e2 = detail::norm3(Vec3(x2 - px));
    const Vec3& xn = e1 <= e2 ? x1 : x2;
    const Vec3& xf = e1 <= e2 ? x2 : x1;
    const Vec3 dd = xf - xn, dn = xn - px;
    return !(dotSeq(dd, dn) < s.splitOrder * detail::norm3(dd) * detail::norm3(dn));
  };
  for (size_t k = 0; k < nv0; ++k) {
    std::vector<int> loc;
    for (int j = 0; j < nf; ++j)
      if (hard[j] == static_cast<int>(k) + 1) loc.push_back(j);
    const int n = static_cast<int>(loc.size());
    if (n < s.splitNmin || n > s.splitNmax) continue;
    std::vector<int> te(n);
    for (int a = 0; a < n; ++a) te[a] = st.fitIdx[loc[a]];
    const Vec3 x0 = V[k].x;
    const VtxFit fa = fitVertex(st.t, te, x0, nullptr, nullptr, JV_CAND_FIT_NITER, JV_CAND_FIT_TOL);
    if (!fa.ok() || JOINTVTX_FAIL_AT(kFsSplitFit)) {
      st.failed = true;
      return;
    }
    if (!allFinite(fa.x)) continue;
    const double chall = npSum(fa.chi2);
    bool have = false;
    double btot = 0., bc1 = 0., bc2 = 0.;
    std::vector<int> bia, bib;
    VtxFit b1, b2;
    for (int na = 1; na < n - 1; ++na) {
      if (na + 1 < 2 || n - na - 1 < 2) continue;
      std::vector<int> comb(na);
      for (int a = 0; a < na; ++a) comb[a] = a + 1;
      while (true) {
        std::vector<int> ia = {0}, ib;
        ia.insert(ia.end(), comb.begin(), comb.end());
        for (int a = 0; a < n; ++a)
          if (std::find(ia.begin(), ia.end(), a) == ia.end()) ib.push_back(a);
        std::vector<int> ta, tb;
        for (int a : ia) ta.push_back(te[a]);
        for (int a : ib) tb.push_back(te[a]);
        const VtxFit f1 = fitVertex(st.t, ta, x0, nullptr, nullptr, JV_CAND_FIT_NITER, JV_CAND_FIT_TOL);
        const VtxFit f2 = fitVertex(st.t, tb, x0, nullptr, nullptr, JV_CAND_FIT_NITER, JV_CAND_FIT_TOL);
        if (!f1.ok() || !f2.ok() || JOINTVTX_FAIL_AT(kFsSplitHalfFit)) {
          st.failed = true;
          return;
        }
        if (allFinite(f1.x) && allFinite(f2.x)) {
          const double s1 = npSum(f1.chi2), s2 = npSum(f2.chi2);
          const double tot = s1 + s2;
          if (!(s1 / ndfOf(ia.size()) >= s.vtxChi2ndf || s2 / ndfOf(ib.size()) >= s.vtxChi2ndf) &&
              geomOk(f1.x, f1.C, f2.x, f2.C) && (!have || tot < btot)) {
            have = true;
            btot = tot;
            bia = ia;
            bib = ib;
            b1 = f1;
            b2 = f2;
            bc1 = s1;
            bc2 = s2;
          }
        }
        int p = na - 1;
        while (p >= 0 && comb[p] == n - na + p) --p;
        if (p < 0) break;
        ++comb[p];
        for (int q = p + 1; q < na; ++q) comb[q] = comb[q - 1] + 1;
      }
    }
    if (!have) continue;
    if (chall - btot <= s.splitDchi2) continue;
    if (bc1 / ndfOf(bia.size()) >= s.vtxChi2ndf || bc2 / ndfOf(bib.size()) >= s.vtxChi2ndf) continue;
    if (!geomOk(b1.x, b1.C, b2.x, b2.C)) continue;
    std::vector<int> ta, tb;
    for (int a : bia) ta.push_back(te[a]);
    for (int a : bib) tb.push_back(te[a]);
    if (std::hypot(b1.x[0], b1.x[1]) >= minAt(st.rmax, ta) || std::hypot(b2.x[0], b2.x[1]) >= minAt(st.rmax, tb))
      continue;
    bool ok = true;
    for (int h = 0; h < 2; ++h) {
      const VtxFit& f = h == 0 ? b1 : b2;
      const std::vector<int>& th = h == 0 ? ta : tb;
      const Vec3 d = f.x - px;
      if (dotSeq(d, massAt(st, th, f.x).P) <= 0) ok = false;
      if (ok && s.splitFinal && s.slFinalMin > JV_SCORE_OFF) {
        Vtx vh;
        vh.x = f.x;
        vh.C = f.C;
        ok = finalOkVertex(st, s, th, vh, pvV, st.s3, tr);
        if (st.failed) return;
        tr.finalOk(th, ok);
      }
    }
    if (!ok) continue;
    Vtx va, vb;
    va.x = b1.x;
    va.C = b1.C;
    va.chi2 = bc1;
    for (int a : bia) va.trk.push_back(loc[a]);
    vb.x = b2.x;
    vb.C = b2.C;
    vb.chi2 = bc2;
    for (int a : bib) vb.trk.push_back(loc[a]);
    V[k] = va;
    V.push_back(vb);
    for (int a : bib) hard[loc[a]] = static_cast<int>(V.size());
  }
  tr.split(V, hard);
}

template <class Tr>
inline void fillResult(const Event&, const Settings& s, EventState& st, Result& r, Tr&) {
  const int n = st.n;
  r.pv = st.px;
  r.pvC = st.pC;
  r.wpv.assign(n, 0.);
  r.assign.assign(n, -1);
  r.wass.assign(n, 0.);
  r.tv0.assign(n, -1);
  const int nf = static_cast<int>(st.fitIdx.size());
  for (int j = 0; j < nf; ++j) {
    const int i = st.fitIdx[j];
    r.wpv[i] = st.pvw.empty() ? st.W(0, j) : st.pvw[j];
    r.assign[i] = st.hard[j];
    r.wass[i] = colMaxNaN(st.W, j);
  }
  r.svs.clear();
  for (size_t k = 0; k < st.V.size(); ++k) {
    SvOut o;
    o.x = st.V[k].x;
    o.C = st.V[k].C;
    for (int j = 0; j < nf; ++j)
      if (st.hard[j] == static_cast<int>(k) + 1) o.trk.push_back(st.fitIdx[j]);
    if (!o.trk.empty()) {
      const InvMass mp = massAt(st, o.trk, o.x);
      o.mass = mp.m;
      o.P = mp.P;
    }
    r.svs.push_back(o);
  }
  if (significances(st, st.px, st.pC, r.s3) != kLinOk || JOINTVTX_FAIL_AT(kFsOutSignif)) st.failed = true;
  r.v0s.clear();
  if (s.v0Class) st.v0.outputs(s.v0OutMin, s.hardMin, n, r.v0s, r.tv0);
  r.status = kResultOk;
}

}  // namespace detail

// one event; tr receives the trace hooks (detail::NoTrace for none)
template <class Tr>
inline Result runEventT(const Event& ev, const Settings& s, Tr& tr) {
  EventState st;
  auto stop = [&tr](int status) {   // a failed step ends the event: no vertices, the status names the step
    Result f;
    f.status = status;
    tr.out(f);
    return f;
  };
  detail::initState(ev, s, st);
  detail::stage0(ev, s, st);
  tr.stage0(st);
  detail::referencePv(ev, s, st);
  if (st.failed) return stop(kResultFailRefPv);
  tr.refPv(st);
  if (detail::significances(st, st.pvx, st.pvC, st.s3) != kLinOk || JOINTVTX_FAIL_AT(kFsSignif))
    return stop(kResultFailSignif);
  tr.s3(st);
  detail::supply(ev, s, st, tr);
  if (st.failed) return stop(kResultFailSupply);
  tr.supply(st);
  std::vector<char> inc(st.n, 0);
  detail::preWindowMask(ev, s, st, inc);
  for (const Vtx& v : st.cands)
    for (int k : v.trk) inc[k] = 1;
  if (s.v0Class) {
    detail::buildV0Set(ev, s, st, inc);
    if (st.failed) return stop(kResultFailV0);
    tr.v0cand(st);
  }
  detail::buildFitList(st, inc);
  if (s.v0Class) st.v0.localize(st.loc);
  tr.fitIdx(st);
  detail::jointFit(ev, s, st, st.cands, st.pvx, st.pvC, tr);
  if (st.failed) return stop(kResultFailJoint);
  detail::secondPass(ev, s, st, tr);
  if (st.failed) return stop(kResultFailSecond);
  detail::hardAssign(s, st);
  if (s.splitCasc && !st.V.empty()) detail::splitCascades(ev, s, st, tr);
  if (st.failed) return stop(kResultFailSplit);
  Result r;
  detail::fillResult(ev, s, st, r, tr);
  if (st.failed) return stop(kResultFailOutput);
  tr.out(r);
  return r;
}

inline Result runEvent(const Event& ev, const Settings& s) {
  detail::NoTrace tr;
  return runEventT(ev, s, tr);
}

}  // namespace AlephJointVtx
}  // namespace FCCAnalyses

#endif
