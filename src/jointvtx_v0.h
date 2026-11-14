#ifndef JOINTVTX_V0_H
#define JOINTVTX_V0_H

/*
  V0 class of the joint vertexer: K0S -> pi+ pi-, Lambda / anti-Lambda -> p pi and photon conversions as two-track
  vertex hypotheses inside the joint fit, competing with the PV, the SV candidates and the outlier term for their
  two tracks.

  v0RefDenominators  denominators of the usable tracks at the reference PV (outlier constant + PV numerator, T = 1;
                     +inf for the other tracks): admission of the candidates
  v0Candidates       admitted candidates of the event (event track indices) in order of decreasing reference log
                     odds: pair-kernel rows (K0S, Lambda, anti-Lambda) then conversions, one Lambda assignment per
                     pair, static terms, prior, heavy-flavour term, admission, component-size limit
  V0Set              the candidates inside the joint fit (fit-list track indices after localize): conflict
                     components, twin-rule denominators, log odds at temperature T, exact marginals; apply scales the
                     non-V0 weights and keeps the state of the call (P, share, last denominators); outputs
  Units cm, GeV.
*/

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <tuple>
#include <vector>

#include <Eigen/Dense>

#include "jointvtx_common.h"
#include "jointvtx_event.h"
#include "jointvtx_track.h"
#include "jointvtx_fit.h"
#include "jointvtx_tables.h"
#include "jointvtx_lookup.h"
#include "jointvtx_conv.h"
#include "jointvtx_v0kin.h"

namespace FCCAnalyses {
namespace AlephJointVtx {

// booked PDG code of a hypothesis code (V0Hyp)
inline int v0PdgOfHyp(int hyp) { return hyp == kV0K0S ? 310 : (hyp == kV0Conversion ? 22 : 3122); }

struct V0Cand {
  int a = -1, b = -1;        // event track indices: positive, negative daughter
  int hyp = 0;               // V0Hyp
  int spec = 0;              // JV_LK_K0S, JV_LK_LAMBDA, JV_LK_CONV
  Vec3 x = Vec3::Zero();     // pair-fit (conversions: constrained-fit) vertex and covariance
  Mat3 C = Mat3::Zero();
  double chi2a = 0., chi2b = 0.;   // track terms (conversions with inner hits of both daughters: half the z-side chi2)
  double m = 0., sm = 0.;    // mass of the hypothesis and its error (conversions: m_ee at the tangent point, 0)
  Vec3 P = Vec3::Zero();     // pair momentum and its covariance (conversions: constrained-fit momentum, covariance zero)
  Mat3 CP = Mat3::Zero();
  double pa = 0., pb = 0., p = 0.;
  bool viola = false, violb = false;   // vertex radius beyond the hit-radius bound of the daughter
  double r = 0., qp = 0., fsig = 0., rv = 0.;   // distance / pointing / flight significance to the reference PV, radius
  double chi2cf = 0., zside = 0., cosp = 0.;    // conversions
  bool fitConv = true;       // convergence flag of the kept pair fit (information only)
  int nExt = 0, nJet = 0;    // heavy-flavour features
  double dl = 0.;
  double z = 0., lm = 0., lp = 0., lh = 0., ha = 1., hb = 1.;   // static terms
  double bp = 0., lx = 0., lrRef = 0.;   // prior, heavy-flavour term, log odds at the reference PV
};

struct V0CandInfo {
  int nPair = 0;             // pair-kernel rows
  int nConv = 0;             // conversions beyond the radius bound
  int nBest = 0;             // after one assignment per pair and species
  int nAdmit = 0;            // above the admission threshold
  int nKept = 0;             // after the component-size limit
  int nNoConvFit = 0;        // kept candidates whose pair fit did not converge
  bool chi2Failed = false;   // a track-point chi2 reported a failure status
  bool fitFailed = false;    // a pair fit reported a failure status
};

namespace detail {

inline LinStatus v0Chi2ToPoint(const std::vector<Trk>& t, const std::vector<int>& idx, const Vec3& x,
                               std::vector<double>& ch) {
  std::vector<double> Lo;
  return chi2ToPoint(t, idx, x, ch, Lo);
}

// n_ext, n_jet and the decay-length density of a candidate
inline bool v0HfFeatures(const Settings& s, const std::vector<Trk>& t, const std::vector<int>& jet,
                         const std::vector<char>& ok, const std::vector<double>& s3, V0Cand& k) {
  bool fail = false;
  const std::vector<int> o = hfOthers(ok, s3, s.v0HfSig, k.a, k.b);
  k.nExt = 0;
  k.nJet = 0;
  if (!o.empty()) {
    std::vector<double> ch;
    if (v0Chi2ToPoint(t, o, k.x, ch) != kLinOk || JOINTVTX_FAIL_AT(kFsV0HfChi2)) fail = true;
    std::vector<int> jo(o.size());
    for (size_t i = 0; i < o.size(); ++i) jo[i] = jet[static_cast<size_t>(o[i])];
    const HfCounts hc = hfCounts(ch, jo, jet[static_cast<size_t>(k.a)], jet[static_cast<size_t>(k.b)], s.v0HfChi2);
    k.nExt = hc.nExt;
    k.nJet = hc.nJet;
  }
  k.dl = hfDecayLength(std::hypot(k.P[0], k.P[1]), k.r, k.spec);
  return !fail;
}

// static terms, prior and heavy-flavour term of a candidate (the track chi2 of a conversion may be replaced)
inline void v0StaticOf(const Settings& s, V0Cand& k) {
  V0StaticIn in;
  in.hyp = k.hyp;
  in.m = k.m;
  in.sm = k.sm;
  in.qp = k.qp;
  in.pa = k.pa;
  in.pb = k.pb;
  in.cosp = k.cosp;
  in.rv = k.rv;
  in.zside = k.zside;
  in.chi2a = k.chi2a;
  in.chi2b = k.chi2b;
  in.viola = k.viola;
  in.violb = k.violb;
  V0StaticCfg cfg;
  cfg.smMin = s.v0SmMin;
  cfg.mass = s.v0Mass;
  cfg.hitScale = s.v0HitScale;
  cfg.convHscale = s.v0ConvHscale;
  const V0StaticOut o = staticTerms(*s.bundle, in, cfg);
  k.z = o.z;
  k.lm = o.lm;
  k.lp = o.lp;
  k.lh = o.lh;
  k.ha = o.ha;
  k.hb = o.hb;
  k.chi2a = o.chi2a;
  k.chi2b = o.chi2b;
  k.bp = priorB(s.bundle->prior, k.spec, k.r, k.fsig);
  k.lx = hfTerm(*s.bundle, k.spec, k.nExt, k.nJet, k.dl);
}

// log odds of a candidate at the reference PV (T = 1) against the reference denominators Aref (event tracks)
inline double v0LogRRef(const Settings& s, const V0Cand& k, const std::vector<double>& Aref) {
  const double lr = k.bp + trackLogF(k.chi2a / powSquare(k.ha), 1., s.v0Nu) +
                    trackLogF(k.chi2b / powSquare(k.hb), 1., s.v0Nu) - std::log(Aref[static_cast<size_t>(k.a)]) -
                    std::log(Aref[static_cast<size_t>(k.b)]) + k.lm + k.lp + k.lh;
  return lr + k.lx;
}

// at most maxComp candidates per connected component, in the given order (decreasing reference log odds)
inline std::vector<V0Cand> v0ComponentLimit(const std::vector<V0Cand>& in, int maxComp) {
  std::vector<V0Cand> keep;
  std::vector<int> cnt;
  std::map<int, int> compOf;
  for (const V0Cand& k : in) {
    const auto ia = compOf.find(k.a), ib = compOf.find(k.b);
    const bool ha = ia != compOf.end(), hb = ib != compOf.end();
    const int ca = ha ? ia->second : -1, cb = hb ? ib->second : -1;
    int cid = ha ? ca : (hb ? cb : static_cast<int>(cnt.size()));
    if (ha && hb && ca != cb) {
      const int n = cnt[static_cast<size_t>(ca)] + cnt[static_cast<size_t>(cb)];
      if (n + 1 > maxComp) continue;
      for (auto& tg : compOf)
        if (tg.second == cb) tg.second = ca;
      cnt[static_cast<size_t>(ca)] = n;
      cnt[static_cast<size_t>(cb)] = 0;
      cid = ca;
    }
    const int cur = cid < static_cast<int>(cnt.size()) ? cnt[static_cast<size_t>(cid)] : 0;
    if (cur + 1 > maxComp) continue;
    if (cid == static_cast<int>(cnt.size())) cnt.push_back(0);
    cnt[static_cast<size_t>(cid)] = cur + 1;
    compOf[k.a] = cid;
    compOf[k.b] = cid;
    keep.push_back(k);
  }
  return keep;
}

}  // namespace detail

// denominators of the usable tracks (ok) at the reference PV pvx with the significances s3 there: outlier constant
// (chi2CutNeg for negative s3) + PV numerator times the lifetime prior, T = 1; +inf for the other tracks
inline LinStatus v0RefDenominators(const Settings& s, const std::vector<Trk>& t, const std::vector<char>& ok,
                                   const Vec3& pvx, const std::vector<double>& s3, std::vector<double>& A) {
  A.assign(t.size(), JV_INF);
  std::vector<int> oi;
  for (size_t i = 0; i < ok.size(); ++i)
    if (ok[i]) oi.push_back(static_cast<int>(i));
  if (oi.empty()) return kLinOk;
  std::vector<double> ch;
  const LinStatus st = detail::v0Chi2ToPoint(t, oi, pvx, ch);
  const double cut = studentF(s.chi2Cut, 1., s.studentNu);
  const double cutNeg = std::isnan(s.chi2CutNeg) ? cut : studentF(s.chi2CutNeg, 1., s.studentNu);
  for (size_t k = 0; k < oi.size(); ++k) {
    const size_t i = static_cast<size_t>(oi[k]);
    const double pr = s.lifetimePrior ? lifetimePrior(s3[i], s.priorS) : 1.;
    A[i] = (s3[i] < 0 ? cutNeg : cut) + studentF(ch[k], 1., s.studentNu) * pr;
  }
  return st;
}

namespace detail {

// K0S / Lambda / anti-Lambda rows of the pair kernel with their quantities at the reference PV (both Lambda assignments)
inline std::vector<V0Cand> v0PairCands(const Settings& s, const std::vector<Trk>& t, const std::vector<int>& q,
                                       const std::vector<char>& ok, const std::vector<double>& rmax, const Vec3& pvx,
                                       const Mat3& pvC, double K, bool* failed = nullptr) {
  std::vector<V0Cand> raw;
  V0PairCuts cut;
  cut.rmax = s.v0Rmax;
  cut.dzSig = s.v0DzSig;
  cut.dzAbs = s.v0DzAbs;
  cut.preKs = s.v0PreKs;
  cut.preLam = s.v0PreLam;
  cut.winKs = s.v0WinKs;
  cut.winLam = s.v0WinLam;
  cut.chi2max = s.v0Chi2Max;
  cut.qtKs = s.v0QtMin[0];
  cut.qtLam = s.v0QtMin[1];
  for (const V0Raw& w : v0RawCandidates(t, q, ok, rmax, pvx, pvC, K, cut, failed)) {
    V0Cand k;
    k.a = w.k.a;
    k.b = w.k.b;
    k.hyp = w.k.hyp;
    k.spec = speciesOfHyp(w.k.hyp);
    k.x = w.k.x;
    k.C = w.k.C;
    k.chi2a = w.k.chi2a;
    k.chi2b = w.k.chi2b;
    k.m = w.k.m;
    k.sm = w.k.sm;
    k.P = w.k.P;
    k.CP = w.k.CP;
    k.pa = w.k.pa;
    k.pb = w.k.pb;
    k.p = w.p;
    k.viola = w.viola;
    k.violb = w.violb;
    k.r = w.r;
    k.qp = w.qp;
    k.fsig = w.fsig;
    k.rv = w.rv;
    k.fitConv = w.k.fitConv;
    raw.push_back(k);
  }
  return raw;
}

// conversion candidates (loose selection, radius from the axis >= rmin) with their quantities at the reference PV
inline std::vector<V0Cand> v0ConvCands(const std::vector<Trk>& t, const std::vector<int>& q, const std::vector<char>& ok,
                                       const std::vector<double>& rmax, const Vec3& pvx, const Mat3& pvC, double K,
                                       double rmin) {
  std::vector<V0Cand> raw;
  std::vector<Vec5> par(t.size());
  std::vector<Mat5> cov(t.size());
  for (size_t i = 0; i < t.size(); ++i) {
    par[i] = t[i].p;
    cov[i] = t[i].C;
  }
  for (const ConvCand& c : convRaw(par, cov, q, ok, rmax, pvx, K, rmin)) {
    V0Cand k;
    k.a = c.a;
    k.b = c.b;
    k.hyp = kV0Conversion;
    k.spec = JV_LK_CONV;
    k.x = c.x;
    k.C = c.C;
    k.chi2a = c.chi2a;
    k.chi2b = c.chi2b;
    k.chi2cf = c.chi2cf;
    k.zside = c.zside;
    k.m = c.m;
    k.P = c.P;
    k.pa = c.pa;
    k.pb = c.pb;
    k.p = c.p;
    k.viola = c.viola;
    k.violb = c.violb;
    k.cosp = c.cosp;
    k.r = c.r;
    k.fsig = v0FlightSig(c.x, c.C, pvx, pvC);
    k.rv = c.rv;
    raw.push_back(k);
  }
  return raw;
}


// one assignment per (pair, species), static terms, prior, heavy-flavour term (features already set), reference log
// odds, admission, stable sort by decreasing log odds, component-size limit
inline std::vector<V0Cand> v0Select(const Settings& s, std::vector<V0Cand>& raw, const std::vector<double>& Aref,
                                    V0CandInfo& inf) {
  // one assignment per (pair, species): the smaller |mass pull| replaces the entry in its first position
  std::vector<V0Cand> best;
  std::map<std::tuple<int, int, int>, size_t> pos;
  for (V0Cand& k : raw) {
    v0StaticOf(s, k);
    const std::tuple<int, int, int> key(k.a, k.b, k.spec);
    const auto it = pos.find(key);
    if (it == pos.end()) {
      pos[key] = best.size();
      best.push_back(k);
    } else if (std::fabs(k.z) < std::fabs(best[it->second].z)) {
      best[it->second] = k;
    }
  }
  inf.nBest = static_cast<int>(best.size());
  std::vector<V0Cand> adm;
  for (V0Cand& k : best) {
    k.lrRef = v0LogRRef(s, k, Aref);
    if (k.lrRef >= s.v0Admit) adm.push_back(k);
  }
  inf.nAdmit = static_cast<int>(adm.size());
  std::stable_sort(adm.begin(), adm.end(), [](const V0Cand& u, const V0Cand& v) { return -u.lrRef < -v.lrRef; });
  std::vector<V0Cand> keep = v0ComponentLimit(adm, s.v0MaxComp);
  inf.nKept = static_cast<int>(keep.size());
  for (const V0Cand& k : keep) inf.nNoConvFit += !k.fitConv;
  return keep;
}

}  // namespace detail

// admitted V0 candidates of the event: t / q / jet / ok / rmax / s3 per event track, (pvx, pvC) = reference PV,
// Aref = v0RefDenominators, K = half the momentum constant of the settings
inline std::vector<V0Cand> v0Candidates(const Settings& s, const std::vector<Trk>& t, const std::vector<int>& q,
                                        const std::vector<int>& jet, const std::vector<char>& ok,
                                        const std::vector<double>& rmax, const Vec3& pvx, const Mat3& pvC,
                                        const std::vector<double>& Aref, double K, const std::vector<double>& s3,
                                        V0CandInfo* info = nullptr) {
  V0CandInfo inf;
  std::vector<V0Cand> raw = detail::v0PairCands(s, t, q, ok, rmax, pvx, pvC, K, &inf.fitFailed);
  inf.nPair = static_cast<int>(raw.size());
  if (s.v0Conv) {
    const double rmin = std::isnan(s.v0ConvRmin) ? -JV_INF : s.v0ConvRmin;
    const std::vector<V0Cand> cr = detail::v0ConvCands(t, q, ok, rmax, pvx, pvC, K, rmin);
    raw.insert(raw.end(), cr.begin(), cr.end());
    inf.nConv = static_cast<int>(cr.size());
  }
  for (V0Cand& k : raw)
    if (!detail::v0HfFeatures(s, t, jet, ok, s3, k)) inf.chi2Failed = true;
  std::vector<V0Cand> keep = detail::v0Select(s, raw, Aref, inf);
  if (info) *info = inf;
  return keep;
}

// connected component of candidates sharing tracks: candidate indices g (ascending) and per member the bit mask of
// the members it conflicts with
struct V0Comp {
  std::vector<int> g;
  std::vector<long long> conf;
};

// The V0 candidates of one event inside the joint fit. Every apply overwrites the state (P, share, lastAa / lastAb /
// lastT); P stays zero and share empty before the first apply.
struct V0Set {
  std::vector<V0Cand> k;            // candidates (event track indices)
  std::vector<int> a, b;            // daughters: event indices until localize, then fit-list indices
  std::vector<double> chi2a, chi2b; // track-term chi2 over the squared hit scale
  std::vector<double> stat;         // static terms: mass + pointing + hits + heavy flavour
  std::vector<double> bprior;
  std::vector<V0Comp> comps;
  double nu = JV_NAN, twinSig = JV_NAN, twinExt = JV_NAN;
  // state of the last apply
  std::vector<double> P;            // marginal weight per candidate
  std::vector<double> share;        // summed marginals per fit track
  std::vector<double> lastAa, lastAb;
  double lastT = JV_NAN;
  bool hasLast = false;

  int n() const { return static_cast<int>(k.size()); }

  void init(const Settings& s, const std::vector<V0Cand>& cands) {
    k = cands;
    const size_t m = k.size();
    a.resize(m);
    b.resize(m);
    chi2a.resize(m);
    chi2b.resize(m);
    stat.resize(m);
    bprior.resize(m);
    for (size_t i = 0; i < m; ++i) {
      a[i] = k[i].a;
      b[i] = k[i].b;
      chi2a[i] = k[i].chi2a / detail::powSquare(k[i].ha);
      chi2b[i] = k[i].chi2b / detail::powSquare(k[i].hb);
      stat[i] = k[i].lm + k[i].lp + k[i].lh;
      bprior[i] = k[i].bp;
    }
    for (size_t i = 0; i < m; ++i) stat[i] = stat[i] + k[i].lx;
    nu = s.v0Nu;
    twinSig = s.v0TwinSig;
    twinExt = s.v0TwinExt;
    P.assign(m, 0.);
    share.clear();
    lastAa.clear();
    lastAb.clear();
    lastT = JV_NAN;
    hasLast = false;
    components();
  }

  void localize(const std::vector<int>& loc) {
    for (size_t i = 0; i < a.size(); ++i) {
      a[i] = loc[static_cast<size_t>(a[i])];
      b[i] = loc[static_cast<size_t>(b[i])];
    }
  }

  // per candidate the denominators of its two tracks: A (per fit track) minus the numerators num(vi, .) of the twin
  // SVs V[vi], vi >= 1 (within twinSig standard deviations of the candidate vertex, V0 + SV covariance, and with at
  // most twinExt summed weight num / A on the other tracks); floored at 1e-300. V: elements with x, C
  template <class VL>
  void denominators(const std::vector<double>& A, const MatX& num, const VL& V, std::vector<double>& Aa,
                    std::vector<double>& Ab) const {
    const int m = n();
    Aa.resize(static_cast<size_t>(m));
    Ab.resize(static_cast<size_t>(m));
    for (int i = 0; i < m; ++i) {
      Aa[i] = A[static_cast<size_t>(a[i])];
      Ab[i] = A[static_cast<size_t>(b[i])];
    }
    const int nv = static_cast<int>(V.size());
    if (std::isnan(twinSig) || nv < 2) return;
    const double ts2 = twinSig * twinSig;
    const bool ext = !std::isnan(twinExt);
    const int nt = static_cast<int>(A.size());
    std::vector<double> wo(static_cast<size_t>(nt));
    for (int i = 0; i < m; ++i) {
      const Vec3& xk = k[static_cast<size_t>(i)].x;
      const Mat3& Ck = k[static_cast<size_t>(i)].C;
      for (int vi = 1; vi < nv; ++vi) {
        const Vec3& vx = V[static_cast<size_t>(vi)].x;
        const Mat3& vC = V[static_cast<size_t>(vi)].C;
        Vec3 d;
        Mat3 S;
        for (int r = 0; r < 3; ++r) {
          d[r] = vx[r] - xk[r];
          for (int c = 0; c < 3; ++c) S(r, c) = vC(r, c) + Ck(r, c);
        }
        Vec3 y;
        if (npSolve<3>(S, d, y) == kLinSingular) continue;
        const double qv = dotSeq(d, y);
        if (qv < ts2 && ext) {
          for (int j = 0; j < nt; ++j) wo[j] = num(vi, j) / npMaximum(A[static_cast<size_t>(j)], 1e-300);
          if (npSum(wo) - wo[static_cast<size_t>(a[i])] - wo[static_cast<size_t>(b[i])] > twinExt) continue;
        }
        if (qv < ts2) {
          Aa[i] = Aa[i] - num(vi, a[i]);
          Ab[i] = Ab[i] - num(vi, b[i]);
        }
      }
    }
    for (int i = 0; i < m; ++i) {
      Aa[i] = npMaximum(Aa[i], 1e-300);
      Ab[i] = npMaximum(Ab[i], 1e-300);
    }
  }

  // log odds per candidate at temperature T (static terms over T, prior not tempered)
  std::vector<double> logR(double T, const std::vector<double>& Aa, const std::vector<double>& Ab) const {
    std::vector<double> lr(static_cast<size_t>(n()));
    for (int i = 0; i < n(); ++i)
      lr[i] = bprior[i] + trackLogF(chi2a[i], T, nu) + trackLogF(chi2b[i], T, nu) - std::log(Aa[i]) - std::log(Ab[i]) +
              stat[i] / T;
    return lr;
  }

  // exact marginals over the sets of candidates without a shared track, per component
  std::vector<double> marginals(const std::vector<double>& lr) const {
    std::vector<double> out(static_cast<size_t>(n()), 0.);
    for (const V0Comp& c : comps) {
      if (c.g.size() == 1) {
        const size_t g0 = static_cast<size_t>(c.g[0]);
        out[g0] = 1. / (1. + std::exp(-npClip(lr[g0], -700., 700.)));
      } else {
        std::vector<double> lg(c.g.size());
        for (size_t u = 0; u < c.g.size(); ++u) lg[u] = lr[static_cast<size_t>(c.g[u])];
        const std::vector<double> pg = componentMarginals(lg, c.conf);
        for (size_t u = 0; u < c.g.size(); ++u) out[static_cast<size_t>(c.g[u])] = pg[u];
      }
    }
    return out;
  }

  // marginals at temperature T against the PV-form denominators A (per fit track; num = numerators per vertex and
  // fit track, V = vertices, 0 = PV, for the twin rule); scales every weight of W (vertices x fit tracks) by
  // 1 - min(share, 1) and stores P, share and the denominators
  template <class VL>
  void apply(MatX& W, const std::vector<double>& A, double T, const MatX& num, const VL& V) {
    std::vector<double> Aa, Ab;
    denominators(A, num, V, Aa, Ab);
    lastAa = Aa;
    lastAb = Ab;
    lastT = T;
    hasLast = true;
    P = marginals(logR(T, Aa, Ab));
    const int nt = static_cast<int>(W.cols());
    share.assign(static_cast<size_t>(nt), 0.);
    for (int i = 0; i < n(); ++i) share[static_cast<size_t>(a[i])] += P[i];
    for (int i = 0; i < n(); ++i) share[static_cast<size_t>(b[i])] += P[i];
    for (int j = 0; j < nt; ++j) {
      const double f = 1. - npMinimum(share[static_cast<size_t>(j)], 1.);
      for (int r = 0; r < W.rows(); ++r) W(r, j) = W(r, j) * f;
    }
  }

  // output V0s (weight above outMin) and per event track the index of its V0 when the weight is above hardMin
  void outputs(double outMin, double hardMin, int nTracks, std::vector<V0Out>& v0s, std::vector<int>& tv0) const {
    v0s.clear();
    tv0.assign(static_cast<size_t>(nTracks), -1);
    for (int i = 0; i < n(); ++i) {
      if (!(P[i] > outMin)) continue;
      const V0Cand& c = k[static_cast<size_t>(i)];
      V0Out o;
      o.t1 = c.a;
      o.t2 = c.b;
      o.pdg = v0PdgOfHyp(c.hyp);
      o.hyp = c.hyp;
      o.x = c.x;
      o.C = c.C;
      o.m = c.m;
      o.w = P[i];
      o.p = c.p;
      o.P = c.P;
      v0s.push_back(o);
      if (P[i] > hardMin) tv0[static_cast<size_t>(c.a)] = tv0[static_cast<size_t>(c.b)] = static_cast<int>(v0s.size()) - 1;
    }
  }

 private:
  void components() {
    const int m = n();
    std::vector<int> parent(static_cast<size_t>(m));
    for (int i = 0; i < m; ++i) parent[i] = i;
    auto find = [&parent](int i) {
      while (parent[i] != i) {
        parent[i] = parent[parent[i]];
        i = parent[i];
      }
      return i;
    };
    std::map<int, int> tr;
    for (int i = 0; i < m; ++i)
      for (int t : {a[i], b[i]}) {
        const auto it = tr.find(t);
        if (it != tr.end()) parent[static_cast<size_t>(find(i))] = find(it->second);
        else tr[t] = i;
      }
    std::vector<int> roots;
    std::vector<std::vector<int>> groups;
    for (int i = 0; i < m; ++i) {
      const int r = find(i);
      size_t gi = 0;
      while (gi < roots.size() && roots[gi] != r) ++gi;
      if (gi == roots.size()) {
        roots.push_back(r);
        groups.emplace_back();
      }
      groups[gi].push_back(i);
    }
    comps.clear();
    for (std::vector<int>& g : groups) {
      std::sort(g.begin(), g.end());
      V0Comp c;
      c.g = g;
      c.conf.assign(g.size(), 0);
      for (size_t u = 0; u < g.size(); ++u)
        for (size_t v = 0; v < g.size(); ++v) {
          const int gu = g[u], gv = g[v];
          if (u != v && (a[gu] == a[gv] || a[gu] == b[gv] || b[gu] == a[gv] || b[gu] == b[gv]))
            c.conf[u] |= 1LL << v;
        }
      comps.push_back(c);
    }
  }
};

}  // namespace AlephJointVtx
}  // namespace FCCAnalyses

#endif
