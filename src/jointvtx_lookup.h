#ifndef JOINTVTX_LOOKUP_H
#define JOINTVTX_LOOKUP_H

/*
  Small lookup and score kernels of the joint vertex fitter: bin lookups, V0 prior / hit / pointing /
  heavy-flavour terms, the exact marginals over sets of V0 candidates, the soft-assignment functions, the
  logistic vertex score and the vertex acceptance window. Pure functions of their arguments; every table
  value comes from jointvtx_tables.h.
*/

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include "jointvtx_common.h"
#include "jointvtx_tables.h"
#include "jointvtx_track.h"

namespace FCCAnalyses {
namespace AlephJointVtx {

// species index of the V0 tables
constexpr int JV_LK_K0S = 0;
constexpr int JV_LK_LAMBDA = 1;
constexpr int JV_LK_CONV = 2;

// species of a V0 hypothesis code (V0Hyp)
inline int speciesOfHyp(int hyp) {
  return hyp == kV0K0S ? JV_LK_K0S : (hyp == kV0Conversion ? JV_LK_CONV : JV_LK_LAMBDA);
}

// ---------------------------------------------------------------- tables
inline const PriorSpecies& priorSpecies(const PriorTable& t, int spec) {
  return spec == JV_LK_K0S ? t.k0s : (spec == JV_LK_LAMBDA ? t.lambda : t.conversion);
}
inline const HfSpecies& hfSpecies(const HfTable& t, int spec) {
  return spec == JV_LK_K0S ? t.k0s : (spec == JV_LK_LAMBDA ? t.lambda : t.conversion);
}

// rate of wrongly attached inner hits per daughter of momentum p [GeV]
inline double hitEps(double p) { return JV_HIT_EPS[static_cast<size_t>(binRight(JV_HIT_EPS_P_BINS, p))]; }

// log prior odds of a V0 candidate at transverse distance r [cm] and 3D flight significance fsig from the
// reference PV
inline double priorB(const PriorTable& t, int spec, double r, double fsig) {
  const PriorSpecies& s = priorSpecies(t, spec);
  const size_t rb = static_cast<size_t>(binRight(t.rBins, r));
  if (r < t.fsRmax) return s.fs[rb][static_cast<size_t>(binRight(t.fsBins, fsig))];
  return s.r[rb];
}

// log L_p = -w log(1 + q / nu): PV pointing significance q of a K0S / Lambda candidate
inline double pointTerm(const BeamSpotBundle& b, double q, int spec) {
  if (spec != JV_LK_K0S && spec != JV_LK_LAMBDA) return JV_NAN;
  const double w = b.pointW[static_cast<size_t>(spec)] * b.pvWscale[static_cast<size_t>(spec)];
  return -w * std::log1p(npMaximum(q, 0.) / JV_POINT_NU_SIGNAL);
}

// V0-against-heavy-flavour term: counts nExt, nJet of tracks, dl log density of the decay length
inline double hfTerm(const BeamSpotBundle& b, int spec, int nExt, int nJet, double dl) {
  const HfSpecies& h = hfSpecies(b.hf, spec);
  double lx = h.b0 + h.wdl * dl;
  if (nExt > 0) lx += h.ext[static_cast<size_t>(std::min<int>(nExt, static_cast<int>(h.ext.size())) - 1)];
  if (nJet > 0) lx += h.jet[static_cast<size_t>(std::min<int>(nJet, static_cast<int>(h.jet.size())) - 1)];
  return lx;
}

// log density of the transverse decay length r [cm] for the lifetime of the species (0 for conversions);
// ptV0 = transverse momentum of the pair
inline double hfDecayLength(double ptV0, double r, int spec) {
  if (spec != JV_LK_K0S && spec != JV_LK_LAMBDA) return 0.;
  const double mass = spec == JV_LK_K0S ? JV_M_KS : JV_M_LAM;
  const double ctau = spec == JV_LK_K0S ? JV_CTAU_K0S : JV_CTAU_LAMBDA;
  const double lam = std::max(ptV0 / mass * ctau, 1e-3);
  return -std::log(lam) - r / lam;
}

// counted features of a V0 candidate: tracks 'others' (not the daughters) displaced from the reference PV;
// chi2 = their chi2 to the V0 vertex, jet = their jet index, ja / jb = jets of the daughters
struct HfCounts {
  int nExt = 0;
  int nJet = 0;
};
inline HfCounts hfCounts(const std::vector<double>& chi2, const std::vector<int>& jet, int ja, int jb, double chi2Max) {
  HfCounts c;
  for (size_t i = 0; i < chi2.size(); ++i) {
    if (chi2[i] < chi2Max) ++c.nExt;
    if ((jet[i] == ja && ja >= 0) || (jet[i] == jb && jb >= 0)) ++c.nJet;
  }
  return c;
}

// ---------------------------------------------------------------- static terms of a V0 candidate
struct V0StaticCfg {
  double smMin = JV_V0_SM_MIN;
  bool mass = JV_V0_MASS;
  double hitScale = JV_V0_HIT_SCALE;
  std::array<double, 2> convHscale = JV_V0_CONV_HSCALE;
};

struct V0StaticIn {
  int hyp = kV0K0S;
  double m = 0., sm = 0., qp = 0., pa = 0., pb = 0., cosp = 0., rv = 0., zside = 0., chi2a = 0., chi2b = 0.;
  bool viola = false, violb = false;
};

// z mass pull, log L_m, log L_p, log L_h, hit scales ha, hb of the two daughter terms; chi2a / chi2b: the track
// chi2 of a conversion with hits in the vertex radius is replaced
struct V0StaticOut {
  double z = 0., lm = 0., lp = 0., lh = 0., ha = 1., hb = 1., chi2a = 0., chi2b = 0.;
};

// the inner radius edges of the conversion hit table are the hit-radius bounds of the track selection; the table
// keeps 13.1, one unit in the last place below the sum 11.3 + 1.8 of JV_HITR_VDET1
static_assert(JV_HIT_TABLE_CONV_R_BINS[1] == JV_HITR_VDET2 && JV_HIT_TABLE_CONV_R_BINS[3] == JV_HITR_ITC,
              "conversion hit table: radius edges differ from the hit-radius bounds");
static_assert(JV_HIT_TABLE_CONV_R_BINS[2] < JV_HITR_VDET1 && JV_HITR_VDET1 - JV_HIT_TABLE_CONV_R_BINS[2] < 1e-14,
              "conversion hit table: radius edge differs from the one-VDET-hit bound");

inline V0StaticOut convTerms(const V0StaticIn& k, const V0StaticCfg& c = V0StaticCfg()) {
  V0StaticOut o;
  o.chi2a = k.chi2a;
  o.chi2b = k.chi2b;
  o.lp = -JV_CONV_POINT_W * std::log1p(std::max(1. - k.cosp, 0.) / JV_CONV_POINT_TAU);
  const int n = static_cast<int>(k.viola) + static_cast<int>(k.violb);
  o.lh = JV_HIT_TABLE_CONV_LOG_LR[static_cast<size_t>(binRight(JV_HIT_TABLE_CONV_R_BINS, k.rv))][static_cast<size_t>(n)];
  if (n == 2) {
    o.chi2a = o.chi2b = 0.5 * k.zside;
    o.ha = o.hb = c.convHscale[1];
  } else if (n == 1) {
    if (k.viola) { o.ha = c.convHscale[0]; o.hb = 1.; }
    else { o.ha = 1.; o.hb = c.convHscale[0]; }
  }
  return o;
}

inline V0StaticOut staticTerms(const BeamSpotBundle& b, const V0StaticIn& k, const V0StaticCfg& c = V0StaticCfg()) {
  if (k.hyp == kV0Conversion) return convTerms(k, c);
  V0StaticOut o;
  o.chi2a = k.chi2a;
  o.chi2b = k.chi2b;
  const double m0 = k.hyp == kV0K0S ? JV_M_KS : JV_M_LAM;
  const double sm = std::max(k.sm, c.smMin);
  o.z = (k.m - m0) / sm;
  o.lm = c.mass ? (-0.5 * o.z * o.z - std::log(sm / 0.001)) : 0.;
  o.lp = pointTerm(b, k.qp, speciesOfHyp(k.hyp));
  if (k.viola) { o.lh += std::log(hitEps(k.pa)); o.ha = c.hitScale; }
  if (k.violb) { o.lh += std::log(hitEps(k.pb)); o.hb = c.hitScale; }
  return o;
}

// ---------------------------------------------------------------- soft assignment
// log of the V0 track-term function (Student-t with nu degrees of freedom) at temperature T
inline double trackLogF(double chi2, double T, double nu = JV_V0_NU) {
  return -0.5 * (nu + 2.) * std::log1p(npMinimum(chi2, 1e6) / (nu * T));
}

// soft-assignment function of the track-vertex compatibility chi2 at temperature T (Student-t form)
inline double studentF(double chi2, double T, double nu = JV_STUDENT_NU) {
  return std::pow(1. + npMinimum(chi2, 1e6) / (nu * T), -0.5 * (nu + 2.));
}

// lifetime prior of a track of signed 3D significance s3 (Lorentz form, scale priorS)
inline double lifetimePrior(double s3, double priorS = JV_PRIOR_S) {
  const double sp = npMaximum(s3, 0.) / priorS;
  return 1. / (1. + sp * sp);
}

// squared significance of the distance of two vertices; +inf if the summed covariance is singular
inline double dsig2(const Vec3& x1, const Mat3& C1, const Vec3& x2, const Mat3& C2) {
  const Vec3 d = x1 - x2;
  Mat3 S = C1 + C2;
  Vec3 y;
  const LinStatus st = npSolve<3>(S, d, y);
  if (st == kLinSingular) return JV_INF;
  return dotSeq(d, y);
}

// ---------------------------------------------------------------- marginals over sets of V0 candidates
// exact marginals of the candidates of one component: all sets without conflicts (conf[k] = bit mask of the
// candidates sharing a track with k) weighted by exp(sum logR); n <= 12
inline std::vector<double> componentMarginals(const std::vector<double>& logR, const std::vector<long long>& conf) {
  const int n = static_cast<int>(logR.size());
  const long long nset = 1LL << n;
  std::vector<double> lw(static_cast<size_t>(nset), -JV_INF);
  double mx = -JV_INF;
  for (long long S = 0; S < nset; ++S) {
    bool ok = true;
    double s = 0.;
    for (int k = 0; k < n; ++k) {
      if (S & (1LL << k)) {
        if (conf[static_cast<size_t>(k)] & S) {
          ok = false;
          break;
        }
        s += logR[static_cast<size_t>(k)];
      }
    }
    if (ok) {
      lw[static_cast<size_t>(S)] = s;
      if (s > mx) mx = s;
    }
  }
  double Z = 0.;
  std::vector<double> P(static_cast<size_t>(n), 0.);
  for (long long S = 0; S < nset; ++S) {
    if (lw[static_cast<size_t>(S)] > -JV_INF) {
      const double w = std::exp(lw[static_cast<size_t>(S)] - mx);
      Z += w;
      for (int k = 0; k < n; ++k)
        if (S & (1LL << k)) P[static_cast<size_t>(k)] += w;
    }
  }
  for (int k = 0; k < n; ++k) P[static_cast<size_t>(k)] /= Z;
  return P;
}


// ---------------------------------------------------------------- vertex score
// features of a vertex candidate (the S0 tracks trk at x with covariance C)
struct CandFeats {
  int n = 0;
  double mass = 0., dis = 0., cos = 0., sigl = 0., fsig2 = 0., chi2 = 0.;
  int n2v = 0, n0v = 0;
  double s0 = 0., s1 = 0., smax = 0., drmax = 0., p = 0., r = 0., pmin = 0.;
  int nsig2 = 0, nneg = 0;
  double sumpos = 0., ptrel = 0., jetcos = 0.;
  int njetsame = 0;
};

// engineered inputs of the logistic score, in the order of JV_SCORE_CAND.names
using ScoreInputs = std::array<double, 30>;

inline ScoreInputs scoreInputs(const CandFeats& f) {
  const double n = static_cast<double>(f.n);
  const double mc = std::sqrt(detail::powSquare(f.mass) + detail::powSquare(f.ptrel)) + f.ptrel;
  const double lf = std::log(f.fsig2 + 1);
  const double lp = std::log(f.p + .1);
  const double ld = std::log(f.dis + 1e-3);
  const double s0c = npClip(f.s0, -5, 3);
  const double is2 = f.n == 2 ? 1. : 0.;
  ScoreInputs x;
  x[0] = static_cast<double>(std::min(f.n, 6));
  x[1] = lp;
  x[2] = std::log(1 - npMinimum(f.cos, 1) + 1e-4);
  x[3] = f.drmax;
  x[4] = npClip(f.s0, -5, 10);
  x[5] = npClip(f.s1, -5, 15);
  x[6] = npClip(f.smax, -3, 20);
  x[7] = static_cast<double>(std::min(f.n2v, 4));
  x[8] = static_cast<double>(std::min(f.n0v, 3));
  x[9] = lf;
  x[10] = std::log(f.sigl + 1e-3);
  x[11] = f.chi2 / static_cast<double>(2 * f.n - 3);
  x[12] = std::log(f.mass + .05);
  x[13] = ld;
  x[14] = std::log(f.pmin + .05);
  x[15] = static_cast<double>(f.nsig2) / n;
  x[16] = static_cast<double>(f.nneg) / n;
  x[17] = f.sumpos / n;
  x[18] = std::log(mc + .05);
  x[19] = std::log(f.ptrel + .02);
  x[20] = std::log(1 - npMinimum(f.jetcos, 1) + 1e-4);
  x[21] = static_cast<double>(f.njetsame) / n;
  x[22] = detail::powSquare(s0c);
  x[23] = detail::powSquare(f.drmax);
  x[24] = detail::powSquare(lf);
  x[25] = detail::powSquare(lp);
  x[26] = ld * lf;
  x[27] = is2;
  x[28] = is2 * lf;
  x[29] = is2 * npClip(f.s1, -5, 15);
  return x;
}

// logistic score: ((x - mu) / sd) . coef + b0
inline double logisticScore(const ScoreCoefs& S, const ScoreInputs& x) {
  double z[30];
  for (size_t i = 0; i < 30; ++i) z[i] = (x[i] - S.mu[i]) / S.sd[i];
  return detail::dotBlasN(z, S.coef.data(), 30) + S.b0;
}

// ---------------------------------------------------------------- candidate features
// eta-phi distance of two tracks from their fit states (index 1 phi, 4 tan lambda)
inline double etaPhiDr(const Vec5& a, const Vec5& b) {
  const double e1 = -std::log(std::tan(0.5 * std::atan2(1., a[4])));
  const double e2 = -std::log(std::tan(0.5 * std::atan2(1., b[4])));
  const double df = wrapPi(a[1] - b[1]);
  return std::hypot(e1 - e2, df);
}

// largest eta-phi distance over the track pairs (Python max over pairs in index order)
inline double maxEtaPhiDr(const std::vector<Vec5>& states) {
  double m = JV_NAN;
  bool first = true;
  for (size_t i = 0; i < states.size(); ++i)
    for (size_t j = i + 1; j < states.size(); ++j) {
      const double d = etaPhiDr(states[i], states[j]);
      if (first) { m = d; first = false; }
      else m = std::max(m, d);
    }
  return m;
}

// momentum sum and invariant mass of tracks with the given momenta (masses: default pion)
struct MassP {
  double m = 0.;
  Vec3 P = Vec3::Zero();
};
inline MassP massOfMomenta(const std::vector<Vec3>& p, const double* masses = nullptr) {
  double E = 0.;
  Vec3 P = Vec3::Zero();
  for (size_t k = 0; k < p.size(); ++k) {
    const double m = masses ? masses[k] : JV_M_PI;
    E += std::sqrt(dotSeq(p[k], p[k]) + m * m);
    P += p[k];
  }
  MassP r;
  r.m = std::sqrt(std::max(E * E - dotSeq(P, P), 0.));
  r.P = P;
  return r;
}

// inputs of candFeats: mass / P from massOfMomenta of the tracks at x; per-track arrays in the order of the track list
struct CandFeatsIn {
  Vec3 x, pvx;
  Mat3 C, pvC;
  bool covFortran = true;         // memory order of C in the reference (see quadForm)
  double chi2 = 0.;
  double mass = 0.;
  Vec3 P = Vec3::Zero();
  std::vector<double> s3, p;
  std::vector<int> nvdet, jet;
  double drmax = 0.;
  std::vector<Vec3> jetAxes;      // one row per jet axis
};

inline CandFeats candFeats(const CandFeatsIn& in) {
  CandFeats f;
  const Vec3 d = in.x - in.pvx;
  const double dis = std::max(detail::norm3(d), 1e-9);
  const double pm = std::max(detail::norm3(in.P), 1e-9);
  const Vec3 u = in.P / pm;
  std::vector<double> ss = in.s3;
  std::sort(ss.begin(), ss.end(), npLess);
  const int n = static_cast<int>(ss.size());
  const Vec3 dn = d / dis;
  const Vec3 cr(in.P[1] * dn[2] - in.P[2] * dn[1], in.P[2] * dn[0] - in.P[0] * dn[2], in.P[0] * dn[1] - in.P[1] * dn[0]);
  f.n = n;
  f.mass = in.mass;
  f.dis = dis;
  f.cos = std::min(dotSeq(d, in.P) / (dis * pm), 1.);
  f.sigl = std::sqrt(std::max(detail::quadForm(u, in.C, in.covFortran), 0.));
  f.fsig2 = dsig2(in.x, in.C, in.pvx, in.pvC);
  f.chi2 = in.chi2;
  f.s0 = ss[0];
  f.s1 = ss[1];
  f.smax = ss[ss.size() - 1];
  f.drmax = in.drmax;
  f.p = pm;
  f.r = std::hypot(in.x[0], in.x[1]);
  double pmin = in.p[0];
  std::vector<double> sp(ss.size());
  for (size_t i = 0; i < ss.size(); ++i) {
    f.n2v += in.nvdet[i] >= 2;
    f.n0v += in.nvdet[i] == 0;
    f.nsig2 += ss[i] > 2;
    f.nneg += ss[i] < -1;
    pmin = npMinimum(pmin, in.p[i]);
    sp[i] = npMinimum(npMaximum(ss[i], 0.), 10.);
  }
  f.pmin = pmin;
  f.sumpos = npSumRange(sp.data(), n);
  f.ptrel = std::sqrt(dotSeq(cr, cr));
  double jc = JV_NAN;
  for (size_t k = 0; k < in.jetAxes.size(); ++k) {
    const Vec3& a = in.jetAxes[k];
    const double v = (a[0] * dn[0] + a[2] * dn[2]) + a[1] * dn[1];
    jc = k == 0 ? v : npMaximum(jc, v);
  }
  f.jetcos = jc;
  int j0 = 0, j1 = 0;
  for (int j : in.jet) { j0 += j == 0; j1 += j == 1; }
  f.njetsame = std::max(j0, j1);
  return f;
}

// ---------------------------------------------------------------- acceptance window and pre-selection
struct AcceptCfg {
  double disLo = JV_ACC_DIS_LO, disHi = JV_ACC_DIS_HI, cosMin = JV_ACC_COS, siglMax = JV_ACC_SIGL, fsig2 = JV_ACC_FSIG2;
  int maxTrk = JV_ACC_MAX_TRK;
};

// window of a vertex at x (covariance C) with ntrk tracks and momentum sum P: distance from the reference PV,
// number of tracks, flight direction, length error and (two tracks) flight significance
inline bool acceptWindow(int ntrk, const Vec3& x, const Mat3& C, const Vec3& pvx, const Mat3& pvC, const Vec3& P,
                         const AcceptCfg& c = AcceptCfg(), bool covFortran = true) {
  const Vec3 d = x - pvx;
  const double dis = detail::norm3(d);
  if (!(c.disLo <= dis && dis <= c.disHi) || ntrk > c.maxTrk || ntrk < 2) return false;
  const double pm = detail::norm3(P);
  if (pm <= 0 || dotSeq(d, P) / (dis * pm) < c.cosMin) return false;
  const Vec3 u = P / pm;
  const double var = detail::quadForm(u, C, covFortran);
  if (!(var > 0) || std::sqrt(var) > c.siglMax) return false;
  if (ntrk == 2 && dsig2(x, C, pvx, pvC) <= detail::powSquare(c.fsig2)) return false;
  return true;
}

// stage-1 PV pre-selection window of a raw track state (d0, phi, omega, z0, tan lambda) w.r.t. the beam spot
inline bool preWindow(const Vec5& par, const Vec3& bs, double d0max = JV_PV_WINDOW_D0, double z0max = JV_PV_WINDOW_Z0) {
  const double s = bs[0] * std::cos(par[1]) + bs[1] * std::sin(par[1]);
  const double d0p = par[0] - bs[0] * std::sin(par[1]) + bs[1] * std::cos(par[1]) - 0.5 * par[2] * s * s;
  const double z0p = par[3] - bs[2] + par[4] * s;
  return std::fabs(d0p) <= d0max && std::fabs(z0p) <= z0max;
}

// ---------------------------------------------------------------- V0-like pair test of the stage-1 veto
inline bool v0Like(int qi, int qj, const Vec3& pi, const Vec3& pj, double ks = JV_KS_WIN, double lam = JV_LAM_WIN,
                   double ee = JV_EE_WIN) {
  if (qi * qj >= 0) return false;
  const std::vector<Vec3> p = {pi, pj};
  if (std::fabs(massOfMomenta(p).m - JV_M_KS) < ks) return true;
  const double m1[2] = {JV_M_P, JV_M_PI};
  const double m2[2] = {JV_M_PI, JV_M_P};
  if (std::fabs(massOfMomenta(p, m1).m - JV_M_LAM) < lam) return true;
  if (std::fabs(massOfMomenta(p, m2).m - JV_M_LAM) < lam) return true;
  const double me[2] = {JV_M_E, JV_M_E};
  if (massOfMomenta(p, me).m < ee) return true;
  return false;
}

// ---------------------------------------------------------------- heavy-flavour features: track selection
// indices of the usable tracks (ok) displaced from the reference PV (s3 > sig) except the two daughters
inline std::vector<int> hfOthers(const std::vector<char>& ok, const std::vector<double>& s3, double sig, int a, int b) {
  std::vector<int> o;
  for (size_t i = 0; i < ok.size(); ++i)
    if (ok[i] && s3[i] > sig && static_cast<int>(i) != a && static_cast<int>(i) != b) o.push_back(static_cast<int>(i));
  return o;
}

}  // namespace AlephJointVtx
}  // namespace FCCAnalyses

#endif
