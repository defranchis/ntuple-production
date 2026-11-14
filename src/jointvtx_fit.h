#ifndef JOINTVTX_FIT_H
#define JOINTVTX_FIT_H

/*
  Vertex fits of the joint vertex fitter (AlephJointVtx) and the jet-projected
  reference primary vertex.

  fitVertexNb: weighted Gauss-Newton fit of one vertex to tracks, relinearised
  at every step (lin1), with an optional Gaussian constraint (information
  Binv, centre c); at most niter undamped steps, stop when max |dx| < tol. The
  per-track start arcs Ls (one per track of the list) persist between calls.
  The convergence flag is returned; the result is that of the last step either
  way.
  teardown: fit, then remove the worst track while its chi2 >= chi2max and more
  than nmin tracks remain (the arcs persist between the refits).
  pairFit: two-track fit started at the perigee midpoint and at each circle
  crossing; lowest chi2 inside the hit-radius bound, else the perigee start.
  jetProjPV: jet-projected primary vertex on a straight-line track model at the
  origin-referenced perigee (hard form): per track the offsets along the theta
  direction, the transverse direction, E (perpendicular to track and jet) and
  U (in the track-jet plane); 6 passes with 4-sigma windows and component-wise
  outlier removal at delta-chi2 >= 9. A track whose in-plane offset is below
  `window` enters with its full information, the others across the jet only.
  Units: cm, GeV.
  Operation order follows the reference implementation (products written as
  loops in the measured BLAS / einsum order). Failures are status values.
*/

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include <Eigen/Dense>

#include "jointvtx_common.h"
#include "jointvtx_track.h"

namespace FCCAnalyses {
namespace AlephJointVtx {

constexpr double JV_JP_MINERR2 = 6.25e-6;  // cm^2, minimum component error (25 um)^2
constexpr double JV_JP_NSIG = 4.0;         // window in standard deviations
constexpr double JV_JP_DCHI2 = 9.0;        // component removal threshold
constexpr int JV_JP_NPASS = 6;             // passes over the track selection
constexpr double JV_JP_PCAP = 45.0;        // GeV, momentum cap of the jet angular error
constexpr int JV_FIT_NITER = 40;           // default iterations of fitVertex / teardown
constexpr double JV_FIT_TOL = 1e-10;       // cm, default step tolerance
constexpr int JV_CAND_FIT_NITER = 15;      // steps of the pair, V0-pair, candidate and secondary-vertex fits and of
                                           // the PV refit of the joint fit
constexpr double JV_CAND_FIT_TOL = 1e-6;   // cm, step tolerance of those fits

using Vec15 = Vec<15>;

struct BeamSpot {
  Vec3 c = Vec3::Zero();    // centre
  Vec3 sig = Vec3::Zero();  // widths x, y, z
};

struct VtxFit {
  Vec3 x = Vec3::Constant(JV_NAN);
  Mat3 C = Mat3::Constant(JV_NAN);
  std::vector<double> chi2;  // per fitted track
  bool conv = false;
  LinStatus status = kLinOk;
  bool ok() const { return status == kLinOk; }
};

struct TeardownFit {
  Vec3 x = Vec3::Constant(JV_NAN);
  Mat3 C = Mat3::Constant(JV_NAN);
  std::vector<int> kept;
  std::vector<double> chi2;  // per kept track
  LinStatus status = kLinOk;
  bool ok() const { return status == kLinOk; }
};

struct PairFit {
  bool found = false;  // false: no start gave a finite solution
  Vec3 x = Vec3::Constant(JV_NAN);
  Mat3 C = Mat3::Constant(JV_NAN);
  double chi2 = JV_NAN;
  LinStatus status = kLinOk;
};

// line-model inputs per track: d0 = -D0, z0, phi0, ct = tanLambda, cov(d0, d0), cov(z0, z0), cov(d0, z0)
struct JpTrk {
  double d0 = 0., z0 = 0., phi = 0., ct = 0., c00 = 0., c11 = 0., c01 = 0.;
};

// per track: directions H[4][3] (theta, transverse, E, U), error^2 E2[5] (theta, transverse, E, U, E without the
// jet term), offsets DT0[4], track-jet angle nu, jet index kj
struct JpPrep {
  std::vector<std::array<double, 12>> H;
  std::vector<std::array<double, 5>> E2;
  std::vector<std::array<double, 4>> DT0;
  std::vector<double> nu;
  std::vector<int> kj;
  int size() const { return static_cast<int>(nu.size()); }
};

enum JpStatus : int {
  kJpOk = 0,
  kJpSingular = 1,   // exactly singular information matrix (the reference raises)
  kJpNoJets = 2,     // no jet axis (the reference raises)
  kJpBadJet = 3,     // forced jet index outside the axes (the reference raises)
  kJpNonFinite = 4   // non-finite position
};

struct JpPV {
  Vec3 x = Vec3::Constant(JV_NAN);
  Mat3 C = Mat3::Constant(JV_NAN);  // inverse information of the last solve
  int nused = 0;                    // tracks with at least one component in the last solve
  JpStatus status = kJpOk;
  bool ok() const { return status == kJpOk; }
};

struct InvMass {
  double m = 0.;
  Vec3 P = Vec3::Zero();
};

// constraint information of a beam spot: diag(1 / sig^2)
inline Mat3 beamInfo(const Vec3& sig) {
  Mat3 B = Mat3::Zero();
  for (int i = 0; i < 3; ++i) B(i, i) = 1. / (sig[i] * sig[i]);
  return B;
}

// weighted fit of the tracks idx (weights w, same length) from x0 with the constraint information Binv (zeros =
// none) around c; Ls: start arcs of all tracks of t, updated
inline VtxFit fitVertexNb(const std::vector<Trk>& t, const std::vector<int>& idx, const std::vector<double>& w,
                          const Vec3& x0, const Mat3& Binv, const Vec3& c, std::vector<double>& Ls, int niter,
                          double tol) {
  VtxFit r;
  const size_t n = idx.size();
  Vec3 x = x0;
  std::vector<Vec3> X;
  std::vector<Mat3> D;
  std::vector<double> chi2;
  auto fail = [&r](LinStatus s) {
    r.status = s;
    r.conv = false;
    setNaN(r.x);
    setNaN(r.C);
    for (double& v : r.chi2) v = JV_NAN;
    return r;
  };
  for (int it = 0; it < niter; ++it) {
    LinStatus s = linAll(t, idx, x, Ls, X, D, chi2);
    if (s != kLinOk) return fail(s);
    Mat3 H = Binv;
    Vec3 g = matVecC<3, 3>(Binv, Vec3(x - c));
    for (size_t k = 0; k < n; ++k) {
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) H(a, b) += w[k] * D[k](a, b);
      const Vec3 v = matVecC<3, 3>(D[k], Vec3(x - X[k]));
      for (int a = 0; a < 3; ++a) g[a] += w[k] * v[a];
    }
    Vec3 dx;
    s = npSolve<3>(H, g, dx);
    if (s != kLinOk) return fail(s);
    dx = -dx;
    if (!allFinite(dx)) break;
    x = x + dx;
    const double m = std::max(std::max(std::fabs(dx[0]), std::fabs(dx[1])), std::fabs(dx[2]));
    if (m < tol) {
      r.conv = true;
      break;
    }
  }
  const LinStatus s = linAll(t, idx, x, Ls, X, D, chi2);
  if (s != kLinOk) return fail(s);
  Mat3 H = Binv;
  for (size_t k = 0; k < n; ++k)
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b) H(a, b) += w[k] * D[k](a, b);
  r.chi2 = chi2;
  const LinStatus si = nbInv<3>(H, r.C);
  if (si != kLinOk) return fail(si);
  r.x = x;
  return r;
}

// fitVertexNb with unit weights (w = nullptr), the beam-spot constraint bs (nullptr = none) and the start arcs Ls
// (nullptr = all zero)
inline VtxFit fitVertex(const std::vector<Trk>& t, const std::vector<int>& idx, const Vec3& x0,
                        const BeamSpot* bs = nullptr, std::vector<double>* Ls = nullptr, int niter = JV_FIT_NITER,
                        double tol = JV_FIT_TOL, const std::vector<double>* w = nullptr) {
  const std::vector<double> w1 = w ? *w : std::vector<double>(idx.size(), 1.);
  const Mat3 Binv = bs ? beamInfo(bs->sig) : Mat3::Zero();
  const Vec3 c = bs ? bs->c : Vec3::Zero();
  std::vector<double> L0;
  if (!Ls) L0.assign(t.size(), 0.);
  return fitVertexNb(t, idx, w1, x0, Binv, c, Ls ? *Ls : L0, niter, tol);
}

// fit, then drop the worst track (np.argmax: first NaN, else first maximum) while its chi2 >= chi2max and more than
// nmin tracks remain
inline TeardownFit teardown(const std::vector<Trk>& t, const std::vector<int>& idx, const Vec3& x0,
                            const BeamSpot* bs, double chi2max, int nmin = 2) {
  TeardownFit r;
  if (idx.empty()) {
    if (!bs) {
      r.x = x0;
      r.C = Mat3::Identity();
    } else {
      r.x = bs->c;
      r.C = Mat3::Zero();
      for (int i = 0; i < 3; ++i) r.C(i, i) = bs->sig[i] * bs->sig[i];
    }
    return r;
  }
  std::vector<int> cur = idx;
  std::vector<double> Ls(t.size(), 0.);
  Vec3 x = x0;
  while (true) {
    const VtxFit f = fitVertex(t, cur, x, bs, &Ls);
    if (!f.ok()) {
      r.status = f.status;
      return r;
    }
    x = f.x;
    const int j = argmaxNaN(f.chi2);
    if (!(f.chi2[j] >= chi2max) || static_cast<int>(cur.size()) - 1 < nmin) {
      r.x = f.x;
      r.C = f.C;
      r.kept = cur;
      r.chi2 = f.chi2;
      return r;
    }
    cur.erase(cur.begin() + j);
  }
}

// two-track vertex of tracks i, j: fits (JV_CAND_FIT_NITER steps) from x0 (default: midpoint of the perigee
// points) and from each circle crossing; the lowest-chi2 solution with transverse radius below min(rmaxI, rmaxJ),
// else the first start's solution
inline PairFit pairFit(const std::vector<Trk>& t, int i, int j, double rmaxI = JV_INF, double rmaxJ = JV_INF,
                       const Vec3* x0 = nullptr, bool useCrossings = true) {
  PairFit res;
  const double rmax = std::min(rmaxI, rmaxJ);
  struct Start {
    Vec3 x;
    double Li, Lj;
  };
  std::vector<Start> starts;
  starts.push_back({x0 ? *x0 : Vec3(0.5 * (pointAt(t[i].p, 0.) + pointAt(t[j].p, 0.))), 0., 0.});
  if (useCrossings) {
    const Crossings cr = circleCrossings(t[i].p, t[j].p);
    for (int k = 0; k < cr.n; ++k) {
      const double Li = arcToXY(t[i].p, cr.x[k], cr.y[k]), Lj = arcToXY(t[j].p, cr.x[k], cr.y[k]);
      const double zi = t[i].p[3] + t[i].p[4] * Li, zj = t[j].p[3] + t[j].p[4] * Lj;
      starts.push_back({Vec3(cr.x[k], cr.y[k], 0.5 * (zi + zj)), Li, Lj});
    }
  }
  const std::vector<int> idx = {i, j};
  bool haveBest = false, haveFirst = false;
  PairFit first;
  for (size_t k = 0; k < starts.size(); ++k) {
    std::vector<double> Ls(t.size(), 0.);
    Ls[i] = starts[k].Li;
    Ls[j] = starts[k].Lj;
    const VtxFit f = fitVertex(t, idx, starts[k].x, nullptr, &Ls, JV_CAND_FIT_NITER, JV_CAND_FIT_TOL);
    if (!f.ok()) {
      PairFit bad;
      bad.status = f.status;
      return bad;
    }
    const double tot = npSum(f.chi2);
    if (!std::isfinite(tot) || !allFinite(f.x)) continue;
    PairFit sol;
    sol.found = true;
    sol.x = f.x;
    sol.C = f.C;
    sol.chi2 = tot;
    if (k == 0) {
      first = sol;
      haveFirst = true;
    }
    if (std::hypot(f.x[0], f.x[1]) < rmax && (!haveBest || tot < res.chi2)) {
      res = sol;
      haveBest = true;
    }
  }
  if (haveBest) return res;
  if (haveFirst) return first;
  return PairFit();
}

// invariant mass [GeV] and momentum sum of the tracks idx with the momenta at the arc of their xy point nearest to
// x (momentum constant ptHalf), scaled by pscale when pscale != 1; masses: one per track, nullptr = pions
inline InvMass invMass(const std::vector<Trk>& t, const std::vector<int>& idx, const Vec3& x,
                       const double* masses = nullptr, double ptHalf = JV_PT_HALF, double pscale = 1.0) {
  double E = 0.;
  Vec3 P = Vec3::Zero();
  for (size_t k = 0; k < idx.size(); ++k) {
    const double m = masses ? masses[k] : JV_M_PI;
    const Vec5& p5 = t[idx[k]].p;
    Vec3 p = momentumAt(p5, arcNear(p5, x), ptHalf);
    if (pscale != 1.0)
      for (int a = 0; a < 3; ++a) p[a] = p[a] * pscale;
    E += std::sqrt(dotSeq(p, p) + m * m);
    for (int a = 0; a < 3; ++a) P[a] += p[a];
  }
  InvMass r;
  r.m = std::sqrt(std::max(E * E - dotSeq(P, P), 0.));
  r.P = P;
  return r;
}

// jet angular error^2 [rad^2] of a jet of momentum p [GeV] (p capped at 45 GeV)
inline double jetSig2Aleph(double p) {
  p = std::min(p, JV_JP_PCAP);
  return 4.0 * detail::powSquare(0.2862 - 0.01205 * p + 0.000135 * p * p);
}

// jet axes from polar and azimuthal angles
inline std::vector<Vec3> jpJetAxes(const std::vector<double>& theta, const std::vector<double>& phi) {
  std::vector<Vec3> ax;
  for (size_t k = 0; k < theta.size() && k < phi.size(); ++k)
    ax.push_back(Vec3(std::sin(theta[k]) * std::cos(phi[k]), std::sin(theta[k]) * std::sin(phi[k]), std::cos(theta[k])));
  return ax;
}

// lifetime-sign axes: the jet axes, or (thrust, -thrust) when any axis is not finite
inline std::vector<Vec3> jetAxes(const std::vector<double>& theta, const std::vector<double>& phi, const Vec3& thrust) {
  std::vector<Vec3> ax = jpJetAxes(theta, phi);
  bool fin = true;
  for (const Vec3& a : ax) fin = fin && allFinite(a);
  if (!fin) ax = {thrust, Vec3(-thrust)};
  return ax;
}

// raw track states (D0, phi, omega, Z0, tanLambda) and raw covariances (15 packed elements) -> line-model inputs
inline std::vector<JpTrk> jpLineParams(const std::vector<Vec5>& parRaw, const std::vector<Vec15>& cov15) {
  std::vector<JpTrk> q(parRaw.size());
  for (size_t k = 0; k < parRaw.size(); ++k) {
    q[k].d0 = -parRaw[k][0];
    q[k].z0 = parRaw[k][3];
    q[k].phi = parRaw[k][1];
    q[k].ct = parRaw[k][4];
    q[k].c00 = cov15[k][0];
    q[k].c11 = cov15[k][9];
    q[k].c01 = -cov15[k][6];
  }
  return q;
}

// per-track directions, errors and offsets; each track goes with the jet of largest t.j unless kjForce (one per
// track, -1 = not forced; empty = none) names a jet the track does not point away from
inline JpStatus jpPrepQfn(const std::vector<JpTrk>& trk, const std::vector<Vec3>& axes, const std::vector<int>& kjForce,
                          JpPrep& Q, double minerr2 = JV_JP_MINERR2, double t3sign = 1.0) {
  const int N = static_cast<int>(trk.size());
  const int nj = static_cast<int>(axes.size());
  Q.H.assign(N, {});
  Q.E2.assign(N, {});
  Q.DT0.assign(N, {});
  Q.nu.assign(N, 0.);
  Q.kj.assign(N, 0);
  if (nj == 0) return kJpNoJets;
  std::vector<std::array<double, 3>> T1(N), T0(N);
  std::vector<double> T_1(N), T_2(N), T_3(N);
  for (int i = 0; i < N; ++i) {
    const JpTrk& q = trk[i];
    const double cosl = 1.0 / std::sqrt(1.0 + q.ct * q.ct);
    const double sinl = q.ct * cosl;
    T1[i] = {cosl * std::cos(q.phi), cosl * std::sin(q.phi), sinl};
    T0[i] = {-q.d0 * std::sin(q.phi), q.d0 * std::cos(q.phi), q.z0};
    T_1[i] = q.c00;
    T_2[i] = q.c11 * (cosl * cosl);
    T_3[i] = t3sign * q.c01 * cosl;
    std::array<double, 12>& H = Q.H[i];
    const double norm = std::sqrt(T1[i][0] * T1[i][0] + T1[i][1] * T1[i][1]);
    H[3] = -T1[i][1] / norm;
    H[4] = T1[i][0] / norm;
    H[5] = 0.;
    H[0] = -T1[i][0] * T1[i][2] / norm;
    H[1] = -T1[i][1] * T1[i][2] / norm;
    H[2] = norm;
    Q.E2[i][1] = npMaximum(T_1[i], minerr2);
    Q.E2[i][0] = npMaximum(T_2[i], minerr2);
  }
  // T1 @ axes.T: gemm left to right; gemv (one row or one axis) even / odd; one row and one axis: dot
  const bool seqOrder = (N >= 2) == (nj >= 2);
  std::vector<double> dots(static_cast<size_t>(N) * nj);
  for (int i = 0; i < N; ++i)
    for (int k = 0; k < nj; ++k) {
      const double p0 = T1[i][0] * axes[k][0], p1 = T1[i][1] * axes[k][1], p2 = T1[i][2] * axes[k][2];
      dots[static_cast<size_t>(i) * nj + k] = seqOrder ? (p0 + p1) + p2 : (p0 + p2) + p1;
    }
  for (int i = 0; i < N; ++i) {
    const double* row = &dots[static_cast<size_t>(i) * nj];
    int kj = argmaxNaN(std::vector<double>(row, row + nj));
    if (!kjForce.empty()) {
      const int kf = kjForce[i];
      const int kc = std::max(kf, 0);
      if (kc >= nj) return kJpBadJet;
      if (kf >= 0 && row[kc] > 0) kj = kf;
    }
    Q.kj[i] = kj;
    Q.nu[i] = std::acos(npMinimum(row[kj], 1.0));
  }
  for (int i = 0; i < N; ++i) {
    std::array<double, 12>& H = Q.H[i];
    std::array<double, 5>& E2 = Q.E2[i];
    const Vec3& J1 = axes[Q.kj[i]];
    if (Q.nu[i] > 0) {
      const double s = std::sin(Q.nu[i]);
      const std::array<double, 3>& tv = T1[i];
      const double E[3] = {(tv[1] * J1[2] - tv[2] * J1[1]) / s, (tv[2] * J1[0] - tv[0] * J1[2]) / s,
                           (tv[0] * J1[1] - tv[1] * J1[0]) / s};
      const double U[3] = {-tv[1] * E[2] + tv[2] * E[1], -tv[2] * E[0] + tv[0] * E[2], -tv[0] * E[1] + tv[1] * E[0]};
      for (int a = 0; a < 3; ++a) {
        H[6 + a] = E[a];
        H[9 + a] = U[a];
      }
      const double cose = (E[0] * H[3] + E[1] * H[4]) + E[2] * H[5];
      const double sine = (E[0] * H[0] + E[1] * H[1]) + E[2] * H[2];
      const double c2 = detail::powSquare(cose), s2 = detail::powSquare(sine);
      const double cr = 2 * T_3[i] * cose * sine;
      E2[2] = std::max(T_1[i] * c2 + T_2[i] * s2 + cr, minerr2);
      E2[3] = std::max(T_2[i] * c2 + T_1[i] * s2 - cr, minerr2);
    } else {
      for (int a = 0; a < 3; ++a) {
        H[6 + a] = H[a];
        H[9 + a] = H[3 + a];
      }
      E2[2] = E2[0];
      E2[3] = E2[1];
    }
    E2[4] = E2[2];
    for (int c = 0; c < 4; ++c)
      Q.DT0[i][c] = (H[3 * c] * T0[i][0] + H[3 * c + 2] * T0[i][2]) + H[3 * c + 1] * T0[i][1];
  }
  return kJpOk;
}

// jet-projected PV from the prepared tracks: per pass the selected components (theta + transverse for a track
// with full information, E only for an across-jet track) are fitted with the constraint bsW around bsc, removing
// the component of largest delta-chi2 while >= 9; stops when the step is within the error ellipse of the last
// solve, or after 6 passes
inline JpPV jpQfndip(const JpPrep& Q, const std::vector<double>& axesSig2, const Mat3& bsW, const Vec3& bsc,
                     double window) {
  JpPV res;
  const int N = Q.size();
  std::vector<std::array<double, 5>> E2 = Q.E2;
  std::vector<int> flag(N, 4);
  const Vec3 bsb = matVecC<3, 3>(bsW, bsc);
  Vec3 old = Vec3::Zero();
  bool haveOld = false;
  Vec3 x = Vec3::Zero();
  Mat3 Ai = Mat3::Zero();
  int nused = 0;
  std::vector<int> ii, tt;
  std::vector<double> pa, pb;
  for (int meta = 0; meta < JV_JP_NPASS; ++meta) {
    ii.clear();
    tt.clear();
    for (int i = 0; i < N; ++i)
      for (int c = 0; c < 3; ++c)
        if ((flag[i] >> c) & 1) {
          ii.push_back(i);
          tt.push_back(c);
        }
    while (true) {
      const int m = static_cast<int>(ii.size());
      std::vector<std::array<double, 3>> hs(m);
      std::vector<double> e2(m), dt(m);
      for (int k = 0; k < m; ++k) {
        for (int a = 0; a < 3; ++a) hs[k][a] = Q.H[ii[k]][3 * tt[k] + a];
        e2[k] = E2[ii[k]][tt[k]];
        dt[k] = Q.DT0[ii[k]][tt[k]];
      }
      Mat3 A;
      Vec3 b;
      pa.resize(m);
      pb.resize(m);
      for (int a = 0; a < 3; ++a) {
        for (int c = 0; c < 3; ++c) {
          for (int k = 0; k < m; ++k) pa[k] = (hs[k][a] / e2[k]) * hs[k][c];
          A(a, c) = bsW(a, c) + detail::gemm3mBlas(pa.data(), m, a, c);
        }
        for (int k = 0; k < m; ++k) pb[k] = hs[k][a] * (dt[k] / e2[k]);
        b[a] = bsb[a] + detail::gemv3mBlas(pb.data(), m);
      }
      const LinStatus s = npInv<3>(A, Ai);
      if (s == kLinSingular) {
        res.status = kJpSingular;
        return res;
      }
      x = matVecC<3, 3>(Ai, b);
      if (m == 0) break;
      int kb = 0;
      std::vector<double> d(m);
      const double xa[3] = {x[0], x[1], x[2]};
      for (int k = 0; k < m; ++k) {
        const double r = dt[k] - detail::rowDot3(hs[k].data(), xa, m == 1);
        double q = 0.;
        for (int a = 0; a < 3; ++a)
          for (int c = 0; c < 3; ++c) q += hs[k][a] * Ai(a, c) * hs[k][c];
        const double den = e2[k] - q;
        d[k] = den > 0 ? r * r / den : 0.0;
      }
      kb = argmaxNaN(d);
      if (d[kb] >= JV_JP_DCHI2) {
        ii.erase(ii.begin() + kb);
        tt.erase(tt.begin() + kb);
      } else {
        break;
      }
    }
    std::vector<int> u = ii;
    std::sort(u.begin(), u.end());
    nused = static_cast<int>(std::unique(u.begin(), u.end()) - u.begin());
    if (haveOld) {
      const Vec3 dx = x - old;
      const double d2 = dotSeq(dx, dx);
      if (d2 == 0 || !(d2 > detail::quadForm(dx, Ai, false) / d2)) break;
    }
    if (meta == JV_JP_NPASS - 1) break;
    old = x;
    haveOld = true;
    for (int i = 0; i < N; ++i) {
      const std::array<double, 12>& H = Q.H[i];
      double L[4];
      for (int c = 0; c < 4; ++c) L[c] = Q.DT0[i][c] - ((H[3 * c] * x[0] + H[3 * c + 2] * x[2]) + H[3 * c + 1] * x[1]);
      double S[4];
      for (int c = 0; c < 4; ++c) S[c] = std::sqrt(E2[i][c]);
      const bool full = (L[3] < window) && (std::fabs(L[0]) < JV_JP_NSIG * S[0]) && (std::fabs(L[1]) < JV_JP_NSIG * S[1]);
      const bool eonly = !full && (std::fabs(L[2]) < JV_JP_NSIG * std::sqrt(E2[i][2]));
      flag[i] = full ? 3 : (eonly ? 4 : 0);
      const double s2 = axesSig2[Q.kj[i]];
      const double nus = Q.nu[i] > 0 ? Q.nu[i] : 1.0;
      const double l2 = L[3] * L[3];
      const double ln = L[3] / nus;
      const double jerr2 = Q.nu[i] > 0 ? npMinimum(l2, s2 * (ln * ln)) : l2;
      if (eonly) E2[i][2] = E2[i][4] + jerr2;
    }
  }
  res.x = x;
  res.C = Ai;
  res.nused = nused;
  if (!allFinite(x)) res.status = kJpNonFinite;
  return res;
}

// jet-projected PV (hard form, window 0) of the tracks given as raw states and covariances; jetOf: jet of each
// track (-1 = none); the beam spot enters as a Gaussian constraint
inline JpPV jetProjPV(const std::vector<Vec5>& parRaw, const std::vector<Vec15>& cov15, const std::vector<int>& jetOf,
                      const std::vector<double>& jetTheta, const std::vector<double>& jetPhi,
                      const std::vector<double>& jetP, const BeamSpot& bs) {
  const std::vector<JpTrk> q = jpLineParams(parRaw, cov15);
  const std::vector<Vec3> axes = jpJetAxes(jetTheta, jetPhi);
  std::vector<double> sg;
  for (double p : jetP) sg.push_back(jetSig2Aleph(p));
  JpPrep Q;
  const JpStatus s = jpPrepQfn(q, axes, jetOf, Q);
  if (s != kJpOk) {
    JpPV r;
    r.status = s;
    return r;
  }
  return jpQfndip(Q, sg, beamInfo(bs.sig), bs.c, 0.0);
}

}  // namespace AlephJointVtx
}  // namespace FCCAnalyses

#endif
