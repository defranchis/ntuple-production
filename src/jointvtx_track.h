#ifndef JOINTVTX_TRACK_H
#define JOINTVTX_TRACK_H

/*
  Track kernels of the joint vertex fitter (AlephJointVtx): helix model, its
  linearisation at a point, track-to-point and track-to-vertex compatibility,
  impact parameters and the two-dimensional circle geometry.

  Trk: helix parameters p = (D, phi0, C, z0, tanLambda) in the vertex-fitter
  convention (D = -d0, C = omega / 2) with their covariance C, plus the VDET and
  ITC hit counts. Event track lists are std::vector<Trk>; kernels over a subset
  take the list and the indices.
  helix1: point X at transverse arc L, direction a and Jacobian A = dX / dp.
  lin1: the track replaced at a point x by its closest approach X in the track
  metric and the rank-2 weight D = W - W a a^T W / (a^T W a), W = (A Cov A^T)^-1.
  Momenta: pT = ptHalf / |C| (ptHalf defaults to the stage1 constant).
  Operation order follows the reference implementation: products are written
  as loops in the measured BLAS order, no Eigen expression products.
  A kernel whose reference raises returns a failure status and NaN outputs.
*/

#include <cmath>
#include <vector>

#include <Eigen/Dense>

#include "jointvtx_common.h"

namespace FCCAnalyses {
namespace AlephJointVtx {

// hit-radius rule: production radius bound = layer radius + tolerance [cm]. Layer radii from the detector axis,
// CERN-PPE/94-170: VDET layers (6.5, 11.3 cm), outermost of the eight ITC wire layers (26 cm). Tolerances: 99 %
// quantile of the fitted minus the true transverse vertex radius of simulated K0S and Lambda decaying at 4-8, 8-14 and
// 22-32 cm; they also cover the flat VDET wafers. To be replaced by the shared hit-reach constants once those are in
// the repository.
constexpr double JV_HITR_VDET2 = 6.5 + 1.1;    // >= 2 VDET hits
constexpr double JV_HITR_VDET1 = 11.3 + 1.8;   // 1 VDET hit
constexpr double JV_HITR_ITC = 26. + 2.9;      // >= 1 ITC hit

struct Trk {
  Vec5 p = Vec5::Zero();
  Mat5 C = Mat5::Zero();
  int nvdet = 0;
  int nitc = 0;
};

struct Lin1 {
  double chi2 = JV_NAN;
  double L = JV_NAN;
  LinStatus status = kLinOk;
};

struct Crossings {
  int n = 0;
  double x[2] = {0., 0.};
  double y[2] = {0., 0.};
};

struct IpInfo {
  std::vector<double> ip3, s3, ip2, s2, sg;
};

namespace detail {

inline double sinc(double u) {
  if (std::fabs(u) < 1e-4) return 1. - u * u / 6. + pow4(u) / 120.;
  return std::sin(u) / u;
}

inline double dsinc(double u) {
  if (std::fabs(u) < 1e-4) return -u / 3. + pow3(u) / 30.;
  return (u * std::cos(u) - std::sin(u)) / (u * u);
}

}  // namespace detail

// helix point X, direction a = dX/dL and Jacobian A = dX/dp at transverse arc L; cp, sp = cos, sin of phi0 = p[1]
inline void helix1(const Vec5& p, double L, double cp, double sp, Vec3& X, Vec3& a, Mat35& A) {
  const double D = p[0], p0 = p[1], C = p[2], z0 = p[3], ct = p[4];
  const double u = C * L;
  const double f = L * detail::sinc(u);
  const double ang = p0 + u;
  const double ca = std::cos(ang), sa = std::sin(ang);
  X[0] = -D * sp + f * ca;
  X[1] = D * cp + f * sa;
  X[2] = z0 + ct * L;
  const double ang2 = p0 + 2 * C * L;
  a[0] = std::cos(ang2);
  a[1] = std::sin(ang2);
  a[2] = ct;
  const double dfdC = L * L * detail::dsinc(u);
  A.setZero();
  A(0, 0) = -sp;
  A(1, 0) = cp;
  A(0, 1) = -D * cp - f * sa;
  A(1, 1) = -D * sp + f * ca;
  A(0, 2) = dfdC * ca - f * L * sa;
  A(1, 2) = dfdC * sa + f * L * ca;
  A(2, 3) = 1.;
  A(2, 4) = L;
}

inline void helix1(const Vec5& p, double L, Vec3& X, Vec3& a, Mat35& A) {
  helix1(p, L, std::cos(p[1]), std::sin(p[1]), X, a, A);
}

inline Vec3 pointAt(const Vec5& p, double L) {
  Vec3 X, a;
  Mat35 A;
  helix1(p, L, X, a, A);
  return X;
}

// inverse of a symmetric 3x3 from the analytic cofactors; eigenvalue-floored inverse when
// ill-conditioned. Non-finite input or a failed eigen-decomposition: status != kLinOk, NaN output
inline LinStatus reginv3(const Mat3& M, Mat3& out) {
  const double a = M(0, 0), b = 0.5 * (M(0, 1) + M(1, 0)), c = 0.5 * (M(0, 2) + M(2, 0));
  const double d = M(1, 1), e = 0.5 * (M(1, 2) + M(2, 1)), f = M(2, 2);
  const double A = d * f - e * e, B = c * e - b * f, Cc = b * e - c * d;
  const double det = a * A + b * B + c * Cc;
  const double tr = std::fabs(a) + std::fabs(d) + std::fabs(f);
  out.setZero();
  if (det > 1e-9 * tr * tr * tr && det > 0) {
    out(0, 0) = A / det;
    out(0, 1) = B / det;
    out(0, 2) = Cc / det;
    out(1, 1) = (a * f - c * c) / det;
    out(1, 2) = (b * c - a * e) / det;
    out(2, 2) = (a * d - b * b) / det;
    out(1, 0) = out(0, 1);
    out(2, 0) = out(0, 2);
    out(2, 1) = out(1, 2);
    return kLinOk;
  }
  Mat3 S;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) S(i, j) = 0.5 * (M(i, j) + M(j, i));
  if (!allFinite(S)) {
    setNaN(out);
    return kLinNonFinite;
  }
  Vec3 w;
  Mat3 U;
  if (detail::eigh3(S, w, U) != 0) {
    setNaN(out);
    return kLinSingular;
  }
  const double wmax = npMaximum(npMaximum(std::fabs(w[0]), std::fabs(w[1])), std::fabs(w[2]));
  const double fl = 1e-12 * wmax;
  for (int k = 0; k < 3; ++k) {
    double wc = w[k] > fl ? w[k] : fl;
    if (wc <= 0) wc = 1e-300;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) out(i, j) += U(i, k) * U(j, k) / wc;
  }
  return kLinOk;
}

// linearise one track at the point x, starting at arc L0 (12 Newton steps at most, step
// clamp 50 cm): fills X (closest approach) and D (3x3 weight); returns chi2 and the arc
inline Lin1 lin1(const Vec5& p, const Mat5& cov, const Vec3& x, double L0, Vec3& X, Mat3& D) {
  Lin1 res;
  Vec3 a = Vec3::Zero();
  Mat35 A = Mat35::Zero();
  Mat3 W = Mat3::Zero();
  Vec3 aw = Vec3::Zero();
  const double cp = std::cos(p[1]), sp = std::sin(p[1]);
  double L = L0;
  for (int it = 0; it < 12; ++it) {
    helix1(p, L, cp, sp, X, a, A);
    res.status = reginv3(detail::sandwich(A, cov), W);
    if (res.status != kLinOk) {
      setNaN(X);
      setNaN(D);
      return res;
    }
    aw = matVecC<3, 3>(W, a);
    const double den = dotSeq(a, aw);
    if (!(den > 0) || !std::isfinite(den)) break;
    const Vec3 dxv = x - X;
    double dL = dotSeq(dxv, aw) / den;
    if (std::fabs(dL) > 50.) dL = dL > 0 ? 50. : -50.;
    L = L + dL;
    if (std::fabs(dL) < 1e-10) break;
  }
  helix1(p, L, cp, sp, X, a, A);
  res.status = reginv3(detail::sandwich(A, cov), W);
  if (res.status != kLinOk) {
    setNaN(X);
    setNaN(D);
    return res;
  }
  aw = matVecC<3, 3>(W, a);
  const double den = dotSeq(a, aw);
  const bool good = den > 0 && std::isfinite(den);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) D(i, j) = W(i, j) - (good ? aw[i] * aw[j] / den : 0.);
  for (int i = 0; i < 3; ++i)
    for (int j = i + 1; j < 3; ++j) {
      const double s = 0.5 * (D(i, j) + D(j, i));
      D(i, j) = s;
      D(j, i) = s;
    }
  const Vec3 r = x - X;
  res.chi2 = detail::quadForm(r, D, false);
  res.L = L;
  return res;
}

// lin1 for the tracks idx at x; Ls (one start arc per track of t) is updated. Xo, Do, chi2o
// have one entry per index. Stops at the first failure and returns its status
inline LinStatus linAll(const std::vector<Trk>& t, const std::vector<int>& idx, const Vec3& x,
                        std::vector<double>& Ls, std::vector<Vec3>& Xo, std::vector<Mat3>& Do,
                        std::vector<double>& chi2o) {
  const size_t n = idx.size();
  Xo.resize(n);
  Do.resize(n);
  chi2o.resize(n);
  for (size_t k = 0; k < n; ++k) {
    const int i = idx[k];
    const Lin1 r = lin1(t[i].p, t[i].C, x, Ls[i], Xo[k], Do[k]);
    if (r.status != kLinOk) return r.status;
    chi2o[k] = r.chi2;
    Ls[i] = r.L;
  }
  return kLinOk;
}

// chi2 of the tracks idx to the point x (track covariance only) and the arc at the closest approach
inline LinStatus chi2ToPoint(const std::vector<Trk>& t, const std::vector<int>& idx, const Vec3& x,
                             std::vector<double>& chi2, std::vector<double>& Lo) {
  const size_t n = idx.size();
  chi2.assign(n, 0.);
  Lo.assign(n, 0.);
  Vec3 X;
  Mat3 D;
  for (size_t k = 0; k < n; ++k) {
    const Trk& tr = t[idx[k]];
    const Lin1 r = lin1(tr.p, tr.C, x, 0., X, D);
    if (r.status != kLinOk) return r.status;
    chi2[k] = r.chi2;
    Lo[k] = r.L;
  }
  return kLinOk;
}

inline LinStatus chi2ToPoint(const std::vector<Trk>& t, const Vec3& x, std::vector<double>& chi2,
                             std::vector<double>& Lo) {
  std::vector<int> idx(t.size());
  for (size_t k = 0; k < idx.size(); ++k) idx[k] = static_cast<int>(k);
  return chi2ToPoint(t, idx, x, chi2, Lo);
}

// momentum vector [GeV] at arc L
inline Vec3 momentumAt(const Vec5& p, double L, double ptHalf = JV_PT_HALF) {
  const double pt = ptHalf / std::max(std::fabs(p[2]), 1e-12);
  const double ang = p[1] + 2 * p[2] * L;
  return Vec3(pt * std::cos(ang), pt * std::sin(ang), pt * p[4]);
}

// per track w.r.t. the vertex pv (covariance pvcov) and a lifetime-sign axis (axes[jetOf[t]], or
// the axis nearest to the track direction when jetOf[t] < 0): 3D impact parameter, signed 3D
// significance, signed transverse impact parameter and significance, lifetime sign.
// ptHalf enters only the choice of the nearest axis
inline LinStatus ipInfo(const std::vector<Trk>& t, const Vec3& pv, const Mat3& pvcov,
                        const std::vector<Vec3>& axes, const std::vector<int>& jetOf, IpInfo& out,
                        double ptHalf = JV_PT_HALF) {
  const size_t n = t.size();
  out.ip3.assign(n, 0.);
  out.s3.assign(n, 0.);
  out.ip2.assign(n, 0.);
  out.s2.assign(n, 0.);
  out.sg.assign(n, 0.);
  const int na = static_cast<int>(axes.size());
  Vec3 X, a;
  Mat3 D;
  Mat35 A;
  for (size_t k = 0; k < n; ++k) {
    const Lin1 l = lin1(t[k].p, t[k].C, pv, 0., X, D);
    if (l.status != kLinOk) return l.status;
    helix1(t[k].p, l.L, X, a, A);
    const Vec3 r = X - pv;
    const Mat3 Ct = detail::sandwich(A, t[k].C);
    Mat3 Cs;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) Cs(i, j) = Ct(i, j) + pvcov(i, j);
    const double d3 = std::sqrt(dotSeq(r, r));
    double var3 = 1.;
    if (d3 > 0) {
      const Vec3 u = r / d3;
      var3 = detail::quadForm(u, Cs, false);
    }
    const Vec3 mom = momentumAt(t[k].p, l.L, ptHalf);
    int j = jetOf[k];
    if (j < 0) {
      double best = -2.;
      for (int q = 0; q < na; ++q) {
        const double dd = dotSeq(axes[q], mom);
        if (dd > best) {
          best = dd;
          j = q;
        }
      }
      if (j < 0) j += na;
    }
    if (j < 0 || j >= na) {
      out.ip3[k] = out.s3[k] = out.ip2[k] = out.s2[k] = out.sg[k] = JV_NAN;
      continue;
    }
    const Vec3 re = r - (dotSeq(r, a) / dotSeq(a, a)) * a;
    const double sgn = dotSeq(re, axes[j]) >= 0 ? 1. : -1.;
    out.sg[k] = sgn;
    out.ip3[k] = d3;
    out.s3[k] = sgn * d3 / std::sqrt(std::max(var3, 1e-14));
    const Vec3 r2(r[0], r[1], 0.);
    const double d2 = std::sqrt(dotSeq(r2, r2));
    double var2 = 1.;
    if (d2 > 0) {
      const Vec3 u2 = r2 / d2;
      var2 = detail::quadForm(u2, Cs, false);
    }
    out.ip2[k] = sgn * d2;
    out.s2[k] = sgn * d2 / std::sqrt(std::max(var2, 1e-14));
  }
  return kLinOk;
}

// per-track compatibility with the vertex at x with information matrix H (inverse covariance):
// mode 0: chi2 = r^T D r; mode 1: chi2 = r^T (D - D (H + D)^-1 D) r; mode 2: leave-one-out
// (the track's weighted contribution w D removed first). dens = ln det(Hi + D) - ln det(Hi);
// cosf = cosine between the track direction at its closest approach and fdir; Lo = arc
inline LinStatus chi2Vertex(const std::vector<Trk>& t, const Vec3& x, const Mat3& H, const std::vector<double>& w,
                            int mode, double eps, const Vec3& fdir, std::vector<double>& chi2,
                            std::vector<double>& dens, std::vector<double>& cosf, std::vector<double>& Lo) {
  const size_t n = t.size();
  chi2.assign(n, 0.);
  dens.assign(n, 0.);
  cosf.assign(n, 0.);
  Lo.assign(n, 0.);
  Vec3 X;
  Mat3 D;
  for (size_t k = 0; k < n; ++k) {
    const Lin1 l = lin1(t[k].p, t[k].C, x, 0., X, D);
    if (l.status != kLinOk) return l.status;
    const double L = l.L;
    Lo[k] = L;
    const double ang = t[k].p[1] + 2 * t[k].p[2] * L;
    const double ct = t[k].p[4];
    const double nrm = std::sqrt(1. + ct * ct);
    cosf[k] = (std::cos(ang) * fdir[0] + std::sin(ang) * fdir[1] + ct * fdir[2]) / nrm;
    if (mode == 0) {
      chi2[k] = l.chi2;
      continue;
    }
    Mat3 Hi;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) Hi(i, j) = H(i, j) + eps * (i == j ? 1. : 0.);
    Vec3 xi = x;
    if (mode == 2 && w[k] > 0) {
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) Hi(i, j) = Hi(i, j) - w[k] * D(i, j);
      Mat3 Ci;
      const LinStatus s = reginv3(Hi, Ci);
      if (s != kLinOk) return s;
      const Vec3 g = w[k] * matVecC<3, 3>(D, Vec3(x - X));
      xi = x + matVecC<3, 3>(Ci, g);
    }
    const Vec3 r = xi - X;
    Mat3 M;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) M(i, j) = Hi(i, j) + D(i, j);
    Mat3 Mi;
    const LinStatus s = reginv3(M, Mi);
    if (s != kLinOk) return s;
    const Vec3 Dr = matVecC<3, 3>(D, r);
    chi2[k] = std::max(dotSeq(r, Dr) - dotSeq(Dr, matVecC<3, 3>(Mi, Dr)), 0.);
    const double dH = npDet<3>(Hi);
    const double dM = npDet<3>(M);
    dens[k] = (dH > 0 && dM > 0) ? std::log(dM / dH) : 0.;
  }
  return kLinOk;
}

// upper bound [cm] on the transverse production radius from the hits: VDET (>= 2, 1), ITC; inf otherwise
inline double maxRadiusFromHits(int nvdet, int nitc) {
  double r = JV_INF;
  if (nvdet >= 2) r = JV_HITR_VDET2;
  else if (nvdet == 1) r = JV_HITR_VDET1;
  if (nitc >= 1) r = npMinimum(r, JV_HITR_ITC);
  return r;
}

inline std::vector<double> maxRadiusFromHits(const std::vector<Trk>& t) {
  std::vector<double> r(t.size());
  for (size_t k = 0; k < t.size(); ++k) r[k] = maxRadiusFromHits(t[k].nvdet, t[k].nitc);
  return r;
}

// xy crossing points of the two track circles: 0, 1 or 2 (disjoint or nested circles: the
// midpoint of the gap)
inline Crossings circleCrossings(const Vec5& p1, const Vec5& p2) {
  Crossings out;
  double xc[2], yc[2], rc[2];
  const Vec5* ps[2] = {&p1, &p2};
  for (int k = 0; k < 2; ++k) {
    const double D = (*ps[k])[0], phi0 = (*ps[k])[1], C = (*ps[k])[2];
    if (C == 0) return out;
    xc[k] = -(D + 0.5 / C) * std::sin(phi0);
    yc[k] = (D + 0.5 / C) * std::cos(phi0);
    rc[k] = 0.5 / std::fabs(C);
  }
  const double dx = xc[1] - xc[0], dy = yc[1] - yc[0];
  const double d = std::hypot(dx, dy);
  if (!std::isfinite(d) || !(d > 0)) return out;
  const double ux = dx / d, uy = dy / d;
  if (d >= rc[0] + rc[1]) {
    const double tt = rc[0] + 0.5 * (d - rc[0] - rc[1]);
    out.n = 1;
    out.x[0] = xc[0] + tt * ux;
    out.y[0] = yc[0] + tt * uy;
    return out;
  }
  if (d <= std::fabs(rc[0] - rc[1])) {
    const double s = rc[0] > rc[1] ? 1. : -1.;
    out.n = 1;
    out.x[0] = 0.5 * (xc[0] + xc[1] + s * (rc[0] + rc[1]) * ux);
    out.y[0] = 0.5 * (yc[0] + yc[1] + s * (rc[0] + rc[1]) * uy);
    return out;
  }
  const double a = (detail::powSquare(rc[0]) - detail::powSquare(rc[1]) + d * d) / (2. * d);
  const double h2 = detail::powSquare(rc[0]) - a * a;
  const double X = xc[0] + a * ux, Y = yc[0] + a * uy;
  if (h2 > 0) {
    const double h = std::sqrt(h2);
    out.n = 2;
    out.x[0] = X - h * uy;
    out.y[0] = Y + h * ux;
    out.x[1] = X + h * uy;
    out.y[1] = Y - h * ux;
    return out;
  }
  out.n = 1;
  out.x[0] = X;
  out.y[0] = Y;
  return out;
}

// transverse arc of the circle point nearest to (x, y)
inline double arcToXY(const Vec5& p, double x, double y) {
  const double D = p[0], phi0 = p[1], C = p[2];
  const double x0 = -D * std::sin(phi0), y0 = D * std::cos(phi0);
  const double xc = -(D + 0.5 / C) * std::sin(phi0), yc = (D + 0.5 / C) * std::cos(phi0);
  const double a0 = std::atan2(y0 - yc, x0 - xc);
  const double a1 = std::atan2(y - yc, x - xc);
  return wrapPi(a1 - a0) / (2 * C);
}

// arc of the track's xy point nearest to x (geometric); 0 for a straight track
inline double arcNear(const Vec5& p, const Vec3& x) {
  if (p[2] == 0) return 0.;
  return arcToXY(p, x[0], x[1]);
}

}  // namespace AlephJointVtx
}  // namespace FCCAnalyses

#endif
