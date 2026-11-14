#ifndef JOINTVTX_CONV_H
#define JOINTVTX_CONV_H

/*
  Photon-conversion finder of the joint vertexer (gamma -> e+ e-).

  TP: tangent point of the two track circles (closed form) with propagated
  errors. CF: constrained conversion fit, 7 parameters (vertex x, y, z, common
  direction phi, common tanLambda, the two curvatures) fitted to the 10 helix
  parameters of the two tracks: Levenberg-Marquardt with a numerical Jacobian,
  started at the TP, with a second-order check and steps out of saddle points.

  Track parameters in the fitter convention p = (D, phi0, C, z0, tanLambda),
  pT = K / |C|. Units: cm, GeV.

  The fit internals use plain arrays, which keep the operation order fixed.
*/

#include <array>
#include <cmath>
#include <vector>

#include <Eigen/Dense>

#include "jointvtx_common.h"

namespace FCCAnalyses {
namespace AlephJointVtx {

constexpr int JV_CONV_NTP = 27;   // TP result vector
constexpr int JV_CONV_NCF = 16;   // CF result vector
constexpr int JV_CONV_NQ = 55;    // quality row: TP, CF, derived columns
constexpr int JV_CF_MAXIT = 500;  // maximum Levenberg-Marquardt iterations
constexpr int JV_CF_NESCAPE = 5;  // maximum restarts from a saddle point
constexpr double JV_CF_NEG = 1e-3;  // negative curvature: Hessian eigenvalue of chi2 / 2 below -JV_CF_NEG
constexpr double JV_CF_TOL = 1e-6;  // convergence: Newton decrement g^T H^-1 g below this

using ConvTp = std::array<double, JV_CONV_NTP>;
using ConvCf = std::array<double, JV_CONV_NCF>;
using ConvQ = std::array<double, JV_CONV_NQ>;

// TP vector: ok, vertex, gap, dz, dct, arcs, photon direction, vertex covariance,
// errors and correlations of (gap, dz, dct), m_ee, photon momentum
enum ConvTpCol : int {
  kTpOk, kTpX, kTpY, kTpZ, kTpGap, kTpDz, kTpDct, kTpL1, kTpL2, kTpPhi, kTpCt,
  kTpSxx, kTpSyy, kTpSzz, kTpSxy, kTpSxz, kTpSyz, kTpSgap, kTpSdz, kTpSdct,
  kTpRgz, kTpRgc, kTpRzc, kTpMee, kTpPx, kTpPy, kTpPz
};
// CF vector: ok (1 local minimum, 2 not converged, 0 undefined), fitted parameters,
// chi2, iterations, photon momentum, chi2 per track, final Newton decrement
enum ConvCfCol : int {
  kCfOk, kCfX, kCfY, kCfZ, kCfPhi, kCfCt, kCfC1, kCfC2, kCfChi2, kCfNiter,
  kCfPx, kCfPy, kCfPz, kCfChi21, kCfChi22, kCfDec
};
// quality row: TP columns, then CF columns (offset kQCf), then these
enum ConvQCol : int {
  kQTpR = JV_CONV_NTP + JV_CONV_NCF, kQCfR, kQTpChi2, kQTpDgam, kQTpCosp, kQCfDgam,
  kQCfCosp, kQCfSl, kQNhx1, kQNhx2, kQP1, kQP2
};
constexpr int kQCf = JV_CONV_NTP;

// prefilter on the TP without errors: max |gap|, max |dz|, max |dct|, min arc, max radius
struct ConvPre {
  double gap = 3.0;
  double dz = 15.0;
  double dct = 0.2;
  double arc = -5.0;
  double r = 200.0;
};

// one cut on a quality column: lo <= value < hi (|value| for the signed columns
// dgam, gap, dz, dct); a NaN value fails
struct ConvCut {
  int col;
  bool hasLo;
  double lo;
  bool hasHi;
  double hi;
};

// loose selection (supply to the joint fit): CF chi2 < 200, TP m_ee < 50 MeV, cos pointing >= 0.9
inline std::vector<ConvCut> convLoose() {
  return {{kQCf + kCfChi2, false, 0.0, true, 200.0},
          {kTpMee, false, 0.0, true, 0.05},
          {kQCfCosp, true, 0.9, false, 0.0}};
}

namespace detail {

using A7 = std::array<double, 7>;
using A10 = std::array<double, 10>;

// symmetric eigen decomposition (cyclic Jacobi): eigenvalues ascending, eigenvectors in the columns of V
template <int N>
inline void symEigen(const double S[N][N], double ev[N], double V[N][N]) {
  double a[N][N], v[N][N];
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j) {
      a[i][j] = S[i][j];
      v[i][j] = i == j ? 1.0 : 0.0;
    }
  for (int sweep = 0; sweep < 100; ++sweep) {
    double off = 0.0, tot = 0.0;
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j) {
        tot += a[i][j] * a[i][j];
        if (i != j) off += a[i][j] * a[i][j];
      }
    if (!(off > 1e-34 * tot)) break;
    for (int p = 0; p < N - 1; ++p)
      for (int q = p + 1; q < N; ++q) {
        if (a[p][q] == 0.0) continue;
        const double theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
        double t;
        if (std::fabs(theta) > 1e150) t = 0.5 / theta;
        else t = (theta >= 0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
        const double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
        for (int k = 0; k < N; ++k) {
          const double akp = a[k][p], akq = a[k][q];
          a[k][p] = c * akp - s * akq;
          a[k][q] = s * akp + c * akq;
        }
        for (int k = 0; k < N; ++k) {
          const double apk = a[p][k], aqk = a[q][k];
          a[p][k] = c * apk - s * aqk;
          a[q][k] = s * apk + c * aqk;
        }
        a[p][q] = 0.0;
        a[q][p] = 0.0;
        for (int k = 0; k < N; ++k) {
          const double vkp = v[k][p], vkq = v[k][q];
          v[k][p] = c * vkp - s * vkq;
          v[k][q] = s * vkp + c * vkq;
        }
      }
  }
  int order[N];
  for (int i = 0; i < N; ++i) order[i] = i;
  for (int i = 1; i < N; ++i)
    for (int k = i; k > 0 && npLess(a[order[k]][order[k]], a[order[k - 1]][order[k - 1]]); --k)
      std::swap(order[k], order[k - 1]);
  for (int c = 0; c < N; ++c) {
    ev[c] = a[order[c]][order[c]];
    for (int k = 0; k < N; ++k) V[k][c] = v[k][order[c]];
  }
}

// ---------------------------------------------------------------- helix geometry
inline void centre(const double* p, double& xc, double& yc, double& rs) {
  rs = 0.5 / p[2];
  xc = -(p[0] + rs) * std::sin(p[1]);
  yc = (p[0] + rs) * std::cos(p[1]);
}

// transverse arc L of the circle point nearest to (x, y), |2 C L| < pi
inline double arcTo(const double* p, double x, double y) {
  double xc, yc, rs;
  centre(p, xc, yc, rs);
  const double x0 = -p[0] * std::sin(p[1]), y0 = p[0] * std::cos(p[1]);
  return wrapPi(std::atan2(y - yc, x - xc) - std::atan2(y0 - yc, x0 - xc)) / (2. * p[2]);
}

inline void momAt(const double* p, double L, double K, double* P) {
  const double pt = K / std::max(std::fabs(p[2]), 1e-12);
  const double ang = p[1] + 2. * p[2] * L;
  P[0] = pt * std::cos(ang);
  P[1] = pt * std::sin(ang);
  P[2] = pt * p[4];
}

// helix (D, phi0, C, z0, ct) through (xv, yv, zv) with direction phi there; the arc from
// the perigee to the vertex is taken on the turn nearest to the reference arc Lref
inline void perigeeFromRef(double xv, double yv, double zv, double phi, double ct, double C, double Lref, double* out) {
  const double rs = 0.5 / C;
  const double xc = xv - rs * std::sin(phi), yc = yv + rs * std::cos(phi);
  const double ra = std::hypot(xc, yc);
  const double rho = C > 0 ? ra : -ra;
  const double phi0 = std::atan2(-xc / rho, yc / rho);
  const double D = rho - rs;
  const double x0 = -D * std::sin(phi0), y0 = D * std::cos(phi0);
  double a = std::atan2(yv - yc, xv - xc) - std::atan2(y0 - yc, x0 - xc);
  const double a0 = 2. * C * Lref;
  a = a0 + wrapPi(a - a0);
  const double L = a / (2. * C);
  out[0] = D;
  out[1] = phi0;
  out[2] = C;
  out[3] = zv - ct * L;
  out[4] = ct;
}

// ---------------------------------------------------------------- tangent point
// o = (x, y, z1, z2, gap, dct, L1, L2); false if undefined (zero or same-sign curvatures, d = 0)
inline bool tpRaw(const double* p1, const double* p2, double* o) {
  if (p1[2] == 0. || p2[2] == 0. || p1[2] * p2[2] > 0.) return false;
  double x1, y1, r1, x2, y2, r2;
  centre(p1, x1, y1, r1);
  centre(p2, x2, y2, r2);
  const double R1 = std::fabs(r1), R2 = std::fabs(r2);
  const double dx = x2 - x1, dy = y2 - y1;
  const double d = std::hypot(dx, dy);
  if (!(d > 0.)) return false;
  const double ux = dx / d, uy = dy / d;
  const double gap = d - R1 - R2;
  const double t = R1 + 0.5 * gap;
  const double x = x1 + t * ux, y = y1 + t * uy;
  const double L1 = arcTo(p1, x, y), L2 = arcTo(p2, x, y);
  o[0] = x;
  o[1] = y;
  o[2] = p1[3] + p1[4] * L1;
  o[3] = p2[3] + p2[4] * L2;
  o[4] = gap;
  o[5] = p1[4] - p2[4];
  o[6] = L1;
  o[7] = L2;
  return true;
}

inline double tpStep(const double* p, int k) {
  if (k == 2) return 1e-6 * std::max(std::fabs(p[2]), 1e-5);
  if (k == 1 || k == 4) return 1e-7;
  return 1e-6;
}

// ---------------------------------------------------------------- products (inner size 5: sum5Blas; longer
// products are summed left to right)
// C (R x M) = A (R x K) @ B (K x M)
template <int R, int K, int M>
inline void mm(const double A[R][K], const double B[][M], double C[R][M]) {
  for (int i = 0; i < R; ++i)
    for (int j = 0; j < M; ++j) {
      if (K == 5) {
        double p[5];
        for (int k = 0; k < 5; ++k) p[k] = A[i][k] * B[k][j];
        C[i][j] = sum5Blas(p, i, j, R, M);
        continue;
      }
      double s = 0.0;
      for (int k = 0; k < K; ++k) s += A[i][k] * B[k][j];
      C[i][j] = s;
    }
}
// C = A @ B^T with B given as (M x K)
template <int R, int K, int M>
inline void mmT(const double A[R][K], const double B[M][K], double C[R][M]) {
  for (int i = 0; i < R; ++i)
    for (int j = 0; j < M; ++j) {
      if (K == 5) {
        double p[5];
        for (int k = 0; k < 5; ++k) p[k] = A[i][k] * B[j][k];
        C[i][j] = sum5Blas(p, i, j, R, M);
        continue;
      }
      double s = 0.0;
      for (int k = 0; k < K; ++k) s += A[i][k] * B[j][k];
      C[i][j] = s;
    }
}
// y = v @ M (gemv on the transpose of a C-ordered M: left to right)
template <int K, int M>
inline void vm(const double* v, const double A[K][M], double* y) {
  for (int j = 0; j < M; ++j) {
    double s = 0.0;
    for (int k = 0; k < K; ++k) s += v[k] * A[k][j];
    y[j] = s;
  }
}
// y = A @ x for a C-ordered A: even and odd columns summed separately, then added
template <int R, int N>
inline void mvEO(const double A[R][N], const double* x, double* y) {
  for (int i = 0; i < R; ++i) {
    double e = A[i][0] * x[0];
    double o = A[i][1] * x[1];
    for (int k = 2; k < N; ++k) {
      if (k % 2 == 0) e += A[i][k] * x[k];
      else o += A[i][k] * x[k];
    }
    y[i] = e + o;
  }
}
template <int N>
inline double dot(const double* a, const double* b) {
  double s = 0.0;
  for (int k = 0; k < N; ++k) s += a[k] * b[k];
  return s;
}

// ---------------------------------------------------------------- constrained fit
struct CfData {
  double m[10];
  double W[10][10];
  double Lr[2];
};

inline void cfRes(const A7& th, const CfData& D, A10& r) {
  double a[5], b[5];
  perigeeFromRef(th[0], th[1], th[2], th[3], th[4], th[5], D.Lr[0], a);
  perigeeFromRef(th[0], th[1], th[2], th[3], th[4], th[6], D.Lr[1], b);
  for (int k = 0; k < 5; ++k) {
    r[k] = a[k] - D.m[k];
    r[5 + k] = b[k] - D.m[5 + k];
  }
  r[1] = wrapPi(r[1]);
  r[6] = wrapPi(r[6]);
}

// r @ W @ r
inline double cfQuad(const A10& r, const double W[10][10]) {
  double y[10];
  vm<10, 10>(r.data(), W, y);
  return dot<10>(y, r.data());
}

inline double cfChi2(const A7& th, const CfData& D) {
  A10 r;
  cfRes(th, D, r);
  return cfQuad(r, D.W);
}

// numerical Jacobian (10, 7) of the residuals, central differences
inline void cfJac(const A7& th, const CfData& D, double J[10][7]) {
  A10 rp, rm;
  for (int k = 0; k < 7; ++k) {
    double h;
    if (k < 3) h = 1e-6;
    else if (k < 5) h = 1e-7;
    else h = 1e-6 * std::max(std::fabs(th[k]), 1e-5);
    A7 tq = th, tm = th;
    tq[k] += h;
    tm[k] -= h;
    cfRes(tq, D, rp);
    cfRes(tm, D, rm);
    for (int a = 0; a < 10; ++a)
      J[a][k] = (a == 1 || a == 6) ? wrapPi(rp[a] - rm[a]) / (2. * h) : (rp[a] - rm[a]) / (2. * h);
  }
}

// H = J^T W J and, if r is given, g = J^T W r
inline void cfNormal(const double J[10][7], const double W[10][10], const double* r, double H[7][7], double* g) {
  double JtW[7][10];
  for (int i = 0; i < 7; ++i)
    for (int j = 0; j < 10; ++j) {
      double s = 0.0;
      for (int k = 0; k < 10; ++k) s += J[k][i] * W[k][j];
      JtW[i][j] = s;
    }
  mm<7, 10, 7>(JtW, J, H);
  if (r) mvEO<7, 10>(JtW, r, g);
}

// Newton decrement g^T H^-1 g; inf if H is singular or the result non-finite
inline double cfDecrement(const double H[7][7], const double* g) {
  double d[7];
  if (npSolveA<7>(H, g, d) != kLinOk) return JV_INF;
  const double v = dot<7>(g, d);
  return std::isfinite(v) ? v : JV_INF;
}

// Levenberg-Marquardt from th (chi2 = its chi2); returns the iterations
inline int cfLm(A7& th, double& chi2, const CfData& D) {
  A10 r, rp;
  double J[10][7], H[7][7], g[7];
  cfRes(th, D, r);
  double lam = 1e-3;
  for (int it = 0; it < JV_CF_MAXIT; ++it) {
    cfJac(th, D, J);
    cfNormal(J, D.W, r.data(), H, g);
    if (cfDecrement(H, g) < JV_CF_TOL) return it + 1;
    bool accepted = false;
    for (int t = 0; t < 40; ++t) {
      double A[7][7];
      for (int i = 0; i < 7; ++i)
        for (int j = 0; j < 7; ++j) A[i][j] = H[i][j];
      for (int k = 0; k < 7; ++k) A[k][k] *= 1. + lam;
      double x[7];
      if (npSolveA<7>(A, g, x) != kLinOk) {
        lam *= 10.;
        continue;
      }
      A7 tn;
      for (int k = 0; k < 7; ++k) tn[k] = th[k] + (-x[k]);
      tn[3] = wrapPi(tn[3]);
      cfRes(tn, D, rp);
      const double cn = cfQuad(rp, D.W);
      if (std::isfinite(cn) && cn < chi2) {
        th = tn;
        chi2 = cn;
        r = rp;
        lam = std::max(lam * 0.2, 1e-9);
        accepted = true;
        break;
      }
      lam *= 10.;
    }
    if (!accepted) return it + 1;
  }
  return JV_CF_MAXIT;
}

// second-order check: full numerical Hessian of chi2 / 2 in units of the Gauss-Newton
// errors; with negative curvature, the lowest-chi2 step along the negative eigenvectors
inline bool cfEscape(const A7& th, double chi2, const CfData& D, A7& best, double& cb) {
  best = th;
  cb = chi2;
  double J[10][7], H[7][7], C[7][7];
  cfJac(th, D, J);
  cfNormal(J, D.W, nullptr, H, nullptr);
  if (nbInvA<7>(H, C) != kLinOk) return false;
  double sg[7];
  for (int k = 0; k < 7; ++k) {
    if (!(C[k][k] > 0)) return false;
    sg[k] = std::sqrt(C[k][k]);
  }
  const double h = 0.01;
  double F[7][7];
  for (int a = 0; a < 7; ++a)
    for (int b = a; b < 7; ++b) {
      double v;
      if (a == b) {
        A7 tq = th, tm = th;
        tq[a] += h * sg[a];
        tm[a] -= h * sg[a];
        v = (cfChi2(tq, D) - 2. * chi2 + cfChi2(tm, D)) / (h * h);
      } else {
        A7 t1 = th, t2 = th, t3 = th, t4 = th;
        t1[a] += h * sg[a];
        t1[b] += h * sg[b];
        t2[a] += h * sg[a];
        t2[b] -= h * sg[b];
        t3[a] -= h * sg[a];
        t3[b] += h * sg[b];
        t4[a] -= h * sg[a];
        t4[b] -= h * sg[b];
        v = (cfChi2(t1, D) - cfChi2(t2, D) - cfChi2(t3, D) + cfChi2(t4, D)) / (4. * h * h);
      }
      F[a][b] = 0.5 * v;
      F[b][a] = 0.5 * v;
    }
  for (int a = 0; a < 7; ++a)
    for (int b = 0; b < 7; ++b)
      if (!std::isfinite(F[a][b])) return false;
  double ev[7], V[7][7];
  symEigen<7>(F, ev, V);
  if (!(ev[0] < -JV_CF_NEG)) return false;
  const double kT[12] = {0.5, -0.5, 1., -1., 2., -2., 4., -4., 8., -8., 16., -16.};
  for (int q = 0; q < 7; ++q) {
    if (!(ev[q] < -JV_CF_NEG)) break;
    for (double t : kT) {
      A7 tn;
      for (int k = 0; k < 7; ++k) tn[k] = th[k] + t * sg[k] * V[k][q];
      tn[3] = wrapPi(tn[3]);
      const double c = cfChi2(tn, D);
      if (std::isfinite(c) && c < cb) {
        best = tn;
        cb = c;
      }
    }
  }
  return true;
}

inline void toArr5(const Mat5& m, double a[5][5]) {
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j) a[i][j] = m(i, j);
}

}  // namespace detail

// TP result vector (ConvTpCol); ok = 0 and NaN entries when undefined
inline ConvTp convPairTp(const Vec5& P1v, const Mat5& C1m, const Vec5& P2v, const Mat5& C2m, double K = JV_PT_HALF) {
  using namespace detail;
  ConvTp r;
  r.fill(JV_NAN);
  r[kTpOk] = 0.;
  const double* p1 = P1v.data();
  const double* p2 = P2v.data();
  double o[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  if (!tpRaw(p1, p2, o)) return r;
  // Jacobian of (x, y, z1, z2, gap) w.r.t. the 10 parameters (dct is linear)
  double J1[5][5] = {}, J2[5][5] = {};
  double q[8] = {0, 0, 0, 0, 0, 0, 0, 0}, qm[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  for (int side = 0; side < 2; ++side)
    for (int k = 0; k < 5; ++k) {
      double pp[5], pm[5];
      for (int i = 0; i < 5; ++i) pp[i] = side == 0 ? p1[i] : p2[i];
      for (int i = 0; i < 5; ++i) pm[i] = pp[i];
      const double h = tpStep(pp, k);
      pp[k] += h;
      pm[k] -= h;
      bool okp, okm;
      if (side == 0) {
        okp = tpRaw(pp, p2, q);
        okm = tpRaw(pm, p2, qm);
      } else {
        okp = tpRaw(p1, pp, q);
        okm = tpRaw(p1, pm, qm);
      }
      if (!(okp && okm)) return r;
      for (int a = 0; a < 5; ++a) {
        const double v = (q[a] - qm[a]) / (2. * h);
        if (side == 0) J1[a][k] = v;
        else J2[a][k] = v;
      }
    }
  double c1[5][5], c2[5][5];
  toArr5(C1m, c1);
  toArr5(C2m, c2);
  double J1c[5][5], J2c[5][5], V1[5][5], V2[5][5], V[5][5];
  mm<5, 5, 5>(J1, c1, J1c);
  mm<5, 5, 5>(J2, c2, J2c);
  mmT<5, 5, 5>(J1c, J1, V1);
  mmT<5, 5, 5>(J2c, J2, V2);
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j) V[i][j] = V1[i][j] + V2[i][j];
  const double g1[5] = {0., 0., 0., 0., 1.}, g2[5] = {0., 0., 0., 0., -1.};
  double Vc1[5], Vc2[5], Vc[5];
  mvEO<5, 5>(J1c, g1, Vc1);
  mvEO<5, 5>(J2c, g2, Vc2);
  for (int i = 0; i < 5; ++i) Vc[i] = Vc1[i] + Vc2[i];
  double y1[5], y2[5];
  vm<5, 5>(g1, c1, y1);
  vm<5, 5>(g2, c2, y2);
  const double vcc = dot<5>(y1, g1) + dot<5>(y2, g2);
  const double s11 = V[2][2], s22 = V[3][3], s12 = V[2][3];
  const double den = s11 + s22 - 2. * s12;
  const double w = den > 0 ? (s22 - s12) / den : 0.5;
  const double T[3][5] = {{1., 0., 0., 0., 0.}, {0., 1., 0., 0., 0.}, {0., 0., w, 1. - w, 0.}};
  double TV[3][5], Cx[3][3];
  mm<3, 5, 5>(T, V, TV);
  mmT<3, 5, 3>(TV, T, Cx);
  const double gz[5] = {0., 0., 1., -1., 0.};
  double yz[5];
  vm<5, 5>(gz, V, yz);
  const double vzz = dot<5>(yz, gz);
  const double vgg = V[4][4];
  const double vgz = dot<5>(V[4], gz);
  const double vgc = Vc[4];
  const double vzc = dot<5>(gz, Vc);
  const double L1 = o[6], L2 = o[7];
  double P1[3], P2[3], P[3];
  momAt(p1, L1, K, P1);
  momAt(p2, L2, K, P2);
  const double E = std::sqrt(dot<3>(P1, P1) + JV_M_E * JV_M_E) + std::sqrt(dot<3>(P2, P2) + JV_M_E * JV_M_E);
  for (int i = 0; i < 3; ++i) P[i] = P1[i] + P2[i];
  r[kTpOk] = 1.;
  r[kTpX] = o[0];
  r[kTpY] = o[1];
  r[kTpZ] = w * o[2] + (1. - w) * o[3];
  r[kTpGap] = o[4];
  r[kTpDz] = o[2] - o[3];
  r[kTpDct] = o[5];
  r[kTpL1] = L1;
  r[kTpL2] = L2;
  r[kTpPhi] = std::atan2(P[1], P[0]);
  r[kTpCt] = P[2] / std::max(std::hypot(P[0], P[1]), 1e-12);
  r[kTpSxx] = Cx[0][0];
  r[kTpSyy] = Cx[1][1];
  r[kTpSzz] = Cx[2][2];
  r[kTpSxy] = Cx[0][1];
  r[kTpSxz] = Cx[0][2];
  r[kTpSyz] = Cx[1][2];
  r[kTpSgap] = std::sqrt(std::max(vgg, 0.));
  r[kTpSdz] = std::sqrt(std::max(vzz, 0.));
  r[kTpSdct] = std::sqrt(std::max(vcc, 0.));
  r[kTpRgz] = vgz / std::max(r[kTpSgap] * r[kTpSdz], 1e-300);
  r[kTpRgc] = vgc / std::max(r[kTpSgap] * r[kTpSdct], 1e-300);
  r[kTpRzc] = vzc / std::max(r[kTpSdz] * r[kTpSdct], 1e-300);
  r[kTpMee] = std::sqrt(std::max(E * E - dot<3>(P, P), 0.));
  r[kTpPx] = P[0];
  r[kTpPy] = P[1];
  r[kTpPz] = P[2];
  return r;
}

// CF result vector (ConvCfCol) and the 3x3 vertex covariance, started at the TP vector tp;
// ok = 0 with NaN entries when undefined or failed
inline void convPairCf(const Vec5& p1, const Mat5& c1, const Vec5& p2, const Mat5& c2, const ConvTp& tp, double K,
                       ConvCf& res, Mat3& Cx) {
  using namespace detail;
  res.fill(JV_NAN);
  Cx.setConstant(JV_NAN);
  res[kCfOk] = 0.;
  if (!(tp[kTpOk] > 0.)) return;
  CfData D;
  for (int i = 0; i < 10; ++i)
    for (int j = 0; j < 10; ++j) D.W[i][j] = 0.0;
  double a1[5][5], a2[5][5], w1[5][5], w2[5][5];
  toArr5(c1, a1);
  toArr5(c2, a2);
  if (nbInvA<5>(a1, w1) != kLinOk) return;
  if (nbInvA<5>(a2, w2) != kLinOk) return;
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j) {
      D.W[i][j] = w1[i][j];
      D.W[5 + i][5 + j] = w2[i][j];
    }
  for (int k = 0; k < 5; ++k) {
    D.m[k] = p1[k];
    D.m[5 + k] = p2[k];
  }
  A7 th = {tp[kTpX], tp[kTpY], tp[kTpZ], tp[kTpPhi], tp[kTpCt], p1[2], p2[2]};
  D.Lr[0] = tp[kTpL1];
  D.Lr[1] = tp[kTpL2];
  A10 r;
  cfRes(th, D, r);
  double chi2 = cfQuad(r, D.W);
  if (!std::isfinite(chi2)) return;
  int ntot = 0;
  bool neg = false;
  for (int e = 0; e < JV_CF_NESCAPE + 1; ++e) {
    ntot += cfLm(th, chi2, D);
    A7 tn;
    double cn;
    neg = cfEscape(th, chi2, D, tn, cn);
    if (!(neg && cn < chi2)) break;
    th = tn;
    chi2 = cn;
  }
  cfRes(th, D, r);
  chi2 = cfQuad(r, D.W);
  double J[10][7], H[7][7], g[7], Cth[7][7];
  cfJac(th, D, J);
  cfNormal(J, D.W, r.data(), H, g);
  const double dec = cfDecrement(H, g);
  if (nbInvA<7>(H, Cth) != kLinOk) return;
  bool fin = std::isfinite(chi2);
  for (int i = 0; i < 7; ++i) {
    fin = fin && std::isfinite(th[i]);
    for (int j = 0; j < 7; ++j) fin = fin && std::isfinite(Cth[i][j]);
  }
  if (!fin) return;
  const double pt1 = K / std::max(std::fabs(th[5]), 1e-12), pt2 = K / std::max(std::fabs(th[6]), 1e-12);
  const double pt = pt1 + pt2;
  res[kCfOk] = (dec < JV_CF_TOL && !neg) ? 1. : 2.;
  for (int k = 0; k < 7; ++k) res[1 + k] = th[k];
  res[kCfChi2] = chi2;
  res[kCfNiter] = ntot;
  res[kCfPx] = pt * std::cos(th[3]);
  res[kCfPy] = pt * std::sin(th[3]);
  res[kCfPz] = pt * th[4];
  double y[5];
  vm<5, 5>(r.data(), w1, y);
  res[kCfChi21] = dot<5>(y, r.data());
  vm<5, 5>(r.data() + 5, w2, y);
  res[kCfChi22] = dot<5>(y, r.data() + 5);
  res[kCfDec] = dec;
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b) Cx(a, b) = Cth[a][b];
}

// prefilter of an opposite-curvature pair on its TP; o receives the raw TP (x, y, z1, z2, gap, dct, L1, L2)
inline bool convPrefilter(const Vec5& p1, const Vec5& p2, const ConvPre& pre, double* o) {
  if (!detail::tpRaw(p1.data(), p2.data(), o)) return false;
  if (std::fabs(o[4]) > pre.gap || std::fabs(o[2] - o[3]) > pre.dz || std::fabs(o[5]) > pre.dct) return false;
  if (o[6] < pre.arc || o[7] < pre.arc) return false;
  if (std::hypot(o[0], o[1]) > pre.r) return false;
  return true;
}

// quality row (ConvQCol) and CF vertex covariance of one pair; rmax1, rmax2 = hit-radius bounds of the tracks
inline void convQuality(const Vec5& p1, const Mat5& c1, const Vec5& p2, const Mat5& c2, const Vec3& pv, double rmax1,
                        double rmax2, double K, ConvQ& out, Mat3& Cx) {
  const ConvTp tp = convPairTp(p1, c1, p2, c2, K);
  ConvCf cf;
  convPairCf(p1, c1, p2, c2, tp, K, cf, Cx);
  for (int k = 0; k < JV_CONV_NTP; ++k) out[k] = tp[k];
  for (int k = 0; k < JV_CONV_NCF; ++k) out[kQCf + k] = cf[k];
  const bool tpok = tp[kTpOk] > 0, cfok = cf[kCfOk] > 0;
  out[kQTpR] = tpok ? std::hypot(tp[kTpX], tp[kTpY]) : JV_NAN;
  out[kQCfR] = cfok ? std::hypot(cf[kCfX], cf[kCfY]) : JV_NAN;
  // TP conversion chi2 of (gap, dz, dct) with their correlations
  if (tpok) {
    const double v[3] = {tp[kTpGap] / npMaximum(tp[kTpSgap], 1e-12), tp[kTpDz] / npMaximum(tp[kTpSdz], 1e-12),
                         tp[kTpDct] / npMaximum(tp[kTpSdct], 1e-12)};
    const double R[3][3] = {{1., tp[kTpRgz], tp[kTpRgc]}, {tp[kTpRgz], 1., tp[kTpRzc]}, {tp[kTpRgc], tp[kTpRzc], 1.}};
    double s[3];
    out[kQTpChi2] = detail::npSolveA<3>(R, v, s) == kLinOk ? detail::dot<3>(v, s) : JV_NAN;
  } else {
    out[kQTpChi2] = JV_NAN;
  }
  // photon pointing: signed transverse distance of the photon line from the PV, cos of the 3D angle
  for (int which = 0; which < 2; ++which) {
    bool ok;
    double x, y, z, px, py, pz;
    if (which == 0) {
      ok = tpok;
      x = tp[kTpX]; y = tp[kTpY]; z = tp[kTpZ]; px = tp[kTpPx]; py = tp[kTpPy]; pz = tp[kTpPz];
    } else {
      ok = cfok;
      x = cf[kCfX]; y = cf[kCfY]; z = cf[kCfZ]; px = cf[kCfPx]; py = cf[kCfPy]; pz = cf[kCfPz];
    }
    double dg = JV_NAN, cp = JV_NAN;
    if (ok) {
      const double dx = x - pv[0], dy = y - pv[1], dz = z - pv[2];
      const double pt = std::max(std::hypot(px, py), 1e-12);
      dg = (dx * py - dy * px) / pt;
      cp = (dx * px + dy * py + dz * pz) /
           std::max(std::sqrt(dx * dx + dy * dy + dz * dz) * std::sqrt(px * px + py * py + pz * pz), 1e-300);
    }
    out[kQTpDgam + 2 * which] = dg;
    out[kQTpCosp + 2 * which] = cp;
  }
  // CF vertex error along the photon momentum
  if (cfok) {
    double u[3] = {cf[kCfPx], cf[kCfPy], cf[kCfPz]};
    const double nu = std::max(std::sqrt(detail::dot<3>(u, u)), 1e-300);
    for (int k = 0; k < 3; ++k) u[k] = u[k] / nu;
    double c[3][3], y[3];
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b) c[a][b] = Cx(a, b);
    detail::vm<3, 3>(u, c, y);
    out[kQCfSl] = std::sqrt(std::max(detail::dot<3>(y, u), 0.));
  } else {
    out[kQCfSl] = JV_NAN;
  }
  // inner hits impossible for the TP radius
  if (tpok) {
    out[kQNhx1] = out[kQTpR] > rmax1 ? 1. : 0.;
    out[kQNhx2] = out[kQTpR] > rmax2 ? 1. : 0.;
  } else {
    out[kQNhx1] = JV_NAN;
    out[kQNhx2] = JV_NAN;
  }
  double P[3];
  detail::momAt(p1.data(), tp[kTpL1], K, P);
  out[kQP1] = std::sqrt(detail::dot<3>(P, P));
  detail::momAt(p2.data(), tp[kTpL2], K, P);
  out[kQP2] = std::sqrt(detail::dot<3>(P, P));
}

// opposite-charge pairs (i < j) of usable tracks (ok) passing the prefilter
inline void convFindPairs(const std::vector<Vec5>& par, const std::vector<int>& q, const std::vector<char>& ok,
                          const ConvPre& pre, std::vector<int>& I, std::vector<int>& J) {
  I.clear();
  J.clear();
  const int n = static_cast<int>(par.size());
  double o[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  for (int i = 0; i < n; ++i) {
    if (!ok[i]) continue;
    for (int j = i + 1; j < n; ++j) {
      if (!ok[j] || q[i] * q[j] >= 0) continue;
      if (convPrefilter(par[i], par[j], pre, o)) {
        I.push_back(i);
        J.push_back(j);
      }
    }
  }
}

// usable tracks: finite parameters and covariance
inline std::vector<char> convUsable(const std::vector<Vec5>& par, const std::vector<Mat5>& cov) {
  std::vector<char> ok(par.size());
  for (size_t i = 0; i < par.size(); ++i) ok[i] = allFinite(par[i]) && allFinite(cov[i]);
  return ok;
}

// quality rows and CF vertex covariances of the pairs (I, J) with one PV
inline void convPairQuality(const std::vector<Vec5>& par, const std::vector<Mat5>& cov, const std::vector<int>& I,
                            const std::vector<int>& J, const Vec3& pv, const std::vector<double>& rmax, double K,
                            std::vector<ConvQ>& Q, std::vector<Mat3>& CX) {
  const size_t n = I.size();
  Q.assign(n, ConvQ{});
  CX.assign(n, Mat3::Zero());
  for (size_t k = 0; k < n; ++k) {
    Q[k].fill(JV_NAN);
    const int i = I[k], j = J[k];
    convQuality(par[i], cov[i], par[j], cov[j], pv, rmax[i], rmax[j], K, Q[k], CX[k]);
  }
}

// selection: TP and CF defined and every cut passed
inline bool convSelect(const ConvQ& q, const std::vector<ConvCut>& sel) {
  if (!(q[kTpOk] > 0) || !(q[kQCf + kCfOk] > 0)) return false;
  for (const ConvCut& c : sel) {
    double v = q[c.col];
    if (c.col == kQTpDgam || c.col == kQCfDgam || c.col == kTpGap || c.col == kTpDz || c.col == kTpDct)
      v = std::fabs(v);
    if (c.hasLo && !(v >= c.lo)) return false;
    if (c.hasHi && !(v < c.hi)) return false;
  }
  return true;
}

// selected conversion candidates of one event
struct ConvFound {
  std::vector<int> i, j;     // track indices (i < j)
  std::vector<Vec3> x;       // CF vertex
  std::vector<Mat3> cov;     // CF vertex covariance
  std::vector<ConvQ> Q;      // quality rows
};

// pairs passing the prefilter and the selection sel (default: loose); rmax = hit-radius bound per track
inline ConvFound convFind(const std::vector<Vec5>& par, const std::vector<Mat5>& cov, const std::vector<int>& q,
                          const std::vector<double>& rmax, const Vec3& pv, double K = JV_PT_HALF,
                          const std::vector<ConvCut>& sel = convLoose(), const ConvPre& pre = ConvPre()) {
  std::vector<int> I, J;
  convFindPairs(par, q, convUsable(par, cov), pre, I, J);
  std::vector<ConvQ> Q;
  std::vector<Mat3> CX;
  convPairQuality(par, cov, I, J, pv, rmax, K, Q, CX);
  ConvFound f;
  for (size_t k = 0; k < I.size(); ++k) {
    if (!convSelect(Q[k], sel)) continue;
    f.i.push_back(I[k]);
    f.j.push_back(J[k]);
    f.x.emplace_back(Q[k][kQCf + kCfX], Q[k][kQCf + kCfY], Q[k][kQCf + kCfZ]);
    f.cov.push_back(CX[k]);
    f.Q.push_back(Q[k]);
  }
  return f;
}

// conversion candidate of the V0 class (a = positive track, b = negative track)
struct ConvCand {
  int a, b;
  Vec3 x;          // CF vertex
  Mat3 C;          // its covariance
  double chi2a, chi2b, chi2cf;   // CF chi2 per track and total
  double zside;    // TP chi2 of (dz, dct) with their correlation
  double m;        // TP m_ee
  Vec3 P;          // CF photon momentum
  double pa, pb, p;
  bool viola, violb;   // CF radius >= hit-radius bound of a / b
  double cosp;     // CF cos pointing to the PV
  double r;        // transverse distance from the PV
  double rv;       // radius from the axis
};

// loose candidates among the tracks with ok set (indices into par / cov / q / rmax), CF radius from the axis >= rmin
inline std::vector<ConvCand> convRaw(const std::vector<Vec5>& par, const std::vector<Mat5>& cov,
                                     const std::vector<int>& q, const std::vector<char>& ok,
                                     const std::vector<double>& rmax, const Vec3& pvx, double K, double rmin) {
  std::vector<ConvCand> out;
  std::vector<int> idx;
  for (size_t k = 0; k < ok.size(); ++k)
    if (ok[k]) idx.push_back(static_cast<int>(k));
  if (idx.size() < 2) return out;
  std::vector<Vec5> sp;
  std::vector<Mat5> sc;
  std::vector<int> sq;
  std::vector<double> sr;
  for (int k : idx) {
    sp.push_back(par[k]);
    sc.push_back(cov[k]);
    sq.push_back(q[k]);
    sr.push_back(rmax[k]);
  }
  const ConvFound R = convFind(sp, sc, sq, sr, pvx, K);
  for (size_t n = 0; n < R.i.size(); ++n) {
    const int i = idx[R.i[n]], j = idx[R.j[n]];
    const bool ai = q[i] > 0;
    const ConvQ& Q = R.Q[n];
    ConvCand c;
    c.a = ai ? i : j;
    c.b = ai ? j : i;
    c.x = R.x[n];
    c.C = R.cov[n];
    c.chi2cf = Q[kQCf + kCfChi2];
    const double c1 = Q[kQCf + kCfChi21], c2 = Q[kQCf + kCfChi22];
    c.chi2a = ai ? c1 : c2;
    c.chi2b = ai ? c2 : c1;
    const double sdz = Q[kTpSdz], sdc = Q[kTpSdct], rzc = Q[kTpRzc];
    const double u = Q[kTpDz] / std::max(sdz, 1e-12), v = Q[kTpDct] / std::max(sdc, 1e-12);
    c.zside = (u * u - 2 * rzc * u * v + v * v) / std::max(1. - rzc * rzc, 1e-6);
    c.m = Q[kTpMee];
    c.P = Vec3(Q[kQCf + kCfPx], Q[kQCf + kCfPy], Q[kQCf + kCfPz]);
    c.rv = std::hypot(c.x[0], c.x[1]);
    c.pa = ai ? Q[kQP1] : Q[kQP2];
    c.pb = ai ? Q[kQP2] : Q[kQP1];
    c.p = std::sqrt(c.P[0] * c.P[0] + c.P[1] * c.P[1] + c.P[2] * c.P[2]);
    c.viola = c.rv >= rmax[c.a];
    c.violb = c.rv >= rmax[c.b];
    c.cosp = Q[kQCfCosp];
    c.r = std::hypot(c.x[0] - pvx[0], c.x[1] - pvx[1]);
    if (c.rv >= rmin) out.push_back(c);
  }
  return out;
}

}  // namespace AlephJointVtx
}  // namespace FCCAnalyses

#endif
