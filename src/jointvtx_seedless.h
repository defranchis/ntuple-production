#ifndef JOINTVTX_SEEDLESS_H
#define JOINTVTX_SEEDLESS_H

/*
  Vertex-density maxima that supply secondary-vertex candidates (AlephJointVtx).

  Each track i gives a tube f_i(r) = exp(-chi2_i(r) / (2 tube)), chi2_i = track chi2 at the point r
  (lin1); the primary vertex gives f_0(r) = k0 exp(-(r - PV)^T C^-1 (r - PV) / 2). The density is
  V(r) = sum f - sum f^2 / sum f. Linearised pair crossings are the grid of starting points; every
  start climbs V by a fixed-point iteration (step limited to 3 mm). A maximum is kept when V and the
  number of member tracks are large enough and it is no duplicate of an earlier one; the dip is the
  lowest V on the straight segment to the primary vertex over the lower of the two end values.

  vdMaxima: the maxima, in the order in which the climbs end.
  vdSupply: the maxima in order of decreasing V, resolved from the primary vertex (dip below R0), each
  with the tracks it may claim; the candidate fit that follows is the caller's.
  Failure statuses (singular or non-finite linear algebra) stand for a raise of the reference, except one:
  for a non-finite primary-vertex covariance the reference's (NumPy) inverse in vdSupply returns without
  raising (non-finite values; finite ones for an infinite diagonal element), while vdSupply fails there.
*/

#include <algorithm>
#include <cmath>
#include <vector>

#include <Eigen/Dense>

#include "jointvtx_common.h"
#include "jointvtx_track.h"

namespace FCCAnalyses {
namespace AlephJointVtx {

// constants of the maxima search that are not settings
constexpr double JV_VD_GRID_R2MAX = 3.5 * 3.5;   // cm^2, largest squared distance of a grid point from the primary vertex
constexpr double JV_VD_STEP_MAX = 0.3;           // cm, largest step of a climb
constexpr double JV_VD_STEP_TOL = 2e-4;          // cm, a climb ends below this step
constexpr double JV_VD_DUP_DIST = 5e-4;          // cm, maxima closer than this are the same
constexpr double JV_VD_WEIGHT_MIN = 1e-8;        // tracks with a smaller climb weight do not enter the step
constexpr double JV_VD_REG = 1e-6;               // cm^-2 added to the diagonal of the climb matrix
constexpr double JV_VD_SUM_MIN = 1e-100;         // density sums below this count as empty
constexpr double JV_VD_DIP_FRAC[5] = {0.15, 0.3, 0.5, 0.7, 0.85};   // points of the dip segment

struct VdMaxima {
  std::vector<Vec3> x;                     // positions
  std::vector<double> v;                   // density at the maxima
  std::vector<double> dip;                 // dip ratio towards the primary vertex
  std::vector<std::vector<double>> chi2;   // chi2 of every track at each maximum
  int nclimb = 0;                          // climbs started
  int ngrid = 0;                           // grid points that passed the cuts
  int nties = 0;                           // equal neighbouring grid scores after sorting
  LinStatus status = kLinOk;
};

namespace detail {

inline double vdValueSum(double s1, double s2) { return s1 > JV_VD_SUM_MIN ? s1 - s2 / s1 : 0.; }

}  // namespace detail

// density from the track chi2 values and the primary vertex term f0
inline double vdValue(const std::vector<double>& chi2, double tube, double f0) {
  double s1 = f0, s2 = f0 * f0;
  for (size_t k = 0; k < chi2.size(); ++k) {
    const double f = std::exp(-0.5 * chi2[k] / tube);
    s1 += f;
    s2 += f * f;
  }
  return detail::vdValueSum(s1, s2);
}

// pvCi: inverse of the primary-vertex covariance in C memory order
inline LinStatus vdMaxima(const std::vector<Trk>& t, const std::vector<double>& rmax, const Vec3& pvx,
                          const Mat3& pvCi, double k0, double tube, double gridChi2, double vmin, int nClimb,
                          int maxStarts, double fwdMargin, double memberChi2, double gridPvChi2, double gridCos,
                          VdMaxima& out) {
  out = VdMaxima();
  const int nt = static_cast<int>(t.size());
  std::vector<double> Ls0(nt, 0.), c0(nt, 0.);
  std::vector<Vec3> X0(nt, Vec3::Zero()), a0(nt, Vec3::Zero());
  std::vector<Mat3> D0(nt, Mat3::Zero());
  Vec3 a, Xt;
  Mat35 A;
  auto fail = [&out](LinStatus s) {
    out.status = s;
    return s;
  };
  for (int k = 0; k < nt; ++k) {
    const Lin1 r = lin1(t[k].p, t[k].C, pvx, 0., X0[k], D0[k]);
    if (r.status != kLinOk) return fail(r.status);
    c0[k] = r.chi2;
    Ls0[k] = r.L;
    helix1(t[k].p, Ls0[k], Xt, a, A);
    a0[k] = a / std::sqrt(dotSeq(a, a));
  }
  // search grid: linearised pair crossings
  const Mat3 eps = Mat3::Identity();
  std::vector<int> si, sj;
  std::vector<Vec3> sx;
  std::vector<double> sc;
  for (int i = 0; i < nt; ++i)
    for (int j = i + 1; j < nt; ++j) {
      if (std::max(c0[i], c0[j]) < gridPvChi2 || dotSeq(a0[i], a0[j]) < gridCos) continue;
      const Mat3 H = D0[i] + D0[j] + eps;
      const Vec3 b = matVecC<3, 3>(D0[i], X0[i] - pvx) + matVecC<3, 3>(D0[j], X0[j] - pvx);
      Vec3 xs;
      const LinStatus st = npSolve<3>(H, b, xs);
      if (st != kLinOk) return fail(st);
      const Vec3 x = xs + pvx;
      const Vec3 ri = x - X0[i], rj = x - X0[j];
      const double ch = detail::quadForm(ri, D0[i], false) + detail::quadForm(rj, D0[j], false);
      if (ch > gridChi2) continue;
      const Vec3 d = x - pvx;
      if (dotSeq(d, d) > JV_VD_GRID_R2MAX) continue;
      if (std::hypot(x[0], x[1]) >= std::min(rmax[i], rmax[j])) continue;
      if (dotSeq(a0[i], d) < -fwdMargin || dotSeq(a0[j], d) < -fwdMargin) continue;
      double s1 = 0., s2 = 0.;
      for (int k = 0; k < nt; ++k) {
        const Vec3 r = x - X0[k];
        const double f = std::exp(-0.5 * detail::quadForm(r, D0[k], false) / tube);
        s1 += f;
        s2 += f * f;
      }
      si.push_back(i);
      sj.push_back(j);
      sx.push_back(x);
      sc.push_back(detail::vdValueSum(s1, s2));
    }
  const int ns = static_cast<int>(sc.size());
  out.ngrid = ns;
  std::vector<double> negsc(ns);
  for (int s = 0; s < ns; ++s) negsc[s] = -sc[s];
  const std::vector<int> order = stableArgsort(negsc);
  for (int s = 1; s < ns; ++s)
    if (negsc[order[s]] == negsc[order[s - 1]]) ++out.nties;

  std::vector<char> covered(static_cast<size_t>(nt) * nt, 0);
  std::vector<Vec3> mx;
  std::vector<double> mv;
  std::vector<std::vector<double>> mc;
  int nclimb = 0;
  std::vector<double> Ls(nt, 0.), ch(nt, 0.), f(nt, 0.);
  std::vector<Vec3> X(nt, Vec3::Zero());
  std::vector<Mat3> D(nt, Mat3::Zero());
  for (int o = 0; o < ns; ++o) {
    const int s_ = order[o];
    const int i = si[s_], j = sj[s_];
    if (covered[static_cast<size_t>(i) * nt + j]) continue;
    if (nclimb >= maxStarts) break;
    ++nclimb;
    Vec3 x = sx[s_];
    Ls = Ls0;
    for (int it = 0; it < nClimb; ++it) {
      double s1 = 0., s2 = 0.;
      Vec3 d = x - pvx;
      double f0 = k0 * std::exp(-0.5 * detail::quadForm(d, pvCi, false));
      for (int k = 0; k < nt; ++k) {
        const Lin1 r = lin1(t[k].p, t[k].C, x, Ls[k], X[k], D[k]);
        if (r.status != kLinOk) return fail(r.status);
        ch[k] = r.chi2;
        Ls[k] = r.L;
        f[k] = std::exp(-0.5 * ch[k] / tube);
        s1 += f[k];
        s2 += f[k] * f[k];
      }
      s1 += f0;
      s2 += f0 * f0;
      if (s1 <= JV_VD_SUM_MIN) break;
      Mat3 H = Mat3::Zero();
      Vec3 g = Vec3::Zero();
      for (int k = 0; k < nt; ++k) {
        const double p = f[k] / s1;
        const double w = f[k] * (1. - 2. * p + s2 / (s1 * s1)) / tube;
        if (w > JV_VD_WEIGHT_MIN) {
          H += w * D[k];
          g += w * matVecC<3, 3>(D[k], X[k] - pvx);
        }
      }
      const double p0 = f0 / s1;
      const double w0 = f0 * (1. - 2. * p0 + s2 / (s1 * s1));
      H += w0 * pvCi;
      H += JV_VD_REG * eps;
      Vec3 sol;
      const LinStatus st = npSolve<3>(H, g, sol);
      if (st != kLinOk) return fail(st);
      Vec3 xn = sol + pvx;
      if (!allFinite(xn)) break;
      const Vec3 dx = xn - x;
      const double nd = std::sqrt(dotSeq(dx, dx));
      if (nd > JV_VD_STEP_MAX) xn = x + dx * JV_VD_STEP_MAX / nd;
      x = xn;
      if (nd < JV_VD_STEP_TOL) break;
    }
    Vec3 d = x - pvx;
    double f0 = k0 * std::exp(-0.5 * detail::quadForm(d, pvCi, false));
    for (int k = 0; k < nt; ++k) {
      const Lin1 r = lin1(t[k].p, t[k].C, x, Ls[k], X[k], D[k]);
      if (r.status != kLinOk) return fail(r.status);
      ch[k] = r.chi2;
      Ls[k] = r.L;
    }
    const double v = vdValue(ch, tube, f0);
    std::vector<int> mem;
    for (int k = 0; k < nt; ++k)
      if (ch[k] < memberChi2) mem.push_back(k);
    for (int a1 : mem)
      for (int a2 : mem) covered[static_cast<size_t>(a1) * nt + a2] = 1;
    covered[static_cast<size_t>(i) * nt + j] = 1;
    if (v < vmin || mem.size() < 2) continue;
    bool dup = false;
    for (size_t m = 0; m < mx.size(); ++m) {
      const Vec3 dm = x - mx[m];
      if (std::sqrt(dotSeq(dm, dm)) < JV_VD_DUP_DIST) {
        dup = true;
        break;
      }
    }
    if (dup) continue;
    mx.push_back(x);
    mv.push_back(v);
    mc.push_back(ch);
    if (static_cast<int>(mx.size()) == maxStarts) break;
  }
  // value at the primary vertex and dip along the straight segment to it
  const double vpv = vdValue(c0, tube, k0);
  std::vector<double> dip(mx.size(), 0.);
  for (size_t m = 0; m < mx.size(); ++m) {
    const double lo = std::min(mv[m], vpv);
    double vm = mv[m];
    for (int s = 0; s < 5; ++s) {
      const Vec3 xs = mx[m] + JV_VD_DIP_FRAC[s] * (pvx - mx[m]);
      const Vec3 d = xs - pvx;
      const double f0 = k0 * std::exp(-0.5 * detail::quadForm(d, pvCi, false));
      for (int k = 0; k < nt; ++k) {
        const Lin1 r = lin1(t[k].p, t[k].C, xs, Ls0[k], X[k], D[k]);
        if (r.status != kLinOk) return fail(r.status);
        ch[k] = r.chi2;
      }
      vm = std::min(vm, vdValue(ch, tube, f0));
    }
    dip[m] = lo > JV_VD_SUM_MIN ? vm / lo : 1.;
  }
  out.x = mx;
  out.v = mv;
  out.dip = dip;
  out.chi2 = mc;
  out.nclimb = nclimb;
  return kLinOk;
}

// settings read by the supply
struct VdParams {
  double k0 = 0., tube = 0., gridChi2 = 0., vmin = 0.;
  int nClimb = 0, maxStarts = 0;
  double fwdMargin = 0., memberChi2 = 0., gridPvChi2 = 0., gridCos = 0.;
  double R0 = 0.;        // a maximum is resolved from the primary vertex if its dip is below this
  double slChi2 = 0.;    // chi2 to the maximum below which a track may claim it
  bool slSymmetric = false;   // also require chi2 to the maximum below the chi2 to the primary vertex
  bool slForward = false;     // tracks must fly away from the primary vertex
};

// a resolved maximum with the tracks (indices into the event track list) that may claim it
struct VdStart {
  Vec3 x = Vec3::Zero();
  std::vector<int> trk;
  double v = 0., dip = 0.;
};

struct VdSupplyInfo {
  int nmax = 0;          // maxima found
  int nties = 0;         // equal neighbouring maxima values after sorting
  int ngridTies = 0;     // equal grid scores (see VdMaxima)
  bool maxima = false;   // the maxima and the track chi2 to the primary vertex were computed
};

// tracks at arc L along the helix move away from the primary vertex: direction . (x - pv) >= 0
inline bool fliesForward(const Vec5& p, double L, const Vec3& d) {
  const double ang = p[1] + 2 * p[2] * L;
  return std::cos(ang) * d[0] + std::sin(ang) * d[1] + p[4] * d[2] >= 0;
}

// idx: the event tracks that enter (the selection mask as indices), rmax: bound on the production
// radius of every event track, pvC: primary-vertex covariance. Starts come in order of decreasing density; on a
// failure out may hold the starts found before it, which the caller discards (the event fails). maxima
// (optional): receives the vertex-density maxima.
inline LinStatus vdSupply(const std::vector<Trk>& t, const std::vector<int>& idx, const std::vector<double>& rmax,
                          const Vec3& pvx, const Mat3& pvC, const VdParams& P, std::vector<VdStart>& out,
                          VdSupplyInfo& info, VdMaxima* maxima = nullptr) {
  out.clear();
  info = VdSupplyInfo();
  const int n = static_cast<int>(idx.size());
  if (n < 2) return kLinOk;
  std::vector<Trk> sub(n);
  std::vector<double> rm(n);
  for (int k = 0; k < n; ++k) {
    sub[k] = t[idx[k]];
    rm[k] = rmax[idx[k]];
  }
  Mat3 Ci;
  LinStatus st = npInv<3>(pvC, Ci);
  if (st != kLinOk) return st;
  VdMaxima Mloc;
  VdMaxima& M = maxima ? *maxima : Mloc;
  st = vdMaxima(sub, rm, pvx, Ci, P.k0, P.tube, P.gridChi2, P.vmin, P.nClimb, P.maxStarts, P.fwdMargin,
                P.memberChi2, P.gridPvChi2, P.gridCos, M);
  info.ngridTies = M.nties;
  if (st != kLinOk) return st;
  std::vector<double> pvchi2, Lpv;
  st = chi2ToPoint(sub, pvx, pvchi2, Lpv);
  if (st != kLinOk) return st;
  info.maxima = true;
  const int nm = static_cast<int>(M.v.size());
  info.nmax = nm;
  std::vector<double> negv(nm);
  for (int m = 0; m < nm; ++m) negv[m] = -M.v[m];
  const std::vector<int> ord = stableArgsort(negv);
  for (int o = 1; o < nm; ++o)
    if (negv[ord[o]] == negv[ord[o - 1]]) ++info.nties;
  for (int o = 0; o < nm; ++o) {
    const int m = ord[o];
    if (M.dip[m] >= P.R0) continue;
    const Vec3& x = M.x[m];
    std::vector<double> ch, L;
    st = chi2ToPoint(sub, x, ch, L);
    if (st != kLinOk) return st;
    const Vec3 d = x - pvx;
    const double rx = std::hypot(x[0], x[1]);
    VdStart s;
    s.x = x;
    s.v = M.v[m];
    s.dip = M.dip[m];
    for (int k = 0; k < n; ++k) {
      bool ok = ch[k] < P.slChi2 && rx < rm[k];
      if (P.slSymmetric) ok = ok && ch[k] < pvchi2[k];
      if (P.slForward) ok = ok && fliesForward(sub[k].p, L[k], d);
      if (ok) s.trk.push_back(idx[k]);
    }
    if (s.trk.size() < 2) continue;
    out.push_back(s);
  }
  return kLinOk;
}

}  // namespace AlephJointVtx
}  // namespace FCCAnalyses

#endif
