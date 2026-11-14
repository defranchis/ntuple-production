#ifndef JOINTVTX_V0KIN_H
#define JOINTVTX_V0KIN_H

/*
  Kinematics of the V0 pair hypotheses of the joint vertex fitter (AlephJointVtx):
  K0S -> pi+ pi-, Lambda -> p pi-, anti-Lambda -> anti-p pi+.

  v0PairKernel: for every opposite-charge pair of the listed tracks, the xy
  crossings of the two track circles inside rmax with compatible z and a mass
  at the crossing inside a pre-window of one species start a two-track vertex
  fit (fitVertexNb, JV_CAND_FIT_NITER undamped steps); the lowest finite chi2
  over the crossings is kept whether or not its fit converged. A pair passing
  chi2 < chi2max and the radius bound gives one row per hypothesis whose mass
  at the fitted vertex is inside the wide window and whose Armenteros qT
  (positive daughter relative to the V0 momentum) is above the species' cut.
  Momenta at fixed arcs; the mass error propagates the two track covariances
  (no vertex error). The positive track is always `a`.
  v0PointingQ: 2D significance^2 of the miss distance of the V0 flight line
  from a vertex; v0FlightSig: 3D flight significance from a vertex.
  v0RawCandidates: the kernel rows of the selected tracks of an event with
  their static quantities relative to the reference primary vertex.
  Units: cm, GeV.
  Operation order follows the reference implementation (products written as
  loops in the measured BLAS order). No function throws.
*/

#include <algorithm>
#include <cmath>
#include <vector>

#include <Eigen/Dense>

#include "jointvtx_common.h"
#include "jointvtx_track.h"
#include "jointvtx_fit.h"

namespace FCCAnalyses {
namespace AlephJointVtx {

// selection of the pair kernel (cm, GeV)
struct V0PairCuts {
  double rmax;     // crossing and vertex radius bound
  double dzSig;    // z compatibility at the crossing, standard deviations
  double dzAbs;    // z compatibility at the crossing, cm
  double preKs;    // mass pre-window at the crossing, K0S
  double preLam;   // mass pre-window at the crossing, Lambda (both assignments)
  double winKs;    // mass window at the fitted vertex, K0S
  double winLam;   // mass window at the fitted vertex, Lambda
  double chi2max;  // pair-fit chi2 bound
  double qtKs;     // Armenteros qT lower bound, K0S
  double qtLam;    // Armenteros qT lower bound, Lambda
};

// one hypothesis of one pair; a = positive track, b = negative track (indices into the track list)
struct V0PairRow {
  int a = -1, b = -1, hyp = 0;
  Vec3 x = Vec3::Zero();
  Mat3 C = Mat3::Zero();
  double chi2a = 0., chi2b = 0.;
  double m = 0., sm = 0.;
  Vec3 P = Vec3::Zero();
  Mat3 CP = Mat3::Zero();
  double pa = 0., pb = 0.;
  bool fitConv = false;  // convergence flag of the kept fit (information only)
};

struct V0MassErr {
  double m, sm;
};

namespace detail {

// circle of the track: centre and radius
inline void v0Circle(const Vec5& p, double& xc, double& yc, double& r) {
  const double D = p[0], phi0 = p[1], C = p[2];
  xc = -(D + 0.5 / C) * std::sin(phi0);
  yc = (D + 0.5 / C) * std::cos(phi0);
  r = 0.5 / std::fabs(C);
}

}  // namespace detail

// xy crossings of the two track circles (disjoint / nested circles: the midpoint of the gap); number of points
inline int v0Crossings(const Vec5& p1, const Vec5& p2, double out[2][2]) {
  if (p1[2] == 0 || p2[2] == 0) return 0;
  double x0, y0, r0, x1, y1, r1;
  detail::v0Circle(p1, x0, y0, r0);
  detail::v0Circle(p2, x1, y1, r1);
  const double dx = x1 - x0, dy = y1 - y0;
  const double d = std::sqrt(dx * dx + dy * dy);
  if (!std::isfinite(d) || !(d > 0)) return 0;
  const double ux = dx / d, uy = dy / d;
  if (d >= r0 + r1) {
    const double t = r0 + 0.5 * (d - r0 - r1);
    out[0][0] = x0 + t * ux;
    out[0][1] = y0 + t * uy;
    return 1;
  }
  if (d <= std::fabs(r0 - r1)) {
    const double s = r0 > r1 ? 1. : -1.;
    out[0][0] = 0.5 * (x0 + x1 + s * (r0 + r1) * ux);
    out[0][1] = 0.5 * (y0 + y1 + s * (r0 + r1) * uy);
    return 1;
  }
  const double a = (r0 * r0 - r1 * r1 + d * d) / (2. * d);
  const double h2 = r0 * r0 - a * a;
  const double X = x0 + a * ux, Y = y0 + a * uy;
  if (h2 > 0) {
    const double h = std::sqrt(h2);
    out[0][0] = X - h * uy;
    out[0][1] = Y + h * ux;
    out[1][0] = X + h * uy;
    out[1][1] = Y - h * ux;
    return 2;
  }
  out[0][0] = X;
  out[0][1] = Y;
  return 1;
}

// transverse arc of the circle point nearest to (x, y); 0 for a straight track
inline double v0Arc(const Vec5& p, double x, double y) {
  if (p[2] == 0) return 0.;
  return arcToXY(p, x, y);
}

// momentum at arc L (fixed arc) and its Jacobian d mom / d par; K = momentum constant / 2 (pT = K / |C|)
inline void v0MomJac(const Vec5& p, double L, double K, Vec3& mom, Mat35& J) {
  const double C = p[2], ct = p[4];
  const double pt = K / std::max(std::fabs(C), 1e-12);
  const double ang = p[1] + 2 * C * L;
  const double ca = std::cos(ang), sa = std::sin(ang);
  mom[0] = pt * ca;
  mom[1] = pt * sa;
  mom[2] = pt * ct;
  J.setZero();
  J(0, 1) = -pt * sa;
  J(1, 1) = pt * ca;
  const double Cs = C != 0 ? C : 1e-12;
  J(0, 2) = -pt / Cs * ca - 2 * L * pt * sa;
  J(1, 2) = -pt / Cs * sa + 2 * L * pt * ca;
  J(2, 2) = -pt / Cs * ct;
  J(2, 4) = pt;
}

// invariant mass of the momenta pa, pb with masses ma, mb and its error from the track covariances ca, cb through
// the momentum Jacobians Ja, Jb; (0, 1) if m^2 <= 0
inline V0MassErr v0MassErr(const Vec3& pa, const Mat35& Ja, const Mat5& ca, const Vec3& pb, const Mat35& Jb,
                           const Mat5& cb, double ma, double mb) {
  const double Ea = std::sqrt(dotSeq(pa, pa) + ma * ma), Eb = std::sqrt(dotSeq(pb, pb) + mb * mb);
  const Vec3 P = pa + pb;
  const double E = Ea + Eb;
  const double m2 = E * E - dotSeq(P, P);
  if (!(m2 > 0)) return {0., 1.};
  const double m = std::sqrt(m2);
  Vec3 ga, gb;
  for (int k = 0; k < 3; ++k) {
    ga[k] = (E * pa[k] / Ea - P[k]) / m;
    gb[k] = (E * pb[k] / Eb - P[k]) / m;
  }
  const Vec5 va = detail::vecMat35(ga, Ja);
  const Vec5 vb = detail::vecMat35(gb, Jb);
  const double s2 = detail::quadForm5C(va, ca) + detail::quadForm5C(vb, cb);
  return {m, std::sqrt(std::max(s2, 1e-18))};
}

// V0 pair rows of the tracks idx (indices into t); q = charges; K = momentum constant / 2. Rows are appended in the
// order pairs (u < v over idx), then hypotheses K0S, Lambda, anti-Lambda. failed (optional): set when a pair fit fails
inline void v0PairKernel(const std::vector<Trk>& t, const std::vector<int>& q, const std::vector<int>& idx, double K,
                         const V0PairCuts& cut, std::vector<V0PairRow>& out, bool* failed = nullptr) {
  const size_t n = idx.size();
  double cr[2][2] = {{0., 0.}, {0., 0.}};
  Vec3 ma, mb;
  Mat35 Ja, Jb;
  const Mat3 Binv = Mat3::Zero();
  const Vec3 c0 = Vec3::Zero();
  std::vector<int> two(2);
  const std::vector<double> w2(2, 1.0);
  std::vector<double> Ls(t.size(), 0.0);
  const double r2max = cut.rmax * cut.rmax;
  for (size_t u = 0; u < n; ++u) {
    for (size_t v = u + 1; v < n; ++v) {
      const int i0 = idx[u], j0 = idx[v];
      if (q[i0] * q[j0] >= 0) continue;
      const int i = q[i0] > 0 ? i0 : j0, j = q[i0] > 0 ? j0 : i0;
      const Vec5& pi = t[i].p;
      const Vec5& pj = t[j].p;
      const Mat5& ci = t[i].C;
      const Mat5& cj = t[j].C;
      const int nc = v0Crossings(pi, pj, cr);
      double best = 1e30;
      Vec3 bx = Vec3::Zero();
      Mat3 bC = Mat3::Zero();
      double bch0 = 0., bch1 = 0.;
      bool bconv = false;
      for (int k = 0; k < nc; ++k) {
        const double cx = cr[k][0], cy = cr[k][1];
        if (cx * cx + cy * cy >= r2max) continue;
        const double Li = v0Arc(pi, cx, cy), Lj = v0Arc(pj, cx, cy);
        const double zi = pi[3] + pi[4] * Li, zj = pj[3] + pj[4] * Lj;
        const double vi = ci(3, 3) + 2 * Li * ci(3, 4) + Li * Li * ci(4, 4);
        const double vj = cj(3, 3) + 2 * Lj * cj(3, 4) + Lj * Lj * cj(4, 4);
        const double dz = std::fabs(zi - zj);
        if (dz > cut.dzAbs && dz > cut.dzSig * std::sqrt(std::max(vi + vj, 1e-18))) continue;
        v0MomJac(pi, Li, K, ma, Ja);
        v0MomJac(pj, Lj, K, mb, Jb);
        const double Ea = std::sqrt(dotSeq(ma, ma) + JV_M_PI * JV_M_PI);
        const double Eb = std::sqrt(dotSeq(mb, mb) + JV_M_PI * JV_M_PI);
        const Vec3 P = ma + mb;
        const double PP = dotSeq(P, P);
        const double mks = std::sqrt(std::max((Ea + Eb) * (Ea + Eb) - PP, 0.));
        const double Eap = std::sqrt(dotSeq(ma, ma) + JV_M_P * JV_M_P);
        const double Ebp = std::sqrt(dotSeq(mb, mb) + JV_M_P * JV_M_P);
        const double ml1 = std::sqrt(std::max((Eap + Eb) * (Eap + Eb) - PP, 0.));
        const double ml2 = std::sqrt(std::max((Ea + Ebp) * (Ea + Ebp) - PP, 0.));
        if (std::fabs(mks - JV_M_KS) > cut.preKs && std::fabs(ml1 - JV_M_LAM) > cut.preLam &&
            std::fabs(ml2 - JV_M_LAM) > cut.preLam)
          continue;
        Ls[i] = Li;
        Ls[j] = Lj;
        two[0] = i;
        two[1] = j;
        const Vec3 x0(cx, cy, 0.5 * (zi + zj));
        const VtxFit f = fitVertexNb(t, two, w2, x0, Binv, c0, Ls, JV_CAND_FIT_NITER, JV_CAND_FIT_TOL);
        if (!f.ok() || JOINTVTX_FAIL_AT(kFsV0PairFit)) {
          if (failed) *failed = true;
          continue;
        }
        const double tot = f.chi2[0] + f.chi2[1];
        if (!(std::isfinite(tot) && allFinite(f.x))) continue;
        if (tot < best) {
          best = tot;
          bx = f.x;
          bC = f.C;
          bch0 = f.chi2[0];
          bch1 = f.chi2[1];
          bconv = f.conv;
        }
      }
      if (!(best < cut.chi2max) || bx[0] * bx[0] + bx[1] * bx[1] >= r2max) continue;
      const double Li = v0Arc(pi, bx[0], bx[1]), Lj = v0Arc(pj, bx[0], bx[1]);
      v0MomJac(pi, Li, K, ma, Ja);
      v0MomJac(pj, Lj, K, mb, Jb);
      const Vec3 P = ma + mb;
      const double pP = std::sqrt(dotSeq(P, P));
      double qt2 = 0.;
      if (pP > 0) {
        const double pl = dotSeq(ma, P) / pP;
        qt2 = dotSeq(ma, ma) - pl * pl;
      }
      const Mat3 CP = detail::sandwich(Ja, ci) + detail::sandwich(Jb, cj);
      for (int hyp = kV0K0S; hyp <= kV0AntiLambda; ++hyp) {
        const double m1 = hyp != kV0Lambda ? JV_M_PI : JV_M_P;
        const double m2 = hyp != kV0AntiLambda ? JV_M_PI : JV_M_P;
        const V0MassErr me = v0MassErr(ma, Ja, ci, mb, Jb, cj, m1, m2);
        if (hyp == kV0K0S && (std::fabs(me.m - JV_M_KS) > cut.winKs || qt2 < cut.qtKs * cut.qtKs)) continue;
        if (hyp > kV0K0S && (std::fabs(me.m - JV_M_LAM) > cut.winLam || qt2 < cut.qtLam * cut.qtLam)) continue;
        V0PairRow r;
        r.a = i;
        r.b = j;
        r.hyp = hyp;
        r.x = bx;
        r.C = bC;
        r.chi2a = bch0;
        r.chi2b = bch1;
        r.m = me.m;
        r.sm = me.sm;
        r.P = P;
        r.CP = CP;
        r.pa = std::sqrt(dotSeq(ma, ma));
        r.pb = std::sqrt(dotSeq(mb, mb));
        r.fitConv = bconv;
        out.push_back(r);
      }
    }
  }
}

// 2D significance^2 of the miss distance of the flight line (x, direction P) from the vertex (pvx, pvC), with the
// direction error CP scaled by the decay length; 0 if |P| = 0 or the 2x2 system is singular
inline double v0PointingQ(const Vec3& x, const Mat3& Cx, const Vec3& P, const Mat3& CP, const Vec3& pvx,
                          const Mat3& pvC) {
  const double pm = std::sqrt(dotSeq(P, P));
  if (!(pm > 0)) return 0.;
  Vec3 u, d, mv;
  for (int k = 0; k < 3; ++k) u[k] = P[k] / pm;
  for (int k = 0; k < 3; ++k) d[k] = x[k] - pvx[k];
  const double L = dotSeq(d, u);
  for (int k = 0; k < 3; ++k) mv[k] = d[k] - L * u[k];
  Vec3 e1 = std::fabs(u[2]) < 0.9 ? detail::npCross(u, Vec3(0., 0., 1.)) : detail::npCross(u, Vec3(1., 0., 0.));
  const double n1 = std::sqrt(dotSeq(e1, e1));
  for (int k = 0; k < 3; ++k) e1[k] /= n1;
  const Vec3 e2 = detail::npCross(u, e1);
  const double s = (L * L) / (pm * pm);
  double M[3][3], B[2][3], T[2][3];
  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 3; ++c) M[r][c] = (Cx(r, c) + pvC(r, c)) + s * CP(r, c);
  for (int c = 0; c < 3; ++c) {
    B[0][c] = e1[c];
    B[1][c] = e2[c];
  }
  for (int r = 0; r < 2; ++r)
    for (int c = 0; c < 3; ++c) T[r][c] = (B[r][0] * M[0][c] + B[r][1] * M[1][c]) + B[r][2] * M[2][c];
  Mat2 Cm;
  for (int r = 0; r < 2; ++r)
    for (int c = 0; c < 2; ++c) Cm(r, c) = (T[r][0] * B[c][0] + T[r][1] * B[c][1]) + T[r][2] * B[c][2];
  Vec2 mm;
  for (int r = 0; r < 2; ++r) mm[r] = (B[r][0] * mv[0] + B[r][2] * mv[2]) + B[r][1] * mv[1];
  Vec2 sol;
  if (npSolve<2>(Cm, mm, sol) == kLinSingular) return 0.;
  return std::max(dotSeq(mm, sol), 0.);
}

// 3D flight significance of x from the vertex (pvx, pvC), covariance C + pvC; 0 if singular
inline double v0FlightSig(const Vec3& x, const Mat3& C, const Vec3& pvx, const Mat3& pvC) {
  Vec3 d;
  for (int k = 0; k < 3; ++k) d[k] = x[k] - pvx[k];
  Mat3 S;
  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 3; ++c) S(r, c) = C(r, c) + pvC(r, c);
  Vec3 sol;
  if (npSolve<3>(S, d, sol) == kLinSingular) return 0.;
  return std::sqrt(std::max(dotSeq(d, sol), 0.));
}

// K0S / Lambda candidate of the event with its static quantities (Lambda: both assignments kept)
struct V0Raw {
  V0PairRow k;            // kernel row (track indices into the event's track list)
  double p = 0.;          // |P|
  bool viola = false;     // vertex radius beyond the hit-radius bound of track a
  bool violb = false;
  double r = 0.;          // transverse distance from the reference PV
  double qp = 0.;         // pointing significance^2 to the reference PV
  double fsig = 0.;       // 3D flight significance from the reference PV
  double rv = 0.;         // vertex radius from the detector axis
};

// kernel candidates of the tracks with ok != 0; rmax = hit-radius bound per track; (pvx, pvC) = reference PV;
// failed (optional): set when a pair fit fails
inline std::vector<V0Raw> v0RawCandidates(const std::vector<Trk>& t, const std::vector<int>& q,
                                          const std::vector<char>& ok, const std::vector<double>& rmax,
                                          const Vec3& pvx, const Mat3& pvC, double K, const V0PairCuts& cut,
                                          bool* failed = nullptr) {
  std::vector<V0Raw> out;
  std::vector<int> idx;
  for (size_t k = 0; k < ok.size(); ++k)
    if (ok[k]) idx.push_back(static_cast<int>(k));
  if (idx.size() < 2) return out;
  std::vector<V0PairRow> rows;
  v0PairKernel(t, q, idx, K, cut, rows, failed);
  out.reserve(rows.size());
  for (const V0PairRow& w : rows) {
    V0Raw c;
    c.k = w;
    c.rv = std::hypot(w.x[0], w.x[1]);
    c.p = std::sqrt(dotSeq(w.P, w.P));
    c.viola = c.rv >= rmax[w.a];
    c.violb = c.rv >= rmax[w.b];
    c.r = std::hypot(w.x[0] - pvx[0], w.x[1] - pvx[1]);
    c.qp = v0PointingQ(w.x, w.C, w.P, w.CP, pvx, pvC);
    c.fsig = v0FlightSig(w.x, w.C, pvx, pvC);
    out.push_back(c);
  }
  return out;
}

}  // namespace AlephJointVtx
}  // namespace FCCAnalyses

#endif
