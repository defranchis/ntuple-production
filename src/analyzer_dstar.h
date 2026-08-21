#ifndef ALEPHDSTAR_H
#define ALEPHDSTAR_H

/*
  D*+ -> D0(K- pi+) pi+_slow reconstruction, a kinematically tagged KAON
  source: dm = m(K pi pi_s) - m(K pi) has its resolution set by the slow pion
  alone, so a narrow dm window plus the D0 mass window isolates a sample in
  which the track given the kaon mass really is a kaon. Output = a dE/dx
  calibration sample, so NO dE/dx quantity enters any selection; the
  daughters' dE/dx is stored, never cut on.

  Two output collections: the stand-alone D0 -> K pi list (one entry per track
  pair AND mass assignment - the kaon hypothesis defines the tag, so the two
  assignments are different objects) and the D* list (one entry per D0
  candidate x third track, carrying d0idx back into the D0 list). Both slow-pion
  charges are kept and flagged (rs); wrong-sign has no signal and measures the
  combinatorial background. No three-track fit: the D0 flies ~0.6 mm while the
  slow pion comes from the PV region, so the D* uses the D0 momenta at the
  fitted D0 vertex and the slow pion's PERIGEE momentum. Promptness is not
  required (D* from b decays are displaced): dpv and its significance are
  stored, only a wide fiducial bounds |vtx-PV|, and cosPoint enters the D*
  tight label alone.

  Input = flipD0_copy'ed trackstates (physical charge = +sign(omega)); lengths
  cm, momenta GeV; the single VertexFitter_Tk call is made with
  rescale_cm_mm=false and its momenta rescaled once by 10, as in AlephV0New.
  All selection values are arguments; the defaults are PROVISIONAL and loose.
*/

#include <algorithm>
#include <array>
#include <cmath>
#include <unordered_map>
#include <vector>

#include <ROOT/RVec.hxx>
#include "TVector3.h"
#include "TLorentzVector.h"

#include "edm4hep/TrackState.h"
#include "edm4hep/MCParticleData.h"
#include "FCCAnalyses/VertexingUtils.h"
#include "FCCAnalyses/VertexFitterSimple.h"

#include "aleph_units.h"

namespace FCCAnalyses {
namespace AlephDstar {

using ROOT::VecOps::RVec;

// PDG 2024 central values [GeV]. The charged-pion constant cannot be called
// M_PI: that name is a <cmath> macro.
constexpr double M_K     = 0.493677;    // charged kaon
constexpr double M_PICH  = 0.13957039;  // charged pion
constexpr double M_D0    = 1.86484;
constexpr double M_DSTAR = 2.01026;
constexpr double DM_NOMINAL = 0.145426; // m(D*+) - m(D0)
constexpr int PDG_DSTAR = 413, PDG_D0 = 421, PDG_K = 321, PDG_PI = 211;

constexpr double E_BEAM = 45.6;   // beam energy [GeV], for xE = E/E_BEAM
constexpr double PRE_MARGIN = 0.05;  // pre-fit mass window margin [GeV]

// Perigee momentum from the helix parameters (pT = kPtPerTeslaCm*Bz/|omega|).
// Used for the pre-fit mass window and for the slow pion, which is not part
// of any fit.
inline TVector3 perigeeMomentum(const edm4hep::TrackState& t, double Bz) {
  double om = std::abs(t.omega);
  if (om <= 0.) return TVector3(0., 0., 0.);
  double pt = AlephUnits::kPtPerTeslaCm * Bz / om;
  return TVector3(pt * std::cos(t.phi), pt * std::sin(t.phi), pt * t.tanLambda);
}

inline TLorentzVector lorentz(const TVector3& p, double m) {
  return TLorentzVector(p, std::sqrt(p.Mag2() + m * m));
}

// Helicity angle: cosine of the kaon direction in the D0 rest frame w.r.t.
// the D0 flight direction in the lab. Flat for a true two-body scalar decay,
// peaked at |cos| = 1 for random combinatorics.
inline double cosThetaStarK(const TVector3& pk, const TVector3& ppi) {
  TLorentzVector k = lorentz(pk, M_K);
  TLorentzVector d = k + lorentz(ppi, M_PICH);
  if (!(d.Vect().Mag() > 0.)) return -99.;
  TLorentzVector kr = k;
  kr.Boost(-d.BoostVector());
  double kk = kr.Vect().Mag();
  if (!(kk > 0.)) return -99.;
  return kr.Vect().Dot(d.Vect().Unit()) / kk;
}

// 3D compatibility significance of two vertex positions: sqrt of the chi2 of
// d = x1 - x2 under the summed position covariances (packed lower triangle
// xx, yx, yy, zx, zy, zz). -1 when the summed covariance is not invertible.
template <typename CovA, typename CovB>
inline float vertexDistSig(const TVector3& d, const CovA& ca, const CovB& cb) {
  double C[3][3] = {
    {double(ca[0]) + cb[0], double(ca[1]) + cb[1], double(ca[3]) + cb[3]},
    {double(ca[1]) + cb[1], double(ca[2]) + cb[2], double(ca[4]) + cb[4]},
    {double(ca[3]) + cb[3], double(ca[4]) + cb[4], double(ca[5]) + cb[5]}};
  double det = C[0][0] * (C[1][1] * C[2][2] - C[1][2] * C[2][1])
             - C[0][1] * (C[1][0] * C[2][2] - C[1][2] * C[2][0])
             + C[0][2] * (C[1][0] * C[2][1] - C[1][1] * C[2][0]);
  if (!(std::abs(det) > 0.) || !std::isfinite(det)) return -1.;
  double inv[3][3];
  inv[0][0] = (C[1][1] * C[2][2] - C[1][2] * C[2][1]) / det;
  inv[0][1] = (C[0][2] * C[2][1] - C[0][1] * C[2][2]) / det;
  inv[0][2] = (C[0][1] * C[1][2] - C[0][2] * C[1][1]) / det;
  inv[1][0] = (C[1][2] * C[2][0] - C[1][0] * C[2][2]) / det;
  inv[1][1] = (C[0][0] * C[2][2] - C[0][2] * C[2][0]) / det;
  inv[1][2] = (C[0][2] * C[1][0] - C[0][0] * C[1][2]) / det;
  inv[2][0] = (C[1][0] * C[2][1] - C[1][1] * C[2][0]) / det;
  inv[2][1] = (C[0][1] * C[2][0] - C[0][0] * C[2][1]) / det;
  inv[2][2] = (C[0][0] * C[1][1] - C[0][1] * C[1][0]) / det;
  double dv[3] = {d.X(), d.Y(), d.Z()};
  double s2 = 0.;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) s2 += dv[i] * inv[i][j] * dv[j];
  return (s2 > 0. && std::isfinite(s2)) ? std::sqrt(s2) : -1.;
}

// ---------------------------------------------------------------------------
// Per-daughter block, one instance per role (K, pi, slow pi). `p`/`costheta`
// are the momentum actually used for that daughter: at the fitted vertex for
// the D0 daughters, at the perigee for the slow pion.
// ---------------------------------------------------------------------------
struct TrkBlock {
  RVec<int>   origIdx;   // index into the original Tracks collection
  RVec<int>   q;         // physical charge, +sign(omega)
  RVec<float> p, costheta;
  RVec<float> d0, z0;    // perigee parameters [cm]
  RVec<float> sigd0;     // sqrt(cov[0]) [cm]
  RVec<int>   nvdet, nitc;
  RVec<float> chi2ndf;   // track fit chi2/ndf
  RVec<int>   isprim;    // 1 = in the fitted primary set (stored, never cut on)
  RVec<int>   pool;      // 0 = primary set, 1 = secondary set, 2 = neither
};

// D0 -> K pi candidates, one entry per (pair, mass assignment).
struct D0Block {
  RVec<float> m_kpi;
  RVec<float> p, px, py, pz, costheta, xE;
  RVec<float> chi2;            // vertex fit chi2, normalised (ndf = 1)
  RVec<float> vx, vy, vz;      // fitted vertex [cm]
  RVec<float> dpv, dpvSig;     // |vtx-PV| [cm] and its 3D significance
  RVec<float> cosPoint;        // cos(angle) between p(D0) and (vtx - PV)
  RVec<float> cosThetaStar;    // kaon helicity angle
  RVec<int>   loose, tight;    // labels, not cuts
  RVec<int>   stage;           // cascade stage 1-3 (0 = cascade off)
  RVec<int>   nsec;            // legs in the secondary pool (0-2)
  TrkBlock    trkK, trkPi;
};

// D* -> D0 pi_slow candidates, one entry per (D0 candidate, third track).
// The D0 quantities are repeated; d0idx points back to the D0 list.
struct DstarBlock {
  RVec<float> m_kpi, dm;
  RVec<float> p, px, py, pz, costheta, xE;
  RVec<float> chi2;
  RVec<float> vx, vy, vz;
  RVec<float> dpv, dpvSig;
  RVec<float> cosPoint, cosThetaStar;
  RVec<int>   rs;              // 1 = right-sign slow pion (charge of the pi)
  RVec<int>   loose, tight;
  RVec<int>   d0idx;           // index into the D0 list, -1 if absent from it
  RVec<int>   stage;           // cascade stage 1-6 (0 = cascade off)
  RVec<int>   nsec;            // legs in the secondary pool (0-3)
  TrkBlock    trkK, trkPi, trkPis;
};

struct DstarCands {
  D0Block    d0;
  DstarBlock ds;
  int        nfits = 0;   // vertex fits actually performed (cache misses)
};

// ---------------------------------------------------------------------------
// The finder.
//   tracks     : flipD0_copy'ed baseline-selected trackstates
//   orig_idx   : original-Tracks index of each entry of `tracks`
//   nvdet/nitc : VDET/ITC hit counts of each entry (-1 = unknown)
//   chi2ndf    : track fit chi2/ndf of each entry (-1 = unknown)
//   isprim     : 1 if the entry is in the fitted primary set (stored, never cut on)
//   pool       : 0 = primary set, 1 = secondary set, 2 = neither (staging class)
//   PV         : fitted primary vertex (positions numerically cm)
//   veto_orig  : original-Tracks indices excluded from the pool (opt-in),
//                empty by default
// The track pool is the FULL baseline-selected list, primary and secondary
// alike, with no masking by the PV split. A cut argument <= 0 disables that
// cut (except the mass window and the fiducial, which are the storage
// definition).
//
// cascade = false is the single all-track pass; true builds the candidates in
// ordered exclusive stages (see below). claim_mode: 0 none, 1 tight, 2 loose.
// ---------------------------------------------------------------------------
inline DstarCands findDstar(
    const RVec<edm4hep::TrackState>& tracks,
    const RVec<int>& orig_idx,
    const RVec<int>& nvdet,
    const RVec<int>& nitc,
    const RVec<float>& chi2ndf,
    const RVec<int>& isprim,
    const RVec<int>& pool,
    const VertexingUtils::FCCAnalysesVertex& PV,
    double solenoidBz = 1.5,
    double m_lo = 1.70, double m_hi = 2.03,  // stored K pi mass window [GeV]
    double chi2_cut = 25.,                   // vertex chi2 (ndf=1); <=0 = off
    double dpv_fid = 10.,                    // |vtx-PV| storage fiducial [cm]; <=0 = off
    double dm_max = 0.20,                    // stored dm ceiling [GeV]
    double p_min = 0.3,                      // K/pi momentum floor [GeV]; <=0 = off
    double ps_min = 0.1,                     // slow-pion momentum floor [GeV]; <=0 = off
    double sigd0_max = -1.,                  // per-track sigma(d0) cap [cm]; <=0 = off
    int    min_vdet_itc = 0,                 // per-track min (nVDET + nITC)
    double trk_chi2ndf_max = -1.,            // per-track chi2/ndf cap; <=0 = off
    double d0_loose_dm = 0.060,              // D0-alone loose label [GeV]
    double d0_tight_dm = 0.030,              // D0-alone tight label [GeV]
    double d0_tight_dpvsig = 3.0,            // D0-alone tight: |vtx-PV| significance
    double d0_tight_cospoint = 0.99,         // D0-alone tight: cosPoint
    double d0_tight_cosstar = 0.8,           // D0-alone tight: |cosThetaStar|
    double ds_loose_dm = 0.060,              // D* loose label: |m - m_D0| [GeV]
    double ds_loose_ddm = 0.0040,            // D* loose label: |dm - dm_nom| [GeV]
    double ds_tight_dm = 0.030,              // D* tight label: |m - m_D0| [GeV]
    double ds_tight_ddm = 0.0015,            // D* tight label: |dm - dm_nom| [GeV]
    double tight_pk = 1.0,                   // tight labels: kaon |p| [GeV]
    double tight_ppi = 1.0,                  // tight labels: pion |p| [GeV]
    double tight_chi2 = 10.,                 // tight labels: vertex chi2
    double tight_ps = -1.,                   // D* tight label: slow-pion |p| [GeV]; <=0 = off
    double tight_cospoint = -1.,             // D* tight label: D0 cosPoint; <=-1 = off
    double pre_margin = PRE_MARGIN,          // pre-fit mass window margin [GeV]
    bool   cascade = false,                  // staged exclusive mode
    int    claim_mode = 1,                   // 0 none, 1 tight, 2 loose
    const RVec<int>& veto_orig = {}) {       // original-Tracks indices to exclude

  DstarCands out;
  const int nTr = tracks.size();
  if (nTr < 2) return out;

  auto nv = [&](const RVec<int>& src, int k) {
    return (k >= 0 && k < (int)src.size()) ? src[k] : -1;
  };
  auto c2 = [&](int k) {
    return (k >= 0 && k < (int)chi2ndf.size()) ? chi2ndf[k] : -1.f;
  };
  auto oidx = [&](int k) {
    return (k >= 0 && k < (int)orig_idx.size()) ? orig_idx[k] : -1;
  };
  auto pushTrk = [&](TrkBlock& b, int k, const TVector3& p) {
    b.origIdx.push_back(oidx(k));
    b.q.push_back((tracks[k].omega > 0) ? 1 : -1);
    b.p.push_back(p.Mag());
    b.costheta.push_back(p.Mag() > 0. ? p.Z() / p.Mag() : -99.);
    b.d0.push_back(tracks[k].D0);
    b.z0.push_back(tracks[k].Z0);
    b.sigd0.push_back(tracks[k].covMatrix[0] > 0. ? std::sqrt(tracks[k].covMatrix[0]) : -1.);
    b.nvdet.push_back(nv(nvdet, k));
    b.nitc.push_back(nv(nitc, k));
    b.chi2ndf.push_back(c2(k));
    b.isprim.push_back(nv(isprim, k));
    b.pool.push_back(nv(pool, k));
  };
  // staging class: SEC = in the secondary set, everything else counts as
  // primary-like (the "neither" tracks are rare and are staged with the
  // primaries, while their true class stays in the stored pool branch)
  auto isSec = [&](int k) { return k < (int)pool.size() && pool[k] == 1; };

  // track-quality prefilter; the momentum floors only decide which ROLES a
  // survivor may take, so a track too soft to be a kaon can still be a slow pion
  const double p_floor = (p_min > 0. && ps_min > 0.) ? std::min(p_min, ps_min)
                                                     : std::max(p_min, ps_min);
  std::vector<int> good;
  good.reserve(nTr);
  std::vector<TVector3> pper(nTr);
  std::vector<char> canKPi(nTr, 0), canPis(nTr, 0);
  for (int i = 0; i < nTr; ++i) {
    const auto& t = tracks[i];
    if (!veto_orig.empty()) {
      int o = oidx(i);
      if (o >= 0 && std::find(veto_orig.begin(), veto_orig.end(), o) != veto_orig.end())
        continue;
    }
    if (sigd0_max > 0.) {
      double c0 = t.covMatrix[0];
      if (!(c0 > 0.) || std::sqrt(c0) > sigd0_max) continue;
    }
    if (min_vdet_itc > 0) {
      int nvd = nv(nvdet, i), nit = nv(nitc, i);
      if (nvd < 0 || nit < 0 || nvd + nit < min_vdet_itc) continue;
    }
    if (trk_chi2ndf_max > 0.) {
      float c = c2(i);
      if (!(c >= 0.) || c > trk_chi2ndf_max) continue;
    }
    pper[i] = perigeeMomentum(t, solenoidBz);
    double pm = pper[i].Mag();
    if (p_floor > 0. && pm < p_floor) continue;
    canKPi[i] = (p_min <= 0. || pm >= p_min) ? 1 : 0;
    canPis[i] = (ps_min <= 0. || pm >= ps_min) ? 1 : 0;
    if (!canKPi[i] && !canPis[i]) continue;
    good.push_back(i);
  }
  if (good.size() < 2) return out;

  TVector3 pv(PV.vertex.position[0], PV.vertex.position[1], PV.vertex.position[2]);
  const double pre_lo = m_lo - pre_margin, pre_hi = m_hi + pre_margin;

  // memoised two-track fit: the same pair is reached by several cascade
  // stages (different slow-pion pools) but is fitted exactly once
  struct FitRes {
    bool ok;
    double chi2, dis;
    float dsig;
    TVector3 x, d, pa, pb;
  };
  std::unordered_map<long long, FitRes> fitcache;
  RVec<edm4hep::TrackState> tr_pair(2);
  // i < j is guaranteed by the callers; pa/pb follow that order
  auto fitPair = [&](int i, int j) -> const FitRes& {
    long long key = (long long)i * nTr + j;
    auto it = fitcache.find(key);
    if (it != fitcache.end()) return it->second;
    FitRes r;
    r.ok = false;
    r.chi2 = 0.; r.dis = 0.; r.dsig = -1.;
    tr_pair[0] = tracks[i];
    tr_pair[1] = tracks[j];
    auto v = VertexFitterSimple::VertexFitter_Tk(
        0, tr_pair, tracks, false, 0., 0., 0., 0., 0., 0., solenoidBz, false);
    ++out.nfits;
    if (v.updated_track_momentum_at_vertex.size() == 2) {
      // cm-as-mm homothety: rescale ONCE at the source (see AlephV0New)
      for (auto& tp : v.updated_track_momentum_at_vertex) tp *= 10.;
      double chi2 = v.vertex.chi2;
      bool pass = (chi2 == chi2) && !(chi2_cut > 0. && chi2 >= chi2_cut);
      if (pass) {
        TVector3 x(v.vertex.position[0], v.vertex.position[1], v.vertex.position[2]);
        TVector3 dvec = x - pv;
        double dis = dvec.Mag();
        if (!(dpv_fid > 0. && dis > dpv_fid)) {
          r.ok = true;
          r.chi2 = chi2; r.x = x; r.d = dvec; r.dis = dis;
          r.dsig = vertexDistSig(dvec, v.vertex.covMatrix, PV.vertex.covMatrix);
          r.pa = v.updated_track_momentum_at_vertex[0];
          r.pb = v.updated_track_momentum_at_vertex[1];
        }
      }
    }
    return fitcache.emplace(key, r).first->second;
  };

  // one accepted (pair, mass assignment) before it is written out
  struct D0Cand {
    int iK, iPi;
    double m, chi2, dis, cosp, cstar, energy;
    float dsig;
    TVector3 p3, x, pK, pPi;
    int loose, tight;
  };

  // D0 candidates whose (K, pi) legs match the staging pattern
  // (-1 = any, 0 = both secondary, 1 = exactly one secondary, 2 = none)
  // and are not yet claimed by an earlier stage.
  auto makeD0 = [&](int kpi_pat, const std::vector<char>& claimed) {
    std::vector<D0Cand> res;
    for (size_t a = 0; a + 1 < good.size(); ++a) {
      const int i = good[a];
      if (!canKPi[i] || claimed[i]) continue;
      for (size_t b = a + 1; b < good.size(); ++b) {
        const int j = good[b];
        if (!canKPi[j] || claimed[j]) continue;
        if (tracks[i].omega * tracks[j].omega > 0) continue;  // opposite charge only
        if (kpi_pat >= 0) {
          int ns = (isSec(i) ? 1 : 0) + (isSec(j) ? 1 : 0);
          if (kpi_pat == 0 && ns != 2) continue;
          if (kpi_pat == 1 && ns != 1) continue;
          if (kpi_pat == 2 && ns != 0) continue;
        }

        // pre-fit window on the perigee momenta: both mass assignments are
        // tried, and the pair is fitted only if at least one is inside
        double m_pre_ij = AlephV0New::invMass(pper[i], M_K, pper[j], M_PICH);
        double m_pre_ji = AlephV0New::invMass(pper[i], M_PICH, pper[j], M_K);
        bool ok_ij = (m_pre_ij >= pre_lo && m_pre_ij <= pre_hi);
        bool ok_ji = (m_pre_ji >= pre_lo && m_pre_ji <= pre_hi);
        if (!ok_ij && !ok_ji) continue;

        const FitRes& fr = fitPair(i, j);
        if (!fr.ok) continue;

        for (int asg = 0; asg < 2; ++asg) {
          if (asg == 0 && !ok_ij) continue;
          if (asg == 1 && !ok_ji) continue;
          D0Cand c;
          c.iK  = (asg == 0) ? i : j;
          c.iPi = (asg == 0) ? j : i;
          c.pK  = (asg == 0) ? fr.pa : fr.pb;
          c.pPi = (asg == 0) ? fr.pb : fr.pa;

          c.m = AlephV0New::invMass(c.pK, M_K, c.pPi, M_PICH);
          if (c.m < m_lo || c.m > m_hi) continue;

          TLorentzVector d0lv = lorentz(c.pK, M_K) + lorentz(c.pPi, M_PICH);
          c.p3 = d0lv.Vect();
          double pmag = c.p3.Mag();
          if (!(pmag > 0.)) continue;
          c.energy = d0lv.E();
          c.chi2 = fr.chi2;
          c.x = fr.x;
          c.dis = fr.dis;
          c.dsig = fr.dsig;
          c.cosp = (fr.dis > 0.) ? c.p3.Dot(fr.d) / (pmag * fr.dis) : -99.;
          c.cstar = cosThetaStarK(c.pK, c.pPi);

          // labels evaluated on the STORED post-fit quantities
          const double dmass = std::abs(c.m - M_D0);
          c.loose = (d0_loose_dm <= 0. || dmass < d0_loose_dm) ? 1 : 0;
          c.tight = (dmass < d0_tight_dm && c.dsig > d0_tight_dpvsig &&
                     c.cosp > d0_tight_cospoint &&
                     std::abs(c.cstar) < d0_tight_cosstar &&
                     c.pK.Mag() > tight_pk && c.pPi.Mag() > tight_ppi &&
                     c.chi2 < tight_chi2) ? 1 : 0;
          res.push_back(c);
        }
      }
    }
    return res;
  };

  auto storeD0 = [&](const D0Cand& c, int stage) {
    out.d0.m_kpi.push_back(c.m);
    out.d0.p.push_back(c.p3.Mag());
    out.d0.px.push_back(c.p3.X());
    out.d0.py.push_back(c.p3.Y());
    out.d0.pz.push_back(c.p3.Z());
    out.d0.costheta.push_back(c.p3.Z() / c.p3.Mag());
    out.d0.xE.push_back(c.energy / E_BEAM);
    out.d0.chi2.push_back(c.chi2);
    out.d0.vx.push_back(c.x.X());
    out.d0.vy.push_back(c.x.Y());
    out.d0.vz.push_back(c.x.Z());
    out.d0.dpv.push_back(c.dis);
    out.d0.dpvSig.push_back(c.dsig);
    out.d0.cosPoint.push_back(c.cosp);
    out.d0.cosThetaStar.push_back(c.cstar);
    out.d0.loose.push_back(c.loose);
    out.d0.tight.push_back(c.tight);
    out.d0.stage.push_back(stage);
    out.d0.nsec.push_back((isSec(c.iK) ? 1 : 0) + (isSec(c.iPi) ? 1 : 0));
    pushTrk(out.d0.trkK, c.iK, c.pK);
    pushTrk(out.d0.trkPi, c.iPi, c.pPi);
  };

  // D* from a D0 candidate plus the slow-pion track k. No three-track fit:
  // the D0 momenta are the at-vertex ones, the slow pion enters at its perigee.
  auto pushDstar = [&](const D0Cand& c, int k, int stage, int d0idx) {
    // the D0 mass enters dm as the STORED (float) value
    const float mf = c.m;
    const double m = mf;
    TLorentzVector d0lv = lorentz(c.pK, M_K) + lorentz(c.pPi, M_PICH);
    TLorentzVector dslv = d0lv + lorentz(pper[k], M_PICH);
    double dm = dslv.M() - m;
    if (!(dm > 0.)) return false;
    if (dm_max > 0. && dm >= dm_max) return false;
    TVector3 p3 = dslv.Vect();
    double pmag = p3.Mag();
    if (!(pmag > 0.)) return false;

    const float chi2f = c.chi2;
    const int qPi = (tracks[c.iPi].omega > 0) ? 1 : -1;
    const int qs = (tracks[k].omega > 0) ? 1 : -1;

    out.ds.m_kpi.push_back(mf);
    out.ds.dm.push_back(dm);
    out.ds.p.push_back(pmag);
    out.ds.px.push_back(p3.X());
    out.ds.py.push_back(p3.Y());
    out.ds.pz.push_back(p3.Z());
    out.ds.costheta.push_back(p3.Z() / pmag);
    out.ds.xE.push_back(dslv.E() / E_BEAM);
    out.ds.chi2.push_back(chi2f);
    out.ds.vx.push_back(c.x.X());
    out.ds.vy.push_back(c.x.Y());
    out.ds.vz.push_back(c.x.Z());
    out.ds.dpv.push_back(c.dis);
    out.ds.dpvSig.push_back(c.dsig);
    out.ds.cosPoint.push_back(c.cosp);
    out.ds.cosThetaStar.push_back(c.cstar);
    out.ds.rs.push_back((qs == qPi) ? 1 : 0);
    out.ds.d0idx.push_back(d0idx);
    out.ds.stage.push_back(stage);
    out.ds.nsec.push_back((isSec(c.iK) ? 1 : 0) + (isSec(c.iPi) ? 1 : 0) +
                          (isSec(k) ? 1 : 0));
    pushTrk(out.ds.trkK, c.iK, c.pK);
    pushTrk(out.ds.trkPi, c.iPi, c.pPi);
    pushTrk(out.ds.trkPis, k, pper[k]);

    const double dmass = std::abs(m - M_D0);
    const double ddm = std::abs(dm - DM_NOMINAL);
    out.ds.loose.push_back((dmass < ds_loose_dm && ddm < ds_loose_ddm) ? 1 : 0);
    out.ds.tight.push_back((dmass < ds_tight_dm && ddm < ds_tight_ddm &&
                            c.pK.Mag() > tight_pk && c.pPi.Mag() > tight_ppi &&
                            chi2f < tight_chi2 &&
                            (tight_ps <= 0. || pper[k].Mag() > tight_ps) &&
                            (tight_cospoint <= -1. || c.cosp > tight_cospoint)) ? 1 : 0);
    return true;
  };

  const std::vector<char> nomask(nTr, 0);
  const bool use_cascade = cascade && !pool.empty();

  if (!use_cascade) {
    // single all-track pass: every D0 candidate, then every third track
    auto cands = makeD0(-1, nomask);
    for (const auto& c : cands) storeD0(c, 0);
    for (size_t c = 0; c < cands.size(); ++c)
      for (size_t s = 0; s < good.size(); ++s) {
        const int k = good[s];
        if (k == cands[c].iK || k == cands[c].iPi) continue;
        if (!canPis[k]) continue;
        pushDstar(cands[c], k, 0, (int)c);
      }
    return out;
  }

  // ordered exclusive stages, most-displaced pattern first; after each stage
  // the legs of the claimed right-sign candidates leave the pool:
  //   1 (sec,sec,sec)  2 (sec,sec,prim)  3 (one D0 leg prim, pis sec)
  //   4 (one D0 leg prim, pis prim)  5 (prim,prim,sec)  6 (prim,prim,prim)
  static const int KPI_PAT[6] = {0, 0, 1, 1, 2, 2};
  static const int PIS_SEC[6] = {1, 0, 1, 0, 1, 0};
  std::vector<char> claimed(nTr, 0);
  for (int st = 0; st < 6; ++st) {
    auto cands = makeD0(KPI_PAT[st], claimed);
    std::vector<std::array<int, 3>> to_claim;
    for (const auto& c : cands)
      for (size_t s = 0; s < good.size(); ++s) {
        const int k = good[s];
        if (k == c.iK || k == c.iPi) continue;
        if (!canPis[k] || claimed[k]) continue;
        if ((isSec(k) ? 1 : 0) != PIS_SEC[st]) continue;
        if (!pushDstar(c, k, st + 1, -1)) continue;
        // claiming: right-sign candidates carrying the claim label only
        if (claim_mode > 0 && out.ds.rs.back() == 1 &&
            ((claim_mode == 1 && out.ds.tight.back()) ||
             (claim_mode == 2 && out.ds.loose.back())))
          to_claim.push_back({c.iK, c.iPi, k});
      }
    for (const auto& t : to_claim)
      for (int q = 0; q < 3; ++q) claimed[t[q]] = 1;
  }

  // D0-alone cascade: its own three stages and its own claiming
  std::vector<char> claimed0(nTr, 0);
  for (int st = 0; st < 3; ++st) {
    auto cands = makeD0(st, claimed0);
    std::vector<std::array<int, 2>> to_claim;
    for (const auto& c : cands) {
      storeD0(c, st + 1);
      if (claim_mode > 0 &&
          ((claim_mode == 1 && out.d0.tight.back()) ||
           (claim_mode == 2 && out.d0.loose.back())))
        to_claim.push_back({c.iK, c.iPi});
    }
    for (const auto& t : to_claim) { claimed0[t[0]] = 1; claimed0[t[1]] = 1; }
  }

  // the two cascades are independent, so link each D* to the entry of the
  // D0 list built from the same two tracks with the same assignment (-1 when
  // that pair was claimed away in the D0 cascade)
  for (size_t k = 0; k < out.ds.m_kpi.size(); ++k) {
    const int wk = out.ds.trkK.origIdx[k], wp = out.ds.trkPi.origIdx[k];
    int found = -1;
    for (size_t c = 0; c < out.d0.m_kpi.size(); ++c)
      if (out.d0.trkK.origIdx[c] == wk && out.d0.trkPi.origIdx[c] == wp) {
        found = (int)c;
        break;
      }
    out.ds.d0idx[k] = found;
  }
  return out;
}

// ---------------------------------------------------------------------------
// Per-track auxiliary quantities, aligned with a selected trackstate
// collection through its original-Tracks index map.
// ---------------------------------------------------------------------------

// Component `which` (0 = VDET, 1 = ITC, 2 = TPC) of the per-track
// subdetectorHitNumbers block. -1 when the block is missing or too short.
inline RVec<int> subdetHits(const RVec<int>& orig_idx,
                            const RVec<unsigned int>& begin,
                            const RVec<unsigned int>& end,
                            const RVec<int>& values, int which) {
  RVec<int> out;
  for (int o : orig_idx) {
    int v = -1;
    if (o >= 0 && o < (int)begin.size() && o < (int)end.size()) {
      unsigned int b = begin[o], e = end[o];
      if (b + which < e && b + which < values.size()) v = values[b + which];
    }
    out.push_back(v);
  }
  return out;
}

// 1 for each entry whose original-Tracks index appears in `set_orig`.
// Stored diagnostic ("was this daughter in the fitted primary set"), never a cut.
inline RVec<int> flagInSet(const RVec<int>& orig_idx, const RVec<int>& set_orig) {
  RVec<int> out;
  for (int o : orig_idx)
    out.push_back((o >= 0 && std::find(set_orig.begin(), set_orig.end(), o) !=
                   set_orig.end()) ? 1 : 0);
  return out;
}

// Staging class of each entry: 0 = in the fitted primary set, 1 = in the
// secondary set, 2 = in neither (baseline-selected only).
inline RVec<int> poolClass(const RVec<int>& orig_idx,
                           const RVec<int>& prim_orig,
                           const RVec<int>& sec_orig) {
  RVec<int> out;
  for (int o : orig_idx) {
    int v = 2;
    if (o >= 0) {
      if (std::find(sec_orig.begin(), sec_orig.end(), o) != sec_orig.end()) v = 1;
      else if (std::find(prim_orig.begin(), prim_orig.end(), o) != prim_orig.end()) v = 0;
    }
    out.push_back(v);
  }
  return out;
}

// Original-Tracks indices of the daughters of the selected V0 candidates
// (select with the tight flag). Only consumed by the OPT-IN pairing veto.
inline RVec<int> claimedOrigIdx(const RVec<int>& d1, const RVec<int>& d2,
                                const RVec<int>& keep) {
  RVec<int> out;
  for (size_t i = 0; i < d1.size() && i < d2.size(); ++i) {
    if (i < keep.size() && !keep[i]) continue;
    if (d1[i] >= 0) out.push_back(d1[i]);
    if (d2[i] >= 0) out.push_back(d2[i]);
  }
  return out;
}

inline RVec<float> trackChi2Ndf(const RVec<int>& orig_idx,
                                const RVec<float>& chi2,
                                const RVec<int>& ndf) {
  RVec<float> out;
  for (int o : orig_idx) {
    float v = -1.f;
    if (o >= 0 && o < (int)chi2.size() && o < (int)ndf.size() && ndf[o] != 0)
      v = chi2[o] / float(ndf[o]);
    out.push_back(v);
  }
  return out;
}

// ---------------------------------------------------------------------------
// MC truth. The MCParticles daughter relations are EMPTY here; provenance
// comes from generatorStatus = 10000*KS + line, with `line` the 1-BASED LUND
// line of the MOTHER (MCParticles entry i is LUND line i+1). Geant
// secondaries carry generatorStatus 0 and have no LUND mother.
// ---------------------------------------------------------------------------
inline int lundMotherIdx(int genStatus) {
  return (genStatus > 0) ? (genStatus % 10000) - 1 : -1;
}

// heaviest quark in a hadron PDG code: 5 (b), 4 (c), else 0
inline int heavyFlavour(int pdg) {
  int a = std::abs(pdg);
  if (a >= 1000000000) return 0;  // nucleus
  int nq3 = (a / 10) % 10, nq2 = (a / 100) % 10, nq1 = (a / 1000) % 10;
  int mx = std::max({nq1, nq2, nq3});
  if (mx >= 5) return (mx == 5) ? 5 : 0;  // top-flavoured codes not expected
  return (mx == 4) ? 4 : 0;
}

// LUND mother index of every MCParticles entry, -1 when there is none.
inline std::vector<int> lundMothers(const RVec<edm4hep::MCParticleData>& mc) {
  const int n = mc.size();
  std::vector<int> moth(n, -1);
  for (int i = 0; i < n; ++i) {
    int m = lundMotherIdx(mc[i].generatorStatus);
    moth[i] = (m >= 0 && m < n) ? m : -1;
  }
  return moth;
}

// 5 if a b hadron is found anywhere in the ancestry, else 4 if a c hadron, else 0.
inline int ancestorFlavour(const std::vector<int>& moth,
                           const RVec<edm4hep::MCParticleData>& mc, int start) {
  int orig = 0;
  for (int a = start, guard = 0; a >= 0 && guard < 200; a = moth[a], ++guard) {
    int hf = heavyFlavour(mc[a].PDG);
    if (hf == 5) return 5;
    if (hf == 4 && orig == 0) orig = 4;  // keep looking: a b ancestor wins
  }
  return orig;
}

// Which reconstructed-track pool the first track linked to MC particle
// `mcidx` ended up in: -1 = no linked track, 0 = fitted PRIMARY set,
// 1 = SECONDARY set, 2 = linked but in neither (e.g. baseline selection).
inline int trackPool(int mcidx, const RVec<RVec<int>>& mcToTracks,
                     const RVec<int>& prim_orig, const RVec<int>& sec_orig) {
  if (mcidx < 0 || mcidx >= (int)mcToTracks.size() || mcToTracks[mcidx].empty())
    return -1;
  int trk = mcToTracks[mcidx][0];
  if (trk < 0) return -1;
  if (std::find(prim_orig.begin(), prim_orig.end(), trk) != prim_orig.end()) return 0;
  if (std::find(sec_orig.begin(), sec_orig.end(), trk) != sec_orig.end()) return 1;
  return 2;
}

// True D0 -> K pi two-body decays, whatever the mother (D* included).
struct TrueD0s {
  RVec<int>   idx;        // MCParticles index of the D0
  RVec<int>   iK, iPi;    // MC indices of the two daughters
  RVec<int>   mothPdg;
  RVec<int>   origin;     // 5 = from b hadron, 4 = from c hadron, 0 = other
  RVec<int>   fromDstar;  // 1 = immediate mother is a D*(2010)
  RVec<float> p, pt, costheta, xE;
  RVec<float> pK, cosK, pPi;
  RVec<float> flight;     // |D0 decay point - D0 production point| [cm]
  RVec<int>   nmatched;   // daughters with >= 1 linked track (0/1/2)
  RVec<int>   K_pool, pi_pool;  // track pool of each daughter (see trackPool)
};

// The two daughters of a two-body D0 -> K-pi+ (charge-conjugate included);
// returns false unless there are exactly those two LUND children.
inline bool d0Daughters(const std::vector<int>& moth,
                        const RVec<edm4hep::MCParticleData>& mc,
                        int i, int& iK, int& iPi) {
  const int n = mc.size();
  iK = -1; iPi = -1;
  int nchild = 0;
  for (int k = 0; k < n; ++k) {
    if (moth[k] != i) continue;
    ++nchild;
    if (std::abs(mc[k].PDG) == PDG_K && iK < 0) iK = k;
    else if (std::abs(mc[k].PDG) == PDG_PI && iPi < 0) iPi = k;
  }
  if (nchild != 2 || iK < 0 || iPi < 0) return false;
  return mc[iK].PDG * mc[iPi].PDG < 0;  // K- pi+ or K+ pi-
}

inline TrueD0s findTrueD0s(const RVec<edm4hep::MCParticleData>& mc,
                           const RVec<RVec<int>>& mcToTracks,
                           const RVec<int>& prim_orig = {},
                           const RVec<int>& sec_orig = {}) {
  TrueD0s out;
  const int n = mc.size();
  if (n == 0) return out;
  std::vector<int> moth = lundMothers(mc);
  for (int i = 0; i < n; ++i) {
    if (std::abs(mc[i].PDG) != PDG_D0) continue;
    int iK = -1, iPi = -1;
    if (!d0Daughters(moth, mc, i, iK, iPi)) continue;

    TVector3 p3(mc[i].momentum.x, mc[i].momentum.y, mc[i].momentum.z);
    TVector3 pk(mc[iK].momentum.x, mc[iK].momentum.y, mc[iK].momentum.z);
    TVector3 pp(mc[iPi].momentum.x, mc[iPi].momentum.y, mc[iPi].momentum.z);
    // the D0 decay point is the production point of its daughters
    TVector3 fl(mc[iK].vertex.x - mc[i].vertex.x,
                mc[iK].vertex.y - mc[i].vertex.y,
                mc[iK].vertex.z - mc[i].vertex.z);
    int mp = (moth[i] >= 0) ? mc[moth[i]].PDG : 0;
    int nm = 0;
    if (iK < (int)mcToTracks.size() && !mcToTracks[iK].empty()) ++nm;
    if (iPi < (int)mcToTracks.size() && !mcToTracks[iPi].empty()) ++nm;

    out.idx.push_back(i);
    out.iK.push_back(iK);
    out.iPi.push_back(iPi);
    out.mothPdg.push_back(mp);
    out.origin.push_back(ancestorFlavour(moth, mc, moth[i]));
    out.fromDstar.push_back((std::abs(mp) == PDG_DSTAR) ? 1 : 0);
    out.p.push_back(p3.Mag());
    out.pt.push_back(p3.Perp());
    out.costheta.push_back(p3.Mag() > 0. ? p3.Z() / p3.Mag() : 0.);
    out.xE.push_back(std::sqrt(p3.Mag2() + double(mc[i].mass) * mc[i].mass) / E_BEAM);
    out.pK.push_back(pk.Mag());
    out.cosK.push_back(pk.Mag() > 0. ? pk.Z() / pk.Mag() : 0.);
    out.pPi.push_back(pp.Mag());
    out.flight.push_back(fl.Mag());
    out.nmatched.push_back(nm);
    out.K_pool.push_back(trackPool(iK, mcToTracks, prim_orig, sec_orig));
    out.pi_pool.push_back(trackPool(iPi, mcToTracks, prim_orig, sec_orig));
  }
  return out;
}

// True D*+ -> D0 pi+ with D0 -> K pi, both two-body.
struct TrueDstars {
  RVec<int>   idx;              // MCParticles index of the D*
  RVec<int>   d0idx;            // MCParticles index of its D0
  RVec<int>   iK, iPi, iPis;    // MC indices of the three charged daughters
  RVec<int>   mothPdg;
  RVec<int>   origin;           // 5 = from b hadron, 4 = from c hadron, 0 = other
  RVec<float> p, pt, costheta, xE, px, py, pz;
  RVec<float> pK, cosK, pPi, pPis, cosPis;
  RVec<float> d0flight;         // D0 flight distance [cm]
  RVec<int>   nmatched;         // daughters with >= 1 linked track (0-3)
  RVec<int>   K_pool, pi_pool, pis_pool;  // track pool of each daughter (see trackPool)
};

inline TrueDstars findTrueDstars(const RVec<edm4hep::MCParticleData>& mc,
                                 const RVec<RVec<int>>& mcToTracks,
                                 const RVec<int>& prim_orig = {},
                                 const RVec<int>& sec_orig = {}) {
  TrueDstars out;
  const int n = mc.size();
  if (n == 0) return out;
  std::vector<int> moth = lundMothers(mc);
  for (int i = 0; i < n; ++i) {
    if (std::abs(mc[i].PDG) != PDG_DSTAR) continue;
    // LUND children of the D*: exactly one D0 and one charged pion
    int id0 = -1, ipis = -1, nchild = 0;
    for (int k = 0; k < n; ++k) {
      if (moth[k] != i) continue;
      ++nchild;
      if (std::abs(mc[k].PDG) == PDG_D0 && id0 < 0) id0 = k;
      else if (std::abs(mc[k].PDG) == PDG_PI && ipis < 0) ipis = k;
    }
    if (nchild != 2 || id0 < 0 || ipis < 0) continue;
    int iK = -1, iPi = -1;
    if (!d0Daughters(moth, mc, id0, iK, iPi)) continue;

    TVector3 p3(mc[i].momentum.x, mc[i].momentum.y, mc[i].momentum.z);
    TVector3 pk(mc[iK].momentum.x, mc[iK].momentum.y, mc[iK].momentum.z);
    TVector3 pp(mc[iPi].momentum.x, mc[iPi].momentum.y, mc[iPi].momentum.z);
    TVector3 ps(mc[ipis].momentum.x, mc[ipis].momentum.y, mc[ipis].momentum.z);
    TVector3 fl(mc[iK].vertex.x - mc[id0].vertex.x,
                mc[iK].vertex.y - mc[id0].vertex.y,
                mc[iK].vertex.z - mc[id0].vertex.z);
    int nm = 0;
    for (int d : {iK, iPi, ipis})
      if (d < (int)mcToTracks.size() && !mcToTracks[d].empty()) ++nm;

    out.idx.push_back(i);
    out.d0idx.push_back(id0);
    out.iK.push_back(iK);
    out.iPi.push_back(iPi);
    out.iPis.push_back(ipis);
    out.mothPdg.push_back((moth[i] >= 0) ? mc[moth[i]].PDG : 0);
    out.origin.push_back(ancestorFlavour(moth, mc, moth[i]));
    out.p.push_back(p3.Mag());
    out.pt.push_back(p3.Perp());
    out.costheta.push_back(p3.Mag() > 0. ? p3.Z() / p3.Mag() : 0.);
    out.xE.push_back(std::sqrt(p3.Mag2() + double(mc[i].mass) * mc[i].mass) / E_BEAM);
    out.px.push_back(p3.X()); out.py.push_back(p3.Y()); out.pz.push_back(p3.Z());
    out.pK.push_back(pk.Mag());
    out.cosK.push_back(pk.Mag() > 0. ? pk.Z() / pk.Mag() : 0.);
    out.pPi.push_back(pp.Mag());
    out.pPis.push_back(ps.Mag());
    out.cosPis.push_back(ps.Mag() > 0. ? ps.Z() / ps.Mag() : 0.);
    out.d0flight.push_back(fl.Mag());
    out.nmatched.push_back(nm);
    out.K_pool.push_back(trackPool(iK, mcToTracks, prim_orig, sec_orig));
    out.pi_pool.push_back(trackPool(iPi, mcToTracks, prim_orig, sec_orig));
    out.pis_pool.push_back(trackPool(ipis, mcToTracks, prim_orig, sec_orig));
  }
  return out;
}

// Per-candidate truth for the D0-alone list.
//   cls 1 = both tracks link to the K and the pi of one true D0 -> K pi with
//           the CORRECT mass assignment, 2 = the same with K and pi swapped,
//           0 = otherwise.
struct D0Truth {
  RVec<int> cls;
  RVec<int> trueidx;                     // index into the TrueD0s lists, -1 if none
  RVec<int> trkK_mcpdg, trkPi_mcpdg;     // PDG of the linked MC particle (0 = unlinked)
  RVec<int> trkK_mothpdg, trkPi_mothpdg; // PDG of its LUND mother (0 = none)
};

// first linked MC particle of an original-Tracks index, -1 when unlinked
inline int firstMCOfTrack(const RVec<RVec<int>>& trackToMCs, int trk, int nmc) {
  if (trk < 0 || trk >= (int)trackToMCs.size() || trackToMCs[trk].empty()) return -1;
  int m = trackToMCs[trk][0];
  return (m >= 0 && m < nmc) ? m : -1;
}

inline D0Truth classifyD0(const D0Block& c,
                          const RVec<RVec<int>>& trackToMCs,
                          const RVec<edm4hep::MCParticleData>& mc,
                          const TrueD0s& td) {
  D0Truth out;
  const int n = mc.size();
  for (size_t k = 0; k < c.m_kpi.size(); ++k) {
    int m1 = firstMCOfTrack(trackToMCs, c.trkK.origIdx[k], n);
    int m2 = firstMCOfTrack(trackToMCs, c.trkPi.origIdx[k], n);
    int mo1 = (m1 >= 0) ? lundMotherIdx(mc[m1].generatorStatus) : -1;
    int mo2 = (m2 >= 0) ? lundMotherIdx(mc[m2].generatorStatus) : -1;
    out.trkK_mcpdg.push_back((m1 >= 0) ? mc[m1].PDG : 0);
    out.trkPi_mcpdg.push_back((m2 >= 0) ? mc[m2].PDG : 0);
    out.trkK_mothpdg.push_back((mo1 >= 0 && mo1 < n) ? mc[mo1].PDG : 0);
    out.trkPi_mothpdg.push_back((mo2 >= 0 && mo2 < n) ? mc[mo2].PDG : 0);
    int cls = 0, ti = -1;
    if (m1 >= 0 && m2 >= 0)
      for (size_t t = 0; t < td.idx.size(); ++t) {
        if (td.iK[t] == m1 && td.iPi[t] == m2) { cls = 1; ti = (int)t; break; }
        if (td.iK[t] == m2 && td.iPi[t] == m1) { cls = 2; ti = (int)t; break; }
      }
    out.cls.push_back(cls);
    out.trueidx.push_back(ti);
  }
  return out;
}

// Per-candidate truth for the D* list.
//   cls 1 = the three tracks link to the K, pi and slow pi of one true D*
//           with the CORRECT K/pi assignment
//       2 = the same true D*, K and pi swapped
//       3 = the K,pi pair is a true D0 (correct assignment) coming from a D*,
//           but the slow-pion track is not that D*'s slow pion
//       4 = the K,pi pair is a true D0 -> K pi (correct assignment) whose
//           mother is not a D*
//       0 = otherwise
// trueidx indexes the TrueDstars lists and is filled for classes 1, 2 and 3.
struct DstarTruth {
  RVec<int> cls;
  RVec<int> trueidx;
  RVec<int> trkK_mcpdg, trkPi_mcpdg, trkPis_mcpdg;
  RVec<int> trkK_mothpdg, trkPi_mothpdg, trkPis_mothpdg;
};

inline DstarTruth classifyDstar(const DstarBlock& c,
                                const RVec<RVec<int>>& trackToMCs,
                                const RVec<edm4hep::MCParticleData>& mc,
                                const TrueDstars& tds,
                                const TrueD0s& td0) {
  DstarTruth out;
  const int n = mc.size();
  for (size_t k = 0; k < c.m_kpi.size(); ++k) {
    int m1 = firstMCOfTrack(trackToMCs, c.trkK.origIdx[k], n);
    int m2 = firstMCOfTrack(trackToMCs, c.trkPi.origIdx[k], n);
    int m3 = firstMCOfTrack(trackToMCs, c.trkPis.origIdx[k], n);
    int mo1 = (m1 >= 0) ? lundMotherIdx(mc[m1].generatorStatus) : -1;
    int mo2 = (m2 >= 0) ? lundMotherIdx(mc[m2].generatorStatus) : -1;
    int mo3 = (m3 >= 0) ? lundMotherIdx(mc[m3].generatorStatus) : -1;
    out.trkK_mcpdg.push_back((m1 >= 0) ? mc[m1].PDG : 0);
    out.trkPi_mcpdg.push_back((m2 >= 0) ? mc[m2].PDG : 0);
    out.trkPis_mcpdg.push_back((m3 >= 0) ? mc[m3].PDG : 0);
    out.trkK_mothpdg.push_back((mo1 >= 0 && mo1 < n) ? mc[mo1].PDG : 0);
    out.trkPi_mothpdg.push_back((mo2 >= 0 && mo2 < n) ? mc[mo2].PDG : 0);
    out.trkPis_mothpdg.push_back((mo3 >= 0 && mo3 < n) ? mc[mo3].PDG : 0);

    int cls = 0, ti = -1;
    // true D* whose D0 daughters are this K,pi pair (either assignment)
    int swapped = 0;
    if (m1 >= 0 && m2 >= 0)
      for (size_t t = 0; t < tds.idx.size(); ++t) {
        if (tds.iK[t] == m1 && tds.iPi[t] == m2) { ti = (int)t; swapped = 0; break; }
        if (tds.iK[t] == m2 && tds.iPi[t] == m1) { ti = (int)t; swapped = 1; break; }
      }
    if (ti >= 0) {
      bool pis_ok = (m3 >= 0 && m3 == tds.iPis[ti]);
      if (!swapped) cls = pis_ok ? 1 : 3;
      else if (pis_ok) cls = 2;
    } else if (m1 >= 0 && m2 >= 0) {
      // no D* parent: a true D0 -> K pi with the correct assignment
      for (size_t t = 0; t < td0.idx.size(); ++t)
        if (td0.iK[t] == m1 && td0.iPi[t] == m2) { cls = 4; break; }
    }
    out.cls.push_back(cls);
    out.trueidx.push_back(ti);
  }
  return out;
}

// 1 if the true D* at position t was reconstructed by at least one candidate
// of class 1 carrying `flag` (the loose or tight label).
inline RVec<int> trueDstarFound(const TrueDstars& tds, const DstarTruth& info,
                                const RVec<int>& flag) {
  RVec<int> out(tds.idx.size(), 0);
  for (size_t k = 0; k < info.cls.size(); ++k) {
    if (info.cls[k] != 1) continue;
    if (k < flag.size() && !flag[k]) continue;
    int t = info.trueidx[k];
    if (t >= 0 && t < (int)out.size()) out[t] = 1;
  }
  return out;
}

// Same for the D0-alone list: class 1 candidate carrying `flag`.
inline RVec<int> trueD0Found(const TrueD0s& td, const D0Truth& info,
                             const RVec<int>& flag) {
  RVec<int> out(td.idx.size(), 0);
  for (size_t k = 0; k < info.cls.size(); ++k) {
    if (info.cls[k] != 1) continue;
    if (k < flag.size() && !flag[k]) continue;
    int t = info.trueidx[k];
    if (t >= 0 && t < (int)out.size()) out[t] = 1;
  }
  return out;
}

} // namespace AlephDstar
} // namespace FCCAnalyses

#endif
