#ifndef ALEPHDSTAR_H
#define ALEPHDSTAR_H

/*
  D*+ -> D0(K- pi+) pi+_slow reconstruction as a kinematically tagged kaon
  source; no dE/dx quantity enters any selection. Input = flipD0_copy'ed
  trackstates (physical charge = +sign(omega)); lengths cm, momenta GeV.
*/

#include <algorithm>
#include <cmath>
#include <vector>

#include <ROOT/RVec.hxx>
#include "TVector3.h"
#include "TLorentzVector.h"

#include "edm4hep/TrackState.h"
#include "edm4hep/MCParticleData.h"
#include "FCCAnalyses/VertexingUtils.h"
#include "FCCAnalyses/VertexFitterSimple.h"

#include "aleph_units.h"
#include "analyzer_trkaux.h"
#include "analyzer_v0new.h"

namespace FCCAnalyses {
namespace AlephDstar {

using ROOT::VecOps::RVec;

// shared track auxiliaries and LUND-truth helpers (see analyzer_trkaux.h)
using AlephTrkAux::perigeeMomentum;
using AlephTrkAux::vertexDistSig;
using AlephTrkAux::fitTracksCm;
using AlephTrkAux::TrkBlock;
using AlephTrkAux::TrkAux;
using AlephTrkAux::pushTrk;
using AlephTrkAux::reserveTrk;
using AlephTrkAux::lundMothers;
using AlephTrkAux::lundChildren;
using AlephTrkAux::ancestorFlavour;
using AlephTrkAux::pushLegTruth;
using AlephTrkAux::setLegTruth;
using AlephTrkAux::mcLinksOfTrack;
using AlephTrkAux::firstOr;
using AlephTrkAux::foundFlags;

constexpr double M_K     = AlephMasses::kK;     // charged kaon
constexpr double M_PICH  = AlephMasses::kPiCh;  // charged pion
constexpr double M_D0    = AlephMasses::kD0;
constexpr double M_DSTAR = AlephMasses::kDstar;
constexpr double DM_NOMINAL = AlephMasses::kDeltaMDstarD0;
constexpr int PDG_DSTAR = 413, PDG_D0 = 421, PDG_K = 321, PDG_PI = 211;

constexpr double E_BEAM = AlephUnits::kEBeam;  // for xE = E/E_BEAM
constexpr double PRE_MARGIN = 0.05;  // pre-fit mass window margin [GeV]

// ---------------------------------------------------------------------------
// Selection: single named source, read directly by findDstar, so stage1.py
// passes no values. Storage stays deliberately loose (any working point is
// re-derivable offline); the labels below are the tuned selection.
// ---------------------------------------------------------------------------
constexpr double M_LO = 1.70, M_HI = 2.03;   // stored K pi mass window [GeV]
constexpr double CHI2_CUT = 25.;             // D0 vertex chi2 (ndf = 1)
constexpr double DPV_FID = 10.;              // |vtx-PV| storage fiducial [cm]
constexpr double DM_MAX = 0.20;              // stored dm ceiling [GeV]
// applied to the PERIGEE momentum, while the stored trk p is the at-vertex one
constexpr double P_MIN = 0.3;                // K/pi |p| floor [GeV]
constexpr double PS_MIN = 0.1;               // slow-pion |p| floor [GeV]

// Working points stored as flags, evaluated on the stored quantities (post-fit
// mass, at-vertex daughter momenta). D0-alone: no dm handle, so the purity has
// to come from displacement, pointing and the helicity angle.
constexpr double D0_LOOSE_DM = 0.060;        // |m(K pi) - m_D0| [GeV]
constexpr double D0_TIGHT_DM = 0.030;
constexpr double D0_TIGHT_DPVSIG = 3.0;      // |vtx-PV| significance
constexpr double D0_TIGHT_COSPOINT = 0.99;
constexpr double D0_TIGHT_COSSTAR = 0.8;     // |cos(theta*)| of the kaon
// D*: the dm window carries the tag, its resolution set by the slow pion alone
constexpr double DS_LOOSE_DM = 0.050, DS_LOOSE_DDM = 0.0030;
constexpr double DS_TIGHT_DM = 0.025, DS_TIGHT_DDM = 0.0015;
constexpr double DS_TIGHT_PS = 0.3;          // slow-pion |p| [GeV]
constexpr double DS_TIGHT_COSPOINT = 0.95;   // D0 vertex pointing
constexpr double TIGHT_PK = 1.0, TIGHT_PPI = 1.0;  // daughter |p| [GeV]
constexpr double TIGHT_CHI2 = 10.;           // D0 vertex chi2

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

// ---------------------------------------------------------------------------
// Kinematics common to the D0 and the D* entry: the vertex is the fitted D0
// vertex in both, and the D* repeats its D0's quantities.
// ---------------------------------------------------------------------------
struct CandKin {
  RVec<float> m_kpi;           // K pi mass at the fitted vertex [GeV]
  RVec<float> p, px, py, pz, costheta, xE;
  RVec<float> chi2;            // vertex fit chi2, normalised (ndf = 1)
  RVec<float> vx, vy, vz;      // fitted vertex [cm]
  RVec<float> dpv, dpvSig;     // |vtx-PV| [cm] and its 3D significance
  RVec<float> cosPoint;        // cos(angle) between p(D0) and (vtx - PV)
  RVec<float> cosThetaStar;    // kaon helicity angle
};

// D0 -> K pi candidates, one entry per (pair, mass assignment).
struct D0Block {
  CandKin     kin;
  RVec<int>   loose, tight;    // labels, not cuts
  RVec<int>   nsec;            // legs in the secondary pool (0-2)
  TrkBlock    trkK, trkPi;
};

// D* -> D0 pi_slow candidates, one entry per (D0 candidate, third track).
// d0idx points back to the D0 list.
struct DstarBlock {
  CandKin     kin;
  RVec<float> dm;
  RVec<int>   rs;              // 1 = right-sign slow pion (charge of the pi)
  RVec<int>   loose, tight;
  RVec<int>   d0idx;           // index into the D0 list, -1 if absent from it
  RVec<int>   nsec;            // legs in the secondary pool (0-3)
  TrkBlock    trkK, trkPi, trkPis;
};

struct DstarCands {
  D0Block    d0;
  DstarBlock ds;
  int        nfits = 0;   // vertex fits actually performed
};

inline void reserveKin(CandKin& k, size_t n) {
  k.m_kpi.reserve(n); k.p.reserve(n); k.px.reserve(n); k.py.reserve(n);
  k.pz.reserve(n); k.costheta.reserve(n); k.xE.reserve(n); k.chi2.reserve(n);
  k.vx.reserve(n); k.vy.reserve(n); k.vz.reserve(n); k.dpv.reserve(n);
  k.dpvSig.reserve(n); k.cosPoint.reserve(n); k.cosThetaStar.reserve(n);
}

// The D0 vertex quantities of one candidate; pmag is |p3|, passed in because
// every caller has it already.
inline void pushKin(CandKin& k, double m, const TVector3& p3, double pmag,
                    double energy, double chi2, const TVector3& x, double dpv,
                    float dpvSig, double cosPoint, double cosThetaStar) {
  k.m_kpi.push_back(m);
  k.p.push_back(pmag);
  k.px.push_back(p3.X());
  k.py.push_back(p3.Y());
  k.pz.push_back(p3.Z());
  k.costheta.push_back(p3.Z() / pmag);
  k.xE.push_back(energy / E_BEAM);
  k.chi2.push_back(chi2);
  k.vx.push_back(x.X());
  k.vy.push_back(x.Y());
  k.vz.push_back(x.Z());
  k.dpv.push_back(dpv);
  k.dpvSig.push_back(dpvSig);
  k.cosPoint.push_back(cosPoint);
  k.cosThetaStar.push_back(cosThetaStar);
}

// One growth-free pass over the output vectors: the D0 list is exactly one
// entry per candidate, the D* list at least that long.
inline void reserveD0(D0Block& d, size_t n) {
  reserveKin(d.kin, n);
  d.loose.reserve(n); d.tight.reserve(n); d.nsec.reserve(n);
  reserveTrk(d.trkK, n); reserveTrk(d.trkPi, n);
}

inline void reserveDstar(DstarBlock& d, size_t n) {
  reserveKin(d.kin, n);
  d.dm.reserve(n); d.rs.reserve(n); d.loose.reserve(n);
  d.tight.reserve(n); d.d0idx.reserve(n); d.nsec.reserve(n);
  reserveTrk(d.trkK, n); reserveTrk(d.trkPi, n); reserveTrk(d.trkPis, n);
}

// ---------------------------------------------------------------------------
// The finder.
//   tracks     : flipD0_copy'ed baseline-selected trackstates
//   orig_idx   : original-Tracks index of each entry of `tracks`
//   nvdet/nitc : VDET/ITC hit counts of each entry (-1 = unknown)
//   chi2ndf    : track fit chi2/ndf of each entry (-1 = unknown)
//   isprim     : 1 if the entry is in the fitted primary set (stored, and read
//                by the dstar_tight primary-pattern veto)
//   pool       : 0 = primary set, 1 = secondary set, 2 = neither (stored)
//   PV         : fitted primary vertex (positions numerically cm)
//   veto_orig  : original-Tracks indices excluded from the pool (empty = none)
// The track pool is the FULL baseline-selected list, primary and secondary
// alike, with no masking by the PV split; every candidate is built in one
// all-track pass, with no exclusive claiming.
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
    const RVec<int>& veto_orig,
    double solenoidBz) {

  DstarCands out;
  const int nTr = tracks.size();
  if (nTr < 2) return out;

  const TrkAux aux{&orig_idx, &nvdet, &nitc, &chi2ndf, &isprim, &pool};
  auto nv = [&](const RVec<int>& src, int k) {
    return (k >= 0 && k < (int)src.size()) ? src[k] : -1;
  };
  auto oidx = [&](int k) {
    return (k >= 0 && k < (int)orig_idx.size()) ? orig_idx[k] : -1;
  };
  // SEC = in the secondary set; the rare "neither" tracks count as
  // primary-like here, while their true class stays in the stored pool branch
  auto isSec = [&](int k) { return k < (int)pool.size() && pool[k] == 1; };

  // momentum prefilter at the softest floor (the slow pion): the K/pi floor
  // only decides which ROLES a survivor may take, so a track too soft to be a
  // kaon can still be a slow pion.
  // eK/ePi = per-track hypothesis energies of the perigee momentum, for the
  // pre-fit mass window.
  std::vector<int> good;
  good.reserve(nTr);
  std::vector<TVector3> pper(nTr);
  std::vector<double> eK(nTr, 0.), ePi(nTr, 0.);
  std::vector<TLorentzVector> pisLv(nTr);  // slow-pion four-momentum
  std::vector<char> canKPi(nTr, 0);
  const std::vector<char> vetoed = AlephTrkAux::memberMask(veto_orig);
  for (int i = 0; i < nTr; ++i) {
    const auto& t = tracks[i];
    if (AlephTrkAux::inMask(vetoed, oidx(i))) continue;
    pper[i] = perigeeMomentum(t, solenoidBz);
    double pm = pper[i].Mag();
    if (pm < PS_MIN) continue;
    canKPi[i] = (pm >= P_MIN) ? 1 : 0;
    eK[i] = std::sqrt(pper[i].Mag2() + M_K * M_K);
    ePi[i] = std::sqrt(pper[i].Mag2() + M_PICH * M_PICH);
    pisLv[i] = TLorentzVector(pper[i], ePi[i]);
    good.push_back(i);
  }
  if (good.size() < 2) return out;

  TVector3 pv(PV.vertex.position[0], PV.vertex.position[1], PV.vertex.position[2]);
  const double pre_lo = M_LO - PRE_MARGIN, pre_hi = M_HI + PRE_MARGIN;

  struct FitRes {
    bool ok;
    double chi2, dis;
    float dsig;
    TVector3 x, d, pa, pb;
  };
  RVec<edm4hep::TrackState> tr_pair(2);
  // i < j is guaranteed by the callers; pa/pb follow that order
  auto fitPair = [&](int i, int j) {
    FitRes r;
    r.ok = false;
    r.chi2 = 0.; r.dis = 0.; r.dsig = -1.;
    tr_pair[0] = tracks[i];
    tr_pair[1] = tracks[j];
    auto v = fitTracksCm(tr_pair, tracks, solenoidBz);
    ++out.nfits;
    if (v.updated_track_momentum_at_vertex.size() == 2) {
      double chi2 = v.vertex.chi2;
      bool pass = (chi2 == chi2) && !(chi2 >= CHI2_CUT);
      if (pass) {
        TVector3 x(v.vertex.position[0], v.vertex.position[1], v.vertex.position[2]);
        TVector3 dvec = x - pv;
        double dis = dvec.Mag();
        if (!(dis > DPV_FID)) {
          r.ok = true;
          r.chi2 = chi2; r.x = x; r.d = dvec; r.dis = dis;
          r.dsig = vertexDistSig(dvec, v.vertex.covMatrix, PV.vertex.covMatrix);
          r.pa = v.updated_track_momentum_at_vertex[0];
          r.pb = v.updated_track_momentum_at_vertex[1];
        }
      }
    }
    return r;
  };

  // one accepted (pair, mass assignment) before it is written out
  struct D0Cand {
    int iK, iPi;
    double m, chi2, dis, cosp, cstar, energy;
    float dsig;
    TVector3 p3, x, pK, pPi;
    int loose, tight;
  };

  // every opposite-charge pair inside the pre-fit window, both assignments
  auto makeD0 = [&]() {
    std::vector<D0Cand> res;
    for (size_t a = 0; a + 1 < good.size(); ++a) {
      const int i = good[a];
      if (!canKPi[i]) continue;
      for (size_t b = a + 1; b < good.size(); ++b) {
        const int j = good[b];
        if (!canKPi[j]) continue;
        if (tracks[i].omega * tracks[j].omega > 0) continue;  // opposite charge only

        // pre-fit window on the perigee momenta: both mass assignments are
        // tried, and the pair is fitted only if at least one is inside
        const TVector3 p_pre = pper[i] + pper[j];
        const double p2_pre = p_pre.Mag2();
        const double e_ij = eK[i] + ePi[j], e_ji = ePi[i] + eK[j];
        double m_pre_ij = std::sqrt(std::max(0., e_ij * e_ij - p2_pre));
        double m_pre_ji = std::sqrt(std::max(0., e_ji * e_ji - p2_pre));
        bool ok_ij = (m_pre_ij >= pre_lo && m_pre_ij <= pre_hi);
        bool ok_ji = (m_pre_ji >= pre_lo && m_pre_ji <= pre_hi);
        if (!ok_ij && !ok_ji) continue;

        const FitRes fr = fitPair(i, j);
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
          if (c.m < M_LO || c.m > M_HI) continue;

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
          c.loose = (dmass < D0_LOOSE_DM) ? 1 : 0;
          c.tight = (dmass < D0_TIGHT_DM && c.dsig > D0_TIGHT_DPVSIG &&
                     c.cosp > D0_TIGHT_COSPOINT &&
                     std::abs(c.cstar) < D0_TIGHT_COSSTAR &&
                     c.pK.Mag() > TIGHT_PK && c.pPi.Mag() > TIGHT_PPI &&
                     c.chi2 < TIGHT_CHI2) ? 1 : 0;
          res.push_back(c);
        }
      }
    }
    return res;
  };

  auto storeD0 = [&](const D0Cand& c) {
    pushKin(out.d0.kin, c.m, c.p3, c.p3.Mag(), c.energy, c.chi2, c.x, c.dis,
            c.dsig, c.cosp, c.cstar);
    out.d0.loose.push_back(c.loose);
    out.d0.tight.push_back(c.tight);
    out.d0.nsec.push_back((isSec(c.iK) ? 1 : 0) + (isSec(c.iPi) ? 1 : 0));
    pushTrk(out.d0.trkK, tracks, c.iK, c.pK, aux);
    pushTrk(out.d0.trkPi, tracks, c.iPi, c.pPi, aux);
  };

  // D* from a D0 candidate plus the slow-pion track k. No three-track fit:
  // the D0 momenta are the at-vertex ones, the slow pion enters at its perigee.
  auto pushDstar = [&](const D0Cand& c, int k, int d0idx) {
    // the D0 mass enters dm as the STORED (float) value
    const float mf = c.m;
    const double m = mf;
    TLorentzVector dslv = TLorentzVector(c.p3, c.energy) + pisLv[k];
    double dm = dslv.M() - m;
    if (!(dm > 0.)) return false;
    if (dm >= DM_MAX) return false;
    TVector3 p3 = dslv.Vect();
    double pmag = p3.Mag();
    if (!(pmag > 0.)) return false;

    const float chi2f = c.chi2;
    const int qPi = (tracks[c.iPi].omega > 0) ? 1 : -1;
    const int qs = (tracks[k].omega > 0) ? 1 : -1;

    pushKin(out.ds.kin, mf, p3, pmag, dslv.E(), chi2f, c.x, c.dis, c.dsig,
            c.cosp, c.cstar);
    out.ds.dm.push_back(dm);
    out.ds.rs.push_back((qs == qPi) ? 1 : 0);
    out.ds.d0idx.push_back(d0idx);
    out.ds.nsec.push_back((isSec(c.iK) ? 1 : 0) + (isSec(c.iPi) ? 1 : 0) +
                          (isSec(k) ? 1 : 0));
    pushTrk(out.ds.trkK, tracks, c.iK, c.pK, aux);
    pushTrk(out.ds.trkPi, tracks, c.iPi, c.pPi, aux);
    pushTrk(out.ds.trkPis, tracks, k, pper[k], aux);

    const double dmass = std::abs(m - M_D0);
    const double ddm = std::abs(dm - DM_NOMINAL);
    // primary-pattern veto: a D0 built from two primary-set tracks with a
    // secondary slow pion is the combinatorial topology, never a real D*
    const bool bad_pattern = (nv(isprim, c.iK) == 1 && nv(isprim, c.iPi) == 1 &&
                              nv(isprim, k) == 0);
    out.ds.loose.push_back((dmass < DS_LOOSE_DM && ddm < DS_LOOSE_DDM) ? 1 : 0);
    out.ds.tight.push_back((dmass < DS_TIGHT_DM && ddm < DS_TIGHT_DDM &&
                            c.pK.Mag() > TIGHT_PK && c.pPi.Mag() > TIGHT_PPI &&
                            chi2f < TIGHT_CHI2 && !bad_pattern &&
                            pper[k].Mag() > DS_TIGHT_PS &&
                            c.cosp > DS_TIGHT_COSPOINT) ? 1 : 0);
    return true;
  };

  // single all-track pass: every D0 candidate, then every third track
  auto cands = makeD0();
  reserveD0(out.d0, cands.size());
  reserveDstar(out.ds, cands.size());
  for (const auto& c : cands) storeD0(c);
  for (size_t c = 0; c < cands.size(); ++c)
    for (size_t s = 0; s < good.size(); ++s) {
      const int k = good[s];
      if (k == cands[c].iK || k == cands[c].iPi) continue;
      pushDstar(cands[c], k, (int)c);
    }
  return out;
}

// ---------------------------------------------------------------------------
// Per-track auxiliary quantities, aligned with a selected trackstate
// collection through its original-Tracks index map.
// ---------------------------------------------------------------------------

// Staging class of each entry: 0 = in the fitted primary set, 1 = in the
// secondary set, 2 = in neither (baseline-selected only).
inline RVec<int> poolClass(const RVec<int>& orig_idx,
                           const RVec<int>& prim_orig,
                           const RVec<int>& sec_orig) {
  RVec<int> out;
  const std::vector<char> sec = AlephTrkAux::memberMask(sec_orig);
  const std::vector<char> prim = AlephTrkAux::memberMask(prim_orig);
  for (int o : orig_idx) {
    int v = 2;
    if (AlephTrkAux::inMask(sec, o)) v = 1;
    else if (AlephTrkAux::inMask(prim, o)) v = 0;
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
// Which reconstructed-track pool the first track linked to MC particle
// `mcidx` ended up in: -1 = no linked track, 0 = fitted PRIMARY set,
// 1 = SECONDARY set, 2 = linked but in neither (e.g. baseline selection).
inline int trackPool(int mcidx, const RVec<RVec<int>>& mcToTracks,
                     const std::vector<char>& prim, const std::vector<char>& sec) {
  if (mcidx < 0 || mcidx >= (int)mcToTracks.size() || mcToTracks[mcidx].empty())
    return -1;
  int trk = mcToTracks[mcidx][0];
  if (trk < 0) return -1;
  if (AlephTrkAux::inMask(prim, trk)) return 0;
  if (AlephTrkAux::inMask(sec, trk)) return 1;
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
inline bool d0Daughters(const std::vector<std::vector<int>>& child,
                        const RVec<edm4hep::MCParticleData>& mc,
                        int i, int& iK, int& iPi) {
  iK = -1; iPi = -1;
  int nchild = 0;
  for (int k : child[i]) {
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
  const std::vector<std::vector<int>> child = lundChildren(moth);
  const std::vector<char> prim = AlephTrkAux::memberMask(prim_orig);
  const std::vector<char> sec = AlephTrkAux::memberMask(sec_orig);
  for (int i = 0; i < n; ++i) {
    if (std::abs(mc[i].PDG) != PDG_D0) continue;
    int iK = -1, iPi = -1;
    if (!d0Daughters(child, mc, i, iK, iPi)) continue;

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
    out.K_pool.push_back(trackPool(iK, mcToTracks, prim, sec));
    out.pi_pool.push_back(trackPool(iPi, mcToTracks, prim, sec));
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
  const std::vector<std::vector<int>> child = lundChildren(moth);
  const std::vector<char> prim = AlephTrkAux::memberMask(prim_orig);
  const std::vector<char> sec = AlephTrkAux::memberMask(sec_orig);
  for (int i = 0; i < n; ++i) {
    if (std::abs(mc[i].PDG) != PDG_DSTAR) continue;
    // LUND children of the D*: exactly one D0 and one charged pion
    int id0 = -1, ipis = -1, nchild = 0;
    for (int k : child[i]) {
      ++nchild;
      if (std::abs(mc[k].PDG) == PDG_D0 && id0 < 0) id0 = k;
      else if (std::abs(mc[k].PDG) == PDG_PI && ipis < 0) ipis = k;
    }
    if (nchild != 2 || id0 < 0 || ipis < 0) continue;
    int iK = -1, iPi = -1;
    if (!d0Daughters(child, mc, id0, iK, iPi)) continue;

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
    out.K_pool.push_back(trackPool(iK, mcToTracks, prim, sec));
    out.pi_pool.push_back(trackPool(iPi, mcToTracks, prim, sec));
    out.pis_pool.push_back(trackPool(ipis, mcToTracks, prim, sec));
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

// First true (K, pi) pair in (iK, iPi) reached by any link of the two legs:
// the direct assignment is tried before the swapped one for each true entry.
// Returns the true index and the matched links.
inline bool matchPair(const RVec<int>& iK, const RVec<int>& iPi,
                      const RVec<int>& L1, const RVec<int>& L2,
                      int& ti, int& m1, int& m2) {
  for (size_t t = 0; t < iK.size(); ++t) {
    for (int a : L1) for (int b : L2)
      if (iK[t] == a && iPi[t] == b) { ti = (int)t; m1 = a; m2 = b; return true; }
    for (int a : L1) for (int b : L2)
      if (iK[t] == b && iPi[t] == a) { ti = (int)t; m1 = a; m2 = b; return true; }
  }
  return false;
}

inline D0Truth classifyD0(const D0Block& c,
                          const RVec<RVec<int>>& trackToMCs,
                          const RVec<edm4hep::MCParticleData>& mc,
                          const TrueD0s& td) {
  D0Truth out;
  const int nmc = mc.size();
  for (size_t k = 0; k < c.kin.m_kpi.size(); ++k) {
    const RVec<int> L1 = mcLinksOfTrack(trackToMCs, c.trkK.origIdx[k], nmc);
    const RVec<int> L2 = mcLinksOfTrack(trackToMCs, c.trkPi.origIdx[k], nmc);
    pushLegTruth(mc, firstOr(L1), out.trkK_mcpdg, out.trkK_mothpdg);
    pushLegTruth(mc, firstOr(L2), out.trkPi_mcpdg, out.trkPi_mothpdg);
    int cls = 0, ti = -1, m1 = -1, m2 = -1;
    if (matchPair(td.iK, td.iPi, L1, L2, ti, m1, m2)) {
      cls = (td.iK[ti] == m1) ? 1 : 2;
      setLegTruth(mc, m1, out.trkK_mcpdg, out.trkK_mothpdg);
      setLegTruth(mc, m2, out.trkPi_mcpdg, out.trkPi_mothpdg);
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
  const int nmc = mc.size();
  for (size_t k = 0; k < c.kin.m_kpi.size(); ++k) {
    const RVec<int> L1 = mcLinksOfTrack(trackToMCs, c.trkK.origIdx[k], nmc);
    const RVec<int> L2 = mcLinksOfTrack(trackToMCs, c.trkPi.origIdx[k], nmc);
    const RVec<int> L3 = mcLinksOfTrack(trackToMCs, c.trkPis.origIdx[k], nmc);
    pushLegTruth(mc, firstOr(L1), out.trkK_mcpdg, out.trkK_mothpdg);
    pushLegTruth(mc, firstOr(L2), out.trkPi_mcpdg, out.trkPi_mothpdg);
    pushLegTruth(mc, firstOr(L3), out.trkPis_mcpdg, out.trkPis_mothpdg);

    int cls = 0, ti = -1, m1 = -1, m2 = -1;
    // true D* whose D0 daughters are this K,pi pair (either assignment)
    if (matchPair(tds.iK, tds.iPi, L1, L2, ti, m1, m2)) {
      const bool swapped = (tds.iK[ti] != m1);
      int m3 = -1;
      for (int cnd : L3) if (cnd == tds.iPis[ti]) { m3 = cnd; break; }
      const bool pis_ok = (m3 >= 0);
      if (!swapped) cls = pis_ok ? 1 : 3;
      else if (pis_ok) cls = 2;
      if (cls) {
        setLegTruth(mc, m1, out.trkK_mcpdg, out.trkK_mothpdg);
        setLegTruth(mc, m2, out.trkPi_mcpdg, out.trkPi_mothpdg);
        if (pis_ok) setLegTruth(mc, m3, out.trkPis_mcpdg, out.trkPis_mothpdg);
      }
    } else {
      // no D* parent: a true D0 -> K pi with the correct assignment
      for (size_t t = 0; t < td0.idx.size() && cls == 0; ++t)
        for (int a : L1) {
          if (td0.iK[t] != a) continue;
          for (int b : L2)
            if (td0.iPi[t] == b) {
              cls = 4;
              setLegTruth(mc, a, out.trkK_mcpdg, out.trkK_mothpdg);
              setLegTruth(mc, b, out.trkPi_mcpdg, out.trkPi_mothpdg);
              break;
            }
          if (cls) break;
        }
    }
    out.cls.push_back(cls);
    // a swapped K,pi pair whose slow pion is wrong stays class 0, so the
    // back-pointer is only meaningful for classes 1-3
    out.trueidx.push_back(cls > 0 ? ti : -1);
  }
  return out;
}

// 1 if the true D* / D0 at position t was reconstructed by at least one
// class-1 candidate carrying `flag` (the loose or tight label).
inline RVec<int> trueDstarFound(const TrueDstars& tds, const DstarTruth& info,
                                const RVec<int>& flag) {
  return foundFlags(tds.idx.size(), info.trueidx, info.cls, flag);
}

inline RVec<int> trueD0Found(const TrueD0s& td, const D0Truth& info,
                             const RVec<int>& flag) {
  return foundFlags(td.idx.size(), info.trueidx, info.cls, flag);
}

} // namespace AlephDstar
} // namespace FCCAnalyses

#endif
