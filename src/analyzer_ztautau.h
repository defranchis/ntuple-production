#ifndef ANALYZER_ZTAUTAU_H
#define ANALYZER_ZTAUTAU_H

// Helpers of the Z -> tau tau (1-prong + 3-prong) selection, stage1_Ztautau.py.
// An argument that is an array of arrays holds one inner array per jet of the
// exclusive two-jet clustering; the flag arrays are AlephSelection::get_isType
// outputs.

#include <cmath>
#include <unordered_map>
#include <vector>

#include "ROOT/RVec.hxx"
#include "TLorentzVector.h"
#include "edm4hep/MCParticleData.h"
#include "edm4hep/ReconstructedParticleData.h"
#include "podio/ObjectID.h"

#include "analyzer.h"
#include "analyzer_alephvars.h"
#include "analyzer_pvnew.h"

namespace FCCAnalyses { namespace MCTruthMatching{

struct RPtoMCLinks{
    ROOT::VecOps::RVec<int> recind;
    ROOT::VecOps::RVec<int> mcind; 
};
inline RPtoMCLinks build_RP2MC_via_track(
    const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& reco,
    const ROOT::VecOps::RVec<podio::ObjectID>& rp_track_idx,
    const ROOT::VecOps::RVec<podio::ObjectID>& link_from,
    const ROOT::VecOps::RVec<podio::ObjectID>& link_to,
    const ROOT::VecOps::RVec<float>& link_weight) {
    
    // 1. Map track index to the highest-weight MC particle index
    std::unordered_map<int, std::pair<int, float>> track_to_mc;
    for (size_t i = 0; i < link_from.size(); ++i) {
        int trk = link_from[i].index;
        int mc = link_to[i].index;
        float w = link_weight[i];
        
        auto it = track_to_mc.find(trk);
        if (it == track_to_mc.end() || w > it->second.second) {
            track_to_mc[trk] = {mc, w};
        }
    }

    // 2. Loop through Reconstructed Particles and perform the matching
    RPtoMCLinks out; // Allocates an empty instance of custom structure (out) to start storing the final successful reco-to-MC matches.
    for (size_t rp = 0; rp < reco.size(); ++rp) {
        int t_begin = reco[rp].tracks_begin;
        int t_end = reco[rp].tracks_end;
        
        for (int t = t_begin; t < t_end; ++t) {
            int trk = rp_track_idx[t].index;
            auto it = track_to_mc.find(trk);
            
            if (it != track_to_mc.end()) {
                out.recind.push_back(rp);
                out.mcind.push_back(it->second.first);
            }
        }
    }
    return out;
}

} 
}

namespace FCCAnalyses { namespace AlephSelection {

inline ROOT::VecOps::RVec<FCCAnalysesJetConstituentsData> get_MCindex(
  const ROOT::VecOps::RVec<int> &recind,
  const ROOT::VecOps::RVec<int> &mcind,
  const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> &RecPart,
  const std::vector<std::vector<int>> &indices
)

{
  ROOT::VecOps::RVec<int> mc_index(RecPart.size(), -1);
  for (size_t k=0; k<recind.size(); ++k)
    mc_index[recind[k]] = mcind[k];

  ROOT::VecOps::RVec<FCCAnalysesJetConstituentsData> out;
    for (const auto &jet_index : indices) {
      FCCAnalysesJetConstituentsData tmp;
      for (const auto &const_index : jet_index)
        tmp.push_back(static_cast<float>(mc_index.at(const_index)));
      out.push_back(tmp);
    }
    return out;
  }
inline ROOT::VecOps::RVec<FCCAnalysesJetConstituentsData> get_true_p_from_index(
    const ROOT::VecOps::RVec<FCCAnalysesJetConstituentsData> &mcIndexNested,
    const ROOT::VecOps::RVec<edm4hep::MCParticleData> &mcParticle)
  
  {
    ROOT::VecOps::RVec<FCCAnalysesJetConstituentsData> out;
    for (const auto &jet_idx : mcIndexNested) {
      FCCAnalysesJetConstituentsData tmp;
      for (float idxf : jet_idx){
        int idx = static_cast<int>(idxf);
        if (idx>=0){
          const auto &p = mcParticle.at(idx).momentum;
          tmp.push_back(std::sqrt(p.x*p.x + p.y*p.y + p.z*p.z));
        } else {
          tmp.push_back(-9999.f);
        }
      }
      out.push_back(tmp);
    }
    return out;
  }
/*
// Charged multiplicity per jet -> tau prong count.
ROOT::VecOps::RVec<int> count_charged(
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<bool>> isChargedHad,
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<bool>> isMu,
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<bool>> isEl) {
  ROOT::VecOps::RVec<int> result;
  for (size_t i = 0; i < isChargedHad.size(); ++i) {
    int n = 0;
    for (size_t j = 0; j < isChargedHad.at(i).size(); ++j) {
      if (isChargedHad.at(i).at(j) || isMu.at(i).at(j) || isEl.at(i).at(j)) ++n;
    }
    result.push_back(n);
  }
  return result;
}
*/
// Charged multiplicity per jet -> tau prong count.
ROOT::VecOps::RVec<int> count_charged(
    const ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> &pfcand_charge) {
  ROOT::VecOps::RVec<int> result;
  result.reserve(pfcand_charge.size());
  for (const auto &jet : pfcand_charge) {
    result.push_back(static_cast<int>(
        std::count_if(jet.begin(), jet.end(), [](int q) { return q != 0; })));
  }
  return result;
}

// True if every charged constituent in the jet is a charged hadron (none
// flagged isMu/isEl) - vetoes tau->e/mu nu-nu; it does not cut on dE/dx
// directly, but ALEPH's electron identification (isEl) uses dE/dx.
ROOT::VecOps::RVec<bool> all_charged_are_hadronic(
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<bool>> isChargedHad,
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<bool>> isMu,
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<bool>> isEl) {
  ROOT::VecOps::RVec<bool> result;
  for (size_t i = 0; i < isChargedHad.size(); ++i) {
    bool ok = true;
    for (size_t j = 0; j < isChargedHad.at(i).size(); ++j) {
      bool charged = isChargedHad.at(i).at(j) || isMu.at(i).at(j) || isEl.at(i).at(j);
      if (charged && !isChargedHad.at(i).at(j)) ok = false;  // charged but not hadronic -> lepton
    }
    result.push_back(ok);
  }
  return result;
}

// Total constituent count per jet (charged + neutral) - compared against
// jet_nCharged to veto extra pi0/photons from rho/a1 decays.
ROOT::VecOps::RVec<int> count_total(
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<bool>> isChargedHad,
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<bool>> isMu,
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<bool>> isEl,
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<bool>> isGamma,
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<bool>> isNeutralHad) {
  ROOT::VecOps::RVec<int> result;
  for (size_t i = 0; i < isChargedHad.size(); ++i) {
    result.push_back(static_cast<int>(isChargedHad.at(i).size()));
  }
  return result;
}
ROOT::VecOps::RVec<float> max_abs_per_jet_charged(
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> values,
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<bool>> isChargedHad,
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<bool>> isMu,
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<bool>> isEl) {
  ROOT::VecOps::RVec<float> result;
  for (size_t i = 0; i < values.size(); ++i) {
    float m = 0;
    for (size_t j = 0; j < values.at(i).size(); ++j) {
      bool charged = isChargedHad.at(i).at(j) || isMu.at(i).at(j) || isEl.at(i).at(j);
      if (!charged) continue;
      float v = values.at(i).at(j);
      if (v == -999) continue;
      if (std::abs(v) > m) m = std::abs(v);
    }
    result.push_back(m);
  }
  return result;
}
// Check that all constituents in a jet (given by their MC indices, as
// float per get_MCindex's convention) share the same production vertex,
// coordinate-by-coordinate, within `tol`. Returns false if fewer than
// 2 valid MC matches exist (nothing to compare).
bool jet_same_vertex(const ROOT::VecOps::RVec<float>& jet_mc_index,
                      const ROOT::VecOps::RVec<float>& MC_vertex_x,
                      const ROOT::VecOps::RVec<float>& MC_vertex_y,
                      const ROOT::VecOps::RVec<float>& MC_vertex_z,
                      float tol=1e-6){
    std::vector<int> idx;
    for (float f: jet_mc_index){
      int i = static_cast<int>(f);
      if (i>=0) idx.push_back(i);
    }
    if (idx.size()<2) return false; //nothing to compare->can't confirm

    for (size_t a = 0; a < idx.size(); ++a){
      for (size_t b = a+1; b < idx.size(); ++b){
        float dx = std::abs(MC_vertex_x[idx[a]]-MC_vertex_x[idx[b]]);
        float dy = std::abs(MC_vertex_y[idx[a]]-MC_vertex_y[idx[b]]);
        float dz = std::abs(MC_vertex_z[idx[a]]-MC_vertex_z[idx[b]]);
        if (dx > tol || dy > tol || dz > tol) return false;
      }
    }
    return true;
}

// Given jet3's matched pion MC indices, find the tau neutrino that shares
// their production vertex (i.e. came from the same tau decay), within tol.
// Returns -1 if none found.
int find_nu_tau_same_vertex(const ROOT::VecOps::RVec<float>& jet_mc_index,
                             const ROOT::VecOps::RVec<edm4hep::MCParticleData>& particles,
                             const ROOT::VecOps::RVec<float>& MC_vertex_x,
                             const ROOT::VecOps::RVec<float>& MC_vertex_y,
                             const ROOT::VecOps::RVec<float>& MC_vertex_z,
                             float tol = 1e-3f) {
    // reference vertex = first valid matched pion's vertex
    int ref = -1;
    for (float f : jet_mc_index) { int i=(int)f; if (i>=0) { ref=i; break; } }
    if (ref < 0) return -1;
    float rx = MC_vertex_x[ref], ry = MC_vertex_y[ref], rz = MC_vertex_z[ref];

    for (size_t i = 0; i < particles.size(); ++i) {
        if (std::abs(particles[i].PDG) != 16) continue;   // nu_tau / anti-nu_tau
        float dx = std::abs(MC_vertex_x[i]-rx);
        float dy = std::abs(MC_vertex_y[i]-ry);
        float dz = std::abs(MC_vertex_z[i]-rz);
        if (dx<tol && dy<tol && dz<tol) return (int)i;
    }
    return -1;
}

// Build the invariant mass of the 3 truth pions + the matched nu_tau,
// i.e. the reconstructed truth-level tau mass. Returns -1 if incomplete.
float truth_tau_mass_3pi_nu(const ROOT::VecOps::RVec<float>& jet_mc_index,
                             const ROOT::VecOps::RVec<edm4hep::MCParticleData>& particles,
                             const ROOT::VecOps::RVec<float>& MC_vertex_x,
                             const ROOT::VecOps::RVec<float>& MC_vertex_y,
                             const ROOT::VecOps::RVec<float>& MC_vertex_z,
                             float tol = 1e-3f) {
    int nu_idx = find_nu_tau_same_vertex(jet_mc_index, particles, MC_vertex_x, MC_vertex_y, MC_vertex_z, tol);
    if (nu_idx < 0) return -1.0f;

    double E=0, px=0, py=0, pz=0;
    int n_pi = 0;
    for (float f : jet_mc_index) {
        int i = (int)f;
        if (i < 0) continue;
        const auto& p = particles[i];
        double m = p.mass;
        double ppx=p.momentum.x, ppy=p.momentum.y, ppz=p.momentum.z;
        double e = std::sqrt(ppx*ppx+ppy*ppy+ppz*ppz+m*m);
        E+=e; px+=ppx; py+=ppy; pz+=ppz;
        n_pi++;
    }
    if (n_pi < 3) return -1.0f;   // require all 3 pions matched

    const auto& nu = particles[nu_idx];
    double npx=nu.momentum.x, npy=nu.momentum.y, npz=nu.momentum.z;
    double ne = std::sqrt(npx*npx+npy*npy+npz*npz);  // neutrino treated massless
    E+=ne; px+=npx; py+=npy; pz+=npz;

    double m2 = E*E - px*px - py*py - pz*pz;
    return m2 > 0 ? (float)std::sqrt(m2) : -1.0f;
}
// Per-jet min TPC hit count, restricted to charged constituents only
// (neutrals have no track and are ignored, not treated as failing the cut).
ROOT::VecOps::RVec<int> min_TPChits_charged(
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> nTrackHits_TPC,
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<bool>> isChargedHad,
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<bool>> isMu,
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<bool>> isEl) {
  ROOT::VecOps::RVec<int> result;
  for (size_t i = 0; i < nTrackHits_TPC.size(); ++i) {
    int m = -1;
    for (size_t j = 0; j < nTrackHits_TPC.at(i).size(); ++j) {
      bool charged = isChargedHad.at(i).at(j) || isMu.at(i).at(j) || isEl.at(i).at(j);
      if (!charged) continue;
      int hits = nTrackHits_TPC.at(i).at(j);
      if (m == -1 || hits < m) m = hits;
    }
    result.push_back(m);
  }
  return result;
}

// Mass of each jet's summed four-vector (the "visible tau mass": includes
// any neutrals - e.g. bremsstrahlung photons - that fell in the jet).
// jet_p4 is the RVec<TLorentzVector>-like object from
// JetConstituentsUtils::compute_tlv_jets; adjust the .M() call if your
// build uses ROOT::Math::PxPyPzEVector instead (same interface).
ROOT::VecOps::RVec<float> get_masses(
    ROOT::VecOps::RVec<TLorentzVector> jet_p4) {
  ROOT::VecOps::RVec<float> result;
  for (size_t i = 0; i < jet_p4.size(); ++i) {
    result.push_back(jet_p4.at(i).M());
  }
  return result;
}

// Invariant mass built from the charged constituents only, per jet - i.e.
// the mass of the pion system itself (trivial = pion mass for 1-prong;
// the 3-pion system mass, e.g. probing the a1 resonance, for 3-prong).
ROOT::VecOps::RVec<float> get_charged_invmass(
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> px,
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> py,
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> pz,
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> e,
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<bool>> isChargedHad,
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<bool>> isMu,
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<bool>> isEl) {
  ROOT::VecOps::RVec<float> result;
  for (size_t i = 0; i < px.size(); ++i) {
    double sumPx = 0, sumPy = 0, sumPz = 0, sumE = 0;
    for (size_t j = 0; j < px.at(i).size(); ++j) {
      bool charged = isChargedHad.at(i).at(j) || isMu.at(i).at(j) || isEl.at(i).at(j);
      if (!charged) continue;
      sumPx += px.at(i).at(j);
      sumPy += py.at(i).at(j);
      sumPz += pz.at(i).at(j);
      sumE  += e.at(i).at(j);
    }
    double m2 = sumE * sumE - sumPx * sumPx - sumPy * sumPy - sumPz * sumPz;
    result.push_back(m2 > 0 ? std::sqrt(m2) : -std::sqrt(-m2));
  }
  return result;
}

// Common vertex of the charged constituents (charge != 0, with a track) of jet
// `jet`: one fit with stage1's standalone vertex fitter (analyzer_pvnew.h), no
// beam-spot constraint and no track removal. `states` = ALEPH-flipped track
// states ordered by the RecoParticle -> Track relation (TrackStateByRP).
// converged == false: the values must not be used (status says why).
inline AlephPVNew::PVFitResult fit_jet_vertex(
    const ROOT::VecOps::RVec<FCCAnalysesJetConstituents> &jcs, int jet,
    const ROOT::VecOps::RVec<edm4hep::TrackState> &states) {
  ROOT::VecOps::RVec<edm4hep::TrackState> trk;
  for (const auto &rp : jcs.at(jet))
    if (rp.charge != 0 && hasOwnTrack(rp, states.size()))
      trk.push_back(states[rp.tracks_begin]);
  return AlephPVNew::detail::fit_core(AlephPVNew::detail::convert(trk), nullptr,
                                      AlephPVNew::FitConfig{});
}

// The three 2-track vertices of the charged constituents (charge != 0, with a
// track) c0, c1, c2 of jet `jet`, in the pair order (c0,c1), (c0,c2), (c1,c2):
// each pair fitted alone with the fitter of fit_jet_vertex (no beam spot).
// leg1/leg2 = the legs' indices in the jet's constituent list (-1: the jet has
// fewer charged constituents). Momenta at the vertex: the helix tangent at each
// track's phase of closest approach to the fitted vertex, pT = kPtPerTeslaCm *
// bz / |omega|; mee = pair mass with the electron mass for both legs
// (AlephVars::pairMass), qt = Armenteros qT (AlephTrkAux::apVars). Undefined:
// charge -999 without a pair; chi2, mee, qt -1 and position -999 when the fit
// did not converge.
struct JetPairVertices {
  ROOT::VecOps::RVec<int> leg1, leg2, converged, status, charge;
  ROOT::VecOps::RVec<float> chi2, x, y, z, mee, qt;
};

inline JetPairVertices fit_jet_pair_vertices(
    const ROOT::VecOps::RVec<FCCAnalysesJetConstituents> &jcs, int jet,
    const ROOT::VecOps::RVec<edm4hep::TrackState> &states, float bz) {
  namespace D = AlephPVNew::detail;
  const auto &jc = jcs.at(jet);
  std::vector<int> charged;
  for (int k = 0; k < (int)jc.size(); ++k)
    if (jc[k].charge != 0 && hasOwnTrack(jc[k], states.size())) charged.push_back(k);
  const AlephPVNew::FitConfig cfg{};
  constexpr int kPairs[3][2] = {{0, 1}, {0, 2}, {1, 2}};
  JetPairVertices out;
  for (const auto &pr : kPairs) {
    const bool have = (int)charged.size() > pr[1];
    const int a = have ? charged[pr[0]] : -1, b = have ? charged[pr[1]] : -1;
    out.leg1.push_back(a);
    out.leg2.push_back(b);
    out.charge.push_back(have ? (int)std::lround(jc[a].charge + jc[b].charge)
                              : AlephVars::kNotAvailableInt);
    AlephPVNew::PVFitResult f;  // default: not converged, kTooFewTracks
    ROOT::VecOps::RVec<edm4hep::TrackState> trk;
    if (have) {
      trk = {states[jc[a].tracks_begin], states[jc[b].tracks_begin]};
      f = D::fit_core(D::convert(trk), nullptr, cfg);
    }
    out.converged.push_back(f.converged ? 1 : 0);
    out.status.push_back(f.status);
    if (!f.converged) {
      out.chi2.push_back(-1.f);
      out.x.push_back(-999.f);
      out.y.push_back(-999.f);
      out.z.push_back(-999.f);
      out.mee.push_back(-1.f);
      out.qt.push_back(-1.f);
      continue;
    }
    out.chi2.push_back(f.chi2);
    out.x.push_back(f.position[0]);
    out.y.push_back(f.position[1]);
    out.z.push_back(f.position[2]);
    // phases of closest approach to the fitted vertex (the fit's own phase refinement)
    const D::TrackSet ts = D::convert(trk);
    const AlephPVNew::Vec3 xv(f.position[0], f.position[1], f.position[2]);
    std::vector<double> L(2, 0.0);
    D::TrackTerms tt;
    for (int r = 0; r < 3; ++r) {
      D::track_terms(ts, xv, L, cfg, tt);
      L = tt.L;
    }
    TVector3 p[2];
    for (int i = 0; i < 2; ++i) {
      const AlephPVNew::Vec3 t = D::helix_dXdL(ts.par[i], L[i]);  // (cos phi, sin phi, tan lambda)
      const double pt = AlephUnits::kPtPerTeslaCm * bz / std::abs(trk[i].omega);
      p[i].SetXYZ(pt * t(0), pt * t(1), pt * t(2));
    }
    out.mee.push_back(AlephVars::pairMass(p[0], AlephMasses::kElectron, p[1], AlephMasses::kElectron));
    double alpha, qt;
    AlephTrkAux::apVars(p[0], p[1], jc[a].charge > 0 ? 1. : -1., alpha, qt);
    out.qt.push_back(qt);
  }
  return out;
}

// values[leg] per pair (undefined for leg -1): the per-constituent value of one
// leg of each pair, from a per-constituent array of the same jet.
template <typename T>
inline ROOT::VecOps::RVec<T> pair_leg_value(const ROOT::VecOps::RVec<T> &values,
                                            const ROOT::VecOps::RVec<int> &legs, T undefined) {
  ROOT::VecOps::RVec<T> out;
  for (int l : legs) out.push_back(l >= 0 && l < (int)values.size() ? values[l] : undefined);
  return out;
}

// 1 = both legs' best-match particles are distinct e+- with the same first
// parent, a photon (pdg 22); per pair, from the pair_leg_value truth arrays.
inline ROOT::VecOps::RVec<int> pair_truth_conversion(
    const ROOT::VecOps::RVec<int> &pdg_a, const ROOT::VecOps::RVec<int> &pdg_b,
    const ROOT::VecOps::RVec<float> &idx_a, const ROOT::VecOps::RVec<float> &idx_b,
    const ROOT::VecOps::RVec<int> &parent_pdg_a,
    const ROOT::VecOps::RVec<float> &parent_idx_a,
    const ROOT::VecOps::RVec<float> &parent_idx_b) {
  ROOT::VecOps::RVec<int> out;
  for (size_t k = 0; k < pdg_a.size(); ++k)
    out.push_back(std::abs(pdg_a[k]) == 11 && std::abs(pdg_b[k]) == 11 && idx_a[k] >= 0 &&
                  idx_b[k] >= 0 && idx_a[k] != idx_b[k] && parent_idx_a[k] >= 0 &&
                  parent_idx_a[k] == parent_idx_b[k] && parent_pdg_a[k] == 22);
  return out;
}

}} // namespace FCCAnalyses::AlephSelection

#endif  // ANALYZER_ZTAUTAU_H
