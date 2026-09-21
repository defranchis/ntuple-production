#ifndef ALEPH_CHARGEDPF_H
#define ALEPH_CHARGEDPF_H

#include "analyzer.h"

namespace FCCAnalyses { namespace AlephChargedPF {

namespace rv = ROOT::VecOps;

/// Value of the truth branches when no MC information is available.
constexpr int kNoTruthInt = -999;
constexpr float kNoTruthFloat = -999.f;

/// Indices of the charged ReconstructedParticles owning a track, as a single group.
inline std::vector<std::vector<int>>
select_indices(const rv::RVec<edm4hep::ReconstructedParticleData> &rps,
               const rv::RVec<int> &rpTrackIndex) {
  std::vector<std::vector<int>> out(1);
  for (size_t i = 0; i < rps.size(); ++i) {
    if (rps[i].charge != 0 &&
        AlephSelection::hasOwnTrack(rps[i], rpTrackIndex.size()))
      out[0].push_back(static_cast<int>(i));
  }
  return out;
}

/// First group of a per-group quantity, as a flat vector.
template <typename T>
inline rv::RVec<T> flatten(const rv::RVec<rv::RVec<T>> &v) {
  return v.empty() ? rv::RVec<T>{} : v[0];
}

/// Index into the Tracks collection of each candidate's own track.
inline rv::RVec<int>
get_trackIndex(const rv::RVec<edm4hep::ReconstructedParticleData> &sel,
               const rv::RVec<int> &rpTrackIndex) {
  rv::RVec<int> out;
  out.reserve(sel.size());
  for (const auto &rp : sel) out.push_back(rpTrackIndex.at(rp.tracks_begin));
  return out;
}

/// Perigee parameter of each candidate's own track state: 0 = phi0, 1 = omega, 2 = tanLambda.
inline rv::RVec<float>
get_trackStateParam(const rv::RVec<edm4hep::ReconstructedParticleData> &sel,
                    const rv::RVec<edm4hep::TrackState> &statesByRP, int which) {
  rv::RVec<float> out;
  out.reserve(sel.size());
  for (const auto &rp : sel) {
    const auto &ts = statesByRP.at(rp.tracks_begin);
    out.push_back(which == 0 ? ts.phi : (which == 1 ? ts.omega : ts.tanLambda));
  }
  return out;
}

/// Vector of n identical entries, for branches that carry no information on data.
inline rv::RVec<float> constant(size_t n, float value) {
  return rv::RVec<float>(n, value);
}

/// Integer-valued counterpart of constant().
inline rv::RVec<int> constant_int(size_t n, int value) {
  return rv::RVec<int>(n, value);
}

/// Truth quantities of the first track->MCParticle link of each candidate's track.
struct MCTruth {
  rv::RVec<int> pdg, parent_pdg, n_links;
  rv::RVec<float> p, vx, vy, vz;
};

/// Resolves trackMCLink (from = Track index, to = MCParticle index); kNoTruth* / 0 if unlinked.
inline MCTruth
get_mcTruth(const rv::RVec<int> &trackIndex, const rv::RVec<int> &linkFrom,
            const rv::RVec<int> &linkTo,
            const rv::RVec<edm4hep::MCParticleData> &mcp,
            const rv::RVec<int> &mcParents) {
  std::unordered_map<int, int> first_link;
  std::unordered_map<int, int> n_link;
  for (size_t i = 0; i < linkFrom.size() && i < linkTo.size(); ++i) {
    const int t = linkFrom[i];
    if (!first_link.count(t)) first_link[t] = static_cast<int>(i);
    ++n_link[t];
  }

  MCTruth out;
  for (const int t : trackIndex) {
    const auto it = first_link.find(t);
    if (it == first_link.end()) {
      out.pdg.push_back(kNoTruthInt); out.parent_pdg.push_back(kNoTruthInt);
      out.p.push_back(kNoTruthFloat);
      out.vx.push_back(kNoTruthFloat); out.vy.push_back(kNoTruthFloat); out.vz.push_back(kNoTruthFloat);
      out.n_links.push_back(0);
      continue;
    }
    const int imc = linkTo.at(it->second);
    if (imc < 0 || imc >= static_cast<int>(mcp.size()))
      throw std::runtime_error("get_mcTruth: trackMCLink target out of range");
    const auto &m = mcp[imc];
    out.pdg.push_back(m.PDG);
    out.p.push_back(static_cast<float>(std::sqrt(m.momentum.x * m.momentum.x +
                                                 m.momentum.y * m.momentum.y +
                                                 m.momentum.z * m.momentum.z)));
    out.vx.push_back(static_cast<float>(m.vertex.x));
    out.vy.push_back(static_cast<float>(m.vertex.y));
    out.vz.push_back(static_cast<float>(m.vertex.z));
    out.n_links.push_back(n_link[t]);

    int parent = kNoTruthInt;
    if (m.parents_begin != m.parents_end &&
        m.parents_begin < mcParents.size()) {
      const int ip = mcParents[m.parents_begin];
      if (ip >= 0 && ip < static_cast<int>(mcp.size()))
        parent = mcp[ip].PDG;
    }
    out.parent_pdg.push_back(parent);
  }
  return out;
}

}} // namespace FCCAnalyses::AlephChargedPF

#endif
