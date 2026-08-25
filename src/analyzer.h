#ifndef SELECTIONUTILS_H
#define SELECTIONUTILS_H

/*
  Selection utilities for filtering particles and events.
  Fully RDataFrame compatible (using ROOT::VecOps::RVec).

  Includes:
    - sel_charged: selects reconstructed particles by absolute charge.
    - sel_class_filter: filters events based on their class bit.
    - sel_runs_filter: filters events by allowed run numbers.
    - get_isEl / get_isMu / get_isChargedHad / get_isNeutralHad / get_isGamma:
      classify jet constituents by particle type.

  Example usage in RDataFrame:

    df = df.Define("charged_particles",
                   "FCCAnalyses::AlephSelection::sel_charged(1)(ReconstructedParticles)")
           .Filter("FCCAnalyses::AlephSelection::sel_class_filter(16)(ClassBitset)")
           .Filter("FCCAnalyses::AlephSelection::sel_runs_filter(allowedRuns)(EventHeader)");

    df = df.Define("isMu", "FCCAnalyses::AlephSelection::get_isMu(JetConstituents)")
           .Define("n_muons_per_jet", "Sum(isMu)");
*/
#include "dedx_valid.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/EventHeaderCollection.h"
#include <set>
#include <bitset>
#include <cmath>
#include <vector>
#include <map>
#include <mutex>
#include <fstream>
#include <string>
#include <cstdlib>
#include <nlohmann/json.hpp>
#include <ROOT/RVec.hxx>

#include "FCCAnalyses/JetConstituentsUtils.h"
#include "FCCAnalyses/ReconstructedParticle.h"
#include "FCCAnalyses/ReconstructedParticle2Track.h"
#include "FCCAnalyses/ReconstructedParticle2MC.h"
#include "FCCAnalyses/JetClusteringUtils.h"
// #include "FCCAnalyses/ExternalRecombiner.h"
#include "FCCAnalyses/MCParticle.h"
#include "FCCAnalyses/VertexingUtils.h"
#include "FCCAnalyses/VertexFinderLCFIPlus.h" 
#include "aleph_units.h"

#include "TVector3.h"

#include "edm4hep/MCParticleData.h"
#include "edm4hep/Track.h"
#include "edm4hep/TrackData.h"
#include "edm4hep/Cluster.h"
#include "edm4hep/ClusterData.h"
#include "edm4hep/CalorimeterHitData.h"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/EDM4hepVersion.h"
#include "edm4hep/RecDqdx.h"

#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"

#include <iostream>
#include <algorithm>
#include <typeinfo>


namespace FCCAnalyses { namespace AlephSelection {

namespace rv = ROOT::VecOps;

// -----------------------------------
// Type aliases for jet constituent data
// -----------------------------------

using FCCAnalysesJetConstituents = rv::RVec<edm4hep::ReconstructedParticleData>;
using FCCAnalysesJetConstituentsData = rv::RVec<float>;

// ----------------------------
// Event and particle selectors
// ----------------------------

//////////////////////////////////////////////////////////////////////////////////////////
// -----------------------------------
// The following will first filter the 
// event to check if it's a qq event
// then it return the PID of qq 
// -----------------------------------

float getJetPID(const ROOT::VecOps::RVec<uint32_t>& ClassBit,
                   const ROOT::VecOps::RVec<edm4hep::MCParticleData>& particles) {
    // Check if bit 15 (16-1) is set in the first element of ClassBit
    if (ClassBit.empty() || !std::bitset<32>(ClassBit[0])[15]) {
        return -1.0f; // Not selected
    }

    // Bit 15 is true — find the first quark (|PDG| in 1..5)
    float result = -1.0f;
    for (const auto& particle : particles) {
        if (std::abs(particle.PDG) > 0 && std::abs(particle.PDG) < 6) {
            result = static_cast<float>(std::abs(particle.PDG));
            break;
        }
    }

    return result;
}

/// Selects charged particles based on their absolute charge.
struct sel_charged {
  const int m_charge;
  sel_charged(int arg_charge) : m_charge(arg_charge) {};

  edm4hep::ReconstructedParticleCollection
  operator()(const edm4hep::ReconstructedParticleCollection& in_coll) const {
    edm4hep::ReconstructedParticleCollection result;
    result.setSubsetCollection();

    for (const auto& i : in_coll) {
      if (std::abs(i.getCharge()) == m_charge) {
        result.push_back(i);
      }
    }
    return result;
  }
};

/// Filters events based on their class bit (RVec-compatible)
struct sel_class_filter {
  const int m_class;
  sel_class_filter(int arg_class) : m_class(arg_class) {};

  bool operator()(const ROOT::VecOps::RVec<uint32_t>& bitset_coll) const {
    if (bitset_coll.empty()) return false;
    std::bitset<32> bits(bitset_coll[0]);
    return bits[m_class - 1];
  }
};

// create a vector of all classes that the event is in (i.e. all bits that are true for the event)
std::vector<int> bitsetToIndices(const ROOT::VecOps::RVec<uint32_t>& bitset_coll) {
    std::vector<int> indices;
    std::bitset<32> bits(bitset_coll[0]);

    for (size_t i = 0; i < bits.size(); ++i) {
        if (bits.test(i)) {  // check if bit i is set
            indices.push_back(static_cast<int>(i) + 1); //class counting starts at 1, but bit indices at 0
        }
    }
    return indices;
}

/// Filters events by run number (RVec-compatible)
struct sel_runs_filter {
  const std::set<int>& m_runs_set;
  sel_runs_filter(const std::set<int>& arg_runs_set) : m_runs_set(arg_runs_set) {};

  bool operator()(const ROOT::VecOps::RVec<edm4hep::EventHeader>& event_header) const {
    if (event_header.empty()) return false;
    return m_runs_set.count(event_header[0].getRunNumber()) > 0;
  }
};

// --------------------------------------
// Jet constituent particle identification
// --------------------------------------

rv::RVec<FCCAnalysesJetConstituentsData>
get_isEl(const rv::RVec<FCCAnalysesJetConstituents>& jcs) {
  rv::RVec<FCCAnalysesJetConstituentsData> out;
  out.reserve(jcs.size());
  for (const auto& jet : jcs) {
    FCCAnalysesJetConstituentsData mask;
    mask.reserve(jet.size());
    for (const auto& c : jet)
      mask.push_back((std::abs(c.charge) > 0 && std::abs(c.mass - 0.000511) < 1e-5) ? 1.f : 0.f);
    out.push_back(std::move(mask));
  }
  return out;
}

rv::RVec<FCCAnalysesJetConstituentsData>
get_isMu(const rv::RVec<FCCAnalysesJetConstituents>& jcs) {
  rv::RVec<FCCAnalysesJetConstituentsData> out;
  out.reserve(jcs.size());
  for (const auto& jet : jcs) {
    FCCAnalysesJetConstituentsData mask;
    mask.reserve(jet.size());
    for (const auto& c : jet)
      mask.push_back((std::abs(c.charge) > 0 && std::abs(c.mass - 0.105658) < 1e-3) ? 1.f : 0.f);
    out.push_back(std::move(mask));
  }
  return out;
}

rv::RVec<FCCAnalysesJetConstituentsData>
get_isChargedHad(const rv::RVec<FCCAnalysesJetConstituents>& jcs) {
  rv::RVec<FCCAnalysesJetConstituentsData> out;
  out.reserve(jcs.size());
  for (const auto& jet : jcs) {
    FCCAnalysesJetConstituentsData mask;
    mask.reserve(jet.size());
    for (const auto& c : jet)
      mask.push_back((std::abs(c.charge) > 0 && std::abs(c.mass - 0.13957) < 1e-3) ? 1.f : 0.f);
    out.push_back(std::move(mask));
  }
  return out;
}

rv::RVec<FCCAnalysesJetConstituentsData>
get_isNeutralHad(const rv::RVec<FCCAnalysesJetConstituents>& jcs) {
  rv::RVec<FCCAnalysesJetConstituentsData> out;
  out.reserve(jcs.size());
  for (const auto& jet : jcs) {
    FCCAnalysesJetConstituentsData mask;
    mask.reserve(jet.size());
    for (const auto& c : jet)
#if edm4hep_VERSION > EDM4HEP_VERSION(0, 10, 5)
      mask.push_back((c.PDG == 130) ? 1.f : 0.f);
#else
      mask.push_back((c.type == 130) ? 1.f : 0.f);
#endif
    out.push_back(std::move(mask));
  }
  return out;
}

rv::RVec<FCCAnalysesJetConstituentsData>
get_isGamma(const rv::RVec<FCCAnalysesJetConstituents>& jcs) {
  rv::RVec<FCCAnalysesJetConstituentsData> out;
  out.reserve(jcs.size());
  for (const auto& jet : jcs) {
    FCCAnalysesJetConstituentsData mask;
    mask.reserve(jet.size());
    for (const auto& c : jet)
#if edm4hep_VERSION > EDM4HEP_VERSION(0, 10, 5)
      mask.push_back((c.PDG == 22) ? 1.f : 0.f);
#else
      mask.push_back((c.type == 22) ? 1.f : 0.f);
#endif
    out.push_back(std::move(mask));
  }
  return out;
}

// --------------------------------------
// Track helpers & selection
// --------------------------------------


// Getters for track fit quality vars
ROOT::VecOps::RVec<float>
get_track_chi2(const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks_in){
  ROOT::VecOps::RVec<float> chi2_out;
  for (const auto &track : tracks_in) {
    chi2_out.push_back(track.chi2);
  }
  return chi2_out;
}

ROOT::VecOps::RVec<float>
get_track_ndf(const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks_in){
  ROOT::VecOps::RVec<float> ndf_out;
  for (const auto &track : tracks_in) {
    ndf_out.push_back(track.ndf);
  }
  return ndf_out;
}

ROOT::VecOps::RVec<float>
get_track_chi2_o_ndf(const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks_in){
  ROOT::VecOps::RVec<float> chi2_out;
  for (const auto &track : tracks_in) {
    chi2_out.push_back(float(track.chi2/track.ndf));
  }
  return chi2_out;
}


// helper for track selection

struct SelectedTracks {
  ROOT::VecOps::RVec<edm4hep::TrackData>  tracks;
  ROOT::VecOps::RVec<edm4hep::TrackState> trackStates;
};


// Base track selection
SelectedTracks
select_tracks_baseline(const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks_in,
              const ROOT::VecOps::RVec<edm4hep::TrackState>& trackstates_in) {
  
  SelectedTracks selected_tracks_and_states;

  // ROOT::VecOps::RVec<edm4hep::TrackData> tracks_out;

  for (const auto &track : tracks_in) {

    // track chi2 selection needs to use track object itself 
    if (track.ndf == 0){
      continue;
    }
    if (track.chi2 / track.ndf > 10.){
      continue;
    }

    // now we need to get the track state to check the other variables:
    auto n_trackstates = track.trackStates_end - track.trackStates_begin;
    // std::cout << n_trackstates << std::endl;

    if (n_trackstates != 1) {
      throw std::runtime_error("Error in track selection: Expected exactly one TrackState per Track");
    }

    // assert(n_trackstates == 1 && "Error in track selection: Expected exactly one TrackState per Track");
    
    for (unsigned int track_state_index = track.trackStates_begin; track_state_index < track.trackStates_end; ++track_state_index) {

      const auto& trackstate = trackstates_in[track_state_index];

      // Make sure covariance Matrix is positive definite
      // Reminder covMatrix convention: https://bib-pubdb1.desy.de/record/81214/files/LC-DET-2006-004%5B1%5D.pdf, sec 5
      const auto& cov_matrix = trackstate.covMatrix;

      if (cov_matrix[0] <= 1e-12 || cov_matrix[2] <= 1e-12 || cov_matrix[9] <= 1e-12) {
        continue;
      }
      if (!std::isfinite(cov_matrix[0]) || !std::isfinite(cov_matrix[2]) || !std::isfinite(cov_matrix[9])) {
        continue;
      }
      
      selected_tracks_and_states.trackStates.push_back(trackstate);
    }

    // if all passed, track is selected
    selected_tracks_and_states.tracks.push_back(track);
    

  }
  return selected_tracks_and_states;
}

//for selecting tracks compatible with primary vertex
SelectedTracks
select_tracks_impactparameters(const SelectedTracks& input,
                               float d0_upper_bound,
                               float z0_upper_bound)
{
    SelectedTracks selected;

    for (size_t i = 0; i < input.tracks.size(); ++i) {

        const auto& track = input.tracks[i];
        const auto& state = input.trackStates[i];

        if (std::abs(state.D0) > d0_upper_bound) continue;
        if (std::abs(state.Z0) > z0_upper_bound) continue;

        selected.tracks.push_back(track);
        selected.trackStates.push_back(state);
    }

    return selected;
}


// Same pre-selection with |D0|, |Z0| re-referenced to the beamspot b [cm]
// instead of the origin. The raw ALEPH/FRFT d0 sign convention is opposite to
// EDM4HEP, hence the +n.b sign. The stored track states are unchanged.
SelectedTracks
select_tracks_impactparameters_bs(const SelectedTracks& input,
                                  float d0_upper_bound,
                                  float z0_upper_bound,
                                  float bsx,
                                  float bsy,
                                  float bsz)
{
    SelectedTracks selected;

    for (size_t i = 0; i < input.tracks.size(); ++i) {

        const auto& track = input.tracks[i];
        const auto& state = input.trackStates[i];

        const double cphi = std::cos(state.phi);
        const double sphi = std::sin(state.phi);

        const double s   = bsx * cphi + bsy * sphi;
        const double d0p = state.D0 - bsx * sphi + bsy * cphi + 0.5 * state.omega * s * s;
        const double z0p = state.Z0 - bsz + state.tanLambda * s;

        if (std::abs(d0p) > d0_upper_bound) continue;
        if (std::abs(z0p) > z0_upper_bound) continue;

        selected.tracks.push_back(track);
        selected.trackStates.push_back(state);
    }

    return selected;
}



// --------------------------------------
// Event primary vertex reconstruction
// --------------------------------------

// Refit reco primary vertex using FCCAna native vertex fitter

// TODO


//Generator level primary vertex based on MC particle info
struct get_EventPrimaryVertexP4 {
  int m_genstatus = 21; // default generator status for incoming hard subprocess
  get_EventPrimaryVertexP4() {};

  TLorentzVector operator()(const ROOT::VecOps::RVec<edm4hep::MCParticleData>& in) const {
    TLorentzVector result(-1e12, -1e12, -1e12, -1e12);
    bool found_py8 = false;

    // First, look for generatorStatus == m_genstatus (e.g., 21 for Pythia8 hard process incoming)
    for (const auto& p : in) {
      //if (p.generatorStatus == m_genstatus) {
      if (1) {
        // vertex.time is in seconds, convert to mm
        TLorentzVector res(p.vertex.x,
                           p.vertex.y,
                           p.vertex.z,
                           p.time * 1.0e3 * 2.99792458e+8);
        result = res;
        found_py8 = true;
        break;
      }
    }

    // Fallback: look for genStatus == 2 with non-zero z vertex
    if (!found_py8) {
      for (const auto& p : in) {
        if (p.generatorStatus == 2 && std::abs(p.vertex.z) > 1.e-12) {
          TLorentzVector res(p.vertex.x,
                             p.vertex.y,
                             p.vertex.z,
                             p.time * 1.0e3 * 2.99792458e+8);
          result = res;
          break;
        }
      }
    }

    return result;
  }
};

// --------------------------------------




float get_EventType(const ROOT::VecOps::RVec<edm4hep::MCParticleData>& in) {
    float result = -1;
    // Look for the first particle with non-zero type

    //std::cout<<"DEBUG: get_EventType called with "<<in.size()<<" particles\n"<<std::endl;

    // for (const auto& p : in) {
    //   std::cout<<"DEBUG: Particle PDG="<<p.PDG<<" px="<<p.momentum.x<<" py="<<p.momentum.y<<" pz="<<p.momentum.z<<" status="<<p.generatorStatus<<std::endl;
    // }

    // first particle encountered with pdg > 0 and < 6
    for (const auto& p : in) {
      if (std::abs(p.PDG) > 0 && std::abs(p.PDG) < 6) {
        result = static_cast<float>(abs(p.PDG));
        break;
      }
    }

    // std::cout<<""<<result<<"\n"<<std::endl;
    return result;
}

rv::RVec<FCCAnalysesJetConstituentsData>
get_isType(const rv::RVec<FCCAnalysesJetConstituentsData>& jcs, float type) {
    rv::RVec<FCCAnalysesJetConstituentsData> out;
    out.reserve(jcs.size());

    for (const auto& jet : jcs) {
        FCCAnalysesJetConstituentsData mask;
        mask.reserve(jet.size());

        for (const auto& c : jet) {
            if (c == type)
                mask.push_back(1);
            else
                mask.push_back(0);
        }

        out.push_back(std::move(mask));
    }

    return out;
}

struct build_constituents_Types {
    rv::RVec<FCCAnalysesJetConstituentsData>
    operator()(const rv::RVec<edm4hep::ParticleIDData> &rpid,
               const std::vector<std::vector<int>> &indices) const
    { 
        rv::RVec<FCCAnalysesJetConstituentsData> jcs;
        for (const auto &jet_index : indices)
        {
            FCCAnalysesJetConstituentsData jc;
            for (const auto &const_index : jet_index)
            {
                jc.push_back(rpid.at(const_index).type);
            }
            jcs.push_back(jc);
        }
        return jcs;
    }
};

// struct build_constituents_Types {
//   // Make the operator static to match usage style
//   static rv::RVec<FCCAnalysesJetConstituentsData>
//   operator()(const rv::RVec<edm4hep::ParticleIDData> &rpid,
//              const std::vector<std::vector<int>> &indices) 
//   {
//     rv::RVec<FCCAnalysesJetConstituentsData> jcs;
//     for (const auto &jet_index : indices)
//     {
//       FCCAnalysesJetConstituentsData jc;
//       for (const auto &const_index : jet_index)
//       {
//         jc.push_back(rpid.at(const_index).type);
//       }
//       jcs.push_back(jc);
//     }
//     return jcs;
//   }
// };

// rv::RVec<FCCAnalysesJetConstituentsData>  build_constituents_Types(const rv::RVec<edm4hep::ParticleIDData> &rpid, const std::vector<std::vector<int>> &indices) 
// { 
//   rv::RVec<FCCAnalysesJetConstituentsData> jcs;
//    for (const auto &jet_index : indices)
//            { FCCAnalysesJetConstituentsData jc;
//               for (const auto &const_index : jet_index) 
//                   { 
//                     jc.push_back(rpid.at(const_index).type);
//                   } 
     
     
//      jcs.push_back(jc); } return jcs; }

// Helper functions for the particle hypothesis p-value using Bethe Bloch fits by Matteo 
// TODO: move to a separate file and load here? 

// Embedded Bethe-Bloch parameters {a, b, c, d, e}
std::unordered_map<std::string, std::unordered_map<std::string, std::vector<double>>> params = {
    {"e", {{"pads", {0.9932659976792287, 1.7094419426691188, 0.07384141776786941, 0.0, -2.0}},
           {"wires", {0.7930482536047483, 2.1221617291314856, 0.048686287341931526, 0.0, -2.0}}}},
    {"mu", {{"pads", {0.8590747, 1.32919203, 0.18728265, 0.01248365, -1.96898002}},
            {"wires", {0.55374855, 1.97190215, 0.26416447, 0.01549576, -1.972309}}}},
    {"pi", {{"pads", {0.536889582, 2.10964282, 0.269949484, 0.00310166058, -3.49991524}},
            {"wires", {0.7922792, 1.30952175, 0.19162276, 0.01523584, -2.4295921}}}},
    {"K", {{"pads", {0.25619823, 3.84049172, 0.53610855, 0.65090594, -2.44656219}},
           {"wires", {0.45920637, 1.82346531, 0.34181439, 0.52306056, -2.21132773}}}},
    {"p", {{"pads", {0.73189858, 1.05891917, 0.256201, 1.34293618, -2.02610177}},
           {"wires", {0.68811606, 1.03354162, 0.24500224, 1.44362792, -2.08656854}}}}
};

// Bethe-Bloch: a * (b + c * log(p) + d * p^e)
double bethe_bloch(double p, const std::vector<double>& par) {
    if (p <= 0.0) return 0.0;
    double a = par[0], b = par[1], c = par[2], d = par[3], e_pow = par[4];
    return a * (b + c * std::log(p) + d * std::pow(p, e_pow));
}

// Single hypothesis/single measurement: signed p-value (-9 if invalid/non-finite), is_wires=true for wires
double signed_p_value(double p, double dedx, double err, const std::string& hypothesis, bool is_wires) {
    std::string sensor = is_wires ? "wires" : "pads";
    
    auto part_it = params.find(hypothesis);
    if (part_it == params.end()) return -9.0;

    auto sens_it = part_it->second.find(sensor);
    if (sens_it == part_it->second.end()) return -9.0;

    const auto& par = sens_it->second;

    if (err <= 0.0) return -9.0;

    double expected = bethe_bloch(p, par);
    double residual = (dedx - expected) / err;
    if (!std::isfinite(residual)) return -9.0;

    double abs_z = std::fabs(residual);
    double sf = 0.5 * std::erfc(abs_z / std::sqrt(2.0));  // norm.sf(|z|)
    return 2.0 * sf * (residual > 0 ? 1.0 : -1.0);
}

// All hypotheses for single measurement: array[5] p-values {e, mu, pi, K, p}
std::array<double, 5> all_hypotheses_pvalues(double p, double dedx, double err, bool is_wires) {
    std::string hypos[5] = {"e", "mu", "pi", "K", "p"};
    std::array<double, 5> results;
    for (int i = 0; i < 5; ++i) {
        results[i] = signed_p_value(p, dedx, err, hypos[i], is_wires);
    }
    return results;
}

// Function to return the dEdx object for each jet constituent and also the array of p-values for PID hypothesis using Bethe-Bloch fits
struct build_constituents_dEdx_PIDhypo{
  struct dEdx_and_PID_result{
    rv::RVec<rv::RVec<edm4hep::RecDqdxData>> dedx_constituents;
    rv::RVec<rv::RVec<std::array<double, 5>>> pid_array_constituents;
    };

    dEdx_and_PID_result
    
    operator()(const rv::RVec<edm4hep::ReconstructedParticleData> &recoParticles,
             const rv::RVec<int> &_recoParticlesIndices,
             const rv::RVec<edm4hep::RecDqdxData> &dEdxCollection,
             const rv::RVec<int> &_dEdxIndicesCollection, 
             const std::vector<std::vector<int>> &jet_indices,
             const rv::RVec<edm4hep::TrackState> &trackStates,
             bool is_wires) const
    { 
        rv::RVec<rv::RVec<edm4hep::RecDqdxData>> dedx_constituents;
        rv::RVec<rv::RVec<std::array<double, 5>>> pid_array_constituents;

        // The links dEdx -> Track and RecoPart -> Track are one-directional, we need a map to store
        // Track.index -> dEdx to not have to loop everytime 
        // in addition, the object itself is stored in <Collection>
        // while the relations (=indices we need for links) are in _<Collection>
        std::unordered_map<int, edm4hep::RecDqdxData> track_index_to_dEdx;
        for (size_t i = 0; i < _dEdxIndicesCollection.size(); ++i) {
          int track_index = _dEdxIndicesCollection[i];
          edm4hep::RecDqdxData dedx_value = dEdxCollection[i];
          track_index_to_dEdx[track_index] = dedx_value;
        }

        //now, for each jet loop over the indices of the jet constituents provided by the JetClusteringUtils
        // retrieve the associated RecoParticle
        // from there get the link to the Track from the corresponding index collection
        for (const auto &jet_const_indices : jet_indices) { //loop over jets
          rv::RVec<edm4hep::RecDqdxData> jet_dEdx;
          rv::RVec<std::array<double, 5>> jet_pid_array;

          for (int constituent_index : jet_const_indices) { // loop over jet constituents
            const auto &recoPart = recoParticles[constituent_index];
            TLorentzVector tlv_recoPart; // needed later to get the momentum total
            tlv_recoPart.SetXYZM(recoPart.momentum.x, recoPart.momentum.y, recoPart.momentum.z, recoPart.mass);

            // Try to find dEdx for this particle, if not found or not good value, use the dummy with defaults
            bool found = false;

            edm4hep::RecDqdxData dEdx_dummy_obj{};
            dEdx_dummy_obj.dQdx.value = -9.0f;
            dEdx_dummy_obj.dQdx.error = -9.0f;
            dEdx_dummy_obj.dQdx.type = -9.0f;

            // also dummy object for the PID hypothesis
            std::array<double, 5> pid_array_dummy{{-9.0f, -9.0f, -9.0f, -9.0f, -9.0f}};

            //loop over tracks associated to the RecoPart (for charged particles should always be exactly one in ALEPH data)
            for (int track = recoPart.tracks_begin; track < recoPart.tracks_end; ++track) {
                 int track_index = _recoParticlesIndices[track]; //this should be the same index used in the link from dEdx to track

                  //find the matching dEdx in the map
                  if (track_index_to_dEdx.count(track_index)) {
                    const auto &dEdx = track_index_to_dEdx[track_index];

                    // A failed leg stores the track's omega as its value;
                    // dQdx.type is the pad-leg status only, so it is not used.
                    const float v = dEdx.dQdx.value;
                    const float omega_sentinel =
                        (track_index >= 0 &&
                         track_index < static_cast<int>(trackStates.size()))
                            ? trackStates[track_index].omega
                            : v; // unknown track: treat as invalid
                    const bool valid = FCCAnalyses::AlephDedx::dEdxValid(
                        v, dEdx.dQdx.error, omega_sentinel);

                    if (valid) {
                      jet_dEdx.push_back(dEdx);
                      jet_pid_array.push_back(all_hypotheses_pvalues(
                          tlv_recoPart.P(), v, dEdx.dQdx.error, is_wires));
                    }
                    else {
                      jet_dEdx.push_back(dEdx_dummy_obj);
                      jet_pid_array.push_back(pid_array_dummy);
                    }

                    found = true;
                    break;
                  }
            }
            // if no track found, i.e. neutral particle, use the dummy 
            if (!found){
              jet_dEdx.push_back(dEdx_dummy_obj);
              jet_pid_array.push_back(pid_array_dummy);
            }
          }
          dedx_constituents.push_back(jet_dEdx); 
          pid_array_constituents.push_back(jet_pid_array);
        }
        return {dedx_constituents, pid_array_constituents};
    }
};

//helpers to read the dEdx objects (to check if can reuse existing FCCAna functions instead):

rv::RVec<rv::RVec<float>> get_dEdx_type(const rv::RVec<rv::RVec<edm4hep::RecDqdxData>> &dedx_vec) {
  rv::RVec<rv::RVec<float>> values;
  for (const auto &inner_vec : dedx_vec) {
    rv::RVec<float> inner_values;
    for (const auto &d : inner_vec) {
      inner_values.push_back(d.dQdx.type);
    }
    values.push_back(inner_values);
  }
  return values;
}

rv::RVec<rv::RVec<float>> get_dEdx_value(const rv::RVec<rv::RVec<edm4hep::RecDqdxData>> &dedx_vec) {
  rv::RVec<rv::RVec<float>> values;
  for (const auto &inner_vec : dedx_vec) {
    rv::RVec<float> inner_values;
    for (const auto &d : inner_vec) {
      inner_values.push_back(d.dQdx.value);
    }
    values.push_back(inner_values);
  }
  return values;
}

rv::RVec<rv::RVec<float>> get_dEdx_error(const rv::RVec<rv::RVec<edm4hep::RecDqdxData>> &dedx_vec) {
  rv::RVec<rv::RVec<float>> values;
  for (const auto &inner_vec : dedx_vec) {
    rv::RVec<float> inner_values;
    for (const auto &d : inner_vec) {
      inner_values.push_back(d.dQdx.error);
    }
    values.push_back(inner_values);
  }
  return values;
}

rv::RVec<rv::RVec<float>> get_PID_pvalue(const rv::RVec<rv::RVec<std::array<double, 5>>> pid_array_vec, int particle_index) {
  rv::RVec<rv::RVec<float>> values;
  for (const auto &inner_vec : pid_array_vec) {
    rv::RVec<float> inner_values;
    for (const auto &d : inner_vec) {
      inner_values.push_back(d[particle_index]);
    }
    values.push_back(inner_values);
  }
  return values;
}

// Return a new collection (same type) with D0 signs flipped.
ROOT::VecOps::RVec<edm4hep::TrackState>
flipD0_copy(const ROOT::VecOps::RVec<edm4hep::TrackState>& tracks) {
  ROOT::VecOps::RVec<edm4hep::TrackState> out;
  out.reserve(tracks.size());
  for (const auto &t : tracks) {
    edm4hep::TrackState tt = t;   // make a copy
    tt.D0 = -tt.D0;               // flip sign
    tt.omega = -tt.omega;

    // Flip the covariance matrix elements that depend on D0 and omega
    // indices to flip 1, 4, 6, 8, 10, 12
    // taken from https://github.com/HEP-FCC/FCCAnalyses/blob/6cccde454007e0ada68a162cc3450ad90c6f65bf/analyzers/dataframe/src/ReconstructedParticle2Track.cc

    tt.covMatrix[1]  = -tt.covMatrix[1];   // cov(D0, phi)
    tt.covMatrix[4]  = -tt.covMatrix[4];   // cov(phi, omega)
    tt.covMatrix[6]  = -tt.covMatrix[6];   // cov(D0, z0)
    tt.covMatrix[8]  = -tt.covMatrix[8];   // cov(omega, z0)
    tt.covMatrix[10] = -tt.covMatrix[10]; // cov(D0, tanLambda)
    tt.covMatrix[12] = -tt.covMatrix[12]; // cov(omega, tanLambda)

    out.push_back(std::move(tt));
  }
  return out;
}



// ===========================================================================
// Beamspot position (data only)
// ===========================================================================
//
// ALEPH's beamspot is offset from the origin by ~0.6 mm in x and ~0.2 mm in y.
// That is 2-3x the transverse beamspot widths used to constrain the primary
// vertex fit, so for DATA the position must be supplied; leaving it at the
// origin biases the constraint. In simulation the beamspot is at the origin by
// construction, so this is not needed there.
//
// One entry point: get_beamspot(run). It loads and caches data/beamspot.json on
// first use and returns the position for that run, or (0,0,0) if the run is not
// listed (same fallback as the reference implementation).
//
// Units: the json stores cm. The FCCAnalyses vertex fitters want the beamspot
// position in the same units as their widths, which we pass as "10 um"
// (res_x_loose/10. etc. in stage1.py), hence the cm -> 10um factor of 1e3.
// Pass `in_10um = false` to get plain cm back instead.
//
// The json path resolves in this order:
//   1. the `path` argument, if non-empty
//   2. $ALEPH_BEAMSPOT_JSON
//   3. <this repo>/data/beamspot.json
//
TVector3 get_beamspot(int run, bool in_10um = true, const std::string &path = "")
{
  // Loaded once on first call and cached. The static initialiser is thread-safe
  // (C++11 magic statics), which matters because RDataFrame runs multi-threaded.
  // Note: only the FIRST call's `path` is used - later calls reuse the cache.
  static const std::map<int, TVector3> coords = [path]() {
    std::map<int, TVector3> m;   // cm; left empty if anything goes wrong -> origin everywhere

    std::string file = path;
    if (file.empty()) {
      if (const char *env = std::getenv("ALEPH_BEAMSPOT_JSON")) file = env;
    }
    if (file.empty()) {
      // default: alongside this header, ../data/beamspot.json
      std::string self = __FILE__;
      size_t slash = self.find_last_of('/');
      file = (slash == std::string::npos ? std::string(".") : self.substr(0, slash))
             + "/../data/beamspot.json";
    }

    std::ifstream in(file);
    if (!in.good()) {
      std::cerr << "WARNING [get_beamspot]: could not open '" << file
                << "' - using a beamspot at the origin for every run. "
                << "That is correct for simulation but WRONG for data." << std::endl;
      return m;
    }
    try {
      nlohmann::json j;
      in >> j;
      for (auto it = j.begin(); it != j.end(); ++it) {
        const auto &v = it.value();
        if (!v.contains("x") || !v.contains("y") || !v.contains("z")) continue;
        m[std::stoi(it.key())] =
            TVector3(v["x"].get<double>(), v["y"].get<double>(), v["z"].get<double>());
      }
    } catch (const std::exception &e) {
      std::cerr << "WARNING [get_beamspot]: failed to parse '" << file
                << "' (" << e.what() << ") - using the origin for every run." << std::endl;
      return std::map<int, TVector3>{};
    }
    std::cout << "INFO [get_beamspot]: loaded " << m.size()
              << " runs from " << file << std::endl;
    return m;
  }();

  TVector3 bs(0., 0., 0.);   // fallback: unknown run, or file missing/unparsable
  auto it = coords.find(run);
  if (it != coords.end()) bs = it->second;
  return in_10um ? bs * 1e3 : bs;   // cm -> 10 um
}

// convenience accessors so stage1.py can Define scalar columns directly
double get_beamspot_x(int run) { return get_beamspot(run).X(); }
double get_beamspot_y(int run) { return get_beamspot(run).Y(); }
double get_beamspot_z(int run) { return get_beamspot(run).Z(); }


auto cast_constituent = [](const auto &jcs, auto &&meth)
    {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      for (const auto &jc : jcs)
        out.emplace_back(meth(jc));
      return out;
    };


rv::RVec<FCCAnalysesJetConstituentsData> get_px(const rv::RVec<FCCAnalysesJetConstituents> &jcs)
    {
      return cast_constituent(jcs, ReconstructedParticle::get_px);
    }


rv::RVec<FCCAnalysesJetConstituentsData> get_py(const rv::RVec<FCCAnalysesJetConstituents> &jcs)
    {
      return cast_constituent(jcs, ReconstructedParticle::get_py);
    }

rv::RVec<FCCAnalysesJetConstituentsData> get_pz(const rv::RVec<FCCAnalysesJetConstituents> &jcs)
    {
      return cast_constituent(jcs, ReconstructedParticle::get_pz);
    }

// ---------------------------------------------------------------------------
// Jet-constituent variables missing from FCCAnalyses' JetConstituentsUtils.
// Definitions follow the reference ntuplizer so the outputs are comparable.
// ---------------------------------------------------------------------------

// ptrel = constituent pT / jet pT  (a ratio, exactly like erel - NOT the
// component of the momentum perpendicular to the jet axis).
rv::RVec<FCCAnalysesJetConstituentsData>
get_ptrel_cluster(const rv::RVec<fastjet::PseudoJet> &jets,
                  const rv::RVec<FCCAnalysesJetConstituents> &jcs)
{
  rv::RVec<FCCAnalysesJetConstituentsData> out;
  for (size_t i = 0; i < jets.size(); ++i) {
    auto &jet_csts = out.emplace_back();
    double pt_jet = jets.at(i).pt();
    auto csts = FCCAnalyses::JetConstituentsUtils::get_jet_constituents(jcs, i);
    for (const auto &jc : csts) {
      TLorentzVector jcvec;
      jcvec.SetXYZM(jc.momentum.x, jc.momentum.y, jc.momentum.z, jc.mass);
      jet_csts.emplace_back(pt_jet > 0. ? jcvec.Pt() / pt_jet : 1.);
    }
  }
  return out;
}

rv::RVec<FCCAnalysesJetConstituentsData>
get_ptrel_log_cluster(const rv::RVec<fastjet::PseudoJet> &jets,
                      const rv::RVec<FCCAnalysesJetConstituents> &jcs)
{
  rv::RVec<FCCAnalysesJetConstituentsData> out;
  for (const auto &jet_csts : get_ptrel_cluster(jets, jcs)) {
    auto &o = out.emplace_back();
    for (const auto &v : jet_csts) o.emplace_back(float(std::log10(v)));
  }
  return out;
}

// --- constituent track parameters w.r.t. the primary vertex ------------------
// Same algebra as ReconstructedParticle2Track::XPtoPar_dxy/dz/phi but with
// every length in cm and every curvature in 1/cm (the upstream helpers mix m and
// mm). `tracks` must be ordered by the RecoParticle->Track relation; neutral
// particles fall outside it and get -9.
struct TrackParamsAtPV {
  rv::RVec<FCCAnalysesJetConstituentsData> dxy, dz, phi0, C, ct;
};

inline TrackParamsAtPV
get_constituent_trackParamsAtPV(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                const rv::RVec<edm4hep::TrackState> &tracks,
                                const TLorentzVector &V, double Bz)
{
  TrackParamsAtPV out;
  for (const auto &jet_csts : jcs) {
    auto &odxy = out.dxy.emplace_back();
    auto &odz  = out.dz.emplace_back();
    auto &ophi = out.phi0.emplace_back();
    auto &oC   = out.C.emplace_back();
    auto &oct  = out.ct.emplace_back();
    for (const auto &rp : jet_csts) {
      // tracks_begin != tracks_end is the "has a track" test: a neutral's
      // tracks_begin is still in range, so the bare range check reads ANOTHER
      // particle's track and propagates NaN into dxy/dz and the btag branches.
      if (!(rp.tracks_begin != rp.tracks_end && rp.tracks_begin < tracks.size())) {
        odxy.push_back(-9.); odz.push_back(-9.); ophi.push_back(-9.);
        oC.push_back(-9.); oct.push_back(-9.);
        continue;
      }
      const auto &ts = tracks.at(rp.tracks_begin);
      const float D0_wrt0 = ts.D0, Z0_wrt0 = ts.Z0, phi0_wrt0 = ts.phi;
      TVector3 X(-D0_wrt0 * TMath::Sin(phi0_wrt0), D0_wrt0 * TMath::Cos(phi0_wrt0), Z0_wrt0);
      TVector3 x = X - V.Vect();
      // The energy-flow momentum is quoted at the track's first point, not at
      // the perigee: keep its magnitude, take the direction from the perigee
      // parameters so that X and p describe the same helix point.
      const double pmag = TVector3(rp.momentum.x, rp.momentum.y, rp.momentum.z).Mag();
      const double tl = ts.tanLambda;
      const double ptp = pmag / TMath::Sqrt(1.0 + tl * tl);
      TVector3 p(ptp * TMath::Cos(phi0_wrt0), ptp * TMath::Sin(phi0_wrt0), ptp * tl);
      const double a = -rp.charge * Bz * AlephUnits::kPtPerTeslaCm;
      const double pt = p.Pt();
      const double r2 = x(0) * x(0) + x(1) * x(1);
      const double cross = x(0) * p(1) - x(1) * p(0);
      const double disc = pt * pt - 2 * a * cross + a * a * r2;
      // T is NaN for disc <= 0; upstream guards the discriminant for dxy only,
      // so dz/phi0 keep using it unguarded and dxy stays at its -9 sentinel.
      const double T = TMath::Sqrt(disc);
      double D = -9.;
      if (disc > 0)
        D = (pt < 10.0) ? (T - pt) / a : (-2 * cross + a * r2) / (T + pt);
      odxy.push_back(D);
      {
        const double C = a / (2 * pt);
        const double Dz = (pt < 10.0) ? (T - pt) / a : (-2 * cross + a * r2) / (T + pt);
        double B = C * TMath::Sqrt(TMath::Max(r2 - Dz * Dz, 0.0) / (1 + 2 * C * Dz));
        if (TMath::Abs(B) > 1.) B = TMath::Sign(1, B);
        const double st = TMath::ASin(B) / C;
        const double ct = p(2) / pt;
        const double dot = x(0) * p(0) + x(1) * p(1);
        odz.push_back((dot > 0.0) ? x(2) - ct * st : x(2) + ct * st);
      }
      ophi.push_back(TMath::ATan2((p(1) - a * x(0)) / T, (p(0) + a * x(1)) / T));
      // Curvature 1/(2R) [1/cm] with the sign of the charge, and cot(theta): both
      // taken straight from the fitted track state rather than from the energy-flow momentum.
      oC.push_back(std::copysign(0.5 * std::abs(ts.omega), rp.charge));
      oct.push_back(ts.tanLambda);
    }
  }
  return out;
}

// --- constituent covariance entry -------------------------------------------
// covMatrix[k] of the constituent's own track state, with the begin!=end test
// that the upstream getters lack: a neutral's tracks_begin can be in range and
// would otherwise return another particle's covariance instead of -9.
inline rv::RVec<FCCAnalysesJetConstituentsData>
get_constituent_trackCov(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                         const rv::RVec<edm4hep::TrackState> &tracks, int k)
{
  rv::RVec<FCCAnalysesJetConstituentsData> out;
  for (const auto &jet_csts : jcs) {
    auto &o = out.emplace_back();
    for (const auto &rp : jet_csts) {
      if (rp.tracks_begin != rp.tracks_end && rp.tracks_begin < tracks.size())
        o.push_back(tracks.at(rp.tracks_begin).covMatrix[k]);
      else
        o.push_back(-9.);
    }
  }
  return out;
}

// --- constituent distance to the jet axis -----------------------------------
// Same algebra as JetConstituentsUtils::get_JetDistVal_clusterV, but the track
// direction is rebuilt from the PV-referenced perigee parameters instead of the
// energy-flow momentum, which is quoted at the track's first point.
inline rv::RVec<FCCAnalysesJetConstituentsData>
get_constituent_jetDistVal(const rv::RVec<fastjet::PseudoJet> &jets,
                           const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                           const rv::RVec<FCCAnalysesJetConstituentsData> &D0,
                           const rv::RVec<FCCAnalysesJetConstituentsData> &Z0,
                           const rv::RVec<FCCAnalysesJetConstituentsData> &phi0,
                           const rv::RVec<FCCAnalysesJetConstituentsData> &ct)
{
  rv::RVec<FCCAnalysesJetConstituentsData> out;
  for (size_t i = 0; i < jets.size(); ++i) {
    auto &o = out.emplace_back();
    TVector3 p_jet(jets[i].px(), jets[i].py(), jets[i].pz());
    const auto &csts = jcs.at(i);
    for (size_t j = 0; j < csts.size(); ++j) {
      const double d0 = D0.at(i).at(j), ph = phi0.at(i).at(j);
      if (d0 == -9) { o.push_back(-9); continue; }
      TVector3 d(-d0 * TMath::Sin(ph), d0 * TMath::Cos(ph), Z0.at(i).at(j));
      TVector3 p_ct(TMath::Cos(ph), TMath::Sin(ph), ct.at(i).at(j));
      TVector3 n = p_ct.Cross(p_jet);
      // the normal is undefined for a track collinear with the jet axis
      if (n.Mag2() <= 0.) { o.push_back(-9); continue; }
      o.push_back(n.Unit().Dot(d));
    }
  }
  return out;
}

// Re-order a track-indexed collection through the ReconstructedParticle->Track
// relation. tracks_begin is an offset into that relation, not a track index, and
// the two are not parallel in the ALEPH files; after this coll.at(p.tracks_begin)
// is the particle's own object.
template <typename T>
rv::RVec<T> reindexByRPLink(const rv::RVec<T> &coll,
                            const rv::RVec<int> &rpTrackIndex) {
  rv::RVec<T> out;
  out.reserve(rpTrackIndex.size());
  for (int idx : rpTrackIndex) {
    // Out-of-range relation entry = corrupt input: fail loudly.
    if (idx < 0 || idx >= static_cast<int>(coll.size()))
      throw std::runtime_error(
          "reindexByRPLink: RP->Track relation entry out of range");
    out.push_back(coll.at(idx));
  }
  return out;
}

// --- per-constituent track quality -----------------------------------------
// A ReconstructedParticle points at its track via tracks_begin; neutral
// constituents have no track, for which we store -1 (as the reference does).

rv::RVec<FCCAnalysesJetConstituentsData>
get_constituent_trackQuality(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                             const rv::RVec<edm4hep::TrackData> &tracks,
                             int mode) // 0 = chi2, 1 = ndof, 2 = chi2/ndof
{
  rv::RVec<FCCAnalysesJetConstituentsData> out;
  for (const auto &jet_csts : jcs) {
    auto &o = out.emplace_back();
    for (const auto &p : jet_csts) {
      float val = -1.;
      // tracks_begin != tracks_end is the actual "has a track" test: for a neutral particle
      // tracks_begin still holds an in-range index, so a bare `tracks_begin < tracks.size()`
      // silently reads ANOTHER particle's track (verified: it mislabels ~99% of neutrals).
      size_t trackIndex = p.tracks_begin;
      if (p.tracks_begin != p.tracks_end && trackIndex < tracks.size()) {
        const edm4hep::TrackData &tr = tracks.at(trackIndex);
        if      (mode == 0) val = tr.chi2;
        else if (mode == 1) val = tr.ndf;
        else                val = (tr.ndf != 0) ? tr.chi2 / float(tr.ndf) : -1.;
      }
      o.emplace_back(val);
    }
  }
  return out;
}

rv::RVec<FCCAnalysesJetConstituentsData>
get_constituent_trackChi2(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                          const rv::RVec<edm4hep::TrackData> &tracks)
{ return get_constituent_trackQuality(jcs, tracks, 0); }

rv::RVec<FCCAnalysesJetConstituentsData>
get_constituent_trackNdof(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                          const rv::RVec<edm4hep::TrackData> &tracks)
{ return get_constituent_trackQuality(jcs, tracks, 1); }

rv::RVec<FCCAnalysesJetConstituentsData>
get_constituent_trackChi2Norm(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                              const rv::RVec<edm4hep::TrackData> &tracks)
{ return get_constituent_trackQuality(jcs, tracks, 2); }

// ORIGINAL-Tracks index of each constituent's own track, -1 when it has none:
// the join key between the pfcand_* block and the finders' *_origIdx branches.
// tracks_begin indexes the RP->Track relation, whose entries are track indices.
rv::RVec<rv::RVec<int>>
get_constituent_trackIdx(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                         const rv::RVec<int> &rpTrackIndex)
{
  rv::RVec<rv::RVec<int>> out;
  for (const auto &jet_csts : jcs) {
    auto &o = out.emplace_back();
    for (const auto &p : jet_csts) {
      int val = -1;
      size_t slot = p.tracks_begin;
      if (p.tracks_begin != p.tracks_end && slot < rpTrackIndex.size())
        val = rpTrackIndex.at(slot);
      o.emplace_back(val);
    }
  }
  return out;
}

// --- per-constituent subdetector hit counts --------------------------------
// subdetectorNumber assumes inside-out ordering: 0 = VDET, 1 = ITC, 2 = TPC.
rv::RVec<FCCAnalysesJetConstituentsData>
get_constituent_nTrackHits(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                           const rv::RVec<edm4hep::TrackData> &tracks,
                           const rv::RVec<int> &subdetectorHitNumbers,
                           int subdetectorNumber)
{
  rv::RVec<FCCAnalysesJetConstituentsData> out;
  for (const auto &jet_csts : jcs) {
    auto &o = out.emplace_back();
    for (const auto &p : jet_csts) {
      float nHits = -1.;
      // see note in get_constituent_trackQuality: neutrals need the begin!=end test
      size_t trackIndex = p.tracks_begin;
      if (p.tracks_begin != p.tracks_end && trackIndex < tracks.size()) {
        const edm4hep::TrackData &tr = tracks.at(trackIndex);
        size_t hitIdx = tr.subdetectorHitNumbers_begin + subdetectorNumber;
        if (hitIdx < subdetectorHitNumbers.size())
          nHits = subdetectorHitNumbers.at(hitIdx);
      }
      o.emplace_back(nHits);
    }
  }
  return out;
}

rv::RVec<FCCAnalysesJetConstituentsData>
get_constituent_nTrackHits_VDET(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                const rv::RVec<edm4hep::TrackData> &tracks,
                                const rv::RVec<int> &shn)
{ return get_constituent_nTrackHits(jcs, tracks, shn, 0); }

rv::RVec<FCCAnalysesJetConstituentsData>
get_constituent_nTrackHits_ITC(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                               const rv::RVec<edm4hep::TrackData> &tracks,
                               const rv::RVec<int> &shn)
{ return get_constituent_nTrackHits(jcs, tracks, shn, 1); }

rv::RVec<FCCAnalysesJetConstituentsData>
get_constituent_nTrackHits_TPC(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                               const rv::RVec<edm4hep::TrackData> &tracks,
                               const rv::RVec<int> &shn)
{ return get_constituent_nTrackHits(jcs, tracks, shn, 2); }


rv::RVec<rv::RVec<int>> mask(const rv::RVec<FCCAnalysesJetConstituentsData> &energies)
    {
    rv::RVec<rv::RVec<int>> out;
    for (const auto &e_vec : energies)  // Iterates over all jets (1, 2, 3, ...)
    {
        rv::RVec<int> jet_mask;
        for (const auto &e : e_vec)  // Iterates over constituents in each jet
        {
            jet_mask.emplace_back(e != 0.0f ? 1 : 0);
        }
        out.emplace_back(jet_mask);
    }
    return out;
}

// -----------------------------------
// Custom helpers for secondary vertexing
// -----------------------------------

using FCCAnalysesVertex = FCCAnalyses::VertexingUtils::FCCAnalysesVertex;


// helper function which assigns secondary vertices as found per event to the jets using a closest dR match 
ROOT::VecOps::RVec<ROOT::VecOps::RVec<FCCAnalysesVertex>>
assign_SV_to_jets(
    const ROOT::VecOps::RVec<FCCAnalysesVertex>& secondary_vertex_objects, 
    const ROOT::VecOps::RVec<fastjet::PseudoJet>& jets)
{   
  // returned object is a vector of vectors where outer index are the jets, and inner index are the secondary vertices assigned to that jet 
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<FCCAnalysesVertex>> result(jets.size());
    if (jets.size() == 0) return result;

    // use helper function from FCCAna to obtain the secondary vertex 3 vector
    ROOT::VecOps::RVec<TVector3> sv_momenta = FCCAnalyses::VertexingUtils::get_p_SV(secondary_vertex_objects);

    for (unsigned int sv_idx = 0; sv_idx < secondary_vertex_objects.size(); sv_idx++) {
        TVector3 sv_p = sv_momenta.at(sv_idx);
        if (sv_p.Mag() < 1e-10) {
            std::cerr << "WARNING: Secondary vertex with zero momentum encountered, skipping jet assignment" << std::endl;
            continue;
        }
        double minDR = 99.;
        unsigned int best_jet = 0;
        for (unsigned int j = 0; j < jets.size(); j++) {
            double dR = sv_p.DeltaR(TVector3(jets[j].px(), jets[j].py(), jets[j].pz()));
            if (dR < minDR) { minDR = dR; best_jet = j; }
        }
        result.at(best_jet).push_back(secondary_vertex_objects.at(sv_idx));
    }
    return result;
}

// V0 rejection with ALEPH-tuned tight constraints, following Luka's implementation in ntuplizer
// inclusive=false (default): a track can only appear in one V0 pair (more physically correct).
// inclusive=true:  matches ntuplizer behaviour — a track can be flagged by multiple V0 pairs
ROOT::VecOps::RVec<edm4hep::TrackState>
V0rejection_ALEPH(
    const ROOT::VecOps::RVec<edm4hep::TrackState>& np_tracks,
    const FCCAnalysesVertex& PV,
    double solenoidBz = 1.5,
    bool inclusive = false)
{
    int nTr = np_tracks.size();
    ROOT::VecOps::RVec<bool> isInV0(nTr, false);
    if (nTr < 2) return np_tracks;

    ROOT::VecOps::RVec<edm4hep::TrackState> tr_pair;
    edm4hep::TrackState tr_i, tr_j;
    tr_pair.push_back(tr_i);
    tr_pair.push_back(tr_j);
    FCCAnalysesVertex V0_vtx;

    for (unsigned int i = 0; i < nTr-1; i++) {
        if (!inclusive && isInV0[i]) continue;
        tr_pair[0] = np_tracks[i];
        for (unsigned int j = i+1; j < nTr; j++) {
            if (!inclusive && isInV0[j]) continue;
            if (tr_pair[0].omega * np_tracks[j].omega > 0) continue;
            tr_pair[1] = np_tracks[j];

            auto cand = FCCAnalyses::VertexFinderLCFIPlus::get_V0candidate(
                V0_vtx, tr_pair, PV, true, 10., solenoidBz);
            if (cand.size() == 0) continue;

            // ALEPH-tuned tight constraints (widened mass windows, reduced distance minimum)
            bool isKs    = cand[0]>0.453 && cand[0]<0.553 && cand[4]>0.1 && cand[5]>0.999;
            bool isLam1  = cand[1]>1.06  && cand[1]<1.16  && cand[4]>0.1 && cand[5]>0.99995;
            bool isLam2  = cand[2]>1.06  && cand[2]<1.16  && cand[4]>0.1 && cand[5]>0.99995;
            bool isGamma = cand[3]<0.005 && cand[4]>0.9   && cand[5]>0.99995;

            if (isKs || isLam1 || isLam2 || isGamma) {
                isInV0[i] = true;
                isInV0[j] = true;
                if (!inclusive) break;
            }
        }
    }

    ROOT::VecOps::RVec<edm4hep::TrackState> result;
    for (unsigned int i = 0; i < nTr; i++)
        if (!isInV0[i]) result.push_back(np_tracks[i]);
    return result;
}

// SV finding with all ALEPH-specific defaults: 1.5 T field, ALEPH-tuned V0 rejection,
// dR prefilter enabled. Set inclusive_v0=true to match ntuplizer behaviour exactly.
ROOT::VecOps::RVec<FCCAnalysesVertex>
get_SV_event_ALEPH(
    const ROOT::VecOps::RVec<edm4hep::TrackState>& np_tracks,
    const ROOT::VecOps::RVec<edm4hep::TrackState>& all_tracks,
    const FCCAnalysesVertex& PV,
    double dR_cut = 0.8,
    bool inclusive_v0 = false)
{
    auto tracks_no_v0 = V0rejection_ALEPH(np_tracks, PV, 1.5, inclusive_v0);
    return FCCAnalyses::VertexFinderLCFIPlus::get_SV_event(
        tracks_no_v0, all_tracks, PV,
        false,         // V0 rejection already done above with ALEPH constraints
        10., 10., 5., // chi2_cut, invM_cut, chi2Tr_cut
        1.5,           // solenoidBz [T]
        dR_cut,       // dR_cut for prefiltering
        true,          // require opposite-charge seed pairs (matches FCCAnalyses@3a4de97 VertexSeed_best)
        false          // LOOSE V0 constraints in per-pair seed screening.
                       // FCCAnalyses@3a4de97 VertexSeed_best does isV0(tr_pair, PV, false) -- explicitly
                       // commented "V0 rejection (loose)" -- while the track-level V0rejection_tight uses
                       // tight. Two different tightnesses; we previously had tight in both.
    );
}

// utils for SV properties:

// SV vertex-fit covariance component ic (packed lower triangle:
// 0=xx 1=yx 2=yy 3=zx 4=zy 5=zz — same order as the Vertex_refit_cov_*
// PV branches), nested per jet like the other sv_* getters.
inline ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>>
svCovComp(
    const ROOT::VecOps::RVec<ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex>>& sv_jets,
    int ic)
{
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> result;
    for (const auto& jet : sv_jets) {
        ROOT::VecOps::RVec<float> row;
        for (const auto& v : jet) row.push_back(v.vertex.covMatrix[ic]);
        result.push_back(row);
    }
    return result;
}

// SV displacement from PV in lab frame x/y/z [mm], per jet
ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>>
get_dx_SV_jets(
    const ROOT::VecOps::RVec<ROOT::VecOps::RVec<FCCAnalysesVertex>>& vertices,
    const TVector3& PV)
{
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>> result;
    for (const auto& jet_vtx : vertices) {
        ROOT::VecOps::RVec<double> temp;
        for (const auto& vtx : jet_vtx)
            temp.push_back(vtx.vertex.position[0] - PV.x());
        result.push_back(temp);
    }
    return result;
}

ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>>
get_dy_SV_jets(
    const ROOT::VecOps::RVec<ROOT::VecOps::RVec<FCCAnalysesVertex>>& vertices,
    const TVector3& PV)
{
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>> result;
    for (const auto& jet_vtx : vertices) {
        ROOT::VecOps::RVec<double> temp;
        for (const auto& vtx : jet_vtx)
            temp.push_back(vtx.vertex.position[1] - PV.y());
        result.push_back(temp);
    }
    return result;
}

ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>>
get_dz_SV_jets(
    const ROOT::VecOps::RVec<ROOT::VecOps::RVec<FCCAnalysesVertex>>& vertices,
    const TVector3& PV)
{
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>> result;
    for (const auto& jet_vtx : vertices) {
        ROOT::VecOps::RVec<double> temp;
        for (const auto& vtx : jet_vtx)
            temp.push_back(vtx.vertex.position[2] - PV.z());
        result.push_back(temp);
    }
    return result;
}

// SV momentum relative to jet momentum (|p_SV| / |p_jet|), per jet
ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>>
get_prel_SV_jets(
    const ROOT::VecOps::RVec<ROOT::VecOps::RVec<FCCAnalysesVertex>>& vertices,
    const ROOT::VecOps::RVec<fastjet::PseudoJet>& jets)
{
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>> result;
    for (unsigned int i = 0; i < jets.size(); i++) {
        double jet_pmag = TVector3(jets[i].px(), jets[i].py(), jets[i].pz()).Mag();
        ROOT::VecOps::RVec<double> sv_pmags = FCCAnalyses::VertexingUtils::get_pMag_SV(vertices[i]);
        ROOT::VecOps::RVec<double> temp;
        for (double pmag : sv_pmags)
            temp.push_back(pmag / jet_pmag);
        result.push_back(temp);
    }
    return result;
}

// Pointing angle of SV wrt PV, per jet.
// Note: FCCAnalyses::VertexingUtils::get_pointingangle_SV has a bug — it uses the absolute
// SV position instead of the displacement vector (SV - PV). This is the corrected version, following Luka's implementation in ntuplizer
ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>>
get_pointingangle_SV(
    const ROOT::VecOps::RVec<ROOT::VecOps::RVec<FCCAnalysesVertex>>& vertices,
    const FCCAnalysesVertex& PV)
{
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>> result;
    edm4hep::Vector3f r_PV = PV.vertex.position;
    for (const auto& jet_vtx : vertices) {
        ROOT::VecOps::RVec<double> temp;
        for (const auto& vtx : jet_vtx) {
            TVector3 p_sum;
            for (const auto& p_tr : vtx.updated_track_momentum_at_vertex)
                p_sum += p_tr;
            edm4hep::Vector3f r_vtx = vtx.vertex.position;
            TVector3 r_vtx_PV(r_vtx[0] - r_PV[0],
                               r_vtx[1] - r_PV[1],
                               r_vtx[2] - r_PV[2]);
            double cosAngle = p_sum.Dot(r_vtx_PV) / (p_sum.Mag() * r_vtx_PV.Mag());
            temp.push_back(cosAngle);
        }
        result.push_back(temp);
    }
    return result;
}

// -----------------------------------
// V0 reconstruction
// -----------------------------------

// Corrected invariant mass of a single SV (CMS-BTV-16-002):
//   sqrt(m^2 + p^2*sin^2(theta)) + p*sin(theta)
// where theta is the angle between the SV momentum sum and the SV-PV displacement vector.
// This accounts for neutral particles that are not reconstructed at the SV.
double get_correctedInvMass_SV_single(
    const FCCAnalysesVertex& sv,
    const FCCAnalysesVertex& PV)
{
    double rawMass = FCCAnalyses::VertexingUtils::get_invM(sv);
    TVector3 p_sum;
    for (const auto& p_tr : sv.updated_track_momentum_at_vertex)
        p_sum += p_tr;
    double p_mag = p_sum.Mag();
    edm4hep::Vector3f r_sv = sv.vertex.position;
    edm4hep::Vector3f r_pv = PV.vertex.position;
    TVector3 flight(r_sv[0]-r_pv[0], r_sv[1]-r_pv[1], r_sv[2]-r_pv[2]);
    double cosTheta = p_sum.Dot(flight) / (p_mag * flight.Mag());
    if (std::abs(cosTheta) > 1.) cosTheta /= std::abs(cosTheta);
    double sin2Theta = 1. - cosTheta*cosTheta;
    double sinTheta  = std::sqrt(sin2Theta);
    return std::sqrt(rawMass*rawMass + p_mag*p_mag*sin2Theta) + p_mag*sinTheta;
}

// Per-jet corrected invariant mass of SVs
ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>>
get_correctedInvMass_SV(
    const ROOT::VecOps::RVec<ROOT::VecOps::RVec<FCCAnalysesVertex>>& vertices,
    const FCCAnalysesVertex& PV)
{
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>> result;
    for (const auto& jet_vtx : vertices) {
        ROOT::VecOps::RVec<double> temp;
        for (const auto& vtx : jet_vtx)
            temp.push_back(get_correctedInvMass_SV_single(vtx, PV));
        result.push_back(temp);
    }
    return result;
}

// Count V0s of a given PDG ID per jet (e.g. 310 for Ks, 3122 for Lambda)
ROOT::VecOps::RVec<int>
count_V0type_jets(
    const ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>>& pdg_per_jet,
    int target_pdg)
{
    ROOT::VecOps::RVec<int> result;
    for (const auto& jet_pdg : pdg_per_jet) {
        int n = 0;
        for (int pdg : jet_pdg) if (pdg == target_pdg) n++;
        result.push_back(n);
    }
    return result;
}

// V0 finding with ALEPH-specific constraints, matching the values used in V0rejection_ALEPH
FCCAnalyses::VertexingUtils::FCCAnalysesV0
get_V0s_ALEPH(
    const ROOT::VecOps::RVec<edm4hep::TrackState>& np_tracks,
    const FCCAnalysesVertex& PV,
    double solenoidBz = 1.5, bool loose_mass_window = false,
    double dR_pair_cut = -1., bool exclusive_tracks = false)
{
  if (loose_mass_window){
      return FCCAnalyses::VertexFinderLCFIPlus::get_V0s(
          np_tracks, PV,
          0.1, 1.4, 0.1, 0.999,    // Ks:     mass window [GeV], dis_min [cm=1mm], cosAng
          0.1, 1.4, 0.1, 0.999,    // Lambda: dis_min 0.1 cm = 1 mm physical
          0.0, -1,  0.9, 0.999,    // Gamma:  invM_max=-1 (never passes, matching ntuplizer loose mode)
          10., solenoidBz, dR_pair_cut, exclusive_tracks
      );
  }

  else{
      return FCCAnalyses::VertexFinderLCFIPlus::get_V0s(
          np_tracks, PV,
          0.453, 0.553, 0.1, 0.999,    // Ks:     mass window [GeV], dis_min [cm=1mm], cosAng
          1.06,  1.16,  0.1, 0.99995,  // Lambda
          0.0,   0.005, 0.9, 0.99995,  // Gamma
          10., solenoidBz, dR_pair_cut, exclusive_tracks
      );
  }
}

// Bundles V0 candidates distributed over jets (vtx, PDG ID, invariant mass together).
// Keeping these three arrays in sync ensures the PDG-to-mass correspondence is never broken.
struct V0sPerJet {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<FCCAnalysesVertex>> vtx;
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>>               pdgAbs;
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>>            invM;
};

// Assign event-level V0 candidates to jets using closest-dR matching
V0sPerJet
assign_V0s_to_jets(
    const FCCAnalyses::VertexingUtils::FCCAnalysesV0& v0s,
    const ROOT::VecOps::RVec<fastjet::PseudoJet>& jets)
{
    V0sPerJet result;
    result.vtx.resize(jets.size());
    result.pdgAbs.resize(jets.size());
    result.invM.resize(jets.size());
    if (jets.size() == 0 || v0s.vtx.size() == 0) return result;

    ROOT::VecOps::RVec<TVector3> v0_momenta = FCCAnalyses::VertexingUtils::get_p_SV(v0s.vtx);
    for (unsigned int i = 0; i < v0s.vtx.size(); i++) {
        TVector3 v0_p = v0_momenta.at(i);
        if (v0_p.Mag() < 1e-10) continue;
        double minDR = 99.;
        unsigned int best_jet = 0;
        for (unsigned int j = 0; j < jets.size(); j++) {
            double dR = v0_p.DeltaR(TVector3(jets[j].px(), jets[j].py(), jets[j].pz()));
            if (dR < minDR) { minDR = dR; best_jet = j; }
        }
        result.vtx[best_jet].push_back(v0s.vtx[i]);
        result.pdgAbs[best_jet].push_back(v0s.pdgAbs[i]);
        result.invM[best_jet].push_back(v0s.invM[i]);
    }
    return result;
}


}} // namespace FCCAnalyses::AlephSelection

#endif

