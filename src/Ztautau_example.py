"""Z -> tau tau example with one charged track per hemisphere (1 + 1 prongs).
Selects lepton-pair events with exactly one selected track in each thrust
hemisphere and writes, per hemisphere, the track's kinematics, fit quality,
dE/dx with ALPHA's particle hypotheses, energy-flow type and (on MC) truth, the
energy-flow neutral energy of the hemisphere and (on MC) the tau decay mode.
Run from this directory (src/), on data and on MC:

    fccanalysis run Ztautau_example.py -i <data.root> -o <out.root> -j 8 -- --doData
    fccanalysis run Ztautau_example.py -i <mc.root> -o <out.root> -j 8

The input is an EDM4hep file of a converter version that writes ALPHA's dE/dx
collections (dEdxWiresPID, dEdxPadsPID) and RunInformation (on MC also
dEdxWiresQDEDXMPID, dEdxWiresQDEDXMMCLink, trackMCLink and the MCParticles
parents and daughters).
An input without them stops with an error naming them, before any event is
read.

Options after `--`: --doData (no truth branches, data beam spot), --year,
--noRunList and --excludeRuns RUN [RUN ...] (the run selection of stage1, from
run_list.py; the last two data only), --beamspotJson PATH (data only). An
unknown option is an error.

Event selection, in this order:
  1. data: the runs of the run list (nominal policy of run_list.py);
  2. ALEPH event class 15 (dilepton) or 24 (extended dilepton; ALPHA guide,
     appendix G), EVENT_CLASSES;
  3. a defined thrust axis: the exact thrust of stage1's analysis particles,
     oriented as in stage1 towards the hemisphere of lower energy;
  4. exactly one track of stage1's baseline selection (chi2/ndf <= 10,
     AlephSelection::kTrackMinTPCHits TPC hits, |z0| <= kTrackMaxAbsZ0,
     finite parameters, finite positive-definite perigee covariance) in each
     hemisphere; a track belongs to hemisphere 0 if its perigee momentum
     has a positive projection on the thrust axis, else to 1.
No charge, momentum, energy or particle-type requirement is applied.

Output branches. Event branches have one value per event (event_class one per
class bit set); every hem_* and trk_* branch is an array of size 2, entry 0 =
hemisphere 0 (the one the thrust axis points into), entry 1 = hemisphere 1;
trk_* describe the one selected track of that hemisphere. h = ele, mu, pi,
kaon, proton. Unit in [], value when undefined in ().
  Event
    event_number, run_number           event and run number
    event_class                        the ALEPH event classes of the event (class N = bit N-1)
    EVT_Ebeam                          LEP beam energy [GeV] (RunInformation)
    EVT_Evis                           summed energy of the analysis particles [GeV]
    EVT_Thrust_Mag, _X, _Y, _Z         thrust and its unit axis
    EVT_Thrust_cosTheta                cosine of the polar angle of the axis
  Hemisphere
    hem_photon_E, hem_photon_n         summed energy [GeV] and number of energy-flow photons (type 4)
                                       with energy >= HEM_NEUTRAL_MIN_E (0 = all; EVT_Evis and the thrust
                                       drop the neutrals below AlephSelection::kNeutralMinE, 1 GeV)
    hem_neutralHad_E, hem_neutralHad_n the same for energy-flow neutral hadrons (type 5)
  Track (perigee of the stored track state, field 1.5 T)
    trk_origIdx                        index of the track in Tracks
    trk_p, trk_pt                      momentum and transverse momentum [GeV]
    trk_costheta, trk_phi              polar-angle cosine and azimuth [rad] (-999)
    trk_q                              charge [e]
    trk_d0, trk_z0                     impact parameters w.r.t. the beam spot (data: per-run position;
                                       MC: origin) [cm], d0 in the LCIO sign convention, as pfcand_d0
    trk_d0Err, trk_z0Err               their uncertainties from the track covariance [cm] (-999)
    trk_nVDET, trk_nITC, trk_nTPC      hits in the vertex detector, inner drift chamber and TPC (-1)
    trk_chi2ndf                        track-fit chi2/ndf (-1)
    trk_efType                         energy-flow type of the track's particle: 0 charged hadron,
                                       1 electron, 2 muon (ALEPH's lepton identification), 3 track of
                                       an ALEPH V0; -1 = no particle
  Track dE/dx, normalised to about 1 for a minimum-ionising particle; -999 = not available
    trk_dEdx_{wires,pads}_...          TPC wires (dEdxWiresPID, ALPHA's QDEDX) / pads (dEdxPadsPID, QPADX):
      _value, _nSamples                measured dE/dx and its number of samples
      _trackLength                     track length used [cm]
      _errorCode                       ALPHA error code, 0 = measured
      _expected_<h>, _sigma_<h>        ALPHA's expected dE/dx and resolution for hypothesis h, momentum
                                       error included (-999 also for a hypothesis ALPHA cannot evaluate,
                                       error code 5: low-momentum kaons and protons)
      _pull_<h>                        (value - expected_<h>) / sigma_<h>
    trk_dEdx_wiresQDEDXM_{value,errorCode,expected_<h>,sigma_<h>,pull_<h>}
                                       MC only: the same from ALPHA's parametrised simulation QDEDXM
  MC truth (PDG codes: +211 = pi+, but +11 = e-)
    trk_truePdg, trk_trueMotherPdg     ALPHA's best-match MC particle of the track (KBESTM) and its first
                                       parent; 0 = no match / no parent; -999 = event without the ALPHA
                                       MC banks
    trk_trueIdx                        its MCParticles index (-1; -999 as truePdg)
    hem_tauMode                        decay mode of the hemisphere's tau: 1 electron, 2 muon, 3 one
                                       charged hadron and no pi0, 4 one charged hadron and one pi0, 5 one
                                       charged hadron and two or more pi0, 6 three charged prongs (any
                                       pi0), 0 other; -1 = not a tau-pair event
    hem_tauNPi0                        pi0 of that decay (-999)
    hem_tauVisP                        magnitude of the summed momentum of its decay products except the
                                       neutrinos [GeV] (-999)

Tau truth. The two taus are MCParticles[0] and [1]; their generator-level
decay chains are followed down to particles without generator daughters, or
to pi0, K0S, K0L, Lambda and photons (a K0S -> pi+ pi- adds no prong).
Charged = |charge| > 0.5. The tau whose visible momentum has the larger
projection on the thrust axis is the tau of hemisphere 0.

Reconstruction chain, every Define written out; the helpers of the hemisphere
split, kinematics and tau truth are in analyzer_ztautau_example.h, the others
are stage1's. Further per-track quantities are reached through trk_origIdx:
arrays over Tracks read with ROOT::VecOps::Take, as the dE/dx below.
"""
import os
import sys
from argparse import ArgumentParser
# fccanalysis loads this file by path, so this directory's modules are made importable here
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_list

EVENT_CLASSES = (15, 24)    # ALEPH classes (ALPHA guide, appendix G), one of them needed: 15 dilepton, 24 extended dilepton
HEM_NEUTRAL_MIN_E = 0.0     # [GeV] lowest energy of a photon / neutral hadron summed per hemisphere (0 = all;
                            # not the kNeutralMinE cut of the analysis particles)
# per-run beam-spot positions of the data
BEAMSPOT_JSON = "/eos/experiment/fcc/ee/analyses/case-studies/aleph/utils/beamspot_position_data/beamspot.json"
# ALPHA dE/dx hypotheses h, in the order of AlephVars in analyzer_alephvars.h (as in stage1)
DEDX_HYPOS = ("ele", "mu", "pi", "kaon", "proton")
# converter collections the branches read; the second set MC only
ALEPH_COLLS = ("dEdxWiresPID", "_dEdxWiresPID_parameters", "dEdxPadsPID", "_dEdxPadsPID_parameters", "RunInformation")
ALEPH_TRUTH_COLLS = ("dEdxWiresQDEDXMPID", "_dEdxWiresQDEDXMPID_parameters", "_dEdxWiresQDEDXMMCLink_from",
                     "_dEdxWiresQDEDXMMCLink_to", "_MCParticles_parents", "_MCParticles_daughters",
                     "trackMCLink", "_trackMCLink_from", "_trackMCLink_to")


class Analysis():

    def __init__(self, cmdline_args):
        parser = ArgumentParser(description='Z -> tau tau example arguments (after --)')
        parser.add_argument('--doData', action='store_true',
                            help='Run on data: run list and beam spot from the JSON, no truth branches.')
        parser.add_argument('--year', default='1994',
                            help='Year of the run list (data).')
        parser.add_argument('--excludeRuns', nargs='+', action='extend', default=[], type=run_list.run_number, metavar='RUN',
                            help='Data only: drop these run numbers in addition to the run list.')
        parser.add_argument('--noRunList', action='store_true',
                            help='Data only: keep every run instead of the run list (--excludeRuns still applies).')
        parser.add_argument('--beamspotJson', default=os.environ.get("ALEPH_BEAMSPOT_JSON", BEAMSPOT_JSON),
                            help='Data only: per-run beam-spot positions (default $ALEPH_BEAMSPOT_JSON, else the EOS file).')
        self.ana_args, unknown = parser.parse_known_args(cmdline_args['remaining'])
        if unknown:
            print(f"----> ERROR: unrecognised arguments: {' '.join(unknown)}")
            sys.exit(1)
        if not self.ana_args.doData and (self.ana_args.excludeRuns or self.ana_args.noRunList):
            print("----> ERROR: --excludeRuns and --noRunList apply to data only (--doData).")
            sys.exit(1)

        # input and output come from the fccanalysis command line: no process list
        if not cmdline_args.get('input') and not cmdline_args.get('input_file_list'):
            sys.exit("Ztautau_example.py: give the input file(s) with -i and the output file with -o")
        self.n_threads = 8      # -j on the command line overrides

        self.include_paths = ["aleph_units.h", "analyzer.h", "analyzer_trkaux.h", "analyzer_alephvars.h",
                              "analyzer_ztautau_example.h"]

    def analyzers(self, df):
        ZT = "FCCAnalyses::AlephZtautau"
        AV = "FCCAnalyses::AlephVars"

        # ---- converter collections --------------------------------------------
        # checked before any event is read: an input of an older converter lacks them
        have = set(str(c) for c in df.GetColumnNames())
        absent = [c for c in ALEPH_COLLS + (() if self.ana_args.doData else ALEPH_TRUTH_COLLS) if c not in have]
        if absent:
            print(f"----> ERROR: the input lacks the converter collections {', '.join(absent)}: it comes from an "
                  "older converter version; convert it again with a current one.")
            sys.exit(1)

        # ---- run list (data) and event class ----------------------------------
        if self.ana_args.doData:
            df = run_list.filter_runs(df, self.ana_args.year, self.ana_args.excludeRuns, self.ana_args.noRunList)
        df = df.Filter(" || ".join(f"AlephSelection::sel_class_filter({c})(ClassBitset)" for c in EVENT_CLASSES))
        df = df.Define("event_number", "EventHeader.eventNumber")
        df = df.Define("run_number", "EventHeader.runNumber")
        df = df.Define("event_class", "AlephSelection::bitsetToIndices(ClassBitset)")
        df = df.Define("EVT_Ebeam", f"{AV}::alephRunInfo(RunInformation, false).beamEnergy")

        # ---- analysis particles and thrust (as stage1) -----------------------
        # RecoParticlesNoLumiCal = energy-flow types 0 to kPFMaxAnalysisType, AnalysisParticles = those without the
        # photons and neutral hadrons below kNeutralMinE; filtered copies with their own ParticleID (.pid)
        df = df.Define("RecoParticlesNoLumiCal_subset", "AlephSelection::filterParticles(RecoParticles, _RecoParticles_tracks.index, ParticleID, _ParticleID_particle.index, AlephSelection::isAnalysisType())")
        df = df.Define("RecoParticlesNoLumiCal",         "RecoParticlesNoLumiCal_subset.particles")
        df = df.Define("ParticleIDNoLumiCal",            "RecoParticlesNoLumiCal_subset.pid")
        df = df.Define("_RecoParticlesNoLumiCal_tracks", "RecoParticlesNoLumiCal_subset.tracks")
        df = df.Define("AnalysisParticles", "AlephSelection::filterParticles(RecoParticlesNoLumiCal_subset, AlephSelection::notSoftNeutral(AlephSelection::kNeutralMinE)).particles")
        df = df.Define("AnalysisParticles_px", "ReconstructedParticle::get_px(AnalysisParticles)")
        df = df.Define("AnalysisParticles_py", "ReconstructedParticle::get_py(AnalysisParticles)")
        df = df.Define("AnalysisParticles_pz", "ReconstructedParticle::get_pz(AnalysisParticles)")
        df = df.Define("AnalysisParticles_e",  "ReconstructedParticle::get_e(AnalysisParticles)")
        df = df.Define("EVT_Evis", "Sum(AnalysisParticles_e)")
        # exact thrust {T, x, y, z} (T = -1 with fewer than two particles), repacked for getThrustPointing,
        # which orients the axis towards the hemisphere of lower energy
        df = df.Define("EVT_thrustExact", "Algorithms::calculate_thrust()(AnalysisParticles_px, AnalysisParticles_py, AnalysisParticles_pz)")
        df = df.Filter("EVT_thrustExact[0] > 0")
        df = df.Define("EVT_thrustNP", "ROOT::VecOps::RVec<float>{EVT_thrustExact[0], EVT_thrustExact[1], 0.f, EVT_thrustExact[2], 0.f, EVT_thrustExact[3], 0.f}")
        df = df.Define("AnalysisParticles_thrustangleNP", "Algorithms::getAxisCosTheta(EVT_thrustNP, AnalysisParticles_px, AnalysisParticles_py, AnalysisParticles_pz)")
        df = df.Define("EVT_thrust", "Algorithms::getThrustPointing(1.)(AnalysisParticles_thrustangleNP, AnalysisParticles_e, EVT_thrustNP)")
        df = df.Define("EVT_Thrust_Mag", "EVT_thrust.at(0)")
        df = df.Define("EVT_Thrust_X",   "EVT_thrust.at(1)")
        df = df.Define("EVT_Thrust_Y",   "EVT_thrust.at(3)")
        df = df.Define("EVT_Thrust_Z",   "EVT_thrust.at(5)")
        df = df.Define("EVT_Thrust_cosTheta", "EVT_Thrust_Z / std::sqrt(EVT_Thrust_X*EVT_Thrust_X + EVT_Thrust_Y*EVT_Thrust_Y + EVT_Thrust_Z*EVT_Thrust_Z)")
        thrust = "EVT_Thrust_X, EVT_Thrust_Y, EVT_Thrust_Z"

        # ---- baseline tracks and the 1 + 1 selection --------------------------
        # baseline: chi2/ndf <= 10, kTrackMinTPCHits TPC hits, |z0| <= kTrackMaxAbsZ0, finite parameters and
        # covariance; .trackStates and .origIdx (index into Tracks) share one order
        df = df.Define("tracks_selected_baseline_result", "AlephSelection::select_tracks_baseline(Tracks, _Tracks_trackStates, _Tracks_subdetectorHitNumbers, AlephSelection::kTrackMinTPCHits, AlephSelection::kTrackMaxAbsZ0)")
        df = df.Define("hemTracks", f"{ZT}::splitTracks(tracks_selected_baseline_result.trackStates, tracks_selected_baseline_result.origIdx, {thrust})")
        # exactly one track per hemisphere: the size-2 output layout assumes it
        df = df.Filter("hemTracks.n[0] == 1 && hemTracks.n[1] == 1")
        df = df.Define("trk_origIdx", "hemTracks.origIdx")

        # ---- energy-flow neutrals per hemisphere ------------------------------
        for name, typ in (("photon", "AlephSelection::kPFPhoton"), ("neutralHad", "AlephSelection::kPFNeutralHad")):
            df = df.Define(f"hemEF_{name}", f"{ZT}::neutralsPerHem(RecoParticlesNoLumiCal, ParticleIDNoLumiCal, {typ}, {HEM_NEUTRAL_MIN_E}, {thrust})")
            df = df.Define(f"hem_{name}_E", f"hemEF_{name}.energy")
            df = df.Define(f"hem_{name}_n", f"hemEF_{name}.count")

        # ---- track kinematics, hits, fit quality, energy-flow type ------------
        # reference of d0, z0: data per-run beam spot from the JSON (10 um units -> cm), MC the origin
        if self.ana_args.doData:
            df = df.Define("BeamspotVec", f'AlephSelection::get_beamspot(run_number[0], true, "{self.ana_args.beamspotJson}")')
            bs = "BeamspotVec.X()*1e-3, BeamspotVec.Y()*1e-3, BeamspotVec.Z()*1e-3"
        else:
            bs = "0., 0., 0."
        df = df.Define("trkKin", f"{ZT}::trackKinematics(trk_origIdx, Tracks, _Tracks_trackStates, FCCAnalyses::AlephUnits::kBz, {bs})")
        kin = ["p", "pt", "costheta", "phi", "q", "d0", "z0", "d0Err", "z0Err"]
        for m in kin:
            df = df.Define(f"trk_{m}", f"trkKin.{m}")
        for k, det in enumerate(("VDET", "ITC", "TPC")):
            df = df.Define(f"trk_n{det}", f"FCCAnalyses::AlephTrkAux::subdetHits(trk_origIdx, Tracks.subdetectorHitNumbers_begin, Tracks.subdetectorHitNumbers_end, _Tracks_subdetectorHitNumbers, {k})")
        df = df.Define("trk_chi2ndf", "FCCAnalyses::AlephTrkAux::trackChi2Ndf(trk_origIdx, Tracks.chi2, Tracks.ndf)")
        # Tracks index -> RecoParticlesNoLumiCal index (-1 = none), then that particle's ParticleID type
        df = df.Define("rpOfTrack", "FCCAnalyses::AlephTrkAux::rpIndexByTrack(RecoParticlesNoLumiCal, _RecoParticlesNoLumiCal_tracks, Tracks.size())")
        df = df.Define("trk_efType", f"{ZT}::trackEFType(trk_origIdx, rpOfTrack, ParticleIDNoLumiCal)")
        self.branches = [f"trk_{m}" for m in kin] + ["trk_nVDET", "trk_nITC", "trk_nTPC", "trk_chi2ndf", "trk_efType"]

        # ---- ALPHA dE/dx (converter collections) -------------------------------
        # AlephVars::trackDedx gives arrays over the whole Tracks collection (stage1's trk_dEdx_*); each
        # hemisphere reads its entry at trk_origIdx. -999 (AlephVars::kNotAvailable) = not available
        df = df.Define("_trk_take", "ROOT::VecOps::RVec<std::size_t>(trk_origIdx.begin(), trk_origIdx.end())")
        # (branch infix; collection; parameter count; empty collection allowed)
        dets = [("wires", "dEdxWiresPID", "kNParDedx", "false"),
                ("pads", "dEdxPadsPID", "kNParDedx", "false")]
        if not self.ana_args.doData:
            dets.append(("wiresQDEDXM", "dEdxWiresQDEDXMPID", "kNParQDEDXM", "true"))
        for name, coll, npar, empty_ok in dets:
            df = df.Define(f"trkDedx_{name}", f'{AV}::trackDedx({coll}, _{coll}_parameters, Tracks.size(), '
                                              f'{AV}::{npar}, "{coll}", {empty_ok})')
            members = [("value", "value", "kNotAvailable")]
            if name != "wiresQDEDXM":
                members += [("nSamples", "nSamples", "kNotAvailable"), ("trackLength", "trackLength", "kNotAvailable")]
            members += [("errorCode", "errorCode", "kNotAvailableInt")]
            members += [(f"{q}_{h}", f"{q}[{k}]", "kNotAvailable") for k, h in enumerate(DEDX_HYPOS) for q in ("expected", "sigma")]
            members += [(f"pull_{h}", f"pull[{k}]", "kNotAvailable") for k, h in enumerate(DEDX_HYPOS)]
            for b, m, na in members:
                self.branches.append(f"trk_dEdx_{name}_{b}")
                df = df.Define(self.branches[-1], f"ROOT::VecOps::Take(trkDedx_{name}.{m}, _trk_take, {AV}::{na})")

        # ---- MC truth -----------------------------------------------------------
        if not self.ana_args.doData:
            # ALPHA's best-match MC particle of each track (KBESTM) and its first parent
            df = df.Define("trkBestMatch", f"{AV}::trackBestMatchTruth(Tracks.size(), dEdxWiresQDEDXMPID.empty(), "
                                           "_dEdxWiresQDEDXMMCLink_from, _dEdxWiresQDEDXMMCLink_to, MCParticles, "
                                           "_MCParticles_parents, _trackMCLink_from, _trackMCLink_to, "
                                           "trackMCLink.weight)")
            for b, m in (("truePdg", "pdg"), ("trueMotherPdg", "motherPdg"), ("trueIdx", "idx")):
                self.branches.append(f"trk_{b}")
                df = df.Define(self.branches[-1], f"ROOT::VecOps::Take(trkBestMatch.{m}, _trk_take, {AV}::kNotAvailableInt)")
            # tau decay mode of each hemisphere, from the MCParticles decay chain
            df = df.Define("tauTruth", f"{ZT}::tauTruthPerHem(MCParticles, _MCParticles_daughters.index, {thrust})")
            for b, m in (("tauMode", "mode"), ("tauNPi0", "nPi0"), ("tauVisP", "visP")):
                self.branches.append(f"hem_{b}")
                df = df.Define(self.branches[-1], f"tauTruth.{m}")
        return df

    def output(self):
        branches = [
            "event_number", "run_number", "event_class", "EVT_Ebeam", "EVT_Evis",
            "EVT_Thrust_Mag", "EVT_Thrust_X", "EVT_Thrust_Y", "EVT_Thrust_Z", "EVT_Thrust_cosTheta",
            "hem_photon_E", "hem_photon_n", "hem_neutralHad_E", "hem_neutralHad_n",
            "trk_origIdx",
        ]
        return branches + self.branches
