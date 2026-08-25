
import os
from argparse import ArgumentParser

BZ = 1.5  # solenoid field [T] — single source for the stage1 Define strings
PVNEW = "FCCAnalyses::AlephPVNew"  # namespace holding the PV selection constants

# Per-daughter dE/dx: the collections to read and the branch suffixes they
# produce, in the order they are written.
DEDX_COLLS = (("pads", "dEdxPads"), ("wires", "dEdxWires"))
DEDX_BRANCHES = tuple(f"dEdx_{_d}_{_q}" for _d, _ in DEDX_COLLS
                      for _q in ("value", "error"))
# particle-flow label joined onto every candidate leg through its origIdx
LEG_PID_BRANCHES = ("isChargedHad",)
# particle-flow type code of a charged hadron (single source: analyzer_trkaux.h)
PF_CHARGED_HAD = "FCCAnalyses::AlephTrkAux::kPFChargedHad"

# V0-module branch lists. Each entry is (branch suffix, Define expression), so
# a name appears once and drives both the Define chain and the output list.
V0N_CAND_DEFINES = (
    ("pdg",         "V0sNew_event.pdgAbs"),
    ("invM",        "V0sNew_event.invM"),
    ("alpha",       "FCCAnalyses::AlephV0New::candAlpha(V0sNew_event, SecondaryTracks_looseBS)"),
    ("qt",          "FCCAnalyses::AlephV0New::candQt(V0sNew_event)"),
    ("chi2",        "FCCAnalyses::AlephTruth::candChi2(V0sNew_event)"),
    ("dxyz",        "FCCAnalyses::AlephTruth::candDxyz(V0sNew_event, VertexObject_looseBS)"),
    ("p",           "FCCAnalyses::AlephTruth::candP(V0sNew_event)"),
    # momentum VECTOR of the same summed vertex momentum as v0n_p
    # (direction-dependent offline studies: pointing at any reference)
    ("px",          "FCCAnalyses::AlephTruth::candPcomp(V0sNew_event, 0)"),
    ("py",          "FCCAnalyses::AlephTruth::candPcomp(V0sNew_event, 1)"),
    ("pz",          "FCCAnalyses::AlephTruth::candPcomp(V0sNew_event, 2)"),
    ("cosPointing", "FCCAnalyses::AlephTruth::candCosPointing(V0sNew_event, VertexObject_looseBS)"),
    ("pointSig",    "FCCAnalyses::AlephV0New::candPointSig(V0sNew_event, VertexObject_looseBS)"),
    # two-tier module: 1 = adopted tight package, 0 = loose training tier.
    # Selecting v0n_tight==1 reproduces the historical tight-only output exactly.
    ("tight",       "FCCAnalyses::AlephV0New::candTight(V0sNew_event, VertexObject_looseBS, SecondaryTracks_looseBS)"),
    # ML-input pulls: cut variables in resolution units (signed; -999 undefined).
    ("bandSig",     "FCCAnalyses::AlephV0New::candBandSig(V0sNew_event, SecondaryTracks_looseBS)"),
    ("massSig",     "FCCAnalyses::AlephV0New::candMassSig(V0sNew_event)"),
    # fitted-vertex position (position-resolution studies)
    ("vx",          "FCCAnalyses::AlephTruth::candVtxPos(V0sNew_event, 0)"),
    ("vy",          "FCCAnalyses::AlephTruth::candVtxPos(V0sNew_event, 1)"),
    ("vz",          "FCCAnalyses::AlephTruth::candVtxPos(V0sNew_event, 2)"),
    # vertex-fit covariance (packed lower triangle, cm^2 - same component
    # order as Vertex_refit_cov_*)
    ("cov_xx",      "FCCAnalyses::AlephV0New::candCovComp(V0sNew_event, 0)"),
    ("cov_yx",      "FCCAnalyses::AlephV0New::candCovComp(V0sNew_event, 1)"),
    ("cov_yy",      "FCCAnalyses::AlephV0New::candCovComp(V0sNew_event, 2)"),
    ("cov_zx",      "FCCAnalyses::AlephV0New::candCovComp(V0sNew_event, 3)"),
    ("cov_zy",      "FCCAnalyses::AlephV0New::candCovComp(V0sNew_event, 4)"),
    ("cov_zz",      "FCCAnalyses::AlephV0New::candCovComp(V0sNew_event, 5)"),
)
V0N_TRKS = ("trk1", "trk2")
V0N_TRUTH_DEFINES = (
    ("v0n_class",                "v0ntruth.cls"),
    ("v0n_trueidx",              "v0ntruth.true_idx"),
    ("v0n_pairmult",             "v0ntruth.pair_mult"),
    ("v0n_trackshared",          "v0ntruth.track_shared"),
    ("v0n_trk1",                 "v0ntruth.trk1"),
    ("v0n_trk2",                 "v0ntruth.trk2"),
    ("truev0_foundnew_any",      "FCCAnalyses::AlephTruth::trueV0FoundAny(trueV0s, v0ntruth)"),
    ("truev0_foundnew_correct",  "FCCAnalyses::AlephTruth::trueV0FoundCorrect(trueV0s, v0ntruth, V0sNew_event)"),
)
# per-jet mirror of the old-finder v0_* block, on the new candidates
V0NJET_DEFINES = (
    ("pdg",           "v0njet_per_jet.pdgAbs"),
    ("invM",          "v0njet_per_jet.invM"),
    ("chi2",          "FCCAnalyses::VertexingUtils::get_chi2_SV(v0njet_jets)"),
    ("chi2_norm",     "FCCAnalyses::VertexingUtils::get_norm_chi2_SV(v0njet_jets)"),
    ("ndof",          "FCCAnalyses::VertexingUtils::get_nDOF_SV(v0njet_jets)"),
    ("ntracks",       "FCCAnalyses::VertexingUtils::get_VertexNtrk(v0njet_jets)"),
    ("p",             "FCCAnalyses::VertexingUtils::get_pMag_SV(v0njet_jets)"),
    ("prel",          "FCCAnalyses::AlephSelection::get_prel_SV_jets(v0njet_jets, jets)"),
    ("thetarel",      "FCCAnalyses::VertexingUtils::get_relTheta_SV(v0njet_jets, jets)"),
    ("phirel",        "FCCAnalyses::VertexingUtils::get_relPhi_SV(v0njet_jets, jets)"),
    ("dxy",           "FCCAnalyses::VertexingUtils::get_dxy_SV(v0njet_jets, VertexObject_looseBS)"),
    ("dxyz",          "FCCAnalyses::VertexingUtils::get_d3d_SV(v0njet_jets, VertexObject_looseBS)"),
    ("cosPointing",   "FCCAnalyses::AlephSelection::get_pointingangle_SV(v0njet_jets, VertexObject_looseBS)"),
    ("correctedMass", "FCCAnalyses::AlephSelection::get_correctedInvMass_SV(v0njet_jets, VertexObject_looseBS)"),
    ("dx",            "FCCAnalyses::AlephSelection::get_dx_SV_jets(v0njet_jets, PrimaryVertexP3)"),
    ("dy",            "FCCAnalyses::AlephSelection::get_dy_SV_jets(v0njet_jets, PrimaryVertexP3)"),
    ("dz",            "FCCAnalyses::AlephSelection::get_dz_SV_jets(v0njet_jets, PrimaryVertexP3)"),
)
# generated V0s and the truth classification of the legacy-finder candidates
TRUEV0_DEFINES = (
    ("pdg",          "trueV0s.pdg"),
    ("p",            "trueV0s.p"),
    ("costheta",     "trueV0s.costheta"),
    ("px",           "trueV0s.px"),
    ("py",           "trueV0s.py"),
    ("pz",           "trueV0s.pz"),
    ("fd",           "trueV0s.fd"),
    ("dpv",          "trueV0s.dpv"),
    ("nmatched",     "trueV0s.nmatched"),
    # daughters surviving into the secondary-track set (0-2): separates
    # PV-claim losses from finder losses
    ("nsec",         "FCCAnalyses::AlephTruth::daughtersInSecondaries(trueV0s, mcToTracks, sec2origIdx)"),
    ("found_any",     "FCCAnalyses::AlephTruth::trueV0FoundAny(trueV0s, v0truth)"),
    ("found_correct", "FCCAnalyses::AlephTruth::trueV0FoundCorrect(trueV0s, v0truth, V0s_event)"),
    # true decay-point components [cm] (position-resolution studies)
    ("x",            "trueV0s.vx"),
    ("y",            "trueV0s.vy"),
    ("z",            "trueV0s.vz"),
)
V0C_DEFINES = (
    ("class",        "v0truth.cls"),
    ("trueidx",      "v0truth.true_idx"),
    ("pairmult",     "v0truth.pair_mult"),
    ("trackshared",  "v0truth.track_shared"),
    ("alpha",        "v0truth.alpha"),
    ("qt",           "v0truth.qt"),
    ("trk1",         "v0truth.trk1"),
    ("trk2",         "v0truth.trk2"),
    # event-order candidate kinematics (independent of jet assignment)
    ("pdg",          "V0s_event.pdgAbs"),
    ("invM",         "V0s_event.invM"),
    ("dxyz",         "FCCAnalyses::AlephTruth::candDxyz(V0s_event, VertexObject_looseBS)"),
    ("p",            "FCCAnalyses::AlephTruth::candP(V0s_event)"),
    ("cosPointing",  "FCCAnalyses::AlephTruth::candCosPointing(V0s_event, VertexObject_looseBS)"),
    # fitted-vertex position (position-resolution studies)
    ("vx",           "FCCAnalyses::AlephTruth::candVtxPos(V0s_event, 0)"),
    ("vy",           "FCCAnalyses::AlephTruth::candVtxPos(V0s_event, 1)"),
    ("vz",           "FCCAnalyses::AlephTruth::candVtxPos(V0s_event, 2)"),
)

# phi->KK branch-name lists: single source for the Define chain and the output
# branch list (per-candidate quantities, and the per-daughter block).
PHIKK_CAND_BRANCHES = ("invM", "p", "px", "py", "pz", "alpha", "qt", "bandEll",
                       "chi2", "vx", "vy", "vz", "dpv", "dpvSig", "same_sign",
                       "wp", "tight")
PHIKK_TRKS = ("trk1", "trk2")
PHIKK_TRK_BRANCHES = ("origIdx", "q", "p", "costheta", "d0", "z0", "sigd0",
                      "nvdet", "nitc", "chi2ndf", "isprim")
TRUEPHI_BRANCHES = ("mothPdg", "origin", "p", "pt", "costheta",
                    "px", "py", "pz", "vx", "vy", "vz", "nmatched",
                    "dauPlus_p", "dauMinus_p")
# per-daughter truth labels, shared by the phi and D* candidate blocks
TRK_TRUTH_BRANCHES = ("mcpdg", "mothpdg")

# D* branch-name lists: single source for the Define chain and the
# output branch list (per-candidate quantities, and the per-daughter block
# instantiated for the K, the pi and the slow pi).
# kinematics the D0 and the D* entry share; they live in the CandKin member
CAND_KIN_BRANCHES = ("m_kpi", "p", "px", "py", "pz", "costheta", "xE", "chi2",
                     "vx", "vy", "vz", "dpv", "dpvSig", "cosPoint",
                     "cosThetaStar")
D0_CAND_BRANCHES = ("m_kpi", "p", "px", "py", "pz", "costheta", "xE", "chi2",
                    "vx", "vy", "vz", "dpv", "dpvSig", "cosPoint",
                    "cosThetaStar", "loose", "tight", "nsec")
DSTAR_CAND_BRANCHES = ("m_kpi", "dm", "p", "px", "py", "pz", "costheta", "xE",
                       "chi2", "vx", "vy", "vz", "dpv", "dpvSig", "cosPoint",
                       "cosThetaStar", "rs", "loose", "tight", "d0idx", "nsec")
DSTAR_TRK_BRANCHES = ("origIdx", "q", "p", "costheta", "d0", "z0", "sigd0",
                      "nvdet", "nitc", "chi2ndf", "isprim", "pool")
TRUED0_BRANCHES = ("p", "pt", "costheta", "xE", "pK", "cosK", "pPi", "origin",
                   "fromDstar", "mothPdg", "nmatched", "flight",
                   "K_pool", "pi_pool")
TRUEDSTAR_BRANCHES = ("p", "pt", "costheta", "xE", "px", "py", "pz", "pK",
                      "cosK", "pPi", "pPis", "cosPis", "origin", "mothPdg",
                      "nmatched", "d0flight", "K_pool", "pi_pool", "pis_pool")
# (branch prefix, member of DstarCands) of every stored daughter leg
DSTAR_TRK_LEGS = (("d0_trkK", "d0.trkK"), ("d0_trkPi", "d0.trkPi"),
                  ("dstar_trkK", "ds.trkK"), ("dstar_trkPi", "ds.trkPi"),
                  ("dstar_trkPis", "ds.trkPis"))
D0_TRKS = ("trkK", "trkPi")
DSTAR_TRKS = ("trkK", "trkPi", "trkPis")


def _cand_member(block, branch):
    """Member path of a D0/D* candidate branch inside DstarCands."""
    return f"{block}.kin.{branch}" if branch in CAND_KIN_BRANCHES \
        else f"{block}.{branch}"

class Analysis():

    def __init__(self, cmdline_args):
        parser = ArgumentParser(
            description='Additional analysis arguments',
            usage='Provide additional arguments after analysis script path')
        parser.add_argument('--tag', required=True, type=str,
                            help='Production tag to indicate version.')
        parser.add_argument('--doData', action='store_true',
                            help='Run on data, instead of MC (which is the default behaviour).')
        parser.add_argument('--year', default='1994',
                            help='MC/data year to run on - currently only 1994 as option.')
        parser.add_argument('--MCtype', default="zqq", type=str,
                            help='Type of MC to run on - currently only zqq as option.')
        parser.add_argument('--MCflavour', default=None, type=str,
                            help='For MC only: filter out events based on truth quark flavours. Default is none. Options: \
                            1 = dd, 2 = uu, 3 = ss, 4 = cc, 5 = bb')
        parser.add_argument('--fraction', default=1.0, type=float,
                            help='Fraction of events to run, default is 1.0 = 100%')
        parser.add_argument('--batch', action='store_true', 
                            help='Submit to HTCondor batch')
        parser.add_argument('--valid', action='store_true', 
                            help='Run tester file only for validation against Lukas ntuples.')
        parser.add_argument('--chunks', default=None, type=int,
                            help='Number of chunks per process/file')
        parser.add_argument('--nthreads', default=None, type=int,
                            help='Override the number of RDataFrame threads (for running several variants concurrently).')
        parser.add_argument('--procfile', default=None, type=str,
                            help='Process a single QQB input file (name without .root, e.g. ZM4212_40_AL) - for condor file-level splitting.')
        # The standalone PV/SV/V0 modules are the DEFAULT chain; the --old*
        # switches below restore the legacy code paths for cross-checks.
        parser.add_argument('--oldV0', action='store_true',
                            help='Legacy V0 only: drop the two-tier V0 module (no v0n_*/v0njet_* branches, no V0 truth branches). Implies --oldSV.')
        parser.add_argument('--oldSV', action='store_true',
                            help='Legacy SV only: drop the standalone SV module (no svn_*/svm_* branches).')
        parser.add_argument('--oldPV', action='store_true',
                            help='Legacy PV chain: get_PrimaryTracks + VertexFitter_Tk and the origin-referenced track pre-selection, instead of the standalone fitter and its beamspot-referenced window (no pv_* flag branches).')
        parser.add_argument('--noPhiKK', action='store_true',
                            help='Skip the phi(1020)->K+K- finder (no phikk_*/truephi_* branches).')
        parser.add_argument('--noDstar', action='store_true',
                            help='Skip the D*+ -> D0(K pi) pi_slow finder (no dstar_*/d0_*/truedstar_*/trued0_* branches).')
        parser.add_argument('--excludeRuns', nargs='+', default=[], type=int, metavar='RUN',
                            help='data only: veto these run numbers before any selection (eventsProcessed still counts the raw input).')
        # Parse additional arguments not known to the FCCAnalyses parsers
        # All command line arguments know to fccanalysis are provided in the
        # `cmdline_arg` dictionary.
        self.ana_args, _ = parser.parse_known_args(cmdline_args['remaining'])

        # Module switches: new PV/SV/V0 by default, --old* opts back out.
        # The SV finder consumes the tight-V0 track veto, so --oldV0 forces --oldSV.
        if self.ana_args.oldV0 and not self.ana_args.oldSV:
            print("----> NOTE: --oldV0 implies --oldSV (the SV finder needs the tight-V0 veto).")
            self.ana_args.oldSV = True
        self.do_v0new = not self.ana_args.oldV0
        self.do_svnew = not self.ana_args.oldSV
        self.do_pvnew = not self.ana_args.oldPV
        self.do_phikk = not self.ana_args.noPhiKK
        self.do_dstar = not self.ana_args.noDstar
        # V0 truth matching needs generator information: MC only.
        self.do_truth = self.do_v0new and not self.ana_args.doData

        #Dictionary for setting output names:
        outnames_dict = {
            # proc: {flavour_id_1:{flavour_name_1}, flavour_id_2:{flavour_name_2}, ..}
            "zqq":{
                "1":"Zdd",
                "2":"Zuu",
                "3":"Zss",
                "4":"Zcc",
                "5":"Zbb",
                }
        }

        # sanity checks for the command line arguments:
        if self.ana_args.doData and self.ana_args.MCtype:
            print("----> WARNING: Incompatible input arguments: --MCtype defined with --doData, will be ignored.")

        if self.ana_args.doData and self.ana_args.MCflavour:
            print("----> WARNING: Incompatible input arguments: --MCflavour defined with --doData, will be ignored.")

        if self.ana_args.MCflavour and not self.ana_args.MCtype:
            print("----> ERROR: Requested truth flavour filter with --MCflavour without specifying --MCtype.")
            exit()
        
        if self.ana_args.MCtype and not self.ana_args.MCtype in outnames_dict:
            print("----> ERROR: Requested unknown --MCtype. Currently only zqq available.")
            exit()
        
        if not self.ana_args.doData and not self.ana_args.MCflavour:
            print(f"----> ERROR: Requested MC run but did not specify --MCflavour. Please pick one..")
            exit()
        
        if self.ana_args.MCflavour and not self.ana_args.MCflavour in outnames_dict[self.ana_args.MCtype]:
            print(f"----> ERROR: Requested unknown --MCflavour for --MCtype {self.ana_args.MCtype}. Check the dictionary.")
            exit()

        #set the input/output directories:
        if self.ana_args.doData:
            self.input_dir = "/eos/experiment/fcc/ee/analyses/case-studies/aleph/LEP1_DATA/"
            _r = os.environ.get("ALEPH_RECLUS_DIR")
            if _r and os.path.isdir(_r):
                print(f"----> INPUT OVERRIDE (ALEPH_RECLUS_DIR): {_r}")
                self.input_dir = _r
            self.output_dir_eos = f"/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedData/{self.ana_args.year}/stage1/{self.ana_args.tag}"

            if self.ana_args.batch:
                self.output_dir = "./data/"
                if self.ana_args.chunks:
                    self.process_list = {
                        "1994" : {"fraction" : self.ana_args.fraction, "chunks":self.ana_args.chunks},           
                    }
                else:
                    self.process_list = {
                        "1994" : {"fraction" : self.ana_args.fraction},           
                    }

                self.n_threads = 8 

            else:
                self.process_list = {
                    "1994" : {"fraction" : self.ana_args.fraction},           
                }

                # ALEPH_OUT_DIR redirects non-batch output to a writable area
                _o = os.environ.get("ALEPH_OUT_DIR")
                self.output_dir = f"{_o}/wp2_data/{self.ana_args.tag}" if _o else "."

                # file-level splitting for condor: one job = one data file
                if self.ana_args.procfile:
                    self.process_list = {
                        f"1994/{self.ana_args.procfile}": {
                            "fraction": self.ana_args.fraction,
                            "output": f"data_{self.ana_args.procfile}",
                        },
                    }

                self.n_threads = 32 

        else:
            self.input_dir = f"/eos/experiment/aleph/EDM4HEP/MC/{self.ana_args.year}/"
            # Optional local re-clustered input copy: the raw files have ~3 TTree
            # clusters, which caps RDataFrame at ~3-4 threads. Opt-in via env var.
            _reclus = os.environ.get("ALEPH_RECLUS_DIR")
            if _reclus and os.path.isdir(_reclus):
                print(f"----> INPUT OVERRIDE (ALEPH_RECLUS_DIR): {_reclus}")
                self.input_dir = _reclus
            # self.output_dir = "."

            #set the output file name depending on resonance flavour 
            output_name = outnames_dict[self.ana_args.MCtype][self.ana_args.MCflavour]

            if self.ana_args.batch:

                self.output_dir_eos = f"/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedMC/{self.ana_args.year}/{self.ana_args.MCtype}/stage1/{self.ana_args.tag}/{output_name}/"
                self.output_dir = f"./{output_name}/"
                
                #split in chunks or not (if works well should make this default)
                if self.ana_args.chunks:
                    self.process_list = {
                        "QQB" : {"fraction" : self.ana_args.fraction, "output":output_name, "chunks":self.ana_args.chunks},        
                    }
                else:
                    self.process_list = {
                        "QQB" : {"fraction" : self.ana_args.fraction, "output":output_name},        
                    }

                self.n_threads = 8
            
            else:
                # ALEPH_OUT_DIR redirects non-batch output to a writable area
                _o = os.environ.get("ALEPH_OUT_DIR")
                self.output_dir = (f"{_o}/wp1_stage1/{self.ana_args.tag}" if _o else
                    f"/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedMC/{self.ana_args.year}/{self.ana_args.MCtype}/stage1/{self.ana_args.tag}")

                #local tester for validation
                if self.ana_args.valid:

                    self.process_list = { 
                        "QQB/ZM4212_39_AL" : {"fraction" : self.ana_args.fraction, "output":"ntuple_valid_tester_{}".format(self.ana_args.MCflavour)},           
                    }
                
                #process full files: 
                else:
                    self.process_list = {
                            "QQB" : {"fraction" : self.ana_args.fraction, "output":output_name},        
                        }

                    # file-level splitting for condor: one job = one input file
                    if self.ana_args.procfile:
                        self.process_list = {
                            f"QQB/{self.ana_args.procfile}": {
                                "fraction": self.ana_args.fraction,
                                "output": f"{output_name}_{self.ana_args.procfile}",
                            },
                        }

            
                self.n_threads = 32 


        #set run options:

        if self.ana_args.nthreads is not None:
            self.n_threads = self.ana_args.nthreads

        # analyzer_truth.h is loaded unconditionally: its truth-FREE helpers
        # (selectedBaselineOriginalIndices / secondaryToOriginalTrack) back the
        # always-written prim2origIdx / sec2origIdx index-map branches, on data too.
        # analyzer_pvnew.h is loaded unconditionally too: both PV chains read the
        # selection constants it defines.
        # analyzer_trkaux.h is unconditional: it carries the vertex-fit glue and
        # the track auxiliaries shared by every finder, plus the always-written
        # per-track membership and particle-flow joins.
        self.include_paths = ["analyzer.h", "analyzer_truth.h", "analyzer_pvnew.h",
                              "analyzer_trkaux.h"]
        for _flag, _hdr in ((self.do_v0new, "analyzer_v0new.h"),
                            (self.do_svnew, "analyzer_svnew.h"),
                            (self.do_phikk, "analyzer_phikk.h"),
                            (self.do_dstar, "analyzer_dstar.h")):
            if _flag and _hdr not in self.include_paths:
                self.include_paths.append(_hdr)

        # #submit to batch if requested:
        # self.run_batch = self.ana_args.batch # no longer supported

    @staticmethod
    def _pv_guard(expr, empty):
        """Empty-return entry guard on the usable-PV predicate: a finder must
        not run on a vertex that is not converged, fully pruned, and
        track-supported (pv_good, goodPV() in analyzer_pvnew.h)."""
        return f"pv_good ? {expr} : {empty}"

    def _define_dedx_join(self, df):
        """Track index -> dE/dx measurement index, once per collection per
        event. A failed leg copies the track omega into dQdx.value, so the
        shared dEdxValid gate is applied here and both branches read -1."""
        for _det, _coll in DEDX_COLLS:
            df = df.Define(f"dedxJoin_{_det}",
                           f"FCCAnalyses::AlephV0New::dedxIndexByTrack({_coll}.dQdx.value, {_coll}.dQdx.error, _{_coll}_track.index, _Tracks_trackStates)")
        return df

    def _define_dedx(self, df, legs):
        """dE/dx value+error per daughter-leg prefix, joined through
        <prefix>_origIdx. STORED for the calibration, never selected on."""
        for _pfx in legs:
            for _det, _coll in DEDX_COLLS:
                for _q in ("value", "error"):
                    df = df.Define(f"{_pfx}_dEdx_{_det}_{_q}",
                                   f"FCCAnalyses::AlephV0New::trackQuantityByIndex({_pfx}_origIdx, {_coll}.dQdx.{_q}, dedxJoin_{_det})")
        return df

    def _define_leg_pid(self, df, legs):
        """Tri-state particle-flow charged-hadron label per daughter-leg prefix,
        joined through <prefix>_origIdx: 1 = PF charged hadron, 0 = another PF
        type, -1 = the track has no linked ReconstructedParticle."""
        for _pfx in legs:
            df = df.Define(f"{_pfx}_isChargedHad",
                           f"FCCAnalyses::AlephTrkAux::legIsChargedHad({_pfx}_origIdx, rpOfTrack, ParticleID)")
        return df

    def analyzers(self, df):

        coll = {
        "GenParticles": "MCParticles",
        "PFParticles": "RecoParticles",
        "PFTracks": "EFlowTrack",
        "PFPhotons": "EFlowPhoton",
        "PFNeutralHadrons": "EFlowNeutralHadron",
        "TrackState": "_Tracks_trackStates",
        "TrackerHits": "TrackerHits",
        "CalorimeterHits": "CalorimeterHits",
        "PathLength": "EFlowTrack_L",
        "Bz": "magFieldBz",
        }

        if self.ana_args.doData:
            if self.ana_args.excludeRuns:
                veto = " && ".join(f"EventHeader.runNumber[0] != {r}" for r in sorted(set(self.ana_args.excludeRuns)))
                df = df.Filter(veto, "excludeRuns")
            #df = df.Filter("AlephSelection::sel_class_filter(16)(ClassBitset)   || AlephSelection::sel_class_filter(17)(ClassBitset) ")
            df = df.Filter("AlephSelection::sel_class_filter(16)(ClassBitset) ")
            df = df.Define("jetPID", "-999")
        else:
            # Using Classbit to filter out QQbar samples and then get a specific flavor of jets
            # d-quark: 1, u-quark:2, s-quark:3, c-quark:4, b-quark: 5
            df = df.Define("jetPID", f"AlephSelection::getJetPID(ClassBitset, {coll['GenParticles']})")
            df = df.Filter(f"jetPID == {self.ana_args.MCflavour}")
        
        # store the classbitset in the output
        df = df.Define("event_class", "AlephSelection::bitsetToIndices(ClassBitset)")
        df = df.Define("event_number", "EventHeader.eventNumber")
        df = df.Define("run_number", "EventHeader.runNumber")

        # Define RP kinematics
        ####################################################################################################
        df = df.Define("RP_px", "ReconstructedParticle::get_px(RecoParticles)")
        df = df.Define("RP_py", "ReconstructedParticle::get_py(RecoParticles)")
        df = df.Define("RP_pz", "ReconstructedParticle::get_pz(RecoParticles)")
        df = df.Define("RP_e", "ReconstructedParticle::get_e(RecoParticles)")
        df = df.Define("RP_m", "ReconstructedParticle::get_mass(RecoParticles)")

        # Define pseudo-jets
        ####################################################################################################
        df = df.Define("pjetc", "JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")

        # Anti-kt clustering and jet constituents
        ####################################################################################################
        df = df.Define("_jet", "JetClustering::clustering_ee_kt(2, 2, 1, 0)(pjetc)")
        df = df.Define("jets","JetClusteringUtils::get_pseudoJets(_jet)" )
        df = df.Define("_jetc", "JetClusteringUtils::get_constituents(_jet)") 
        df = df.Define("jetc", "JetConstituentsUtils::build_constituents_cluster(RecoParticles, _jetc)")
        df = df.Define("jetConstitutentsTypes", f"AlephSelection::build_constituents_Types()(ParticleID, _jetc)")
        df = df.Define("JetClustering_d23", "std::sqrt(JetClusteringUtils::get_exclusive_dmerge(_jet, 2))")
        df = df.Define("JetClustering_d34", "std::sqrt(JetClusteringUtils::get_exclusive_dmerge(_jet, 3))")

        ############################################# Event Level Variables #######################################################
        df = df.Define("jet_p4", "JetConstituentsUtils::compute_tlv_jets(jets)" )
        df = df.Define("event_invariant_mass", "JetConstituentsUtils::InvariantMass(jet_p4[0], jet_p4[1])")


        # Beamspot CENTRE, per run, in 10um units (see AlephSelection::get_beamspot).
        # At the origin in simulation, offset by ~0.6/0.2 mm in x/y in data.
        # The json path is passed explicitly (condor workers may not read AFS);
        # Override the path with $ALEPH_BEAMSPOT_JSON.
        if self.ana_args.doData:
            beamspot_json = os.environ.get(
                "ALEPH_BEAMSPOT_JSON",
                "/eos/experiment/fcc/ee/analyses/case-studies/aleph/utils/beamspot_position_data/beamspot.json")
            df = df.Define("BeamspotVec", 'AlephSelection::get_beamspot(run_number[0], true, "{}")'.format(beamspot_json))
            df = df.Define("Beamspot_x", "BeamspotVec.X()")
            df = df.Define("Beamspot_y", "BeamspotVec.Y()")
            df = df.Define("Beamspot_z", "BeamspotVec.Z()")
        else:
            df = df.Define("Beamspot_x", "0.0")
            df = df.Define("Beamspot_y", "0.0")
            df = df.Define("Beamspot_z", "0.0")

        # ==== Track selection (to harmonize with Luka's code)
        # Note: The selection strategy here only works if there is one trackstate stored pre track.
        # The code includes an assertion for that, if it is somehow not the case it will fail. 
        # df = df.Define("n_tracks_all", f"AlephSelection::select_tracks( {coll['PFTracks']} )")
        df = df.Define("n_tracks_all", "Tracks.size()")
        df = df.Define("chi2_tracks_all","AlephSelection::get_track_chi2( Tracks )") #TODO: use collection here
        df = df.Define("ndf_tracks_all","AlephSelection::get_track_ndf( Tracks )") #TODO: use collection here
        df = df.Define("chi2_o_ndf_tracks_all","AlephSelection::get_track_chi2_o_ndf( Tracks )") #TODO: use collection here
        
        # baseline track selection: positive definite cov matrix & chi2 < 10 
        df = df.Define("tracks_selected_baseline_result","AlephSelection::select_tracks_baseline( Tracks, _Tracks_trackStates )") #TODO: use collection here  0.75, 2.0
        df = df.Define("tracks_selected_baseline","tracks_selected_baseline_result.tracks") 
        df = df.Define("trackstates_selected_baseline","tracks_selected_baseline_result.trackStates") 

        # Upper bounds on the impact parameters, pre-selecting tracks for the PV fit.
        # Referenced to the run beamspot (Beamspot_* are in 10um units -> cm); with
        # --oldPV to the origin, off-centre by the beamspot offset in data.
        ip_window = "{0}::PVN_D0_MAX, {0}::PVN_Z0_MAX".format(PVNEW)
        if self.do_pvnew:
            df = df.Define("tracks_selected_for_vertexfit_result","AlephSelection::select_tracks_impactparameters_bs( tracks_selected_baseline_result, {}, Beamspot_x*1e-3, Beamspot_y*1e-3, Beamspot_z*1e-3 )".format(ip_window)) 
        else:
            df = df.Define("tracks_selected_for_vertexfit_result","AlephSelection::select_tracks_impactparameters( tracks_selected_baseline_result, {} )".format(ip_window)) 
        df = df.Define("tracks_selected_for_vertexfit","tracks_selected_for_vertexfit_result.tracks") 
        df = df.Define("trackstates_selected_for_vertexfit","tracks_selected_for_vertexfit_result.trackStates") 

        df = df.Define("n_tracks_sel", "tracks_selected_baseline.size()")
        df = df.Define("n_trackstates_sel", "trackstates_selected_baseline.size()") #for debug

        df = df.Define("n_tracks_sel_vertexfit", "tracks_selected_for_vertexfit.size()")

        # need to flip the sign of d0 and omega (why??)
        df = df.Define("trackstates_selected_for_vertexfit_flipped","AlephSelection::flipD0_copy(trackstates_selected_for_vertexfit )")
        df = df.Define("trackstates_selected_baseline_flipped","AlephSelection::flipD0_copy(trackstates_selected_baseline )")

        # ===== VERTEX

        # run primary vertex fit using FCCAna native fitter

        # Beam-spot constraint widths and the track-compatibility chi2max: both PV
        # chains read them from analyzer_pvnew.h, in cm.
        bs_sig_cm = "{0}::PVN_BS_SIGMA_X, {0}::PVN_BS_SIGMA_Y, {0}::PVN_BS_SIGMA_Z".format(PVNEW)
        chi2max = "{}::PVN_CHI2_MAX".format(PVNEW)
        # the legacy chain reads lengths in the raw Beamspot_* unit (10 um = 1e-3 cm): cm x 1e3 for all three
        bs_sig_legacy = "{0}::PVN_BS_SIGMA_X*1e3, {0}::PVN_BS_SIGMA_Y*1e3, {0}::PVN_BS_SIGMA_Z*1e3".format(PVNEW)

        if self.do_pvnew:
            # Standalone PV fitter (analyzer_pvnew.h), all lengths in cm: the
            # selection fit and the final fit share the same beam-spot constraint.
            bs_cm = ("FCCAnalyses::AlephPVNew::BeamSpot{{Beamspot_x*1e-3, Beamspot_y*1e-3, "
                     "Beamspot_z*1e-3, {}}}").format(bs_sig_cm)
            df = df.Define("PVSelNew", "FCCAnalyses::AlephPVNew::select_primary_tracks(trackstates_selected_for_vertexfit_flipped, {}, {})".format(bs_cm, chi2max))
            # Two flags: the selection fit can fail independently of the position
            # fit, so one flag cannot cover both. int-typed.
            df = df.Define("pv_converged",       "int(PVSelNew.fit.converged)")
            df = df.Define("pv_split_converged", "int(PVSelNew.split_converged)")
            # fewer than 2 IP-preselected tracks entered the pruning: the fit
            # "converges" at/near the beam spot with no track information, so
            # both flags above can still read 1. int-typed.
            df = df.Define("pv_trivial",         "int(PVSelNew.trivial)")
            # the three flags combined into the "usable PV" predicate, from the
            # single named source in analyzer_pvnew.h. Also the consumer guard
            # (finder entry ternaries, Vertex_refit_tlv); the three raw flags
            # stay stored as diagnostics.
            df = df.Define("pv_good",            "int(FCCAnalyses::AlephPVNew::goodPV(PVSelNew))")
            # split from the pruning when it converged, else the
            # beamspot-as-fixed-PV fallback (never the unpruned return)
            df = df.Define("RecoedPrimaryTracks_looseBS", "FCCAnalyses::AlephPVNew::primaryTracksFromSel(trackstates_selected_for_vertexfit_flipped, PVSelNew, Beamspot_x*1e-3, Beamspot_y*1e-3, Beamspot_z*1e-3, {})".format(chi2max))
            # position always written (the garbage IS the diagnostic),
            # covariance zeroed on non-convergence
            df = df.Define("VertexObject_looseBS", "FCCAnalyses::AlephPVNew::toFCCVertex(PVSelNew)")
            df = df.Define("Vertex_refit_looseBS", "VertexObject_looseBS.vertex")
            # jet-level IP variables fail open with huge finite values under
            # a garbage PV -> substitute the beam-spot position on the flag
            df = df.Define("Vertex_refit_tlv", "pv_good ? TLorentzVector(Vertex_refit_looseBS.position.x, Vertex_refit_looseBS.position.y, Vertex_refit_looseBS.position.z, 0.) : TLorentzVector(Beamspot_x*1e-3, Beamspot_y*1e-3, Beamspot_z*1e-3, 0.)")
        else:
            # Guard: with fewer than 2 IP-preselected tracks there is no meaningful primary vertex,
            # so return NO primary tracks (the PV fit then falls back to the dummy beamspot vertex).
            # FCCAnalyses' get_PrimaryTracks instead returns `seltracks` unchanged, i.e. the single
            # track - that is what the reference wrapper (getPrimaryTracks in analyzer_pvtools.cxx,
            # `if(tracksToUse.size() < 2){ return primaryTracks; }`) guards against. Without this we
            # get nPrim=1 where the reference has nPrim=0 (~1400 events / 1.05M in the full sweep).
            # note: the {{}} is an escaped literal {} for str.format - it is the empty RVec, not a placeholder
            df = df.Define("RecoedPrimaryTracks_looseBS", "trackstates_selected_for_vertexfit_flipped.size() < 2 ? ROOT::VecOps::RVec<edm4hep::TrackState>{{}} : VertexFitterSimple::get_PrimaryTracks(trackstates_selected_for_vertexfit_flipped, true, {}, Beamspot_x, Beamspot_y, Beamspot_z, {})".format(bs_sig_legacy, chi2max))
            df = df.Define("VertexObject_looseBS", "VertexFitterSimple::VertexFitter_Tk(1, RecoedPrimaryTracks_looseBS, true, {}, Beamspot_x, Beamspot_y, Beamspot_z)".format(bs_sig_legacy))
            df = df.Define("Vertex_refit_looseBS", "VertexingUtils::get_VertexData(VertexObject_looseBS)")
            df = df.Define("Vertex_refit_tlv", "TLorentzVector(Vertex_refit_looseBS.position.x, Vertex_refit_looseBS.position.y, Vertex_refit_looseBS.position.z, 0.)")
        # for retrieving secondary tracks, use the full list of selected tracks 
        df = df.Define("SecondaryTracks_looseBS", "VertexFitterSimple::get_NonPrimaryTracks(trackstates_selected_baseline_flipped, RecoedPrimaryTracks_looseBS)")

        # original-Tracks index maps for the primary/secondary splits (truth-free
        # track-state matching, written for data too). sec2origIdx is the join
        # that v0n reco_ind / svn_trk_idx need to reach the original tracks.
        df = df.Define("selBaselineOrigIdx", "FCCAnalyses::AlephTruth::selectedBaselineOriginalIndices(Tracks, _Tracks_trackStates, trackstates_selected_baseline)")
        df = df.Define("sec2origIdx",        "FCCAnalyses::AlephTruth::secondaryToOriginalTrack(SecondaryTracks_looseBS, trackstates_selected_baseline_flipped, selBaselineOrigIdx)")
        df = df.Define("prim2origIdx",       "FCCAnalyses::AlephTruth::secondaryToOriginalTrack(RecoedPrimaryTracks_looseBS, trackstates_selected_baseline_flipped, selBaselineOrigIdx)")

        # track<->pfcand join: the RecoParticles->Tracks relation, flattened
        # (recopart_tracks_index[recopart_tracks_begin[i]] = track index of
        # reco particle i when it has a track; begin==end for neutrals)
        df = df.Define("recopart_tracks_index", "_RecoParticles_tracks.index")
        df = df.Define("recopart_tracks_begin", "RecoParticles.tracks_begin")
        df = df.Define("recopart_tracks_end",   "RecoParticles.tracks_end")

        df = df.Define("Vertex_refit_x", "Vertex_refit_looseBS.position.x")
        df = df.Define("Vertex_refit_y", "Vertex_refit_looseBS.position.y")
        df = df.Define("Vertex_refit_z", "Vertex_refit_looseBS.position.z")

        # PV fit covariance (lower-triangular xx, yx, yy, zx, zy, zz)
        df = df.Define("Vertex_refit_cov_xx", "Vertex_refit_looseBS.covMatrix.values[0]")
        df = df.Define("Vertex_refit_cov_yx", "Vertex_refit_looseBS.covMatrix.values[1]")
        df = df.Define("Vertex_refit_cov_yy", "Vertex_refit_looseBS.covMatrix.values[2]")
        df = df.Define("Vertex_refit_cov_zx", "Vertex_refit_looseBS.covMatrix.values[3]")
        df = df.Define("Vertex_refit_cov_zy", "Vertex_refit_looseBS.covMatrix.values[4]")
        df = df.Define("Vertex_refit_cov_zz", "Vertex_refit_looseBS.covMatrix.values[5]")

        # PV fit quality (chi2/ndf as the fitter stores it): a silently
        # non-converged fit sits orders of magnitude above any genuine vertex.
        df = df.Define("Vertex_refit_chi2", "Vertex_refit_looseBS.chi2")

        df = df.Define("n_primary_tracks", "ReconstructedParticle2Track::getTK_n(RecoedPrimaryTracks_looseBS)")
        df = df.Define("n_secondary_tracks", "ReconstructedParticle2Track::getTK_n(SecondaryTracks_looseBS)")

        # for reference: vertex as stored - can be removed?
        # guarded: the Vertices.size()>0 filter below is disabled (Luka does not apply it), so this
        # must not index an empty collection. Currently 'pv' is not snapshotted and RDF never
        # evaluates it, but keep the guard so adding it to the output later cannot crash.
        df = df.Define(
            "pv",
            "Vertices.size() > 0 ? TLorentzVector(Vertices[0].position.x, Vertices[0].position.y, Vertices[0].position.z, 0.0) : TLorentzVector(0., 0., 0., 0.)",
        )
        df = df.Define("VertexX", "Vertices.position.x")
        df = df.Define("VertexY", "Vertices.position.y")
        df = df.Define("VertexZ", "Vertices.position.z")

        # TEST FILTER         
        # df = df.Filter("Vertices.size() > 0")  # to remove eventually


        # gen level vertex for checks, fill dummies for data
        if self.ana_args.doData:
            df = df.Define("gen_vertex_x", "-999")
            df = df.Define("gen_vertex_y", "-999")
            df = df.Define("gen_vertex_z", "-999")

            # refit vertex resolution:
            df = df.Define("res_vertex_x", "-999")
            df = df.Define("res_vertex_y", "-999")
            df = df.Define("res_vertex_z", "-999")
        
        else:
            df = df.Define("pv_gen_level", f'AlephSelection::get_EventPrimaryVertexP4()({coll["GenParticles"]})')
            df = df.Define("gen_vertex_x", "pv_gen_level.X()")
            df = df.Define("gen_vertex_y", "pv_gen_level.Y()")
            df = df.Define("gen_vertex_z", "pv_gen_level.Z()")

            # refit vertex resolution:
            df = df.Define("res_vertex_x", "Vertex_refit_x - gen_vertex_x")
            df = df.Define("res_vertex_y", "Vertex_refit_y - gen_vertex_y")
            df = df.Define("res_vertex_z", "Vertex_refit_z - gen_vertex_z")

        # check without track selection
        # df = df.Define("res_vertex_x_all_tracks", "Vertex_refit_x_all_tracks - gen_vertex_x")
        # df = df.Define("res_vertex_y_all_tracks", "Vertex_refit_y_all_tracks - gen_vertex_y")
        # df = df.Define("res_vertex_z_all_tracks", "Vertex_refit_z_all_tracks - gen_vertex_z")

        ############################################# Secondary Vertices #######################################################
        # first we find the secondary vertices per event ...        
        old_sv_expr = ("FCCAnalyses::AlephSelection::get_SV_event_ALEPH("
            "SecondaryTracks_looseBS, "               # non-primary tracks
            "trackstates_selected_baseline_flipped, " # all tracks
            "VertexObject_looseBS, "                  # primary vertex
            "0.8, "                                   # dR prefilter cut
            "false)"                                  # exclusive V0 rejection (skip+break), matching FCCAnalyses@3a4de97 isV0 - the code that produced ntuples-withks
        )
        if self.do_pvnew:
            # the old LCFIPlus finder fails OPEN under a garbage PV (its only
            # PV-dependent cut is an angle<0 rejection) -> hard skip on the flag
            old_sv_expr = self._pv_guard(
                old_sv_expr,
                "ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex>{}")
        df = df.Define("SVs_looseBS", old_sv_expr)

        #.. then we assign them to the closest jet based on dR (also tracks to be moved between jets, in contrast to using get_SV_jet ! )
        df = df.Define("sv_jets", "FCCAnalyses::AlephSelection::assign_SV_to_jets(SVs_looseBS, jets)")

        # secondary vertex multiplicities
        df = df.Define("n_sv_event", "int(SVs_looseBS.size())")
        df = df.Define("n_sv_jets",  "FCCAnalyses::VertexingUtils::get_n_SV_jets(sv_jets)")

        # secondary vertex  properties
        df = df.Define("sv_chi2",        "FCCAnalyses::VertexingUtils::get_chi2_SV(sv_jets)")
        df = df.Define("sv_chi2_norm",   "FCCAnalyses::VertexingUtils::get_norm_chi2_SV(sv_jets)")
        df = df.Define("sv_ndof",        "FCCAnalyses::VertexingUtils::get_nDOF_SV(sv_jets)")
        df = df.Define("sv_ntracks",     "FCCAnalyses::VertexingUtils::get_VertexNtrk(sv_jets)")
        df = df.Define("sv_mass",        "FCCAnalyses::VertexingUtils::get_invM(sv_jets)")
        df = df.Define("sv_p",           "FCCAnalyses::VertexingUtils::get_pMag_SV(sv_jets)")
        df = df.Define("sv_thetarel",    "FCCAnalyses::VertexingUtils::get_relTheta_SV(sv_jets, jets)")
        df = df.Define("sv_phirel",      "FCCAnalyses::VertexingUtils::get_relPhi_SV(sv_jets, jets)")
        df = df.Define("sv_dxy",         "FCCAnalyses::VertexingUtils::get_dxy_SV(sv_jets, VertexObject_looseBS)")
        df = df.Define("sv_dxyz",        "FCCAnalyses::VertexingUtils::get_d3d_SV(sv_jets, VertexObject_looseBS)")
        # for pointing angle, use custom defined function following luka's code
        df = df.Define("sv_cosPointing",    "FCCAnalyses::AlephSelection::get_pointingangle_SV(sv_jets, VertexObject_looseBS)")
        df = df.Define("sv_prel",           "FCCAnalyses::AlephSelection::get_prel_SV_jets(sv_jets, jets)")
        df = df.Define("sv_correctedMass",  "FCCAnalyses::AlephSelection::get_correctedInvMass_SV(sv_jets, VertexObject_looseBS)")

        # displacement of SVs wrt to primary vertex
        df = df.Define("PrimaryVertexP3",
             "TVector3(VertexObject_looseBS.vertex.position[0], "
             "VertexObject_looseBS.vertex.position[1], "
             "VertexObject_looseBS.vertex.position[2])")
        df = df.Define("sv_dx", "FCCAnalyses::AlephSelection::get_dx_SV_jets(sv_jets, PrimaryVertexP3)")
        df = df.Define("sv_dy", "FCCAnalyses::AlephSelection::get_dy_SV_jets(sv_jets, PrimaryVertexP3)")
        df = df.Define("sv_dz", "FCCAnalyses::AlephSelection::get_dz_SV_jets(sv_jets, PrimaryVertexP3)")
        # legacy-SV vertex-fit covariance (nested per jet like the other sv_*
        # branches; packed lower triangle, same component order as the
        # Vertex_refit_cov_* branches)
        for ic, cc in enumerate(("xx", "yx", "yy", "zx", "zy", "zz")):
            df = df.Define(f"sv_cov_{cc}", f"FCCAnalyses::AlephSelection::svCovComp(sv_jets, {ic})")

        ############################################# V0 Reconstruction #######################################################
        df = df.Define("V0s_event",
            "FCCAnalyses::AlephSelection::get_V0s_ALEPH("
            "SecondaryTracks_looseBS, "
            "VertexObject_looseBS,"
            f"{BZ}," #solenoidBz
            "true," #loose_mass_window
            "-1.," #dR preselection on track pairs (<=0 disables) - 0.4 tested, made it much worse
            "true)" #exclusive tracks (each track in at most one V0) - TESTING against ntuples-withks
        )
        df = df.Define("v0s_per_jet", "FCCAnalyses::AlephSelection::assign_V0s_to_jets(V0s_event, jets)")
        df = df.Define("v0_jets",  "v0s_per_jet.vtx")
        df = df.Define("v0_pdg",   "v0s_per_jet.pdgAbs")
        df = df.Define("v0_invM",  "v0s_per_jet.invM")
        df = df.Define("n_v0_event",   "int(V0s_event.vtx.size())")
        df = df.Define("n_v0_jets",    "FCCAnalyses::VertexingUtils::get_n_SV_jets(v0_jets)")
        df = df.Define("n_v0_ks",      "FCCAnalyses::AlephSelection::count_V0type_jets(v0_pdg, 310)")
        df = df.Define("n_v0_lambda",  "FCCAnalyses::AlephSelection::count_V0type_jets(v0_pdg, 3122)")
        df = df.Define("v0_chi2",          "FCCAnalyses::VertexingUtils::get_chi2_SV(v0_jets)")
        df = df.Define("v0_chi2_norm",     "FCCAnalyses::VertexingUtils::get_norm_chi2_SV(v0_jets)")
        df = df.Define("v0_ndof",          "FCCAnalyses::VertexingUtils::get_nDOF_SV(v0_jets)")
        df = df.Define("v0_ntracks",       "FCCAnalyses::VertexingUtils::get_VertexNtrk(v0_jets)")
        df = df.Define("v0_p",             "FCCAnalyses::VertexingUtils::get_pMag_SV(v0_jets)")
        df = df.Define("v0_prel",          "FCCAnalyses::AlephSelection::get_prel_SV_jets(v0_jets, jets)")
        df = df.Define("v0_thetarel",      "FCCAnalyses::VertexingUtils::get_relTheta_SV(v0_jets, jets)")
        df = df.Define("v0_phirel",        "FCCAnalyses::VertexingUtils::get_relPhi_SV(v0_jets, jets)")
        df = df.Define("v0_dxy",           "FCCAnalyses::VertexingUtils::get_dxy_SV(v0_jets, VertexObject_looseBS)")
        df = df.Define("v0_dxyz",          "FCCAnalyses::VertexingUtils::get_d3d_SV(v0_jets, VertexObject_looseBS)")
        df = df.Define("v0_cosPointing",   "FCCAnalyses::AlephSelection::get_pointingangle_SV(v0_jets, VertexObject_looseBS)")
        df = df.Define("v0_correctedMass", "FCCAnalyses::AlephSelection::get_correctedInvMass_SV(v0_jets, VertexObject_looseBS)")
        df = df.Define("v0_dx",  "FCCAnalyses::AlephSelection::get_dx_SV_jets(v0_jets, PrimaryVertexP3)")
        df = df.Define("v0_dy",  "FCCAnalyses::AlephSelection::get_dy_SV_jets(v0_jets, PrimaryVertexP3)")
        df = df.Define("v0_dz",  "FCCAnalyses::AlephSelection::get_dz_SV_jets(v0_jets, PrimaryVertexP3)")

        ############################################# V0 truth matching (MC only) #############################################
        if self.do_truth:
            # many-to-many track<->MC maps from the (non-empty) trackMCLink ObjectIDs
            df = df.Define("mcToTracks",  f"FCCAnalyses::AlephTruth::buildMCToTracks({coll['GenParticles']}.size(), _trackMCLink_from, _trackMCLink_to)")
            df = df.Define("trackToMCs",  "FCCAnalyses::AlephTruth::buildTrackToMCs(Tracks.size(), _trackMCLink_from, _trackMCLink_to)")
            # mother-anchored true V0s (geometric daughter recovery, cm units)
            df = df.Define("trueV0s",     f"FCCAnalyses::AlephTruth::findTrueV0s({coll['GenParticles']}, mcToTracks)")
            # (selBaselineOrigIdx / sec2origIdx are defined unconditionally in the
            # PV block above — truth-free track-state matching, available on data)
            # recover which track pair each candidate came from (compiled get_V0s leaves reco_ind empty);
            # classifyV0s cross-checks this replica against V0s_event pdg/invM and throws on mismatch
            df = df.Define("v0pairs",       f"FCCAnalyses::AlephTruth::rerunV0Pairing(SecondaryTracks_looseBS, VertexObject_looseBS, {BZ})")
            # truth classification of the reco V0 candidates (event order = V0s_event order)
            df = df.Define("v0truth",       f"FCCAnalyses::AlephTruth::classifyV0s(V0s_event, v0pairs, SecondaryTracks_looseBS, sec2origIdx, trackToMCs, {coll['GenParticles']}, trueV0s)")
            for _b, _e in TRUEV0_DEFINES:
                df = df.Define(f"truev0_{_b}", _e)
            for _b, _e in V0C_DEFINES:
                df = df.Define(f"v0c_{_b}", _e)

        # track -> ReconstructedParticle join, for the per-leg PF label and the
        # pfcand join key (begin != end is the "has a track" test)
        df = df.Define("rpOfTrack",
                       "FCCAnalyses::AlephTrkAux::rpIndexByTrack(RecoParticles.tracks_begin, RecoParticles.tracks_end, _RecoParticles_tracks.index, Tracks.size())")
        # dE/dx track->measurement join, shared by every daughter-leg block
        if self.do_v0new or self.do_phikk or self.do_dstar:
            df = self._define_dedx_join(df)

        ############################################# Standalone two-tier V0 module ###########################################
        if self.do_v0new:
            v0n_expr = f"FCCAnalyses::AlephV0New::findV0s(SecondaryTracks_looseBS, VertexObject_looseBS, {BZ})"
            if self.do_pvnew:
                # explicit empty-return entry guard on the flag (the window
                # cuts would empty it anyway, a silent efficiency loss; the
                # guard makes the failure explicit and empties pointSig too)
                v0n_expr = self._pv_guard(
                    v0n_expr, "FCCAnalyses::VertexingUtils::FCCAnalysesV0{}")
            df = df.Define("V0sNew_event", v0n_expr)
            df = df.Define("n_v0n_event",  "int(V0sNew_event.vtx.size())")
            # truth-free kinematic branches (available on data)
            for _b, _e in V0N_CAND_DEFINES:
                df = df.Define(f"v0n_{_b}", _e)
            # per-daughter joins + dE/dx (truth-free: reco_ind -> sec2origIdx).
            # dE/dx validity: value != omega(track) (failed-leg sentinel),
            # finite positive value and error; invalid -> -1 in both branches.
            for _i, _t in enumerate(V0N_TRKS):
                df = df.Define(f"v0n_{_t}_origIdx",
                               f"FCCAnalyses::AlephV0New::candDaughterOrigIdx(V0sNew_event, sec2origIdx, {_i})")
            df = self._define_dedx(df, [f"v0n_{_t}" for _t in V0N_TRKS])
            df = self._define_leg_pid(df, [f"v0n_{_t}" for _t in V0N_TRKS])
            # per-jet new-module V0s: mirror of the old-finder v0_* block on
            # V0sNew_event, so v0_* vs v0njet_* is an apples-to-apples comparison
            # at jet level (prel = pT wrt jet axis, thetarel/phirel wrt the jet).
            df = df.Define("v0njet_per_jet", "FCCAnalyses::AlephSelection::assign_V0s_to_jets(V0sNew_event, jets)")
            df = df.Define("v0njet_jets",  "v0njet_per_jet.vtx")
            for _b, _e in V0NJET_DEFINES:
                df = df.Define(f"v0njet_{_b}", _e)
            df = df.Define("n_v0njet_jets",    "FCCAnalyses::VertexingUtils::get_n_SV_jets(v0njet_jets)")
            df = df.Define("n_v0njet_ks",      "FCCAnalyses::AlephSelection::count_V0type_jets(v0njet_pdg, 310)")
            df = df.Define("n_v0njet_lambda",  "FCCAnalyses::AlephSelection::count_V0type_jets(v0njet_pdg, 3122)")
            # truth classification (MC only; reco_ind is filled by the new module)
            if self.do_truth:
                df = df.Define("v0npairs",     "FCCAnalyses::AlephTruth::pairsFromRecoInd(V0sNew_event)")
                df = df.Define("v0ntruth",     f"FCCAnalyses::AlephTruth::classifyV0s(V0sNew_event, v0npairs, SecondaryTracks_looseBS, sec2origIdx, trackToMCs, {coll['GenParticles']}, trueV0s)")
                for _n, _e in V0N_TRUTH_DEFINES:
                    df = df.Define(_n, _e)

        ############################################# Standalone SV module ####################################################
        if self.do_svnew:
            # V0-first: svn_* = SV finding after masking the tight-claimed V0 tracks;
            # svm_* = unmasked control twin from the SAME event, for the interplay study.
            SVNEW = "FCCAnalyses::AlephSVNew"
            # two-track seed pass, shared by both masking modes
            seed_expr = f"{SVNEW}::svSeedPass(SecondaryTracks_looseBS, VertexObject_looseBS, {BZ})"
            if self.do_pvnew:
                seed_expr = self._pv_guard(seed_expr, f"{SVNEW}::SVSeeds{{}}")
            df = df.Define("SVSeeds_event", seed_expr)
            for pfx, mode in (("svn", f"{SVNEW}::SVN_MASK_MODE"), ("svm", f"{SVNEW}::SVN_MASK_NONE")):
                svn_expr = f"{SVNEW}::findSVs(SecondaryTracks_looseBS, VertexObject_looseBS, V0sNew_event, v0n_tight, {mode}, {BZ}, SVSeeds_event)"
                if self.do_pvnew:
                    # explicit entry guard (see V0sNew_event above)
                    svn_expr = self._pv_guard(
                        svn_expr, "FCCAnalyses::VertexingUtils::FCCAnalysesV0{}")
                df = df.Define(f"SVs_{pfx}", svn_expr)
                df = df.Define(f"n_{pfx}_event",    f"int(SVs_{pfx}.vtx.size())")
                df = df.Define(f"{pfx}_mass",        f"SVs_{pfx}.invM")
                df = df.Define(f"{pfx}_chi2",        f"FCCAnalyses::AlephTruth::candChi2(SVs_{pfx})")
                df = df.Define(f"{pfx}_dxyz",        f"FCCAnalyses::AlephTruth::candDxyz(SVs_{pfx}, VertexObject_looseBS)")
                df = df.Define(f"{pfx}_p",           f"FCCAnalyses::AlephTruth::candP(SVs_{pfx})")
                df = df.Define(f"{pfx}_cosPointing", f"FCCAnalyses::AlephTruth::candCosPointing(SVs_{pfx}, VertexObject_looseBS)")
                df = df.Define(f"{pfx}_pointSig",    f"FCCAnalyses::AlephV0New::candPointSig(SVs_{pfx}, VertexObject_looseBS)")
                df = df.Define(f"{pfx}_ntracks",     f"FCCAnalyses::AlephSVNew::candNtracks(SVs_{pfx})")
                df = df.Define(f"{pfx}_sigL",        f"FCCAnalyses::AlephSVNew::candSigL(SVs_{pfx})")
                df = df.Define(f"{pfx}_trk_sv",      f"FCCAnalyses::AlephSVNew::candTrkSV(SVs_{pfx})")
                df = df.Define(f"{pfx}_trk_idx",     f"FCCAnalyses::AlephSVNew::candTrkIdx(SVs_{pfx})")
                # displacement VECTOR wrt PV (dxyz is only the magnitude) -> offline
                # position matching against the true SV positions
                for ic, cc in enumerate("xyz"):
                    df = df.Define(f"{pfx}_d{cc}", f"FCCAnalyses::AlephSVNew::candDcomp(SVs_{pfx}, VertexObject_looseBS, {ic})")
                # SV vertex-fit covariance (packed lower triangle, cm^2, same
                # component order as Vertex_refit_cov_*)
                for ic, cc in enumerate(("xx", "yx", "yy", "zx", "zy", "zz")):
                    df = df.Define(f"{pfx}_cov_{cc}", f"FCCAnalyses::AlephV0New::candCovComp(SVs_{pfx}, {ic})")

            # V0-candidate pointing at the nearest svn vertex (largest cosine).
            # Feature only; SVs sharing a daughter track are excluded
            # (self-pointing). Sentinels cos=-2, sig=-1, idx=-1.
            df = df.Define("v0n_svnpoint",    "FCCAnalyses::AlephV0New::candSVPointing(V0sNew_event, SVs_svn, sec2origIdx)")
            df = df.Define("v0n_svnCosPoint", "v0n_svnpoint.cosPoint")
            df = df.Define("v0n_svnPointSig", "v0n_svnpoint.pointSig")
            df = df.Define("v0n_svnIdx",      "v0n_svnpoint.svIdx")

        ############################################# exclusive-finder track auxiliaries ######################################
        if self.do_phikk or self.do_dstar:
            # Both finders build candidates from the FULL baseline-selected
            # track list: primary and secondary tracks alike, no masking by the
            # PV split. selBaselineOrigIdx maps that collection to the original
            # Tracks, which is how these per-track auxiliaries are joined.
            TRKAUX = "FCCAnalyses::AlephTrkAux"
            df = df.Define("trkaux_nvdet", f"{TRKAUX}::subdetHits(selBaselineOrigIdx, Tracks.subdetectorHitNumbers_begin, Tracks.subdetectorHitNumbers_end, _Tracks_subdetectorHitNumbers, 0)")
            df = df.Define("trkaux_nitc",  f"{TRKAUX}::subdetHits(selBaselineOrigIdx, Tracks.subdetectorHitNumbers_begin, Tracks.subdetectorHitNumbers_end, _Tracks_subdetectorHitNumbers, 1)")
            df = df.Define("trkaux_chi2ndf", f"{TRKAUX}::trackChi2Ndf(selBaselineOrigIdx, Tracks.chi2, Tracks.ndf)")
            df = df.Define("trkaux_isprim",  f"{TRKAUX}::flagInSet(selBaselineOrigIdx, prim2origIdx)")
            # tracks claimed as daughters of a tight Ks/Lambda leave BOTH pools;
            # the claim list is a V0-module product, so under --oldV0 it is empty
            if self.do_v0new:
                df = df.Define("v0n_claimed_orig", f"{TRKAUX}::claimedOrigIdx(v0n_trk1_origIdx, v0n_trk2_origIdx, v0n_tight)")
            else:
                df = df.Define("v0n_claimed_orig", "ROOT::VecOps::RVec<int>{}")

        ############################################# phi(1020) -> K+K- module ################################################
        if self.do_phikk:
            # every selection value is a constant in analyzer_phikk.h
            phikk_expr = ("FCCAnalyses::AlephPhiKK::findPhiKK(trackstates_selected_baseline_flipped, "
                          "selBaselineOrigIdx, trkaux_nvdet, trkaux_nitc, trkaux_chi2ndf, "
                          f"trkaux_isprim, VertexObject_looseBS, {BZ}, v0n_claimed_orig, "
                          "Beamspot_x*1e-3, Beamspot_y*1e-3, Beamspot_z*1e-3)")
            if self.do_pvnew:
                # explicit entry guard (see V0sNew_event above)
                phikk_expr = self._pv_guard(
                    phikk_expr, "FCCAnalyses::AlephPhiKK::PhiKKCands{}")
            df = df.Define("PhiKKCands_event", phikk_expr)
            df = df.Define("n_phikk_event", "int(PhiKKCands_event.invM.size())")
            for _b in PHIKK_CAND_BRANCHES:
                df = df.Define(f"phikk_{_b}", f"PhiKKCands_event.{_b}")
            for _t in PHIKK_TRKS:
                for _b in PHIKK_TRK_BRANCHES:
                    df = df.Define(f"phikk_{_t}_{_b}", f"PhiKKCands_event.{_t}.{_b}")
            df = self._define_dedx(df, [f"phikk_{_t}" for _t in PHIKK_TRKS])
            df = self._define_leg_pid(df, [f"phikk_{_t}" for _t in PHIKK_TRKS])
            if self.do_truth:
                df = df.Define("truePhis", f"FCCAnalyses::AlephPhiKK::findTruePhis({coll['GenParticles']}, mcToTracks)")
                for _b in TRUEPHI_BRANCHES:
                    df = df.Define(f"truephi_{_b}", f"truePhis.{_b}")
                df = df.Define("n_truephi_event",    "int(truePhis.idx.size())")
                df = df.Define("phikktruth", f"FCCAnalyses::AlephPhiKK::classifyPhiKK(PhiKKCands_event, trackToMCs, {coll['GenParticles']}, truePhis)")
                df = df.Define("phikk_class",        "phikktruth.cls")
                df = df.Define("phikk_trueidx",      "phikktruth.truephi_idx")
                for _t in PHIKK_TRKS:
                    for _b in TRK_TRUTH_BRANCHES:
                        df = df.Define(f"phikk_{_t}_{_b}", f"phikktruth.{_t}_{_b}")
                df = df.Define("truephi_found", "FCCAnalyses::AlephPhiKK::truePhiFound(truePhis, phikktruth)")

        ############################################# D*->D0(K pi) pi_slow module #############################################
        if self.do_dstar:
            # primary/secondary class of every pool track (0 prim / 1 sec / 2 neither)
            df = df.Define("dstar_pool_all",    "FCCAnalyses::AlephDstar::poolClass(selBaselineOrigIdx, prim2origIdx, sec2origIdx)")
            # every selection value is a constant in analyzer_dstar.h
            dstar_expr = ("FCCAnalyses::AlephDstar::findDstar(trackstates_selected_baseline_flipped, "
                          "selBaselineOrigIdx, trkaux_nvdet, trkaux_nitc, trkaux_chi2ndf, "
                          "trkaux_isprim, dstar_pool_all, VertexObject_looseBS, v0n_claimed_orig, "
                          f"{BZ}, Beamspot_x*1e-3, Beamspot_y*1e-3, Beamspot_z*1e-3)")
            if self.do_pvnew:
                # explicit entry guard (see V0sNew_event above)
                dstar_expr = self._pv_guard(
                    dstar_expr, "FCCAnalyses::AlephDstar::DstarCands{}")
            df = df.Define("DstarCands_event", dstar_expr)
            df = df.Define("n_d0_event",    "int(DstarCands_event.d0.kin.m_kpi.size())")
            df = df.Define("n_dstar_event", "int(DstarCands_event.ds.kin.m_kpi.size())")
            # two-track fits actually performed (cache misses): the combinatorial cost
            df = df.Define("n_d0fits_event", "DstarCands_event.nfits")
            for _b in D0_CAND_BRANCHES:
                df = df.Define(f"d0_{_b}", f"DstarCands_event.{_cand_member('d0', _b)}")
            for _b in DSTAR_CAND_BRANCHES:
                df = df.Define(f"dstar_{_b}", f"DstarCands_event.{_cand_member('ds', _b)}")
            for _pfx, _mem in DSTAR_TRK_LEGS:
                for _b in DSTAR_TRK_BRANCHES:
                    df = df.Define(f"{_pfx}_{_b}", f"DstarCands_event.{_mem}.{_b}")
            df = self._define_dedx(df, [_pfx for _pfx, _ in DSTAR_TRK_LEGS])
            df = self._define_leg_pid(df, [_pfx for _pfx, _ in DSTAR_TRK_LEGS])
            if self.do_truth:
                # the primary/secondary pool of each true daughter, measured
                # through its linked track (prim2origIdx / sec2origIdx)
                df = df.Define("trueD0s",    f"FCCAnalyses::AlephDstar::findTrueD0s({coll['GenParticles']}, mcToTracks, prim2origIdx, sec2origIdx)")
                df = df.Define("trueDstars", f"FCCAnalyses::AlephDstar::findTrueDstars({coll['GenParticles']}, mcToTracks, prim2origIdx, sec2origIdx)")
                for _b in TRUED0_BRANCHES:
                    df = df.Define(f"trued0_{_b}", f"trueD0s.{_b}")
                for _b in TRUEDSTAR_BRANCHES:
                    df = df.Define(f"truedstar_{_b}", f"trueDstars.{_b}")
                df = df.Define("n_trued0_event",    "int(trueD0s.idx.size())")
                df = df.Define("n_truedstar_event", "int(trueDstars.idx.size())")
                df = df.Define("d0truth",    f"FCCAnalyses::AlephDstar::classifyD0(DstarCands_event.d0, trackToMCs, {coll['GenParticles']}, trueD0s)")
                df = df.Define("dstartruth", f"FCCAnalyses::AlephDstar::classifyDstar(DstarCands_event.ds, trackToMCs, {coll['GenParticles']}, trueDstars, trueD0s)")
                df = df.Define("d0_class",      "d0truth.cls")
                df = df.Define("d0_trueidx",    "d0truth.trueidx")
                df = df.Define("dstar_class",   "dstartruth.cls")
                df = df.Define("dstar_trueidx", "dstartruth.trueidx")
                for _t in D0_TRKS:
                    for _b in TRK_TRUTH_BRANCHES:
                        df = df.Define(f"d0_{_t}_{_b}", f"d0truth.{_t}_{_b}")
                for _t in DSTAR_TRKS:
                    for _b in TRK_TRUTH_BRANCHES:
                        df = df.Define(f"dstar_{_t}_{_b}", f"dstartruth.{_t}_{_b}")
                # per-true-particle efficiency flags (class-1 candidate carrying the label)
                df = df.Define("trued0_found_loose",    "FCCAnalyses::AlephDstar::trueD0Found(trueD0s, d0truth, d0_loose)")
                df = df.Define("trued0_found_tight",    "FCCAnalyses::AlephDstar::trueD0Found(trueD0s, d0truth, d0_tight)")
                df = df.Define("truedstar_found_loose", "FCCAnalyses::AlephDstar::trueDstarFound(trueDstars, dstartruth, dstar_loose)")
                df = df.Define("truedstar_found_tight", "FCCAnalyses::AlephDstar::trueDstarFound(trueDstars, dstartruth, dstar_tight)")

        ############################################# per-track membership ####################################################
        # One pass over the finished candidate lists: which sets and stored
        # candidates each ORIGINAL track belongs to. trk_nCand is the
        # multiplicity the offline 1/n de-duplication weight needs.
        _EMPTY = "ROOT::VecOps::RVec<int>{}"
        _v0 = ("v0n_trk1_origIdx, v0n_trk2_origIdx, v0n_tight" if self.do_v0new
               else f"{_EMPTY}, {_EMPTY}, {_EMPTY}")
        _phi = ("phikk_trk1_origIdx, phikk_trk2_origIdx, phikk_wp" if self.do_phikk
                else f"{_EMPTY}, {_EMPTY}, {_EMPTY}")
        _ds = ("d0_trkK_origIdx, d0_trkPi_origIdx, dstar_trkK_origIdx, "
               "dstar_trkPi_origIdx, dstar_trkPis_origIdx, dstar_tight"
               if self.do_dstar else ", ".join([_EMPTY] * 6))
        _sv = "svn_trk_idx" if self.do_svnew else _EMPTY
        df = df.Define("trkTags",
                       "FCCAnalyses::AlephTrkAux::trackTags(Tracks.size(), "
                       f"selBaselineOrigIdx, prim2origIdx, {_sv}, sec2origIdx, "
                       f"{_v0}, {_phi}, {_ds})")
        df = df.Define("trk_member", "trkTags.member")
        df = df.Define("trk_nCand",  "trkTags.nCand")

        ############################################# Particle Flow Level Variables #######################################################
        df = df.Define("pfcand_isMu",     "AlephSelection::get_isType(jetConstitutentsTypes,2)")
        df = df.Define("pfcand_isEl",     "AlephSelection::get_isType(jetConstitutentsTypes,1)")
        df = df.Define("pfcand_isGamma",  "AlephSelection::get_isType(jetConstitutentsTypes,4)")
        df = df.Define("pfcand_isChargedHad", f"AlephSelection::get_isType(jetConstitutentsTypes,{PF_CHARGED_HAD})")
        df = df.Define("pfcand_isNeutralHad", "AlephSelection::get_isType(jetConstitutentsTypes,5)")


        ############################################# Kinematics and PID #######################################################

        df = df.Define("pfcand_e",        "JetConstituentsUtils::get_e(jetc)") 
        df = df.Define("pfcand_p",        "JetConstituentsUtils::get_p(jetc)") 
        df = df.Define("pfcand_px",        "AlephSelection::get_px(jetc)")
        df = df.Define("pfcand_py",        "AlephSelection::get_py(jetc)")
        df = df.Define("pfcand_pz",        "AlephSelection::get_pz(jetc)")
        df = df.Define("pfcand_mask",        "AlephSelection::mask(pfcand_e)")


        df = df.Define("pfcand_theta",    "JetConstituentsUtils::get_theta(jetc)") 
        df = df.Define("pfcand_phi",      "JetConstituentsUtils::get_phi(jetc)") 
        df = df.Define("pfcand_charge",   "JetConstituentsUtils::get_charge(jetc)") 
        df = df.Define("pfcand_erel",     "JetConstituentsUtils::get_erel_cluster(jets, jetc)")
        df = df.Define("pfcand_erel_log", "JetConstituentsUtils::get_erel_log_cluster(jets, jetc)")
        df = df.Define("pfcand_thetarel", "JetConstituentsUtils::get_thetarel_cluster(jets, jetc)")
        df = df.Define("pfcand_phirel",   "JetConstituentsUtils::get_phirel_cluster(jets, jetc)")

        # transverse momentum: ptrel is the ratio pT_constituent / pT_jet (same convention as erel)
        df = df.Define("pfcand_pt",        "JetConstituentsUtils::get_pt(jetc)")
        df = df.Define("pfcand_ptrel",     "AlephSelection::get_ptrel_cluster(jets, jetc)")
        df = df.Define("pfcand_ptrel_log", "AlephSelection::get_ptrel_log_cluster(jets, jetc)")

        # re-index through the RecoParticle->Track relation so that .at(tracks_begin) picks the particle's own track
        df = df.Define("TracksByRP", "AlephSelection::reindexByRPLink(Tracks, _RecoParticles_tracks.index)")

        # track fit quality per constituent (-1 for neutrals, which have no track)
        df = df.Define("pfcand_trackChi2",     "AlephSelection::get_constituent_trackChi2(jetc, TracksByRP)")
        df = df.Define("pfcand_trackNdof",     "AlephSelection::get_constituent_trackNdof(jetc, TracksByRP)")
        df = df.Define("pfcand_trackChi2Norm", "AlephSelection::get_constituent_trackChi2Norm(jetc, TracksByRP)")

        # original-Tracks index of each constituent's own track (-1 = no track):
        # the join key between pfcand_* and the finders' *_origIdx branches
        df = df.Define("pfcand_trackIdx", "AlephSelection::get_constituent_trackIdx(jetc, _RecoParticles_tracks.index)")

        # subdetector hit counts per constituent (inside-out: VDET, ITC, TPC)
        df = df.Define("pfcand_nTrackHits_VDET", "AlephSelection::get_constituent_nTrackHits_VDET(jetc, TracksByRP, _Tracks_subdetectorHitNumbers)")
        df = df.Define("pfcand_nTrackHits_ITC",  "AlephSelection::get_constituent_nTrackHits_ITC(jetc, TracksByRP, _Tracks_subdetectorHitNumbers)")
        df = df.Define("pfcand_nTrackHits_TPC",  "AlephSelection::get_constituent_nTrackHits_TPC(jetc, TracksByRP, _Tracks_subdetectorHitNumbers)")

        df = df.Define("Bz", f"{BZ}")

        ############################################# Track Parameters and Covariance #######################################################

        df = df.Define("TrackStateFlipped",f"AlephSelection::flipD0_copy( {coll['TrackState']} )")
        # re-index through the RecoParticle->Track relation so that .at(tracks_begin) picks the particle's own track
        df = df.Define("TrackStateByRP", "AlephSelection::reindexByRPLink(TrackStateFlipped, _RecoParticles_tracks.index)")

        # dxy/dz/phi0 w.r.t. the PV, computed once with every curvature term in cm; C and ct
        # are the fitted curvature and dip angle (see analyzer.h)
        df = df.Define("pfcand_trkparPV",   "AlephSelection::get_constituent_trackParamsAtPV(jetc, TrackStateByRP, Vertex_refit_tlv, Bz)")
        df = df.Define("pfcand_dxy",        "pfcand_trkparPV.dxy")
        df = df.Define("pfcand_dz",         "pfcand_trkparPV.dz")
        df = df.Define("pfcand_phi0",       "pfcand_trkparPV.phi0")
        df = df.Define("pfcand_C",          "pfcand_trkparPV.C")
        df = df.Define("pfcand_ct",         "pfcand_trkparPV.ct")
        # track covariance, lower-triangular packing in the order (d0, phi0, omega, z0, tanLambda):
        # cov(a,b) with a >= b sits at index a*(a+1)/2 + b
        df = df.Define("pfcand_dptdpt",     "AlephSelection::get_constituent_trackCov(jetc, TrackStateByRP, 5)")
        df = df.Define("pfcand_dxydxy",     "AlephSelection::get_constituent_trackCov(jetc, TrackStateByRP, 0)")
        df = df.Define("pfcand_dzdz",       "AlephSelection::get_constituent_trackCov(jetc, TrackStateByRP, 9)")
        df = df.Define("pfcand_dphidphi",   "AlephSelection::get_constituent_trackCov(jetc, TrackStateByRP, 2)")
        df = df.Define("pfcand_detadeta",   "AlephSelection::get_constituent_trackCov(jetc, TrackStateByRP, 14)")
        df = df.Define("pfcand_dxydz",      "AlephSelection::get_constituent_trackCov(jetc, TrackStateByRP, 6)") # do we not need to recalculate this?
        df = df.Define("pfcand_dphidxy",    "AlephSelection::get_constituent_trackCov(jetc, TrackStateByRP, 1)")
        df = df.Define("pfcand_phidz",      "AlephSelection::get_constituent_trackCov(jetc, TrackStateByRP, 7)")
        df = df.Define("pfcand_phictgtheta","AlephSelection::get_constituent_trackCov(jetc, TrackStateByRP, 11)")
        df = df.Define("pfcand_dxyctgtheta","AlephSelection::get_constituent_trackCov(jetc, TrackStateByRP, 10)")
        df = df.Define("pfcand_dlambdadz",  "AlephSelection::get_constituent_trackCov(jetc, TrackStateByRP, 13)")
        df = df.Define("pfcand_cctgtheta",  "AlephSelection::get_constituent_trackCov(jetc, TrackStateByRP, 12)")
        df = df.Define("pfcand_phic",       "AlephSelection::get_constituent_trackCov(jetc, TrackStateByRP, 4)")
        df = df.Define("pfcand_dxyc",       "AlephSelection::get_constituent_trackCov(jetc, TrackStateByRP, 3)")
        df = df.Define("pfcand_cdz",        "AlephSelection::get_constituent_trackCov(jetc, TrackStateByRP, 8)")

        ############################################# Btag Variables #######################################################

        df = df.Define("pfcand_btagSip2dVal",   "JetConstituentsUtils::get_Sip2dVal_clusterV(jets, pfcand_dxy, pfcand_phi0, Bz)") 
        df = df.Define("pfcand_btagSip2dSig",   "JetConstituentsUtils::get_Sip2dSig(pfcand_btagSip2dVal, pfcand_dxydxy)") 
        df = df.Define("pfcand_btagSip3dVal",   "JetConstituentsUtils::get_Sip3dVal_clusterV(jets, pfcand_dxy, pfcand_dz, pfcand_phi0, Bz)") 
        df = df.Define("pfcand_btagSip3dSig",   "JetConstituentsUtils::get_Sip3dSig(pfcand_btagSip3dVal, pfcand_dxydxy, pfcand_dzdz)") 
        df = df.Define("pfcand_btagJetDistVal","AlephSelection::get_constituent_jetDistVal(jets, jetc, pfcand_dxy, pfcand_dz, pfcand_phi0, pfcand_ct)")
        df = df.Define("pfcand_btagJetDistSig","JetConstituentsUtils::get_JetDistSig(pfcand_btagJetDistVal, pfcand_dxydxy, pfcand_dzdz)")


        ############################################# Jet Level Variables and selection #######################################################
        
        df=df.Define("event_njet",   "JetConstituentsUtils::count_jets(jetc)")
        df = df.Filter("event_njet > 1")

        ##############################################################################################################
        df = df.Define("sumTLVs", "JetConstituentsUtils::sum_tlv_constituents(jetc)")

        df = df.Define("jet_p", "ROOT::VecOps::RVec<Double_t>({sumTLVs[0].P(), sumTLVs[1].P()})")
        df = df.Define("jet_e", "ROOT::VecOps::RVec<Double_t>({sumTLVs[0].E(), sumTLVs[1].E()})")
        df = df.Define("jet_mass", "ROOT::VecOps::RVec<Double_t>({sumTLVs[0].M(), sumTLVs[1].M()})")
        df = df.Define("jet_phi", "ROOT::VecOps::RVec<Double_t>({sumTLVs[0].Phi(), sumTLVs[1].Phi()})")
        df = df.Define("jet_theta", "ROOT::VecOps::RVec<Double_t>({sumTLVs[0].Theta(), sumTLVs[1].Theta()})")
        df = df.Define("jet_pT", "ROOT::VecOps::RVec<Double_t>({sumTLVs[0].Pt(), sumTLVs[1].Pt()})")
        df = df.Define("jet_eta", "ROOT::VecOps::RVec<Double_t>({sumTLVs[0].Eta(), sumTLVs[1].Eta()})")
        # Leading jet
        df = df.Define("jet_p_leading",      "sumTLVs[0].P()")
        df = df.Define("jet_e_leading",      "sumTLVs[0].E()")
        df = df.Define("jet_mass_leading",   "sumTLVs[0].M()")
        df = df.Define("jet_phi_leading",    "sumTLVs[0].Phi()")
        df = df.Define("jet_theta_leading",  "sumTLVs[0].Theta()")
        df = df.Define("jet_pT_leading",     "sumTLVs[0].Pt()")
        df = df.Define("jet_eta_leading",    "sumTLVs[0].Eta()")
        
        # Subleading jet
        df = df.Define("jet_p_subleading",      "sumTLVs[1].P()")
        df = df.Define("jet_e_subleading",      "sumTLVs[1].E()")
        df = df.Define("jet_mass_subleading",   "sumTLVs[1].M()")
        df = df.Define("jet_phi_subleading",    "sumTLVs[1].Phi()")
        df = df.Define("jet_theta_subleading",  "sumTLVs[1].Theta()")
        df = df.Define("jet_pT_subleading",     "sumTLVs[1].Pt()")
        df = df.Define("jet_eta_subleading",    "sumTLVs[1].Eta()")


        df = df.Define("jet_nconst", "JetConstituentsUtils::count_consts(jetc)") 
        ##
        df = df.Define(f"jet_nmu",    f"JetConstituentsUtils::count_type(pfcand_isMu)") 
        df = df.Define(f"jet_nel",    f"JetConstituentsUtils::count_type(pfcand_isEl)") 
        df = df.Define(f"jet_nchad",  f"JetConstituentsUtils::count_type(pfcand_isChargedHad)") 
        df = df.Define(f"jet_ngamma", f"JetConstituentsUtils::count_type(pfcand_isGamma)") 
        df = df.Define(f"jet_nnhad",  f"JetConstituentsUtils::count_type(pfcand_isNeutralHad)")

        # df = df.Define("dEdxPadsValue" , "dEdxPads.dQdx.value")
        # df = df.Define("dEdxPadsError" , "dEdxPads.dQdx.error")
        # df = df.Define("dEdxWiresValue" , "dEdxWires.dQdx.value")
        # df = df.Define("dEdxWiresError" , "dEdxPads.dQdx.error")

        # df = df.Define("jet_constituents_dEdx_pads_objs", "AlephSelection::build_constituents_dEdx()(RecoParticles, _RecoParticles_tracks.index, dEdxPads, _dEdxPads_track.index, _jetc)" )
        # df = df.Define("pfcand_dEdx_pads_type", "AlephSelection::get_dEdx_type(jet_constituents_dEdx_pads_objs)")
        # df = df.Define("pfcand_dEdx_pads_value", "AlephSelection::get_dEdx_value(jet_constituents_dEdx_pads_objs)")
        # df = df.Define("pfcand_dEdx_pads_error", "AlephSelection::get_dEdx_error(jet_constituents_dEdx_pads_objs)")

        # df = df.Define("jet_constituents_dEdx_wires_objs", "AlephSelection::build_constituents_dEdx()(RecoParticles, _RecoParticles_tracks.index, dEdxWires, _dEdxWires_track.index, _jetc)" )
        # df = df.Define("pfcand_dEdx_wires_type", "AlephSelection::get_dEdx_type(jet_constituents_dEdx_wires_objs)")
        # df = df.Define("pfcand_dEdx_wires_value", "AlephSelection::get_dEdx_value(jet_constituents_dEdx_wires_objs)")
        # df = df.Define("pfcand_dEdx_wires_error", "AlephSelection::get_dEdx_error(jet_constituents_dEdx_wires_objs)")

        # Get the dE/dx value and matching PID hypothesis pvalue from Bethe-Bloch fits for the jet constituents

        ## Pads
        df = df.Define("jet_constituents_dEdx_PIDhypo_pads_result", "AlephSelection::build_constituents_dEdx_PIDhypo()(RecoParticles, _RecoParticles_tracks.index, dEdxPads, _dEdxPads_track.index, _jetc, _Tracks_trackStates, false)" )
        df = df.Define("jet_constituents_dEdx_pads_objs", "jet_constituents_dEdx_PIDhypo_pads_result.dedx_constituents")
        df = df.Define("pfcand_dEdx_pads_type", "AlephSelection::get_dEdx_type(jet_constituents_dEdx_pads_objs)")
        df = df.Define("pfcand_dEdx_pads_value", "AlephSelection::get_dEdx_value(jet_constituents_dEdx_pads_objs)")
        df = df.Define("pfcand_dEdx_pads_error", "AlephSelection::get_dEdx_error(jet_constituents_dEdx_pads_objs)")

        ## extract the pvalues for different PID hyptheses based on Bethe-Bloch dE/dx vs p fits - order: ["e", "mu", "pi", "K", "p"]
        df = df.Define("jet_constituents_PID_pvals_pads", "jet_constituents_dEdx_PIDhypo_pads_result.pid_array_constituents")
        df = df.Define("pfcand_PID_pval_pads_ele", "AlephSelection::get_PID_pvalue(jet_constituents_PID_pvals_pads, 0)")
        df = df.Define("pfcand_PID_pval_pads_mu", "AlephSelection::get_PID_pvalue(jet_constituents_PID_pvals_pads, 1)")
        df = df.Define("pfcand_PID_pval_pads_pi", "AlephSelection::get_PID_pvalue(jet_constituents_PID_pvals_pads, 2)")
        df = df.Define("pfcand_PID_pval_pads_kaon", "AlephSelection::get_PID_pvalue(jet_constituents_PID_pvals_pads, 3)")
        df = df.Define("pfcand_PID_pval_pads_proton", "AlephSelection::get_PID_pvalue(jet_constituents_PID_pvals_pads, 4)")

        ## Wires
        df = df.Define("jet_constituents_dEdx_PIDhypo_wires_result", "AlephSelection::build_constituents_dEdx_PIDhypo()(RecoParticles, _RecoParticles_tracks.index, dEdxWires, _dEdxWires_track.index, _jetc, _Tracks_trackStates, true)" )
        df = df.Define("jet_constituents_dEdx_wires_objs", "jet_constituents_dEdx_PIDhypo_wires_result.dedx_constituents")
        df = df.Define("pfcand_dEdx_wires_type", "AlephSelection::get_dEdx_type(jet_constituents_dEdx_wires_objs)")
        df = df.Define("pfcand_dEdx_wires_value", "AlephSelection::get_dEdx_value(jet_constituents_dEdx_wires_objs)")
        df = df.Define("pfcand_dEdx_wires_error", "AlephSelection::get_dEdx_error(jet_constituents_dEdx_wires_objs)")

        ## extract the pvalues for different PID hyptheses based on Bethe-Bloch dE/dx vs p fits - order: ["e", "mu", "pi", "K", "p"]
        df = df.Define("jet_constituents_PID_pvals_wires", "jet_constituents_dEdx_PIDhypo_wires_result.pid_array_constituents")
        df = df.Define("pfcand_PID_pval_wires_ele", "AlephSelection::get_PID_pvalue(jet_constituents_PID_pvals_wires, 0)")
        df = df.Define("pfcand_PID_pval_wires_mu", "AlephSelection::get_PID_pvalue(jet_constituents_PID_pvals_wires, 1)")
        df = df.Define("pfcand_PID_pval_wires_pi", "AlephSelection::get_PID_pvalue(jet_constituents_PID_pvals_wires, 2)")
        df = df.Define("pfcand_PID_pval_wires_kaon", "AlephSelection::get_PID_pvalue(jet_constituents_PID_pvals_wires, 3)")
        df = df.Define("pfcand_PID_pval_wires_proton", "AlephSelection::get_PID_pvalue(jet_constituents_PID_pvals_wires, 4)")

        #for debug:
        df = df.Define("pfcand_dEdx_len", "pfcand_dEdx_wires_value[0].size()")
        df = df.Define("pfcand_pval_ele_len", "pfcand_PID_pval_wires_ele[0].size()")
        df = df.Define("pfcand_E_len", "pfcand_e[0].size()")



        ### Thrust variables
        df = df.Define("EVT_thrustNP",      'Algorithms::minimize_thrust("Minuit2","Migrad")(RP_px, RP_py, RP_pz)')
        df = df.Define("RP_thrustangleNP",  'Algorithms::getAxisCosTheta(EVT_thrustNP, RP_px, RP_py, RP_pz)')
        df = df.Define("EVT_thrust",        'Algorithms::getThrustPointing(1.)(RP_thrustangleNP, RP_e, EVT_thrustNP)')
        df = df.Define("EVT_Thrust_Mag",    "EVT_thrust.at(0)")  # thrust magnitude T (keep if you want it)
        df = df.Define("EVT_Thrust_X",      "EVT_thrust.at(1)")
        df = df.Define("EVT_Thrust_Y",      "EVT_thrust.at(3)")
        df = df.Define("EVT_Thrust_Z",      "EVT_thrust.at(5)")
        df = df.Define("EVT_Thrust_cosTheta", "EVT_Thrust_Z / sqrt(EVT_Thrust_X*EVT_Thrust_X + EVT_Thrust_Y*EVT_Thrust_Y + EVT_Thrust_Z*EVT_Thrust_Z)")
        

        return df

    def output(self):

        module_branches = []
        if self.do_truth:
            module_branches = [
                f"truev0_{b}" for b, _ in TRUEV0_DEFINES
            ] + [
                f"v0c_{b}" for b, _ in V0C_DEFINES
            ]
        if self.do_v0new:
            module_branches += ["n_v0n_event"] + [
                f"v0n_{b}" for b, _ in V0N_CAND_DEFINES
            ] + [
                f"v0n_{t}_origIdx" for t in V0N_TRKS
            ] + [
                f"v0n_{t}_{b}" for t in V0N_TRKS
                for b in DEDX_BRANCHES + LEG_PID_BRANCHES
            ] + [
                # per-jet new-module V0s (mirror of the old v0_* block) -> jet-level apples-to-apples
                "n_v0njet_jets", "n_v0njet_ks", "n_v0njet_lambda",
            ] + [
                f"v0njet_{b}" for b, _ in V0NJET_DEFINES
            ]
            if self.do_truth:
                module_branches += [n for n, _ in V0N_TRUTH_DEFINES]
        if self.do_svnew:
            for pfx in ("svn", "svm"):
                module_branches += [
                    f"n_{pfx}_event", f"{pfx}_mass", f"{pfx}_chi2", f"{pfx}_dxyz",
                    f"{pfx}_p", f"{pfx}_cosPointing", f"{pfx}_pointSig",
                    f"{pfx}_ntracks", f"{pfx}_sigL", f"{pfx}_trk_sv", f"{pfx}_trk_idx",
                    f"{pfx}_dx", f"{pfx}_dy", f"{pfx}_dz",
                    # SV vertex-fit covariance
                    f"{pfx}_cov_xx", f"{pfx}_cov_yx", f"{pfx}_cov_yy",
                    f"{pfx}_cov_zx", f"{pfx}_cov_zy", f"{pfx}_cov_zz",
                ]
            # V0 -> nearest-svn pointing feature
            module_branches += ["v0n_svnCosPoint", "v0n_svnPointSig", "v0n_svnIdx"]
            # (sec2origIdx lives in the always-written list — it is truth-free)
        if self.do_pvnew:
            # the flag surface of the standalone PV fitter
            module_branches += ["pv_converged", "pv_split_converged", "pv_trivial",
                               "pv_good"]
        if self.do_phikk:
            module_branches += ["n_phikk_event"] + [
                f"phikk_{b}" for b in PHIKK_CAND_BRANCHES
            ] + [
                f"phikk_{t}_{b}" for t in PHIKK_TRKS
                for b in PHIKK_TRK_BRANCHES + DEDX_BRANCHES + LEG_PID_BRANCHES
            ]
            if self.do_truth:
                module_branches += ["n_truephi_event", "truephi_found"] + [
                    f"truephi_{b}" for b in TRUEPHI_BRANCHES
                ] + ["phikk_class", "phikk_trueidx"] + [
                    f"phikk_{t}_{b}" for t in PHIKK_TRKS for b in TRK_TRUTH_BRANCHES
                ]
        if self.do_dstar:
            module_branches += ["n_d0_event", "n_dstar_event", "n_d0fits_event"] + [
                f"d0_{b}" for b in D0_CAND_BRANCHES
            ] + [
                f"dstar_{b}" for b in DSTAR_CAND_BRANCHES
            ] + [
                f"{pfx}_{b}" for pfx, _ in DSTAR_TRK_LEGS
                for b in DSTAR_TRK_BRANCHES + DEDX_BRANCHES + LEG_PID_BRANCHES
            ]
            if self.do_truth:
                module_branches += ["n_trued0_event", "n_truedstar_event",
                                   "trued0_found_loose", "trued0_found_tight",
                                   "truedstar_found_loose", "truedstar_found_tight",
                                   "d0_class", "d0_trueidx",
                                   "dstar_class", "dstar_trueidx"] + [
                    f"trued0_{b}" for b in TRUED0_BRANCHES
                ] + [
                    f"truedstar_{b}" for b in TRUEDSTAR_BRANCHES
                ] + [
                    f"d0_{t}_{b}" for t in D0_TRKS for b in TRK_TRUTH_BRANCHES
                ] + [
                    f"dstar_{t}_{b}" for t in DSTAR_TRKS for b in TRK_TRUTH_BRANCHES
                ]

        return module_branches + [
            #DEBUG
            "pfcand_dEdx_len", "pfcand_E_len", "pfcand_pval_ele_len",

            # Event variables
            "event_class",
            "event_number",
            "run_number",
            #"event_type",
            "event_invariant_mass",
            "event_njet",  
            "VertexX", 
            "VertexY", 
            "VertexZ",

            #refitted vertices
            "n_primary_tracks",
            "n_secondary_tracks",
            "Beamspot_x",
            "Beamspot_y",
            "Beamspot_z",
            "Vertex_refit_x",
            "Vertex_refit_y",
            "Vertex_refit_z",
            "Vertex_refit_cov_xx",
            "Vertex_refit_cov_yx",
            "Vertex_refit_cov_yy",
            "Vertex_refit_cov_zx",
            "Vertex_refit_cov_zy",
            "Vertex_refit_cov_zz",
            # PV fit quality + index maps + track<->pfcand join
            "Vertex_refit_chi2",
            "prim2origIdx",
            "sec2origIdx",
            "recopart_tracks_index",
            "recopart_tracks_begin",
            "recopart_tracks_end",
            "pfcand_trackIdx",
            # per-original-track membership bitmask + stored-candidate multiplicity
            "trk_member",
            "trk_nCand",

            # gen level vertex & resolutions
            "gen_vertex_x",
            "gen_vertex_y",
            "gen_vertex_z",

            # vertex resolution
            "res_vertex_x",
            "res_vertex_y",
            "res_vertex_z",
            
            # "res_vertex_x_all_tracks",
            # "res_vertex_y_all_tracks",
            # "res_vertex_z_all_tracks",

            # secondary vertices:
            "n_sv_event",
            "n_sv_jets",
            "sv_chi2",
            "sv_chi2_norm",
            "sv_ndof",
            "sv_ntracks",
            "sv_mass",
            "sv_p",
            "sv_thetarel",
            "sv_phirel",
            "sv_dxy",
            "sv_dxyz",
            "sv_cosPointing",
            "sv_prel",
            "sv_correctedMass",
            "sv_dx",
            "sv_dy",
            "sv_dz",
            # legacy-SV vertex-fit covariance
            "sv_cov_xx",
            "sv_cov_yx",
            "sv_cov_yy",
            "sv_cov_zx",
            "sv_cov_zy",
            "sv_cov_zz",

            # V0 candidates:
            "n_v0_event",
            "n_v0_jets",
            "n_v0_ks",
            "n_v0_lambda",
            "v0_pdg",
            "v0_invM",
            "v0_chi2",
            "v0_chi2_norm",
            "v0_ndof",
            "v0_ntracks",
            "v0_p",
            "v0_prel",
            "v0_thetarel",
            "v0_phirel",
            "v0_dxy",
            "v0_dxyz",
            "v0_cosPointing",
            "v0_correctedMass",
            "v0_dx",
            "v0_dy",
            "v0_dz",

            # Track variables
            "n_tracks_all",
            "n_tracks_sel",
            "n_trackstates_sel",
            "n_tracks_sel_vertexfit",
            "chi2_tracks_all",
            "ndf_tracks_all",
            "chi2_o_ndf_tracks_all",

            # Jet variables
            "JetClustering_d23",
            "JetClustering_d34", 
            "jet_mass",
            "jet_p",
            "jet_e", 
            "jet_phi", 
            "jet_theta", 
            "jet_pT",
            "jet_eta",
            "jet_p_leading",
            "jet_e_leading",
            "jet_mass_leading",
            "jet_phi_leading",
            "jet_theta_leading",
            "jet_pT_leading",
            "jet_eta_leading",
            "jet_p_subleading",
            "jet_e_subleading",
            "jet_mass_subleading",
            "jet_phi_subleading",
            "jet_theta_subleading",
            "jet_pT_subleading",
            "jet_eta_subleading", 
            "jet_nnhad",
            "jet_ngamma",
            "jet_nchad",
            "jet_nel", 
            "jet_nmu", 
            "jet_nconst",  

            "jetPID",

            # Pfcand/jet constituent variables
            "pfcand_isMu", 
            "pfcand_isEl", 
            "pfcand_isChargedHad", 
            "pfcand_isGamma", 
            "pfcand_isNeutralHad",
            "pfcand_e", 
            "pfcand_p", 
            "pfcand_px",
            "pfcand_py",
            "pfcand_pz",
            "pfcand_mask",
            "pfcand_theta", 
            "pfcand_phi", 
            "pfcand_charge", 
            "pfcand_erel",
            "pfcand_erel_log",
            "pfcand_thetarel",
            "pfcand_phirel",

            "pfcand_pt",
            "pfcand_ptrel",
            "pfcand_ptrel_log",
            "pfcand_trackChi2",
            "pfcand_trackNdof",
            "pfcand_trackChi2Norm",
            "pfcand_nTrackHits_VDET",
            "pfcand_nTrackHits_ITC",
            "pfcand_nTrackHits_TPC", 
            "pfcand_dxy", 
            "pfcand_dz", 
            "pfcand_phi0", 
            "pfcand_C", 
            "pfcand_ct",
            "pfcand_dptdpt", 
            "pfcand_dxydxy", 
            "pfcand_dzdz", 
            "pfcand_dphidphi", 
            "pfcand_detadeta",
            "pfcand_dxydz", 
            "pfcand_dphidxy", 
            "pfcand_phidz", 
            "pfcand_phictgtheta", 
            "pfcand_dxyctgtheta",
            "pfcand_dlambdadz", 
            "pfcand_cctgtheta", 
            "pfcand_phic", 
            "pfcand_dxyc", 
            "pfcand_cdz",
            "pfcand_btagSip2dVal", 
            "pfcand_btagSip2dSig",
            "pfcand_btagSip3dVal", 
            "pfcand_btagSip3dSig", 
            "pfcand_btagJetDistVal", 
            "pfcand_btagJetDistSig",

            # jet constituent PID 
            "pfcand_dEdx_pads_type", 
            "pfcand_dEdx_pads_value", 
            "pfcand_dEdx_pads_error",
            "pfcand_PID_pval_pads_ele",
            "pfcand_PID_pval_pads_mu",
            "pfcand_PID_pval_pads_pi",
            "pfcand_PID_pval_pads_kaon",
            "pfcand_PID_pval_pads_proton",

            "pfcand_dEdx_wires_type", 
            "pfcand_dEdx_wires_value", 
            "pfcand_dEdx_wires_error",
            "pfcand_PID_pval_wires_ele",
            "pfcand_PID_pval_wires_mu",
            "pfcand_PID_pval_wires_pi",
            "pfcand_PID_pval_wires_kaon",
            "pfcand_PID_pval_wires_proton",


            "EVT_Thrust_Mag",
            "EVT_Thrust_X",
            "EVT_Thrust_Y",
            "EVT_Thrust_Z",
            "EVT_Thrust_cosTheta",
            
            # to check if needed still? 
            # "dEdxPadsValue", "dEdxPadsError", "dEdxWiresValue", "dEdxWiresError",
            # #"Bz",

                
            ]
