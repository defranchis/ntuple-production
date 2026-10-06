
import os
import sys
from argparse import ArgumentParser
# fccanalysis loads this file by path and its batch workers do not inherit PYTHONPATH,
# so this directory's modules are made importable here.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_list

BZ = "FCCAnalyses::AlephUnits::kBz"  # solenoid field [T]

# Legacy V0 finder call options; its mass windows live in analyzer_trkaux.h.
V0_LEGACY_LOOSE_MASS_WINDOW = "true"
V0_LEGACY_DR_PAIR_CUT = "-1."   # dR preselection on track pairs (<= 0 disables)
V0_LEGACY_EXCLUSIVE_TRACKS = "true"  # skips used tracks, except that a booked pair's first track keeps pairing: tracks can be shared

# Per-daughter dE/dx, joined through <prefix>_origIdx and STORED for the
# calibration, never selected on: the collections to read, in written order.
DEDX_COLLS = (("pads", "dEdxPads"), ("wires", "dEdxWires"))
DEDX_LEG_DEFINES = tuple(
    (f"dEdx_{_d}_{_q}",
     "FCCAnalyses::AlephV0New::trackQuantityByIndex({pfx}_origIdx, "
     f"{_c}.dQdx.{_q}, dedxJoin_{_d})")
    for _d, _c in DEDX_COLLS for _q in ("value", "error"))
# tri-state particle-flow label (see legIsChargedHad) joined onto every
# candidate leg through its origIdx
LEG_PID_DEFINES = (
    ("isChargedHad",
     "FCCAnalyses::AlephTrkAux::legIsChargedHad({pfx}_origIdx, rpOfTrack, ParticleID)"),
)
PF_CHARGED_HAD = "FCCAnalyses::AlephTrkAux::kPFChargedHad"

# V0-module branch lists: (branch suffix, Define expression), one entry per branch.
V0N_CAND_DEFINES = (
    ("pdg",         "V0sNew_event.pdgAbs"),
    ("invM",        "V0sNew_event.invM"),
    ("alpha",       "FCCAnalyses::AlephV0New::candAlpha(V0sNew_event, SecondaryTracks_looseBS)"),
    ("qt",          "FCCAnalyses::AlephV0New::candQt(V0sNew_event)"),
    ("chi2",        "FCCAnalyses::AlephTruth::candChi2(V0sNew_event)"),
    ("dxyz",        "FCCAnalyses::AlephTruth::candDxyz(V0sNew_event, VertexObject_looseBS)"),
    # summed daughter momentum at the fitted vertex [GeV]
    ("px",          "FCCAnalyses::AlephTruth::candPcomp(V0sNew_event, 0)"),
    ("py",          "FCCAnalyses::AlephTruth::candPcomp(V0sNew_event, 1)"),
    ("pz",          "FCCAnalyses::AlephTruth::candPcomp(V0sNew_event, 2)"),
    ("cosPointing", "FCCAnalyses::AlephTruth::candCosPointing(V0sNew_event, VertexObject_looseBS)"),
    ("pointSig",    "FCCAnalyses::AlephV0New::candPointSig(V0sNew_event, VertexObject_looseBS)"),
    # 1 = tight tier, 0 = loose training tier
    ("tight",       "V0sNew_event.tight"),
    # ML-input pulls: cut variables in resolution units (signed; -999 undefined).
    ("bandSig",     "FCCAnalyses::AlephV0New::candBandSig(V0sNew_event, SecondaryTracks_looseBS)"),
    ("massSig",     "FCCAnalyses::AlephV0New::candMassSig(V0sNew_event)"),
    # fitted-vertex position [cm]
    ("vx",          "FCCAnalyses::AlephTruth::candVtxPos(V0sNew_event, 0)"),
    ("vy",          "FCCAnalyses::AlephTruth::candVtxPos(V0sNew_event, 1)"),
    ("vz",          "FCCAnalyses::AlephTruth::candVtxPos(V0sNew_event, 2)"),
    # vertex-fit covariance (packed lower triangle xx,yx,yy,zx,zy,zz, cm^2)
    ("cov_xx",      "FCCAnalyses::AlephV0New::candCovComp(V0sNew_event, 0)"),
    ("cov_yx",      "FCCAnalyses::AlephV0New::candCovComp(V0sNew_event, 1)"),
    ("cov_yy",      "FCCAnalyses::AlephV0New::candCovComp(V0sNew_event, 2)"),
    ("cov_zx",      "FCCAnalyses::AlephV0New::candCovComp(V0sNew_event, 3)"),
    ("cov_zy",      "FCCAnalyses::AlephV0New::candCovComp(V0sNew_event, 4)"),
    ("cov_zz",      "FCCAnalyses::AlephV0New::candCovComp(V0sNew_event, 5)"),
)
# daughter legs in momentum order: trk1 = higher momentum at the fitted vertex
V0N_TRKS = ("trk1", "trk2")
# Jet-relative tagger inputs, one entry per stored candidate (v0n_pdg order).
# Float sentinel -1; a candidate without a jet has v0n_jetIdx == -1.
V0N_TAG_DEFINES = (
    ("jetIdx",    "FCCAnalyses::AlephV0New::candJetIdx(V0sNew_event, jets)"),
    ("z",         "FCCAnalyses::AlephV0New::candJetVar(V0sNew_event, jets, v0n_jetIdx, 0)"),
    ("zL",        "FCCAnalyses::AlephV0New::candJetVar(V0sNew_event, jets, v0n_jetIdx, 1)"),
    ("ptRel",     "FCCAnalyses::AlephV0New::candJetVar(V0sNew_event, jets, v0n_jetIdx, 2)"),
    ("dRjet",     "FCCAnalyses::AlephV0New::candJetVar(V0sNew_event, jets, v0n_jetIdx, 3)"),
    ("rankInJet", "FCCAnalyses::AlephV0New::candRankInJet(V0sNew_event, v0n_jetIdx)"),
    # flight length from the PV; the 3D one is v0n_dxyz
    ("Lxy",       "FCCAnalyses::AlephV0New::candLxy(V0sNew_event, VertexObject_looseBS)"),
    ("LxySig",    "FCCAnalyses::AlephV0New::candFlightSig(V0sNew_event, VertexObject_looseBS, 0)"),
    ("LxyzSig",   "FCCAnalyses::AlephV0New::candFlightSig(V0sNew_event, VertexObject_looseBS, 1)"),
    ("baryon",    "FCCAnalyses::AlephV0New::candBaryon(V0sNew_event, SecondaryTracks_looseBS)"),
    ("nShared",   "FCCAnalyses::AlephV0New::candNShared(v0n_trk1_origIdx, v0n_trk2_origIdx)"),
)
# Per-daughter tagger inputs: {i} = leg index (0/1), {pfx} = the leg's prefix.
V0N_LEG_TAG_DEFINES = (
    ("q",    "FCCAnalyses::AlephV0New::candDaughterCharge(V0sNew_event, SecondaryTracks_looseBS, {i})"),
    ("p",    "FCCAnalyses::AlephV0New::candDaughterP(V0sNew_event, {i})"),
    ("nTPC", "FCCAnalyses::AlephTrkAux::subdetHits({pfx}_origIdx, Tracks.subdetectorHitNumbers_begin, Tracks.subdetectorHitNumbers_end, _Tracks_subdetectorHitNumbers, 2)"),
)

# SV-module branches (branch suffix, Define expression); per candidate, in event order:
SVN_CAND_DEFINES = (
    ("mass",        "SVs_svn.invM"),
    ("chi2",        "FCCAnalyses::AlephTruth::candChi2(SVs_svn)"),
    ("dxyz",        "FCCAnalyses::AlephTruth::candDxyz(SVs_svn, VertexObject_looseBS)"),
    ("dx",          "FCCAnalyses::AlephSVNew::candDcomp(SVs_svn, VertexObject_looseBS, 0)"),
    ("dy",          "FCCAnalyses::AlephSVNew::candDcomp(SVs_svn, VertexObject_looseBS, 1)"),
    ("dz",          "FCCAnalyses::AlephSVNew::candDcomp(SVs_svn, VertexObject_looseBS, 2)"),
    ("p",           "FCCAnalyses::AlephTruth::candP(SVs_svn)"),
    ("cosPointing", "FCCAnalyses::AlephTruth::candCosPointing(SVs_svn, VertexObject_looseBS)"),
    ("pointSig",    "FCCAnalyses::AlephV0New::candPointSig(SVs_svn, VertexObject_looseBS)"),
    ("ntracks",     "FCCAnalyses::AlephSVNew::candNtracks(SVs_svn)"),
    ("sigL",        "FCCAnalyses::AlephSVNew::candSigL(SVs_svn)"),
    ("cov_xx",      "FCCAnalyses::AlephV0New::candCovComp(SVs_svn, 0)"),
    ("cov_yx",      "FCCAnalyses::AlephV0New::candCovComp(SVs_svn, 1)"),
    ("cov_yy",      "FCCAnalyses::AlephV0New::candCovComp(SVs_svn, 2)"),
    ("cov_zx",      "FCCAnalyses::AlephV0New::candCovComp(SVs_svn, 3)"),
    ("cov_zy",      "FCCAnalyses::AlephV0New::candCovComp(SVs_svn, 4)"),
    ("cov_zz",      "FCCAnalyses::AlephV0New::candCovComp(SVs_svn, 5)"),
    # per constituent track, flat over the candidates:
    ("trk_sv",      "FCCAnalyses::AlephSVNew::candTrkSV(SVs_svn)"),
    ("trk_idx",     "FCCAnalyses::AlephSVNew::candTrkIdx(SVs_svn)"),
    ("trk_origIdx", "FCCAnalyses::AlephSVNew::candTrkOrigIdx(SVs_svn, sec2origIdx)"),
)
# V0 pointing at the svn vertex nearest in angle, per V0 candidate (v0n_pdg order):
V0N_SVN_DEFINES = (("svnCosPoint", "cosPoint"), ("svnPointSig", "pointSig"),
                   ("svnIdx", "svIdx"))

# phi->KK branch names: single source for the Define chain and the output list
PHIKK_CAND_BRANCHES = ("invM", "p", "px", "py", "pz", "alpha", "qt", "bandEll",
                       "chi2", "vx", "vy", "vz", "dpv", "dpvSig", "same_sign",
                       "wp", "tight")
PHIKK_TRKS = ("trk1", "trk2")
PHIKK_TRK_BRANCHES = ("origIdx", "q", "p", "costheta", "d0", "z0", "sigd0",
                      "nvdet", "nitc", "chi2ndf", "isprim")

# D* branch names; those read from the CandKin member (.kin): p/px/py/pz/costheta/xE are the D*'s, the rest the D0's
CAND_KIN_BRANCHES = ("m_kpi", "p", "px", "py", "pz", "costheta", "xE", "chi2",
                     "vx", "vy", "vz", "dpv", "dpvSig", "cosPoint",
                     "cosThetaStar")
DSTAR_CAND_BRANCHES = ("m_kpi", "dm", "p", "px", "py", "pz", "costheta", "xE",
                       "chi2", "vx", "vy", "vz", "dpv", "dpvSig", "cosPoint",
                       "cosThetaStar", "rs", "loose", "tight", "d0idx", "nsec")
DSTAR_TRK_BRANCHES = ("origIdx", "q", "p", "costheta", "d0", "z0", "sigd0",
                      "nvdet", "nitc", "chi2ndf", "isprim", "pool")
# (branch prefix, member of DstarCands) of every stored daughter leg
DSTAR_TRK_LEGS = (("dstar_trkK", "ds.trkK"), ("dstar_trkPi", "ds.trkPi"),
                  ("dstar_trkPis", "ds.trkPis"))
# the internal D0 legs: no output block, only their track indices are used
D0_TRK_LEGS = (("d0_trkK", "d0.trkK"), ("d0_trkPi", "d0.trkPi"))


def _cand_member(block, branch):
    """Member path of a D0/D* candidate branch inside DstarCands."""
    return f"{block}.kin.{branch}" if branch in CAND_KIN_BRANCHES \
        else f"{block}.{branch}"


PVNEW = "FCCAnalyses::AlephPVNew"  # namespace holding the PV selection constants

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
                            help='Fraction of events to run, default is 1.0 = 100%%')
        parser.add_argument('--batch', action='store_true', 
                            help='Submit to HTCondor batch')
        parser.add_argument('--valid', action='store_true', 
                            help='Run tester file only for validation against Lukas ntuples.')
        parser.add_argument('--chunks', default=None, type=int,
                            help='Number of chunks per process/file')
        parser.add_argument('--excludeRuns', nargs='+', action='extend', default=[], type=run_list.run_number, metavar='RUN',
                            help='data only: drop these run numbers in addition to the run list (eventsProcessed still counts the raw input).')
        parser.add_argument('--noRunList', action='store_true',
                            help='data only: keep every run instead of the data/lumi run list (--excludeRuns still applies).')
        parser.add_argument('--oldV0', action='store_true',
                            help='Legacy V0 only: drop the two-tier V0 module (no v0n_* branches).')
        parser.add_argument('--oldSV', action='store_true',
                            help='Drop the secondary-vertex module (no svn_* branches). The legacy sv_* block is unaffected.')
        parser.add_argument('--noV0TagVars', action='store_true',
                            help='Drop the jet-relative V0 tagger inputs (v0n_jetIdx/z/zL/ptRel/... and the per-leg q/p/nTPC). Implied by --oldV0.')
        parser.add_argument('--noPhiKK', action='store_true',
                            help='Skip the phi(1020)->K+K- finder (no phikk_* branches).')
        parser.add_argument('--noDstar', action='store_true',
                            help='Skip the D*+ -> D0(K pi) pi_slow finder (no dstar_* branches).')
        parser.add_argument('--noDedxGate', action='store_true',
                            help='accept every linked dE/dx measurement as valid, i.e. switch off the failed-leg omega sentinel gate; for converters that no longer copy omega into a failed leg.')
        parser.add_argument('--oldTrackSel', action='store_true',
                            help='baseline track selection without the minimum-TPC-hits and |z0| requirements, for the vertex fit and every finder (V0, secondary vertex, phi->KK, D*); also changes trk_member bit 9.')
        parser.add_argument('--oldPV', action='store_true',
                            help='Legacy PV chain: get_PrimaryTracks + VertexFitter_Tk and the origin-referenced track pre-selection, instead of the standalone fitter and its beamspot-referenced window (no pv_* flag branches).')
        # Parse additional arguments not known to the FCCAnalyses parsers
        # All command line arguments know to fccanalysis are provided in the
        # `cmdline_arg` dictionary.
        self.ana_args, unknown = parser.parse_known_args(cmdline_args['remaining'])
        if unknown:
            print(f"----> ERROR: unrecognised arguments: {' '.join(unknown)}")
            sys.exit(1)
        if not self.ana_args.doData and (self.ana_args.excludeRuns or self.ana_args.noRunList):
            print("----> ERROR: --excludeRuns and --noRunList apply to data only (--doData); Monte Carlo has no run list.")
            sys.exit(1)

        self.do_v0new = not self.ana_args.oldV0
        self.do_v0tagvars = self.do_v0new and not self.ana_args.noV0TagVars
        self.do_svnew = not self.ana_args.oldSV
        self.do_phikk = not self.ana_args.noPhiKK
        self.do_dstar = not self.ana_args.noDstar

        self.do_pvnew = not self.ana_args.oldPV

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
            sys.exit(1)
        
        if self.ana_args.MCtype and not self.ana_args.MCtype in outnames_dict:
            print("----> ERROR: Requested unknown --MCtype. Currently only zqq available.")
            sys.exit(1)
        
        if not self.ana_args.doData and not self.ana_args.MCflavour:
            print(f"----> ERROR: Requested MC run but did not specify --MCflavour. Please pick one..")
            sys.exit(1)
        
        if self.ana_args.MCflavour and not self.ana_args.MCflavour in outnames_dict[self.ana_args.MCtype]:
            print(f"----> ERROR: Requested unknown --MCflavour for --MCtype {self.ana_args.MCtype}. Check the dictionary.")
            sys.exit(1)

        #set the input/output directories:
        if self.ana_args.doData:
            self.input_dir = "/eos/experiment/fcc/ee/analyses/case-studies/aleph/LEP1_DATA/"
            self.output_dir_eos = f"/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedData/{self.ana_args.year}/stage1/{self.ana_args.tag}"
            self.output_dir = "."
            
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

                self.n_threads = 32 

        else:
            self.input_dir = f"/eos/experiment/aleph/EDM4HEP/MC/{self.ana_args.year}/"
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
                self.output_dir = f"/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedMC/{self.ana_args.year}/{self.ana_args.MCtype}/stage1/{self.ana_args.tag}"

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

            
                self.n_threads = 32 


        #set run options:
        
        # analyzer_truth.h and analyzer_trkaux.h are unconditional: the helpers the
        # V0 daughter branches join through (sec2origIdx index map, candidate
        # getters, vertex-fit glue) carry no truth and run on data too.
        self.include_paths = ["aleph_units.h", "aleph_reco_config.h", "analyzer.h", "analyzer_pvnew.h", "analyzer_truth.h", "analyzer_trkaux.h"]
        # the SV module reuses the V0 module's candidate getters
        if self.do_v0new or self.do_svnew:
            self.include_paths.append("analyzer_v0new.h")
        if self.do_svnew:
            self.include_paths.append("analyzer_svnew.h")
        if self.do_phikk:
            self.include_paths.append("analyzer_phikk.h")
        if self.do_dstar:
            self.include_paths.append("analyzer_dstar.h")

        # dE/dx validity gate, shared by the pfcand block and the candidate legs
        self.dedx_gate = "false" if self.ana_args.noDedxGate else "true"

        # #submit to batch if requested:
        # self.run_batch = self.ana_args.batch # no longer supported

    def _define_legs(self, df, legs, table):
        """One branch per daughter-leg prefix and (suffix, expression) table
        entry, named <prefix>_<suffix>; {pfx} = the prefix, {i} = the leg index."""
        for _i, _pfx in enumerate(legs):
            for _b, _e in table:
                df = df.Define(f"{_pfx}_{_b}", _e.format(i=_i, pfx=_pfx))
        return df

    @staticmethod
    def _pv_guard(expr, empty):
        """Empty-return entry guard on the usable-PV predicate: a finder runs only
        on a vertex that converged, whose kept tracks are all compatible with it,
        and that is track-supported (pv_good, goodPV() in analyzer_pvnew.h)."""
        return f"pv_good ? {expr} : {empty}"

    @staticmethod
    def _oldpv_guard(expr, empty):
        """The same under --oldPV: a PV of fewer than kPVMinTracks tracks is the
        default vertex at the origin (analyzer_trkaux.h)."""
        return f"VertexObject_looseBS.ntracks >= FCCAnalyses::AlephTrkAux::kPVMinTracks ? {expr} : {empty}"

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
            # Run selection from the data/lumi list minus --excludeRuns (--noRunList: every run).
            # The list is read where the graph is built; $ALEPH_RUN_LIST_<year> overrides its path.
            df = run_list.filter_runs(df, self.ana_args.year, self.ana_args.excludeRuns, self.ana_args.noRunList)
            #df = df.Filter("AlephSelection::sel_class_filter(16)(ClassBitset)   || AlephSelection::sel_class_filter(17)(ClassBitset) ")
            df = df.Filter("AlephSelection::sel_class_filter(16)(ClassBitset) ")
            df = df.Define("jetPID", "-999.f")
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

        # Analysis particles: RecoParticles without the photons and neutral hadrons below
        # kNeutralMinE; the jets, EVT_Evis and the thrust are built from them
        df = df.Define("RP_sel_idx", "AlephSelection::analysisParticleIndices(ParticleID, RP_e, AlephSelection::kNeutralMinE)")
        df = df.Define("RP_sel_px", "ROOT::VecOps::Take(RP_px, RP_sel_idx)")
        df = df.Define("RP_sel_py", "ROOT::VecOps::Take(RP_py, RP_sel_idx)")
        df = df.Define("RP_sel_pz", "ROOT::VecOps::Take(RP_pz, RP_sel_idx)")
        df = df.Define("RP_sel_e", "ROOT::VecOps::Take(RP_e, RP_sel_idx)")

        # Define pseudo-jets
        ####################################################################################################
        df = df.Define("pjetc", "JetClusteringUtils::set_pseudoJets(RP_sel_px, RP_sel_py, RP_sel_pz, RP_sel_e)")

        # Exclusive ee_kt (Durham) clustering to exactly 2 jets, E-ordered, E-scheme; jet constituents
        ####################################################################################################
        df = df.Define("_jet", "JetClustering::clustering_ee_kt(2, 2, 1, 0)(pjetc)")
        df = df.Define("jets","JetClusteringUtils::get_pseudoJets(_jet)" )
        # constituent indices into RecoParticles (and the index-parallel ParticleID)
        df = df.Define("_jetc", "AlephSelection::constituentsToRecoIndices(JetClusteringUtils::get_constituents(_jet), RP_sel_idx)")
        df = df.Define("jetc", "JetConstituentsUtils::build_constituents_cluster(RecoParticles, _jetc)")
        df = df.Define("jetConstitutentsTypes", f"AlephSelection::build_constituents_Types()(ParticleID, _jetc)")
        df = df.Define("JetClustering_d23", "std::sqrt(JetClusteringUtils::get_exclusive_dmerge(_jet, 2))")
        df = df.Define("JetClustering_d34", "std::sqrt(JetClusteringUtils::get_exclusive_dmerge(_jet, 3))")

        ############################################# Event Level Variables #######################################################
        df = df.Define("jet_p4", "JetConstituentsUtils::compute_tlv_jets(jets)" )
        df = df.Define("event_invariant_mass", "JetConstituentsUtils::InvariantMass(jet_p4[0], jet_p4[1])")

        ### Thrust variables
        # exact thrust {T, x, y, z} of the analysis particles, repacked into the {T, x, ex, y, ey, z, ez} layout of getAxisCosTheta/getThrustPointing
        df = df.Define("EVT_thrustExact",   "Algorithms::calculate_thrust()(RP_sel_px, RP_sel_py, RP_sel_pz)")
        df = df.Define("EVT_thrustNP",      "EVT_thrustExact[0] < 0 ? ROOT::VecOps::RVec<float>{-1.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f} : "
                                            "ROOT::VecOps::RVec<float>{EVT_thrustExact[0], EVT_thrustExact[1], 0.f, EVT_thrustExact[2], 0.f, EVT_thrustExact[3], 0.f}")
        df = df.Define("RP_sel_thrustangleNP", 'Algorithms::getAxisCosTheta(EVT_thrustNP, RP_sel_px, RP_sel_py, RP_sel_pz)')
        df = df.Define("EVT_thrust",        'Algorithms::getThrustPointing(1.)(RP_sel_thrustangleNP, RP_sel_e, EVT_thrustNP)')
        df = df.Define("EVT_Thrust_Mag",    "EVT_thrust.at(0)")
        df = df.Define("EVT_Thrust_X",      "EVT_thrust.at(1)")
        df = df.Define("EVT_Thrust_Y",      "EVT_thrust.at(3)")
        df = df.Define("EVT_Thrust_Z",      "EVT_thrust.at(5)")
        df = df.Define("EVT_Thrust_cosTheta", "EVT_Thrust_Mag < 0 ? -2.f : EVT_Thrust_Z / sqrt(EVT_Thrust_X*EVT_Thrust_X + EVT_Thrust_Y*EVT_Thrust_Y + EVT_Thrust_Z*EVT_Thrust_Z)")

        # per-run beamspot centre, 10 um units; override the json with $ALEPH_BEAMSPOT_JSON
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

        df = df.Define("Beamspot_x_cm", "Beamspot_x*1e-3")
        df = df.Define("Beamspot_y_cm", "Beamspot_y*1e-3")
        df = df.Define("Beamspot_z_cm", "Beamspot_z*1e-3")

        # ==== Track selection (to harmonize with Luka's code)
        # Note: The selection strategy here only works if there is one trackstate stored pre track.
        # The code includes an assertion for that, if it is somehow not the case it will fail. 
        # df = df.Define("n_tracks_all", f"AlephSelection::select_tracks( {coll['PFTracks']} )")
        df = df.Define("n_tracks_all", "Tracks.size()")
        df = df.Define("chi2_tracks_all","AlephSelection::get_track_chi2( Tracks )") #TODO: use collection here
        df = df.Define("ndf_tracks_all","AlephSelection::get_track_ndf( Tracks )") #TODO: use collection here
        df = df.Define("chi2_o_ndf_tracks_all","AlephSelection::get_track_chi2_o_ndf( Tracks )") #TODO: use collection here
        
        # baseline track selection
        if self.ana_args.oldTrackSel:
            min_tpc_hits, max_abs_z0 = "0", "std::numeric_limits<double>::infinity()"
        else:
            min_tpc_hits, max_abs_z0 = "AlephSelection::kTrackMinTPCHits", "AlephSelection::kTrackMaxAbsZ0"
        df = df.Define("tracks_selected_baseline_result",f"AlephSelection::select_tracks_baseline( Tracks, _Tracks_trackStates, _Tracks_subdetectorHitNumbers, {min_tpc_hits}, {max_abs_z0} )") #TODO: use collection here
        df = df.Define("tracks_selected_baseline","tracks_selected_baseline_result.tracks") 
        df = df.Define("trackstates_selected_baseline","tracks_selected_baseline_result.trackStates") 

        # impose upper bounds on impact parameters to pre-select compatible tracks for the primary vertex fit 
        if self.do_pvnew:
            df = df.Define("tracks_selected_for_vertexfit_result","AlephSelection::select_tracks_impactparameters_bs( tracks_selected_baseline_result, {0}::PVN_D0_MAX, {0}::PVN_Z0_MAX, Beamspot_x_cm, Beamspot_y_cm, Beamspot_z_cm )".format(PVNEW)) 
        else:
            df = df.Define("tracks_selected_for_vertexfit_result","AlephSelection::select_tracks_impactparameters( tracks_selected_baseline_result, FCCAnalyses::AlephReco::kPVTrackD0Max, FCCAnalyses::AlephReco::kPVTrackZ0Max )") 
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

        # Luka's loose BS constraints from looking at data, plus the chi2 below which
        # a track stays in the fit; values in aleph_reco_config.h, in the fitter's
        # unit of 10 um
        res_x_loose = "FCCAnalyses::AlephReco::kBeamSigmaXFit"
        res_y_loose = "FCCAnalyses::AlephReco::kBeamSigmaYFit"
        res_z_loose = "FCCAnalyses::AlephReco::kBeamSigmaZFit"
        chi2max = "FCCAnalyses::AlephReco::kPVChi2Max"

        if self.do_pvnew:
            bs_cm = "{}::beamSpot(Beamspot_x_cm, Beamspot_y_cm, Beamspot_z_cm)".format(PVNEW)
            df = df.Define("PVSelNew", "{}::select_primary_tracks(trackstates_selected_for_vertexfit_flipped, {})".format(PVNEW, bs_cm))
            df = df.Define("pv_converged",       "int(PVSelNew.fit.converged)")
            df = df.Define("pv_split_converged", "int(PVSelNew.split_converged)")
            df = df.Define("pv_trivial",         "int(PVSelNew.trivial)")
            df = df.Define("pv_good",            "int({}::goodPV(PVSelNew))".format(PVNEW))
            df = df.Define("RecoedPrimaryTracks_looseBS", "{}::primaryTracksFromSel(trackstates_selected_for_vertexfit_flipped, PVSelNew, Beamspot_x_cm, Beamspot_y_cm, Beamspot_z_cm)".format(PVNEW))
            df = df.Define("VertexObject_looseBS", "{}::toFCCVertex(PVSelNew)".format(PVNEW))
            df = df.Define("Vertex_refit_looseBS", "VertexObject_looseBS.vertex")
            df = df.Define("Vertex_refit_tlv", "pv_good ? TLorentzVector(Vertex_refit_looseBS.position.x, Vertex_refit_looseBS.position.y, Vertex_refit_looseBS.position.z, 0.) : TLorentzVector(Beamspot_x_cm, Beamspot_y_cm, Beamspot_z_cm, 0.)")
        else:
            # no primary tracks with fewer than 2 pre-selected tracks (get_PrimaryTracks would keep the one);
            # with fewer than 2 kept primary tracks, pruning included, VertexFitter_Tk returns the default vertex
            df = df.Define("RecoedPrimaryTracks_looseBS", "trackstates_selected_for_vertexfit_flipped.size() < 2 ? ROOT::VecOps::RVec<edm4hep::TrackState>{{}} : VertexFitterSimple::get_PrimaryTracks(trackstates_selected_for_vertexfit_flipped, true, {},{},{}, Beamspot_x, Beamspot_y, Beamspot_z, {})".format(res_x_loose, res_y_loose, res_z_loose, chi2max))
            df = df.Define("VertexObject_looseBS", "VertexFitterSimple::VertexFitter_Tk(1, RecoedPrimaryTracks_looseBS, true, {},{},{}, Beamspot_x, Beamspot_y, Beamspot_z)".format(res_x_loose, res_y_loose, res_z_loose))
            df = df.Define("Vertex_refit_looseBS", "VertexingUtils::get_VertexData(VertexObject_looseBS)")
            df = df.Define("Vertex_refit_tlv", "TLorentzVector(Vertex_refit_looseBS.position.x, Vertex_refit_looseBS.position.y, Vertex_refit_looseBS.position.z, 0.)")
        # for retrieving secondary tracks, use the full list of selected tracks 
        df = df.Define("SecondaryTracks_looseBS", "VertexFitterSimple::get_NonPrimaryTracks(trackstates_selected_baseline_flipped, RecoedPrimaryTracks_looseBS)")

        # original-Tracks index map of the secondary-track split, the join the V0
        # daughter branches reach the original tracks through; truth-free, so it
        # runs on data too
        df = df.Define("selBaselineOrigIdx", "tracks_selected_baseline_result.origIdx")
        df = df.Define("sec2origIdx",        "FCCAnalyses::AlephTruth::secondaryToOriginalTrack(SecondaryTracks_looseBS, trackstates_selected_baseline_flipped, selBaselineOrigIdx)")
        df = df.Define("prim2origIdx",       "FCCAnalyses::AlephTruth::secondaryToOriginalTrack(RecoedPrimaryTracks_looseBS, trackstates_selected_baseline_flipped, selBaselineOrigIdx)")

        df = df.Define("Vertex_refit_x", "Vertex_refit_looseBS.position.x")
        df = df.Define("Vertex_refit_y", "Vertex_refit_looseBS.position.y")
        df = df.Define("Vertex_refit_z", "Vertex_refit_looseBS.position.z")

        df = df.Define("Vertex_refit_cov_xx", "Vertex_refit_looseBS.covMatrix.values[0]")
        df = df.Define("Vertex_refit_cov_yx", "Vertex_refit_looseBS.covMatrix.values[1]")
        df = df.Define("Vertex_refit_cov_yy", "Vertex_refit_looseBS.covMatrix.values[2]")
        df = df.Define("Vertex_refit_cov_zx", "Vertex_refit_looseBS.covMatrix.values[3]")
        df = df.Define("Vertex_refit_cov_zy", "Vertex_refit_looseBS.covMatrix.values[4]")
        df = df.Define("Vertex_refit_cov_zz", "Vertex_refit_looseBS.covMatrix.values[5]")

        df = df.Define("Vertex_refit_chi2", "Vertex_refit_looseBS.chi2")

        df = df.Define("n_primary_tracks", "ReconstructedParticle2Track::getTK_n(RecoedPrimaryTracks_looseBS)")
        df = df.Define("n_secondary_tracks", "ReconstructedParticle2Track::getTK_n(SecondaryTracks_looseBS)")

        # for comparison test, fit vertex with tracks all tracks:
        # df = df.Define("RecoedPrimaryTracks_looseBS_all_tracks", "VertexFitterSimple::get_PrimaryTracks(_Tracks_trackStates, true, {},{},{},0.,0.,0., {})".format(res_x_loose, res_y_loose, res_z_loose, chi2max))
        # df = df.Define("VertexObject_looseBS_all_tracks", "VertexFitterSimple::VertexFitter_Tk(1, RecoedPrimaryTracks_looseBS_all_tracks, true, {},{},{},0.,0.,0.)".format(res_x_loose, res_y_loose, res_z_loose))
        # df = df.Define("Vertex_refit_looseBS_all_tracks", "VertexingUtils::get_VertexData(VertexObject_looseBS_all_tracks)")
        # df = df.Define("Vertex_refit_tlv_all_tracks", "TLorentzVector(Vertex_refit_looseBS_all_tracks.position.x, Vertex_refit_looseBS_all_tracks.position.y, Vertex_refit_looseBS_all_tracks.position.z, 0.)")

        # df = df.Define("Vertex_refit_x_all_tracks", "Vertex_refit_looseBS_all_tracks.position.x")
        # df = df.Define("Vertex_refit_y_all_tracks", "Vertex_refit_looseBS_all_tracks.position.y")
        # df = df.Define("Vertex_refit_z_all_tracks", "Vertex_refit_looseBS_all_tracks.position.z")

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
            df = df.Define("gen_vertex_x", "-999.")
            df = df.Define("gen_vertex_y", "-999.")
            df = df.Define("gen_vertex_z", "-999.")

            # refit vertex resolution:
            df = df.Define("res_vertex_x", "-999.")
            df = df.Define("res_vertex_y", "-999.")
            df = df.Define("res_vertex_z", "-999.")
        
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
        sv_expr = ("FCCAnalyses::AlephSelection::get_SV_event_ALEPH("
            "SecondaryTracks_looseBS, "               # non-primary tracks
            "trackstates_selected_baseline_flipped, " # all tracks
            "VertexObject_looseBS, "                  # primary vertex
            "0.8, "                                   # dR prefilter cut
            "false)"                                  # exclusive V0 rejection (skip+break), matching FCCAnalyses@3a4de97 isV0 - the code that produced ntuples-withks
        )
        if self.do_pvnew:
            sv_expr = self._pv_guard(
                sv_expr,
                "ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex>{}")
        df = df.Define("SVs_looseBS", sv_expr)

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
        # vertex-fit covariance [cm^2]
        for ic, cc in enumerate(("xx", "yx", "yy", "zx", "zy", "zz")):
            df = df.Define(f"sv_cov_{cc}", f"FCCAnalyses::AlephSelection::svCovComp(sv_jets, {ic})")
        df = df.Define("sv_trk_sv",      "FCCAnalyses::AlephSelection::svTrkSV(sv_jets)")
        df = df.Define("sv_trk_origIdx", "FCCAnalyses::AlephSelection::svTrkOrigIdx(sv_jets, selBaselineOrigIdx)")

        ############################################# V0 Reconstruction #######################################################
        v0_expr = ("FCCAnalyses::AlephSelection::get_V0s_ALEPH("
            "SecondaryTracks_looseBS, "
            "VertexObject_looseBS,"
            f"{BZ},"
            f"{V0_LEGACY_LOOSE_MASS_WINDOW},"
            f"{V0_LEGACY_DR_PAIR_CUT},"
            f"{V0_LEGACY_EXCLUSIVE_TRACKS})"
        )
        if self.do_pvnew:
            v0_expr = self._pv_guard(
                v0_expr, "FCCAnalyses::VertexingUtils::FCCAnalysesV0{}")
        df = df.Define("V0s_event", v0_expr)
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

        # joins feeding the candidate-leg branches below, built once per event:
        # track -> ReconstructedParticle for the per-leg PF label, and track ->
        # dE/dx measurement index, where the shared validity gate is applied so
        # that a failed leg reads -1 in both the value and the error branch
        if self.do_v0new or self.do_phikk or self.do_dstar:
            df = df.Define("rpOfTrack",
                           "FCCAnalyses::AlephTrkAux::rpIndexByTrack(RecoParticles.tracks_begin, RecoParticles.tracks_end, _RecoParticles_tracks.index, Tracks.size())")
            for _det, _coll in DEDX_COLLS:
                df = df.Define(f"dedxJoin_{_det}",
                               f"FCCAnalyses::AlephV0New::dedxIndexByTrack({_coll}.dQdx.value, {_coll}.dQdx.error, _{_coll}_track.index, Tracks, _Tracks_trackStates, {self.dedx_gate})")

        ############################################# Standalone two-tier V0 module ###########################################
        if self.do_v0new:
            v0n_expr = f"FCCAnalyses::AlephV0New::findV0s(SecondaryTracks_looseBS, VertexObject_looseBS, {BZ})"
            if self.do_pvnew:
                v0n_expr = self._pv_guard(v0n_expr, "FCCAnalyses::AlephV0New::V0Collection{}")
            df = df.Define("V0sNew_event", v0n_expr)
            df = df.Define("n_v0n_event",  "int(V0sNew_event.vtx.size())")
            for _b, _e in V0N_CAND_DEFINES:
                df = df.Define(f"v0n_{_b}", _e)
            # per-daughter joins: candidate reco_ind -> sec2origIdx -> Tracks index
            for _i, _t in enumerate(V0N_TRKS):
                df = df.Define(f"v0n_{_t}_origIdx",
                               f"FCCAnalyses::AlephV0New::candDaughterOrigIdx(V0sNew_event, sec2origIdx, {_i})")
            _legs = [f"v0n_{_t}" for _t in V0N_TRKS]
            df = self._define_legs(df, _legs, DEDX_LEG_DEFINES)
            df = self._define_legs(df, _legs, LEG_PID_DEFINES)
            if self.do_v0tagvars:
                for _b, _e in V0N_TAG_DEFINES:
                    df = df.Define(f"v0n_{_b}", _e)
                df = self._define_legs(df, _legs, V0N_LEG_TAG_DEFINES)

        ############################################# Secondary-vertex module #################################################
        if self.do_svnew:
            SVNEW = "FCCAnalyses::AlephSVNew"
            # tight Ks/Lambda daughters are masked; under --oldV0 nothing is
            if self.do_v0new:
                v0_mask = f"V0sNew_event, v0n_tight, {SVNEW}::SVN_MASK_MODE"
            else:
                v0_mask = f"FCCAnalyses::VertexingUtils::FCCAnalysesV0{{}}, ROOT::VecOps::RVec<int>{{}}, {SVNEW}::SVN_MASK_NONE"
            seed_expr = f"{SVNEW}::svSeedPass(SecondaryTracks_looseBS, VertexObject_looseBS, {BZ})"
            svn_expr = (f"{SVNEW}::findSVs(SecondaryTracks_looseBS, VertexObject_looseBS, {v0_mask}, {BZ}, SVSeeds_event, "
                        "TVector3(EVT_Thrust_X, EVT_Thrust_Y, EVT_Thrust_Z))")
            if self.do_pvnew:
                seed_expr = self._pv_guard(seed_expr, f"{SVNEW}::SVSeeds{{}}")
                svn_expr = self._pv_guard(svn_expr, "FCCAnalyses::VertexingUtils::FCCAnalysesV0{}")
            else:
                seed_expr = self._oldpv_guard(seed_expr, f"{SVNEW}::SVSeeds{{}}")
                svn_expr = self._oldpv_guard(svn_expr, "FCCAnalyses::VertexingUtils::FCCAnalysesV0{}")
            df = df.Define("SVSeeds_event", seed_expr)
            df = df.Define("SVs_svn", svn_expr)
            df = df.Define("n_svn_event", "int(SVs_svn.vtx.size())")
            for _b, _e in SVN_CAND_DEFINES:
                df = df.Define(f"svn_{_b}", _e)
            if self.do_v0new:
                df = df.Define("v0n_svnpoint", "FCCAnalyses::AlephV0New::candSVPointing(V0sNew_event, SVs_svn, sec2origIdx)")
                for _b, _m in V0N_SVN_DEFINES:
                    df = df.Define(f"v0n_{_b}", f"v0n_svnpoint.{_m}")

        ############################################# exclusive-finder track auxiliaries ######################################
        if self.do_phikk or self.do_dstar:
            # both finders run on the full baseline-selected track list, joined to the original Tracks by selBaselineOrigIdx
            TRKAUX = "FCCAnalyses::AlephTrkAux"
            df = df.Define("trkaux_nvdet", f"{TRKAUX}::subdetHits(selBaselineOrigIdx, Tracks.subdetectorHitNumbers_begin, Tracks.subdetectorHitNumbers_end, _Tracks_subdetectorHitNumbers, 0)")
            df = df.Define("trkaux_nitc",  f"{TRKAUX}::subdetHits(selBaselineOrigIdx, Tracks.subdetectorHitNumbers_begin, Tracks.subdetectorHitNumbers_end, _Tracks_subdetectorHitNumbers, 1)")
            df = df.Define("trkaux_chi2ndf", f"{TRKAUX}::trackChi2Ndf(selBaselineOrigIdx, Tracks.chi2, Tracks.ndf)")
            df = df.Define("trkaux_isprim",  f"{TRKAUX}::flagInSet(selBaselineOrigIdx, prim2origIdx)")
            # tight Ks/Lambda daughters leave both pools; empty under --oldV0
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
                          "Beamspot_x_cm, Beamspot_y_cm, Beamspot_z_cm)")
            if self.do_pvnew:
                phikk_expr = self._pv_guard(phikk_expr, "FCCAnalyses::AlephPhiKK::PhiKKCands{}")
            else:
                phikk_expr = self._oldpv_guard(phikk_expr, "FCCAnalyses::AlephPhiKK::PhiKKCands{}")
            df = df.Define("PhiKKCands_event", phikk_expr)
            df = df.Define("n_phikk_event", "int(PhiKKCands_event.invM.size())")
            for _b in PHIKK_CAND_BRANCHES:
                df = df.Define(f"phikk_{_b}", f"PhiKKCands_event.{_b}")
            for _t in PHIKK_TRKS:
                for _b in PHIKK_TRK_BRANCHES:
                    df = df.Define(f"phikk_{_t}_{_b}", f"PhiKKCands_event.{_t}.{_b}")
            _legs = [f"phikk_{_t}" for _t in PHIKK_TRKS]
            df = self._define_legs(df, _legs, DEDX_LEG_DEFINES)
            df = self._define_legs(df, _legs, LEG_PID_DEFINES)

        ############################################# D*->D0(K pi) pi_slow module #############################################
        if self.do_dstar:
            # primary/secondary class of every pool track (0 prim / 1 sec)
            df = df.Define("dstar_pool_all", "FCCAnalyses::AlephDstar::poolClass(selBaselineOrigIdx, prim2origIdx, sec2origIdx)")
            # every selection value is a constant in analyzer_dstar.h
            dstar_expr = ("FCCAnalyses::AlephDstar::findDstar(trackstates_selected_baseline_flipped, "
                          "selBaselineOrigIdx, trkaux_nvdet, trkaux_nitc, trkaux_chi2ndf, "
                          "trkaux_isprim, dstar_pool_all, VertexObject_looseBS, v0n_claimed_orig, "
                          f"{BZ}, Beamspot_x_cm, Beamspot_y_cm, Beamspot_z_cm)")
            if self.do_pvnew:
                dstar_expr = self._pv_guard(dstar_expr, "FCCAnalyses::AlephDstar::DstarCands{}")
            else:
                dstar_expr = self._oldpv_guard(dstar_expr, "FCCAnalyses::AlephDstar::DstarCands{}")
            df = df.Define("DstarCands_event", dstar_expr)
            df = df.Define("n_dstar_event", "int(DstarCands_event.ds.kin.m_kpi.size())")
            # two-track fits actually performed: the combinatorial cost
            df = df.Define("n_d0fits_event", "DstarCands_event.nfits")
            for _b in DSTAR_CAND_BRANCHES:
                df = df.Define(f"dstar_{_b}", f"DstarCands_event.{_cand_member('ds', _b)}")
            for _pfx, _mem in DSTAR_TRK_LEGS:
                for _b in DSTAR_TRK_BRANCHES:
                    df = df.Define(f"{_pfx}_{_b}", f"DstarCands_event.{_mem}.{_b}")
            _legs = [_pfx for _pfx, _ in DSTAR_TRK_LEGS]
            df = self._define_legs(df, _legs, DEDX_LEG_DEFINES)
            df = self._define_legs(df, _legs, LEG_PID_DEFINES)
            # D0 leg indices for the membership pass; the D0 list is not written
            for _pfx, _mem in D0_TRK_LEGS:
                df = df.Define(f"{_pfx}_origIdx", f"DstarCands_event.{_mem}.origIdx")

        ############################################# per-track membership ####################################################
        # one pass over the finished candidate lists, per original track
        _EMPTY = "ROOT::VecOps::RVec<int>{}"
        # SV constituents come back in the baseline-selected frame the finder was given
        df = df.Define("svTrkIdx", "FCCAnalyses::AlephTrkAux::svTrackIdx(SVs_looseBS)")
        _v0 = ("v0n_trk1_origIdx, v0n_trk2_origIdx, v0n_tight" if self.do_v0new
               else f"{_EMPTY}, {_EMPTY}, {_EMPTY}")
        _phi = ("phikk_trk1_origIdx, phikk_trk2_origIdx, phikk_wp" if self.do_phikk
                else f"{_EMPTY}, {_EMPTY}, {_EMPTY}")
        _ds = ("d0_trkK_origIdx, d0_trkPi_origIdx, dstar_trkK_origIdx, "
               "dstar_trkPi_origIdx, dstar_trkPis_origIdx, dstar_tight"
               if self.do_dstar else ", ".join([_EMPTY] * 6))
        _svn = "svn_trk_origIdx" if self.do_svnew else _EMPTY
        df = df.Define("trkTags",
                       "FCCAnalyses::AlephTrkAux::trackTags(Tracks.size(), "
                       "selBaselineOrigIdx, prim2origIdx, svTrkIdx, selBaselineOrigIdx, "
                       f"{_svn}, {_v0}, {_phi}, {_ds})")
        df = df.Define("trk_member", "trkTags.member")
        df = df.Define("trk_nCand",  "trkTags.nCand")

        ############################################# Particle Flow Level Variables #######################################################
        df = df.Define("pfcand_isMu",     "AlephSelection::get_isType(jetConstitutentsTypes,2)")
        df = df.Define("pfcand_isEl",     "AlephSelection::get_isType(jetConstitutentsTypes,1)")
        df = df.Define("pfcand_isGamma",  "AlephSelection::get_isType(jetConstitutentsTypes,AlephSelection::kPFPhoton)")
        df = df.Define("pfcand_isChargedHad", f"AlephSelection::get_isType(jetConstitutentsTypes,{PF_CHARGED_HAD})")
        df = df.Define("pfcand_isNeutralHad", "AlephSelection::get_isType(jetConstitutentsTypes,AlephSelection::kPFNeutralHad)")


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

        # tracks re-indexed through the RecoParticle->Track relation (see analyzer.h)
        df = df.Define("TracksByRP", "AlephSelection::reindexByRPLink(Tracks, _RecoParticles_tracks.index)")

        # track fit quality per constituent (-1 for neutrals, which have no track)
        df = df.Define("pfcand_trackChi2",     "AlephSelection::get_constituent_trackChi2(jetc, TracksByRP)")
        df = df.Define("pfcand_trackNdof",     "AlephSelection::get_constituent_trackNdof(jetc, TracksByRP)")
        df = df.Define("pfcand_trackChi2Norm", "AlephSelection::get_constituent_trackChi2Norm(jetc, TracksByRP)")

        # original-Tracks index of each constituent's track (-1 = none), the join key to the finders' *_origIdx branches
        df = df.Define("pfcand_trackIdx", "AlephSelection::get_constituent_trackIdx(jetc, _RecoParticles_tracks.index)")

        # subdetector hit counts per constituent (inside-out: VDET, ITC, TPC)
        # the offsets inside TrackData index the flat hit-number array, so it is passed as is
        df = df.Define("pfcand_nTrackHits_VDET", "AlephSelection::get_constituent_nTrackHits_VDET(jetc, TracksByRP, _Tracks_subdetectorHitNumbers)")
        df = df.Define("pfcand_nTrackHits_ITC",  "AlephSelection::get_constituent_nTrackHits_ITC(jetc, TracksByRP, _Tracks_subdetectorHitNumbers)")
        df = df.Define("pfcand_nTrackHits_TPC",  "AlephSelection::get_constituent_nTrackHits_TPC(jetc, TracksByRP, _Tracks_subdetectorHitNumbers)")

        df = df.Define("Bz", f'{BZ}') # luka reads this from the event ? 

        ############################################# Track Parameters and Covariance #######################################################

        df = df.Define("TrackStateFlipped",f"AlephSelection::flipD0_copy( {coll['TrackState']} )")
        df = df.Define("TrackStateByRP", "AlephSelection::trackStatesByRPLink(TrackStateFlipped, Tracks, _RecoParticles_tracks.index)")

        # dxy/dz/phi0 at the PV; C and ct are the fitted curvature and dip angle
        df = df.Define("pfcand_trkparPV",   "AlephSelection::get_constituent_trackParamsAtPV(jetc, TrackStateByRP, Vertex_refit_tlv, Bz)")
        df = df.Define("pfcand_dxy",        "pfcand_trkparPV.dxy")
        df = df.Define("pfcand_dz",         "pfcand_trkparPV.dz")
        df = df.Define("pfcand_phi0",       "pfcand_trkparPV.phi0")
        df = df.Define("pfcand_C",          "pfcand_trkparPV.C")
        df = df.Define("pfcand_ct",         "pfcand_trkparPV.ct")
        # raw perigee d0/z0 of the constituent's track: origin-referenced, cm, -9 for neutrals
        df = df.Define("pfcand_d0",         "AlephSelection::get_constituent_D0(jetc, TrackStateByRP)")
        df = df.Define("pfcand_z0",         "AlephSelection::get_constituent_Z0(jetc, TrackStateByRP)")
        # covariance, lower-triangular in (d0, phi0, omega, z0, tanLambda): cov(a,b) at a*(a+1)/2 + b
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
        df = df.Define("jet_constituents_dEdx_PIDhypo_pads_result", f"AlephSelection::build_constituents_dEdx_PIDhypo()(RecoParticles, _RecoParticles_tracks.index, dEdxPads, _dEdxPads_track.index, _jetc, Tracks, _Tracks_trackStates, false, {self.dedx_gate})" )
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
        df = df.Define("jet_constituents_dEdx_PIDhypo_wires_result", f"AlephSelection::build_constituents_dEdx_PIDhypo()(RecoParticles, _RecoParticles_tracks.index, dEdxWires, _dEdxWires_track.index, _jetc, Tracks, _Tracks_trackStates, true, {self.dedx_gate})" )
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



        df = df.Define("EVT_Evis",          "Sum(RP_sel_e)")  # total visible energy: sum over the analysis particles [GeV]
        df = df.Define("jetPartonPDG", f"AlephSelection::getJetPartonPDG({coll['GenParticles']}, jets, 0.8)")
        df = df.Define("jetPartonCosTheta", f"AlephSelection::getJetPartonCosTheta({coll['GenParticles']}, jets)")
        

        return df

    def output(self):

        module_branches = []
        if self.do_v0new:
            module_branches += ["n_v0n_event"] + [
                f"v0n_{b}" for b, _ in V0N_CAND_DEFINES
            ] + [
                f"v0n_{t}_origIdx" for t in V0N_TRKS
            ] + [
                f"v0n_{t}_{b}" for t in V0N_TRKS
                for b, _ in DEDX_LEG_DEFINES + LEG_PID_DEFINES
            ]
            if self.do_v0tagvars:
                module_branches += [
                    f"v0n_{b}" for b, _ in V0N_TAG_DEFINES
                ] + [
                    f"v0n_{t}_{b}" for t in V0N_TRKS
                    for b, _ in V0N_LEG_TAG_DEFINES
                ]
        if self.do_svnew:
            module_branches += ["n_svn_event"] + [
                f"svn_{b}" for b, _ in SVN_CAND_DEFINES
            ]
            if self.do_v0new:
                module_branches += [f"v0n_{b}" for b, _ in V0N_SVN_DEFINES]
        if self.do_phikk:
            module_branches += ["n_phikk_event"] + [
                f"phikk_{b}" for b in PHIKK_CAND_BRANCHES
            ] + [
                f"phikk_{t}_{b}" for t in PHIKK_TRKS
                for b in PHIKK_TRK_BRANCHES
                + tuple(b for b, _ in DEDX_LEG_DEFINES + LEG_PID_DEFINES)
            ]
        if self.do_dstar:
            module_branches += ["n_dstar_event", "n_d0fits_event"] + [
                f"dstar_{b}" for b in DSTAR_CAND_BRANCHES
            ] + [
                f"{pfx}_{b}" for pfx, _ in DSTAR_TRK_LEGS
                for b in DSTAR_TRK_BRANCHES
                + tuple(b for b, _ in DEDX_LEG_DEFINES + LEG_PID_DEFINES)
            ]

        pv_branches = ["pv_converged", "pv_split_converged", "pv_trivial",
                       "pv_good"] if self.do_pvnew else []

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
            "Beamspot_x_cm",
            "Beamspot_y_cm",
            "Beamspot_z_cm",
            "Vertex_refit_x",
            "Vertex_refit_y",
            "Vertex_refit_z",
            "Vertex_refit_cov_xx",
            "Vertex_refit_cov_yx",
            "Vertex_refit_cov_yy",
            "Vertex_refit_cov_zx",
            "Vertex_refit_cov_zy",
            "Vertex_refit_cov_zz",
            "Vertex_refit_chi2",
            *pv_branches,

            # track <-> pfcand join key
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
            "sv_cov_xx",
            "sv_cov_yx",
            "sv_cov_yy",
            "sv_cov_zx",
            "sv_cov_zy",
            "sv_cov_zz",
            "sv_trk_sv",
            "sv_trk_origIdx",

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
            "pfcand_d0",
            "pfcand_z0",
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
            "EVT_Evis",

            "jetPartonPDG",
            "jetPartonCosTheta",
            # to check if needed still? 
            # "dEdxPadsValue", "dEdxPadsError", "dEdxWiresValue", "dEdxWiresError",
            # #"Bz",

                
            ]
