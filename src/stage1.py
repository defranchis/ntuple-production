
import os
from argparse import ArgumentParser, BooleanOptionalAction

BZ = 1.5  # solenoid field [T] — single source for the stage1 Define strings
PVNEW = "FCCAnalyses::AlephPVNew"  # namespace holding the PV selection constants

# --dstar branch-name lists: single source for the Define chain and the
# output branch list (per-candidate quantities, and the per-daughter block
# instantiated for the K, the pi and the slow pi).
D0_CAND_BRANCHES = ("m_kpi", "p", "px", "py", "pz", "costheta", "xE", "chi2",
                    "vx", "vy", "vz", "dpv", "dpvSig", "cosPoint",
                    "cosThetaStar", "loose", "tight", "stage", "nsec")
DSTAR_CAND_BRANCHES = ("m_kpi", "dm", "p", "px", "py", "pz", "costheta", "xE",
                       "chi2", "vx", "vy", "vz", "dpv", "dpvSig", "cosPoint",
                       "cosThetaStar", "rs", "loose", "tight", "d0idx",
                       "stage", "nsec")
DSTAR_TRK_BRANCHES = ("origIdx", "q", "p", "costheta", "d0", "z0", "sigd0",
                      "nvdet", "nitc", "chi2ndf", "isprim", "pool")
TRUED0_BRANCHES = ("p", "pt", "costheta", "xE", "pK", "cosK", "pPi", "origin",
                   "fromDstar", "mothPdg", "nmatched", "flight",
                   "K_pool", "pi_pool")
TRUEDSTAR_BRANCHES = ("p", "pt", "costheta", "xE", "px", "py", "pz", "pK",
                      "cosK", "pPi", "pPis", "cosPis", "origin", "mothPdg",
                      "nmatched", "d0flight", "K_pool", "pi_pool", "pis_pool")

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
        parser.add_argument('--v0nWideLamLoose', action='store_true',
                            help='TAIL-MEASUREMENT VARIANT: loose Lambda AP band ramp edges doubled (0.40/0.80; stored acceptance 0.8x that) to measure the band tail. Not for standard productions.')
        parser.add_argument('--v0nLamPointKsTiers', action='store_true',
                            help='SIZING VARIANT: tight-tier Lambda pointing aligned to the Ks p-tiers. Not for standard productions; candTight still encodes the adopted package.')
        parser.add_argument('--phiKK', action='store_true',
                            help='Run the standalone phi(1020)->K+K- finder (phikk_* branches); '
                                 'opt-in extension of the standalone V0 machinery. Pairs are formed '
                                 'from the FULL selected track list; truth branches added on MC.')
        parser.add_argument('--phiKKmLo', default=None, type=float,
                            help='--phiKK: low edge of the STORED K+K- mass window [GeV]; default AlephPhiKK::M_LO.')
        parser.add_argument('--phiKKmHi', default=None, type=float,
                            help='--phiKK: high edge of the STORED K+K- mass window [GeV]; default AlephPhiKK::M_HI.')
        parser.add_argument('--phiKKchi2', default=None, type=float,
                            help='--phiKK: vertex-fit chi2 cut (ndf=1); loose sanity value, <=0 disables; '
                                 'default AlephPhiKK::CHI2_CUT.')
        parser.add_argument('--phiKKapBand', default=-1., type=float,
                            help='--phiKK: |bandEll-1| Armenteros band cut; <=0 = off (default). '
                                 'For equal-mass daughters the band is a reparametrisation of the mass window.')
        parser.add_argument('--phiKKdpv', default=-1., type=float,
                            help='--phiKK: PV-compatibility cut |vtx-PV| [cm]; <=0 = off (default). '
                                 'Promptness is NOT required: phi from b/c decays are displaced.')
        parser.add_argument('--phiKKdpvFid', default=None, type=float,
                            help='--phiKK: |vtx-PV| storage fiducial [cm]; <=0 = off; default AlephPhiKK::DPV_FID. '
                                 'Sanity bound only: near-collinear pairs give a vertex '
                                 'unconstrained along the flight direction, so a real phi '
                                 'can reach ~17 cm.')
        parser.add_argument('--phiKKdpvSig', default=-1., type=float,
                            help='--phiKK: PV-compatibility significance cut; <=0 = off (default).')
        parser.add_argument('--phiKKsigd0', default=-1., type=float,
                            help='--phiKK: per-track sigma(d0) cap [cm]; <=0 = off (default).')
        parser.add_argument('--phiKKminHits', default=0, type=int,
                            help='--phiKK: per-track minimum (nVDET + nITC) hits; 0 = off (default).')
        parser.add_argument('--phiKKtrkChi2', default=-1., type=float,
                            help='--phiKK: per-track fit chi2/ndf cap; <=0 = off (default).')
        parser.add_argument('--phiKKpMin', default=None, type=float,
                            help='--phiKK: per-track momentum floor [GeV]; <=0 = off; default AlephPhiKK::P_MIN_DEF. '
                                 'Applied to the perigee momentum, while the stored trk p is the '
                                 'at-vertex one (~1%% tail difference).')
        parser.add_argument('--phiKKnoSameSign', action='store_true',
                            help='--phiKK: do NOT reconstruct same-charge pairs (they are the '
                                 'data-driven combinatorial control and are stored by default).')
        parser.add_argument('--phiKKvetoV0', action='store_true',
                            help='--phiKK: exclude tracks already claimed by a tight Ks/Lambda '
                                 'candidate from the pairing (needs the V0 module, i.e. not --oldV0). '
                                 'Off by default.')
        parser.add_argument('--dstar', action='store_true',
                            help='Run the standalone D*+ -> D0(K pi) pi_slow finder (dstar_* branches '
                                 'plus the stand-alone D0 list d0_*); opt-in extension of the standalone '
                                 'V0 machinery. Candidates are built from the FULL selected track list; '
                                 'truth branches added on MC.')
        parser.add_argument('--dstarMLo', default=1.70, type=float,
                            help='--dstar: low edge of the STORED K pi mass window [GeV] (provisional).')
        parser.add_argument('--dstarMHi', default=2.03, type=float,
                            help='--dstar: high edge of the STORED K pi mass window [GeV] (provisional).')
        parser.add_argument('--dstarChi2', default=25., type=float,
                            help='--dstar: D0 vertex-fit chi2 cut (ndf=1); loose sanity value, <=0 disables (provisional).')
        parser.add_argument('--dstarDpvMax', default=10., type=float,
                            help='--dstar: |vtx-PV| storage fiducial [cm], default 10; <=0 = off. Sanity '
                                 'bound only - promptness is NOT required (D* from b decays are displaced), '
                                 'dpv/dpvSig are stored and never cut on (provisional).')
        parser.add_argument('--dstarDmMax', default=0.20, type=float,
                            help='--dstar: ceiling on the stored dm = m(K pi pi_s) - m(K pi) [GeV] (provisional).')
        parser.add_argument('--dstarPMin', default=0.3, type=float,
                            help='--dstar: momentum floor [GeV] for the K and pi candidates; <=0 = off. '
                                 'Applied to the perigee momentum, while the stored trk p is the at-vertex one (provisional).')
        parser.add_argument('--dstarPsMin', default=0.1, type=float,
                            help='--dstar: momentum floor [GeV] for the slow pion; <=0 = off (provisional).')
        parser.add_argument('--dstarSigd0', default=-1., type=float,
                            help='--dstar: per-track sigma(d0) cap [cm]; <=0 = off (default).')
        parser.add_argument('--dstarMinHits', default=0, type=int,
                            help='--dstar: per-track minimum (nVDET + nITC) hits; 0 = off (default).')
        parser.add_argument('--dstarTrkChi2', default=-1., type=float,
                            help='--dstar: per-track fit chi2/ndf cap; <=0 = off (default).')
        parser.add_argument('--d0LooseDm', default=0.060, type=float,
                            help='--dstar: d0_loose label, |m(K pi) - m_D0| [GeV] (provisional).')
        parser.add_argument('--d0TightDm', default=0.030, type=float,
                            help='--dstar: d0_tight label, |m(K pi) - m_D0| [GeV] (provisional).')
        parser.add_argument('--d0TightDpvSig', default=3.0, type=float,
                            help='--dstar: d0_tight label, minimum 3D |vtx-PV| significance (provisional).')
        parser.add_argument('--d0TightCosPoint', default=0.99, type=float,
                            help='--dstar: d0_tight label, minimum cos(angle) between p(D0) and (vtx-PV) (provisional).')
        parser.add_argument('--d0TightCosStar', default=0.8, type=float,
                            help='--dstar: d0_tight label, maximum |cos(theta*)| of the kaon in the D0 frame (provisional).')
        parser.add_argument('--dstarLooseDm', default=0.050, type=float,
                            help='--dstar: dstar_loose label, |m(K pi) - m_D0| [GeV] (provisional).')
        parser.add_argument('--dstarLooseDdm', default=0.0030, type=float,
                            help='--dstar: dstar_loose label, |dm - 0.145426| [GeV] (provisional).')
        parser.add_argument('--dstarTightDm', default=0.025, type=float,
                            help='--dstar: dstar_tight label, |m(K pi) - m_D0| [GeV] (provisional).')
        parser.add_argument('--dstarTightDdm', default=0.0015, type=float,
                            help='--dstar: dstar_tight label, |dm - 0.145426| [GeV] (provisional).')
        parser.add_argument('--dstarTightPK', default=1.0, type=float,
                            help='--dstar: tight labels, minimum kaon momentum [GeV] (provisional).')
        parser.add_argument('--dstarTightPPi', default=1.0, type=float,
                            help='--dstar: tight labels, minimum pion momentum [GeV] (provisional).')
        parser.add_argument('--dstarTightChi2', default=10., type=float,
                            help='--dstar: tight labels, maximum D0 vertex chi2 (provisional).')
        parser.add_argument('--dstarTightPs', default=0.3, type=float,
                            help='--dstar: dstar_tight label, minimum slow-pion momentum [GeV]; '
                                 '<=0 = off. Not applied to dstar_loose (provisional).')
        parser.add_argument('--dstarTightCosPoint', default=0.95, type=float,
                            help='--dstar: dstar_tight label, minimum cosPoint of the D0 vertex; '
                                 '<=-1 = off. Not applied to dstar_loose (provisional).')
        parser.add_argument('--dstarCascade', action='store_true',
                            help='--dstar: STAGED EXCLUSIVE mode. Instead of one all-track pass, '
                                 'candidates are built in six ordered stages by the primary/secondary '
                                 'pool pattern of the (K, pi, pi_s) legs, most displaced first; after '
                                 'each stage the legs of the claimed candidates leave the pool, so the '
                                 'later stages see far fewer combinations. The D0-alone list gets the '
                                 'same treatment in three stages. Off by default (single pass).')
        parser.add_argument('--dstarClaim', default='tight', choices=['tight', 'loose', 'none'],
                            help='--dstarCascade: which candidates claim their tracks between stages '
                                 '(right-sign only for the D* cascade); none = stage without claiming.')
        parser.add_argument('--dstarVetoV0', action=BooleanOptionalAction, default=True,
                            help='--dstar: exclude tracks already claimed by a tight Ks/Lambda '
                                 'candidate from the D* and D0 pools (needs the V0 module, i.e. not '
                                 '--oldV0). On by default; --no-dstarVetoV0 turns it off.')
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
            import os as _os
            _r = _os.environ.get("ALEPH_RECLUS_DIR")
            if _r and _os.path.isdir(_r):
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
                _o = _os.environ.get("ALEPH_OUT_DIR")
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
            import os
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
                #local tester for validation
                if self.ana_args.valid:

                    # local tester: node-local disk (/eos/experiment is read-only
                    # from the batch nodes), single input file
                    self.output_dir = f"/tmp/aleph_valid_runs/{self.ana_args.tag}"

                    # re-clustered input copy (30 TTree clusters instead of 3):
                    # lifts the RDataFrame thread cap from 3 to ~30 cores
                    import os
                    if os.path.isdir("/tmp/reclus_input/QQB"):
                        self.input_dir = "/tmp/reclus_input/"

                    self.process_list = {
                        "QQB/ZM4212_40_AL" : {"fraction" : self.ana_args.fraction, "output":"ntuple_valid_tester_{}".format(self.ana_args.MCflavour)},
                    }
                
                #process full files:
                else:
                    self.process_list = {
                            "QQB" : {"fraction" : self.ana_args.fraction, "output":output_name},
                        }

                    # ALEPH_OUT_DIR redirects non-batch output to a writable area
                    _o = os.environ.get("ALEPH_OUT_DIR")
                    self.output_dir = (f"{_o}/wp1_stage1/{self.ana_args.tag}" if _o else
                        f"/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedMC/{self.ana_args.year}/{self.ana_args.MCtype}/stage1/{self.ana_args.tag}")

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
        self.include_paths = ["analyzer.h", "analyzer_truth.h", "analyzer_pvnew.h"]
        if self.do_v0new:
            self.include_paths.append("analyzer_v0new.h")
        if self.do_svnew:
            self.include_paths.append("analyzer_svnew.h")
        if self.ana_args.phiKK:
            if self.ana_args.phiKKvetoV0 and not self.do_v0new:
                print("----> ERROR: --phiKKvetoV0 needs the V0 module (the veto list comes from "
                      "the tight Ks/Lambda claims); it is incompatible with --oldV0.")
                exit()
            # analyzer_v0new.h carries trackQuantityByIndex, the shared dE/dx
            # lookup; appended only when it is not already in the list
            if "analyzer_v0new.h" not in self.include_paths:
                self.include_paths.append("analyzer_v0new.h")
            self.include_paths.append("analyzer_phikk.h")
        if self.ana_args.dstar:
            if self.ana_args.dstarVetoV0 and not self.do_v0new:
                print("----> ERROR: --dstarVetoV0 (on by default) needs the V0 module (the veto list "
                      "comes from the tight Ks/Lambda claims); use --no-dstarVetoV0 to disable it.")
                exit()
            # analyzer_v0new.h carries invMass and trackQuantityByIndex, the
            # shared helpers; appended only when it is not already in the list
            if "analyzer_v0new.h" not in self.include_paths:
                self.include_paths.append("analyzer_v0new.h")
            self.include_paths.append("analyzer_dstar.h")

        # #submit to batch if requested:
        # self.run_batch = self.ana_args.batch # no longer supported

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
        # a reference copy lives in data/. Override with $ALEPH_BEAMSPOT_JSON.
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
        # the legacy fitter wants x,y in 10um and z in mm: cm x 1e3 for all three
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
            # single named source in analyzer_pvnew.h. Stored only; the existing
            # consumer ternaries keep their pv_converged-only policy.
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
            df = df.Define("Vertex_refit_tlv", "pv_converged ? TLorentzVector(Vertex_refit_looseBS.position.x, Vertex_refit_looseBS.position.y, Vertex_refit_looseBS.position.z, 0.) : TLorentzVector(Beamspot_x*1e-3, Beamspot_y*1e-3, Beamspot_z*1e-3, 0.)")
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

        # for comparison test, fit vertex with tracks all tracks:
        # df = df.Define("RecoedPrimaryTracks_looseBS_all_tracks", "VertexFitterSimple::get_PrimaryTracks(_Tracks_trackStates, true, {},0.,0.,0., {})".format(bs_sig_legacy, chi2max))
        # df = df.Define("VertexObject_looseBS_all_tracks", "VertexFitterSimple::VertexFitter_Tk(1, RecoedPrimaryTracks_looseBS_all_tracks, true, {},0.,0.,0.)".format(bs_sig_legacy))
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
            old_sv_expr = ("pv_converged ? " + old_sv_expr +
                           " : ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex>{}")
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
            "1.5," #solenoidBz
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
            df = df.Define("truev0_pdg",      "trueV0s.pdg")
            df = df.Define("truev0_p",        "trueV0s.p")
            df = df.Define("truev0_costheta", "trueV0s.costheta")
            df = df.Define("truev0_px",       "trueV0s.px")
            df = df.Define("truev0_py",       "trueV0s.py")
            df = df.Define("truev0_pz",       "trueV0s.pz")
            df = df.Define("truev0_fd",       "trueV0s.fd")
            df = df.Define("truev0_dpv",      "trueV0s.dpv")
            df = df.Define("truev0_nmatched", "trueV0s.nmatched")
            # true decay-point components [cm] (position-resolution studies)
            df = df.Define("truev0_x",        "trueV0s.vx")
            df = df.Define("truev0_y",        "trueV0s.vy")
            df = df.Define("truev0_z",        "trueV0s.vz")
            # (selBaselineOrigIdx / sec2origIdx are defined unconditionally in the
            # PV block above — truth-free track-state matching, available on data)
            # daughters surviving into the secondary-track set (0-2): separates PV-claim losses from finder losses
            df = df.Define("truev0_nsec",        "FCCAnalyses::AlephTruth::daughtersInSecondaries(trueV0s, mcToTracks, sec2origIdx)")
            # recover which track pair each candidate came from (compiled get_V0s leaves reco_ind empty);
            # classifyV0s cross-checks this replica against V0s_event pdg/invM and throws on mismatch
            df = df.Define("v0pairs",       "FCCAnalyses::AlephTruth::rerunV0Pairing(SecondaryTracks_looseBS, VertexObject_looseBS)")
            # truth classification of the reco V0 candidates (event order = V0s_event order)
            df = df.Define("v0truth",       f"FCCAnalyses::AlephTruth::classifyV0s(V0s_event, v0pairs, SecondaryTracks_looseBS, sec2origIdx, trackToMCs, {coll['GenParticles']}, trueV0s)")
            df = df.Define("v0c_class",       "v0truth.cls")
            df = df.Define("v0c_trueidx",     "v0truth.true_idx")
            df = df.Define("v0c_pairmult",    "v0truth.pair_mult")
            df = df.Define("v0c_trackshared", "v0truth.track_shared")
            df = df.Define("v0c_alpha",       "v0truth.alpha")
            df = df.Define("v0c_qt",          "v0truth.qt")
            df = df.Define("v0c_trk1",        "v0truth.trk1")
            df = df.Define("v0c_trk2",        "v0truth.trk2")
            # event-order candidate kinematics (independent of jet assignment)
            df = df.Define("v0c_pdg",         "V0s_event.pdgAbs")
            df = df.Define("v0c_invM",        "V0s_event.invM")
            df = df.Define("v0c_dxyz",        "FCCAnalyses::AlephTruth::candDxyz(V0s_event, VertexObject_looseBS)")
            df = df.Define("v0c_p",           "FCCAnalyses::AlephTruth::candP(V0s_event)")
            df = df.Define("v0c_cosPointing", "FCCAnalyses::AlephTruth::candCosPointing(V0s_event, VertexObject_looseBS)")
            # fitted-vertex position (position-resolution studies)
            df = df.Define("v0c_vx",          "FCCAnalyses::AlephTruth::candVtxPos(V0s_event, 0)")
            df = df.Define("v0c_vy",          "FCCAnalyses::AlephTruth::candVtxPos(V0s_event, 1)")
            df = df.Define("v0c_vz",          "FCCAnalyses::AlephTruth::candVtxPos(V0s_event, 2)")
            # per-true-V0 found flags
            df = df.Define("truev0_found_any",     "FCCAnalyses::AlephTruth::trueV0FoundAny(trueV0s, v0truth)")
            df = df.Define("truev0_found_correct", "FCCAnalyses::AlephTruth::trueV0FoundCorrect(trueV0s, v0truth, V0s_event)")

        ############################################# Standalone two-tier V0 module ###########################################
        if self.do_v0new:
            if self.ana_args.v0nLamPointKsTiers and self.ana_args.v0nWideLamLoose:
                print("----> ERROR: --v0nLamPointKsTiers and --v0nWideLamLoose are separate "
                      "single-purpose variants; combining them is not supported.")
                exit()
            v0n_finder = ("findV0sLamKsPointing" if self.ana_args.v0nLamPointKsTiers
                          else "findV0sWideLamLoose" if self.ana_args.v0nWideLamLoose
                          else "findV0s")
            v0n_expr = f"FCCAnalyses::AlephV0New::{v0n_finder}(SecondaryTracks_looseBS, VertexObject_looseBS, {BZ})"
            if self.do_pvnew:
                # explicit empty-return entry guard on the flag (the window
                # cuts would empty it anyway, a silent efficiency loss; the
                # guard makes the failure explicit and empties pointSig too)
                v0n_expr = ("pv_converged ? " + v0n_expr +
                            " : FCCAnalyses::VertexingUtils::FCCAnalysesV0{}")
            df = df.Define("V0sNew_event", v0n_expr)
            df = df.Define("n_v0n_event",  "int(V0sNew_event.vtx.size())")
            df = df.Define("v0n_pdg",      "V0sNew_event.pdgAbs")
            df = df.Define("v0n_invM",     "V0sNew_event.invM")
            # truth-free kinematic branches (available on data)
            df = df.Define("v0n_alpha",       "FCCAnalyses::AlephV0New::candAlpha(V0sNew_event, SecondaryTracks_looseBS)")
            df = df.Define("v0n_qt",          "FCCAnalyses::AlephV0New::candQt(V0sNew_event)")
            df = df.Define("v0n_chi2",        "FCCAnalyses::AlephTruth::candChi2(V0sNew_event)")
            df = df.Define("v0n_dxyz",        "FCCAnalyses::AlephTruth::candDxyz(V0sNew_event, VertexObject_looseBS)")
            df = df.Define("v0n_p",           "FCCAnalyses::AlephTruth::candP(V0sNew_event)")
            # momentum VECTOR of the same summed vertex momentum as v0n_p
            # (direction-dependent offline studies: pointing at any reference)
            for ic, cc in enumerate("xyz"):
                df = df.Define(f"v0n_p{cc}",  f"FCCAnalyses::AlephTruth::candPcomp(V0sNew_event, {ic})")
            df = df.Define("v0n_cosPointing", "FCCAnalyses::AlephTruth::candCosPointing(V0sNew_event, VertexObject_looseBS)")
            df = df.Define("v0n_pointSig",    "FCCAnalyses::AlephV0New::candPointSig(V0sNew_event, VertexObject_looseBS)")
            # two-tier module: 1 = adopted tight package, 0 = loose training tier.
            # Selecting v0n_tight==1 reproduces the historical tight-only output exactly.
            df = df.Define("v0n_tight",       "FCCAnalyses::AlephV0New::candTight(V0sNew_event, VertexObject_looseBS, SecondaryTracks_looseBS)")
            # ML-input pulls: cut variables in resolution units (signed; -999 undefined).
            df = df.Define("v0n_bandSig",     "FCCAnalyses::AlephV0New::candBandSig(V0sNew_event, SecondaryTracks_looseBS)")
            df = df.Define("v0n_massSig",     "FCCAnalyses::AlephV0New::candMassSig(V0sNew_event)")
            # V0 vertex-fit covariance (packed lower triangle, cm^2 — same
            # component order as Vertex_refit_cov_*)
            for ic, cc in enumerate(("xx", "yx", "yy", "zx", "zy", "zz")):
                df = df.Define(f"v0n_cov_{cc}", f"FCCAnalyses::AlephV0New::candCovComp(V0sNew_event, {ic})")
            # per-daughter joins + dE/dx (truth-free: reco_ind -> sec2origIdx).
            # dE/dx validity: value != omega(track) (failed-leg sentinel),
            # finite positive value and error; invalid -> -1 in both branches.
            df = df.Define("v0n_trk1_origIdx", "FCCAnalyses::AlephV0New::candDaughterOrigIdx(V0sNew_event, sec2origIdx, 0)")
            df = df.Define("v0n_trk2_origIdx", "FCCAnalyses::AlephV0New::candDaughterOrigIdx(V0sNew_event, sec2origIdx, 1)")
            for trk in ("trk1", "trk2"):
                for det, dedx_coll in (("pads", "dEdxPads"), ("wires", "dEdxWires")):
                    df = df.Define(f"v0n_{trk}_dEdx_{det}_value",
                                   f"FCCAnalyses::AlephV0New::trackQuantityByIndex(v0n_{trk}_origIdx, {dedx_coll}.dQdx.value, {dedx_coll}.dQdx.value, {dedx_coll}.dQdx.error, _{dedx_coll}_track.index, _Tracks_trackStates)")
                    df = df.Define(f"v0n_{trk}_dEdx_{det}_error",
                                   f"FCCAnalyses::AlephV0New::trackQuantityByIndex(v0n_{trk}_origIdx, {dedx_coll}.dQdx.error, {dedx_coll}.dQdx.value, {dedx_coll}.dQdx.error, _{dedx_coll}_track.index, _Tracks_trackStates)")
            # fitted-vertex position (position-resolution studies)
            df = df.Define("v0n_vx",          "FCCAnalyses::AlephTruth::candVtxPos(V0sNew_event, 0)")
            df = df.Define("v0n_vy",          "FCCAnalyses::AlephTruth::candVtxPos(V0sNew_event, 1)")
            df = df.Define("v0n_vz",          "FCCAnalyses::AlephTruth::candVtxPos(V0sNew_event, 2)")
            # per-jet new-module V0s: mirror of the old-finder v0_* block on
            # V0sNew_event, so v0_* vs v0njet_* is an apples-to-apples comparison
            # at jet level (prel = pT wrt jet axis, thetarel/phirel wrt the jet).
            df = df.Define("v0njet_per_jet", "FCCAnalyses::AlephSelection::assign_V0s_to_jets(V0sNew_event, jets)")
            df = df.Define("v0njet_jets",  "v0njet_per_jet.vtx")
            df = df.Define("v0njet_pdg",   "v0njet_per_jet.pdgAbs")
            df = df.Define("v0njet_invM",  "v0njet_per_jet.invM")
            df = df.Define("n_v0njet_jets",    "FCCAnalyses::VertexingUtils::get_n_SV_jets(v0njet_jets)")
            df = df.Define("n_v0njet_ks",      "FCCAnalyses::AlephSelection::count_V0type_jets(v0njet_pdg, 310)")
            df = df.Define("n_v0njet_lambda",  "FCCAnalyses::AlephSelection::count_V0type_jets(v0njet_pdg, 3122)")
            df = df.Define("v0njet_chi2",          "FCCAnalyses::VertexingUtils::get_chi2_SV(v0njet_jets)")
            df = df.Define("v0njet_chi2_norm",     "FCCAnalyses::VertexingUtils::get_norm_chi2_SV(v0njet_jets)")
            df = df.Define("v0njet_ndof",          "FCCAnalyses::VertexingUtils::get_nDOF_SV(v0njet_jets)")
            df = df.Define("v0njet_ntracks",       "FCCAnalyses::VertexingUtils::get_VertexNtrk(v0njet_jets)")
            df = df.Define("v0njet_p",             "FCCAnalyses::VertexingUtils::get_pMag_SV(v0njet_jets)")
            df = df.Define("v0njet_prel",          "FCCAnalyses::AlephSelection::get_prel_SV_jets(v0njet_jets, jets)")
            df = df.Define("v0njet_thetarel",      "FCCAnalyses::VertexingUtils::get_relTheta_SV(v0njet_jets, jets)")
            df = df.Define("v0njet_phirel",        "FCCAnalyses::VertexingUtils::get_relPhi_SV(v0njet_jets, jets)")
            df = df.Define("v0njet_dxy",           "FCCAnalyses::VertexingUtils::get_dxy_SV(v0njet_jets, VertexObject_looseBS)")
            df = df.Define("v0njet_dxyz",          "FCCAnalyses::VertexingUtils::get_d3d_SV(v0njet_jets, VertexObject_looseBS)")
            df = df.Define("v0njet_cosPointing",   "FCCAnalyses::AlephSelection::get_pointingangle_SV(v0njet_jets, VertexObject_looseBS)")
            df = df.Define("v0njet_correctedMass", "FCCAnalyses::AlephSelection::get_correctedInvMass_SV(v0njet_jets, VertexObject_looseBS)")
            df = df.Define("v0njet_dx",  "FCCAnalyses::AlephSelection::get_dx_SV_jets(v0njet_jets, PrimaryVertexP3)")
            df = df.Define("v0njet_dy",  "FCCAnalyses::AlephSelection::get_dy_SV_jets(v0njet_jets, PrimaryVertexP3)")
            df = df.Define("v0njet_dz",  "FCCAnalyses::AlephSelection::get_dz_SV_jets(v0njet_jets, PrimaryVertexP3)")
            # truth classification (MC only; reco_ind is filled by the new module)
            if self.do_truth:
                df = df.Define("v0npairs",     "FCCAnalyses::AlephTruth::pairsFromRecoInd(V0sNew_event)")
                df = df.Define("v0ntruth",     f"FCCAnalyses::AlephTruth::classifyV0s(V0sNew_event, v0npairs, SecondaryTracks_looseBS, sec2origIdx, trackToMCs, {coll['GenParticles']}, trueV0s)")
                df = df.Define("v0n_class",       "v0ntruth.cls")
                df = df.Define("v0n_trueidx",     "v0ntruth.true_idx")
                df = df.Define("v0n_pairmult",    "v0ntruth.pair_mult")
                df = df.Define("v0n_trackshared", "v0ntruth.track_shared")
                df = df.Define("v0n_trk1",        "v0ntruth.trk1")
                df = df.Define("v0n_trk2",        "v0ntruth.trk2")
                df = df.Define("truev0_foundnew_any",     "FCCAnalyses::AlephTruth::trueV0FoundAny(trueV0s, v0ntruth)")
                df = df.Define("truev0_foundnew_correct", "FCCAnalyses::AlephTruth::trueV0FoundCorrect(trueV0s, v0ntruth, V0sNew_event)")

        ############################################# Standalone SV module ####################################################
        if self.do_svnew:
            # V0-first: svn_* = SV finding after masking the tight-claimed V0 tracks;
            # svm_* = unmasked control twin from the SAME event, for the interplay study.
            SVNEW = "FCCAnalyses::AlephSVNew"
            for pfx, mode in (("svn", f"{SVNEW}::SVN_MASK_MODE"), ("svm", f"{SVNEW}::SVN_MASK_NONE")):
                svn_expr = f"{SVNEW}::findSVs(SecondaryTracks_looseBS, VertexObject_looseBS, V0sNew_event, v0n_tight, {mode}, {BZ})"
                if self.do_pvnew:
                    # explicit entry guard (see V0sNew_event above)
                    svn_expr = ("pv_converged ? " + svn_expr +
                                " : FCCAnalyses::VertexingUtils::FCCAnalysesV0{}")
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

        ############################################# phi(1020) -> K+K- module (--phiKK) ######################################
        if self.ana_args.phiKK:
            a = self.ana_args
            # Pairs are formed from the FULL baseline-selected track list:
            # primary and secondary tracks alike, no masking by the PV split
            # or by other finders' claims (the Ks/Lambda veto is opt-in).
            # selBaselineOrigIdx maps that collection to the original Tracks,
            # which is how the per-track auxiliaries below are joined.
            df = df.Define("phikk_nvdet_all", "FCCAnalyses::AlephPhiKK::subdetHits(selBaselineOrigIdx, Tracks.subdetectorHitNumbers_begin, Tracks.subdetectorHitNumbers_end, _Tracks_subdetectorHitNumbers, 0)")
            df = df.Define("phikk_nitc_all",  "FCCAnalyses::AlephPhiKK::subdetHits(selBaselineOrigIdx, Tracks.subdetectorHitNumbers_begin, Tracks.subdetectorHitNumbers_end, _Tracks_subdetectorHitNumbers, 1)")
            df = df.Define("phikk_chi2ndf_all", "FCCAnalyses::AlephPhiKK::trackChi2Ndf(selBaselineOrigIdx, Tracks.chi2, Tracks.ndf)")
            df = df.Define("phikk_isprim_all",  "FCCAnalyses::AlephPhiKK::flagInSet(selBaselineOrigIdx, prim2origIdx)")
            if a.phiKKvetoV0:
                df = df.Define("phikk_veto_orig", "FCCAnalyses::AlephPhiKK::claimedOrigIdx(v0n_trk1_origIdx, v0n_trk2_origIdx, v0n_tight)")
            else:
                df = df.Define("phikk_veto_orig", "ROOT::VecOps::RVec<int>{}")
            # unset knobs are passed by constant NAME so the header stays the single source
            def knob(v, cname):
                return f"FCCAnalyses::AlephPhiKK::{cname}" if v is None else str(v)
            phikk_args = ", ".join([
                str(BZ), knob(a.phiKKmLo, "M_LO"), knob(a.phiKKmHi, "M_HI"), knob(a.phiKKchi2, "CHI2_CUT"),
                str(a.phiKKapBand), str(a.phiKKdpv), knob(a.phiKKdpvFid, "DPV_FID"),
                str(a.phiKKdpvSig),
                str(a.phiKKsigd0), str(a.phiKKminHits), str(a.phiKKtrkChi2),
                knob(a.phiKKpMin, "P_MIN_DEF"),
                "false" if a.phiKKnoSameSign else "true",
                "FCCAnalyses::AlephPhiKK::PRE_MARGIN",
                "phikk_veto_orig",
            ])
            df = df.Define("PhiKKCands_event",
                           "FCCAnalyses::AlephPhiKK::findPhiKK(trackstates_selected_baseline_flipped, "
                           "selBaselineOrigIdx, phikk_nvdet_all, phikk_nitc_all, phikk_chi2ndf_all, "
                           f"phikk_isprim_all, VertexObject_looseBS, {phikk_args})")
            df = df.Define("n_phikk_event", "int(PhiKKCands_event.invM.size())")
            for _b in ("invM", "p", "px", "py", "pz", "alpha", "qt", "bandEll",
                       "chi2", "vx", "vy", "vz", "dpv", "dpvSig", "same_sign",
                       "wp", "tight"):
                df = df.Define(f"phikk_{_b}", f"PhiKKCands_event.{_b}")
            for _t in ("trk1", "trk2"):
                for _b in ("origIdx", "q", "p", "costheta", "d0", "z0", "sigd0",
                           "nvdet", "nitc", "chi2ndf", "isprim"):
                    df = df.Define(f"phikk_{_t}_{_b}", f"PhiKKCands_event.{_t}_{_b}")
                # daughter dE/dx: STORED for the calibration, never selected on.
                # A failed leg copies the track omega into dQdx.value, so the
                # shared dEdxValid gate gives -1 in both value and error.
                for _det, _coll in (("pads", "dEdxPads"), ("wires", "dEdxWires")):
                    df = df.Define(f"phikk_{_t}_dEdx_{_det}_value",
                                   f"FCCAnalyses::AlephV0New::trackQuantityByIndex(phikk_{_t}_origIdx, {_coll}.dQdx.value, {_coll}.dQdx.value, {_coll}.dQdx.error, _{_coll}_track.index, _Tracks_trackStates)")
                    df = df.Define(f"phikk_{_t}_dEdx_{_det}_error",
                                   f"FCCAnalyses::AlephV0New::trackQuantityByIndex(phikk_{_t}_origIdx, {_coll}.dQdx.error, {_coll}.dQdx.value, {_coll}.dQdx.error, _{_coll}_track.index, _Tracks_trackStates)")
            if self.do_truth:
                df = df.Define("truePhis", f"FCCAnalyses::AlephPhiKK::findTruePhis({coll['GenParticles']}, mcToTracks)")
                for _b in ("mothPdg", "origin", "p", "pt", "costheta",
                           "px", "py", "pz", "vx", "vy", "vz", "nmatched"):
                    df = df.Define(f"truephi_{_b}", f"truePhis.{_b}")
                df = df.Define("truephi_dauPlus_p",  "truePhis.dauPlus_p")
                df = df.Define("truephi_dauMinus_p", "truePhis.dauMinus_p")
                df = df.Define("n_truephi_event",    "int(truePhis.idx.size())")
                df = df.Define("phikktruth", f"FCCAnalyses::AlephPhiKK::classifyPhiKK(PhiKKCands_event, trackToMCs, {coll['GenParticles']}, truePhis)")
                df = df.Define("phikk_class",        "phikktruth.cls")
                df = df.Define("phikk_trueidx",      "phikktruth.truephi_idx")
                for _t in ("trk1", "trk2"):
                    df = df.Define(f"phikk_{_t}_mcpdg",   f"phikktruth.{_t}_mcpdg")
                    df = df.Define(f"phikk_{_t}_mothpdg", f"phikktruth.{_t}_mothpdg")
                df = df.Define("truephi_found", "FCCAnalyses::AlephPhiKK::truePhiFound(truePhis, phikktruth)")

        ############################################# D*->D0(K pi) pi_slow module (--dstar) ###################################
        if self.ana_args.dstar:
            a = self.ana_args
            # Candidates are built from the FULL baseline-selected track list:
            # primary and secondary tracks alike, no masking by the PV split.
            # selBaselineOrigIdx maps that collection to the original Tracks,
            # which is how the per-track auxiliaries below are joined.
            df = df.Define("dstar_nvdet_all", "FCCAnalyses::AlephDstar::subdetHits(selBaselineOrigIdx, Tracks.subdetectorHitNumbers_begin, Tracks.subdetectorHitNumbers_end, _Tracks_subdetectorHitNumbers, 0)")
            df = df.Define("dstar_nitc_all",  "FCCAnalyses::AlephDstar::subdetHits(selBaselineOrigIdx, Tracks.subdetectorHitNumbers_begin, Tracks.subdetectorHitNumbers_end, _Tracks_subdetectorHitNumbers, 1)")
            df = df.Define("dstar_chi2ndf_all", "FCCAnalyses::AlephDstar::trackChi2Ndf(selBaselineOrigIdx, Tracks.chi2, Tracks.ndf)")
            df = df.Define("dstar_isprim_all",  "FCCAnalyses::AlephDstar::flagInSet(selBaselineOrigIdx, prim2origIdx)")
            # staging class of every pool track (0 prim / 1 sec / 2 neither)
            df = df.Define("dstar_pool_all",    "FCCAnalyses::AlephDstar::poolClass(selBaselineOrigIdx, prim2origIdx, sec2origIdx)")
            if a.dstarVetoV0:
                df = df.Define("dstar_veto_orig", "FCCAnalyses::AlephDstar::claimedOrigIdx(v0n_trk1_origIdx, v0n_trk2_origIdx, v0n_tight)")
            else:
                df = df.Define("dstar_veto_orig", "ROOT::VecOps::RVec<int>{}")
            dstar_claim = {"tight": 1, "loose": 2, "none": 0}[a.dstarClaim]
            dstar_args = ", ".join([
                str(BZ), str(a.dstarMLo), str(a.dstarMHi), str(a.dstarChi2),
                str(a.dstarDpvMax), str(a.dstarDmMax),
                str(a.dstarPMin), str(a.dstarPsMin),
                str(a.dstarSigd0), str(a.dstarMinHits), str(a.dstarTrkChi2),
                str(a.d0LooseDm), str(a.d0TightDm), str(a.d0TightDpvSig),
                str(a.d0TightCosPoint), str(a.d0TightCosStar),
                str(a.dstarLooseDm), str(a.dstarLooseDdm),
                str(a.dstarTightDm), str(a.dstarTightDdm),
                str(a.dstarTightPK), str(a.dstarTightPPi), str(a.dstarTightChi2),
                str(a.dstarTightPs), str(a.dstarTightCosPoint),
                "FCCAnalyses::AlephDstar::PRE_MARGIN",
                "true" if a.dstarCascade else "false", str(dstar_claim),
                "dstar_veto_orig",
            ])
            df = df.Define("DstarCands_event",
                           "FCCAnalyses::AlephDstar::findDstar(trackstates_selected_baseline_flipped, "
                           "selBaselineOrigIdx, dstar_nvdet_all, dstar_nitc_all, dstar_chi2ndf_all, "
                           f"dstar_isprim_all, dstar_pool_all, VertexObject_looseBS, {dstar_args})")
            df = df.Define("n_d0_event",    "int(DstarCands_event.d0.m_kpi.size())")
            df = df.Define("n_dstar_event", "int(DstarCands_event.ds.m_kpi.size())")
            # two-track fits actually performed (cache misses): the combinatorial cost
            df = df.Define("n_d0fits_event", "DstarCands_event.nfits")
            for _b in D0_CAND_BRANCHES:
                df = df.Define(f"d0_{_b}", f"DstarCands_event.d0.{_b}")
            for _b in DSTAR_CAND_BRANCHES:
                df = df.Define(f"dstar_{_b}", f"DstarCands_event.ds.{_b}")
            for _pfx, _mem in (("d0_trkK", "d0.trkK"), ("d0_trkPi", "d0.trkPi"),
                               ("dstar_trkK", "ds.trkK"), ("dstar_trkPi", "ds.trkPi"),
                               ("dstar_trkPis", "ds.trkPis")):
                for _b in DSTAR_TRK_BRANCHES:
                    df = df.Define(f"{_pfx}_{_b}", f"DstarCands_event.{_mem}.{_b}")
                # daughter dE/dx: STORED for the calibration, never selected on.
                # A failed leg copies the track omega into dQdx.value, so the
                # shared dEdxValid gate gives -1 in both value and error.
                for _det, _coll in (("pads", "dEdxPads"), ("wires", "dEdxWires")):
                    df = df.Define(f"{_pfx}_dEdx_{_det}_value",
                                   f"FCCAnalyses::AlephV0New::trackQuantityByIndex({_pfx}_origIdx, {_coll}.dQdx.value, {_coll}.dQdx.value, {_coll}.dQdx.error, _{_coll}_track.index, _Tracks_trackStates)")
                    df = df.Define(f"{_pfx}_dEdx_{_det}_error",
                                   f"FCCAnalyses::AlephV0New::trackQuantityByIndex({_pfx}_origIdx, {_coll}.dQdx.error, {_coll}.dQdx.value, {_coll}.dQdx.error, _{_coll}_track.index, _Tracks_trackStates)")
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
                for _t in ("trkK", "trkPi"):
                    df = df.Define(f"d0_{_t}_mcpdg",   f"d0truth.{_t}_mcpdg")
                    df = df.Define(f"d0_{_t}_mothpdg", f"d0truth.{_t}_mothpdg")
                for _t in ("trkK", "trkPi", "trkPis"):
                    df = df.Define(f"dstar_{_t}_mcpdg",   f"dstartruth.{_t}_mcpdg")
                    df = df.Define(f"dstar_{_t}_mothpdg", f"dstartruth.{_t}_mothpdg")
                # per-true-particle efficiency flags (class-1 candidate carrying the label)
                df = df.Define("trued0_found_loose",    "FCCAnalyses::AlephDstar::trueD0Found(trueD0s, d0truth, d0_loose)")
                df = df.Define("trued0_found_tight",    "FCCAnalyses::AlephDstar::trueD0Found(trueD0s, d0truth, d0_tight)")
                df = df.Define("truedstar_found_loose", "FCCAnalyses::AlephDstar::trueDstarFound(trueDstars, dstartruth, dstar_loose)")
                df = df.Define("truedstar_found_tight", "FCCAnalyses::AlephDstar::trueDstarFound(trueDstars, dstartruth, dstar_tight)")

        ############################################# Particle Flow Level Variables #######################################################
        df = df.Define("pfcand_isMu",     "AlephSelection::get_isType(jetConstitutentsTypes,2)")
        df = df.Define("pfcand_isEl",     "AlephSelection::get_isType(jetConstitutentsTypes,1)")
        df = df.Define("pfcand_isGamma",  "AlephSelection::get_isType(jetConstitutentsTypes,4)")
        df = df.Define("pfcand_isChargedHad", "AlephSelection::get_isType(jetConstitutentsTypes,0)")
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

        # subdetector hit counts per constituent (inside-out: VDET, ITC, TPC)
        df = df.Define("pfcand_nTrackHits_VDET", "AlephSelection::get_constituent_nTrackHits_VDET(jetc, TracksByRP, _Tracks_subdetectorHitNumbers)")
        df = df.Define("pfcand_nTrackHits_ITC",  "AlephSelection::get_constituent_nTrackHits_ITC(jetc, TracksByRP, _Tracks_subdetectorHitNumbers)")
        df = df.Define("pfcand_nTrackHits_TPC",  "AlephSelection::get_constituent_nTrackHits_TPC(jetc, TracksByRP, _Tracks_subdetectorHitNumbers)")

        df = df.Define("Bz", '1.5') # luka reads this from the event ? 

        ############################################# Track Parameters and Covariance #######################################################

        df = df.Define("TrackStateFlipped",f"AlephSelection::flipD0_copy( {coll['TrackState']} )")
        # re-index through the RecoParticle->Track relation so that .at(tracks_begin) picks the particle's own track
        df = df.Define("TrackStateByRP", "AlephSelection::reindexByRPLink(TrackStateFlipped, _RecoParticles_tracks.index)")

        # dxy/dz/phi0/C/ct w.r.t. the PV, computed once with every curvature term in cm (see analyzer.h)
        df = df.Define("pfcand_trkparPV",   "AlephSelection::get_constituent_trackParamsAtPV(jetc, TrackStateByRP, Vertex_refit_tlv, Bz)")
        df = df.Define("pfcand_dxy",        "pfcand_trkparPV.dxy")
        df = df.Define("pfcand_dz",         "pfcand_trkparPV.dz")
        df = df.Define("pfcand_phi0",       "pfcand_trkparPV.phi0")
        df = df.Define("pfcand_C",          "pfcand_trkparPV.C")
        df = df.Define("pfcand_ct",         "pfcand_trkparPV.ct")
        df = df.Define("pfcand_dptdpt",     f'JetConstituentsUtils::get_omega_cov(jetc, TrackStateByRP)') 
        df = df.Define("pfcand_dxydxy",     f'JetConstituentsUtils::get_d0_cov(jetc, TrackStateByRP)') 
        df = df.Define("pfcand_dzdz",       f'JetConstituentsUtils::get_z0_cov(jetc, TrackStateByRP)') 
        df = df.Define("pfcand_dphidphi",   f'JetConstituentsUtils::get_phi0_cov(jetc, TrackStateByRP)') 
        df = df.Define("pfcand_detadeta",   f'JetConstituentsUtils::get_tanlambda_cov(jetc, TrackStateByRP)') 
        df = df.Define("pfcand_dxydz",      f'JetConstituentsUtils::get_d0_z0_cov(jetc, TrackStateByRP)') # do we not need to recalculate this? 
        df = df.Define("pfcand_dphidxy",    f'JetConstituentsUtils::get_phi0_d0_cov(jetc, TrackStateByRP)') 
        df = df.Define("pfcand_phidz",      f'JetConstituentsUtils::get_phi0_z0_cov(jetc, TrackStateByRP)') 
        df = df.Define("pfcand_phictgtheta",f'JetConstituentsUtils::get_tanlambda_phi0_cov(jetc, TrackStateByRP)') 
        df = df.Define("pfcand_dxyctgtheta",f'JetConstituentsUtils::get_tanlambda_d0_cov(jetc, TrackStateByRP)') 
        df = df.Define("pfcand_dlambdadz",  f'JetConstituentsUtils::get_tanlambda_z0_cov(jetc, TrackStateByRP)') 
        df = df.Define("pfcand_cctgtheta",  f'JetConstituentsUtils::get_omega_tanlambda_cov(jetc, TrackStateByRP)') 
        df = df.Define("pfcand_phic",       f'JetConstituentsUtils::get_omega_phi0_cov(jetc, TrackStateByRP)') 
        df = df.Define("pfcand_dxyc",       f'JetConstituentsUtils::get_omega_d0_cov(jetc, TrackStateByRP)') 
        df = df.Define("pfcand_cdz",        f'JetConstituentsUtils::get_omega_z0_cov(jetc, TrackStateByRP)')

        ############################################# Btag Variables #######################################################

        df = df.Define("pfcand_btagSip2dVal",   "JetConstituentsUtils::get_Sip2dVal_clusterV(jets, pfcand_dxy, pfcand_phi0, Bz)") 
        df = df.Define("pfcand_btagSip2dSig",   "JetConstituentsUtils::get_Sip2dSig(pfcand_btagSip2dVal, pfcand_dxydxy)") 
        df = df.Define("pfcand_btagSip3dVal",   "JetConstituentsUtils::get_Sip3dVal_clusterV(jets, pfcand_dxy, pfcand_dz, pfcand_phi0, Bz)") 
        df = df.Define("pfcand_btagSip3dSig",   "JetConstituentsUtils::get_Sip3dSig(pfcand_btagSip3dVal, pfcand_dxydxy, pfcand_dzdz)") 
        df = df.Define("pfcand_btagJetDistVal","JetConstituentsUtils::get_JetDistVal_clusterV(jets, jetc, pfcand_dxy, pfcand_dz, pfcand_phi0, Bz)") 
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

        truth_branches = []
        if self.do_truth:
            truth_branches = [
                "truev0_pdg", "truev0_p", "truev0_costheta",
                "truev0_px", "truev0_py", "truev0_pz",
                "truev0_fd", "truev0_dpv",
                "truev0_nmatched", "truev0_nsec", "truev0_found_any", "truev0_found_correct",
                "truev0_x", "truev0_y", "truev0_z",
                "v0c_class", "v0c_trueidx", "v0c_pairmult", "v0c_trackshared",
                "v0c_alpha", "v0c_qt", "v0c_trk1", "v0c_trk2",
                "v0c_pdg", "v0c_invM", "v0c_dxyz", "v0c_p", "v0c_cosPointing",
                "v0c_vx", "v0c_vy", "v0c_vz",
            ]
        if self.do_v0new:
            truth_branches += [
                "n_v0n_event", "v0n_pdg", "v0n_invM", "v0n_alpha", "v0n_qt",
                "v0n_chi2", "v0n_dxyz", "v0n_p", "v0n_px", "v0n_py", "v0n_pz",
                "v0n_cosPointing", "v0n_pointSig",
                "v0n_tight", "v0n_bandSig", "v0n_massSig", "v0n_vx", "v0n_vy", "v0n_vz",
                # vertex-fit covariance + truth-free per-daughter joins/dE/dx
                "v0n_cov_xx", "v0n_cov_yx", "v0n_cov_yy",
                "v0n_cov_zx", "v0n_cov_zy", "v0n_cov_zz",
                "v0n_trk1_origIdx", "v0n_trk2_origIdx",
                "v0n_trk1_dEdx_pads_value", "v0n_trk1_dEdx_pads_error",
                "v0n_trk1_dEdx_wires_value", "v0n_trk1_dEdx_wires_error",
                "v0n_trk2_dEdx_pads_value", "v0n_trk2_dEdx_pads_error",
                "v0n_trk2_dEdx_wires_value", "v0n_trk2_dEdx_wires_error",
                # per-jet new-module V0s (mirror of the old v0_* block) -> jet-level apples-to-apples
                "n_v0njet_jets", "n_v0njet_ks", "n_v0njet_lambda",
                "v0njet_pdg", "v0njet_invM", "v0njet_chi2", "v0njet_chi2_norm",
                "v0njet_ndof", "v0njet_ntracks", "v0njet_p", "v0njet_prel",
                "v0njet_thetarel", "v0njet_phirel", "v0njet_dxy", "v0njet_dxyz",
                "v0njet_cosPointing", "v0njet_correctedMass",
                "v0njet_dx", "v0njet_dy", "v0njet_dz",
            ]
            if self.do_truth:
                truth_branches += [
                    "v0n_class", "v0n_trueidx", "v0n_pairmult", "v0n_trackshared",
                    "v0n_trk1", "v0n_trk2",
                    "truev0_foundnew_any", "truev0_foundnew_correct",
                ]
        if self.do_svnew:
            for pfx in ("svn", "svm"):
                truth_branches += [
                    f"n_{pfx}_event", f"{pfx}_mass", f"{pfx}_chi2", f"{pfx}_dxyz",
                    f"{pfx}_p", f"{pfx}_cosPointing", f"{pfx}_pointSig",
                    f"{pfx}_ntracks", f"{pfx}_sigL", f"{pfx}_trk_sv", f"{pfx}_trk_idx",
                    f"{pfx}_dx", f"{pfx}_dy", f"{pfx}_dz",
                    # SV vertex-fit covariance
                    f"{pfx}_cov_xx", f"{pfx}_cov_yx", f"{pfx}_cov_yy",
                    f"{pfx}_cov_zx", f"{pfx}_cov_zy", f"{pfx}_cov_zz",
                ]
            # V0 -> nearest-svn pointing feature
            truth_branches += ["v0n_svnCosPoint", "v0n_svnPointSig", "v0n_svnIdx"]
            # (sec2origIdx lives in the always-written list — it is truth-free)
        if self.do_pvnew:
            # the flag surface of the standalone PV fitter
            truth_branches += ["pv_converged", "pv_split_converged", "pv_trivial",
                               "pv_good"]
        if self.ana_args.phiKK:
            truth_branches += ["n_phikk_event"] + [
                f"phikk_{b}" for b in (
                    "invM", "p", "px", "py", "pz", "alpha", "qt", "bandEll",
                    "chi2", "vx", "vy", "vz", "dpv", "dpvSig", "same_sign",
                    "wp", "tight")
            ] + [
                f"phikk_{t}_{b}" for t in ("trk1", "trk2")
                for b in ("origIdx", "q", "p", "costheta", "d0", "z0", "sigd0",
                          "nvdet", "nitc", "chi2ndf", "isprim",
                          "dEdx_pads_value", "dEdx_pads_error",
                          "dEdx_wires_value", "dEdx_wires_error")
            ]
            if self.do_truth:
                truth_branches += ["n_truephi_event", "truephi_found"] + [
                    f"truephi_{b}" for b in (
                        "mothPdg", "origin", "p", "pt", "costheta",
                        "px", "py", "pz", "vx", "vy", "vz", "nmatched",
                        "dauPlus_p", "dauMinus_p")
                ] + ["phikk_class", "phikk_trueidx"] + [
                    f"phikk_{t}_{b}" for t in ("trk1", "trk2")
                    for b in ("mcpdg", "mothpdg")
                ]
        if self.ana_args.dstar:
            truth_branches += ["n_d0_event", "n_dstar_event", "n_d0fits_event"] + [
                f"d0_{b}" for b in D0_CAND_BRANCHES
            ] + [
                f"dstar_{b}" for b in DSTAR_CAND_BRANCHES
            ] + [
                f"{pfx}_{b}"
                for pfx in ("d0_trkK", "d0_trkPi", "dstar_trkK", "dstar_trkPi", "dstar_trkPis")
                for b in DSTAR_TRK_BRANCHES + (
                    "dEdx_pads_value", "dEdx_pads_error",
                    "dEdx_wires_value", "dEdx_wires_error")
            ]
            if self.do_truth:
                truth_branches += ["n_trued0_event", "n_truedstar_event",
                                   "trued0_found_loose", "trued0_found_tight",
                                   "truedstar_found_loose", "truedstar_found_tight",
                                   "d0_class", "d0_trueidx",
                                   "dstar_class", "dstar_trueidx"] + [
                    f"trued0_{b}" for b in TRUED0_BRANCHES
                ] + [
                    f"truedstar_{b}" for b in TRUEDSTAR_BRANCHES
                ] + [
                    f"d0_{t}_{b}" for t in ("trkK", "trkPi") for b in ("mcpdg", "mothpdg")
                ] + [
                    f"dstar_{t}_{b}" for t in ("trkK", "trkPi", "trkPis")
                    for b in ("mcpdg", "mothpdg")
                ]

        return truth_branches + [
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
