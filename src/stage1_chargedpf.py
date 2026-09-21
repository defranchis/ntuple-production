from argparse import ArgumentParser

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
        parser.add_argument('--inputDir', default=None, type=str,
                            help='Override the default EDM4HEP input directory.')
        parser.add_argument('--outdir', default=None, type=str,
                            help='Override the default output directory.')
        parser.add_argument('--nthreads', default=16, type=int,
                            help='Number of threads to run on, default is 16.')
        # Parse additional arguments not known to the FCCAnalyses parsers
        # All command line arguments know to fccanalysis are provided in the
        # `cmdline_arg` dictionary.
        self.ana_args, _ = parser.parse_known_args(cmdline_args['remaining'])

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
        if self.ana_args.doData and self.ana_args.MCflavour:
            print("----> WARNING: Incompatible input arguments: --MCflavour defined with --doData, will be ignored.")

        if self.ana_args.MCtype and not self.ana_args.MCtype in outnames_dict:
            print("----> ERROR: Requested unknown --MCtype. Currently only zqq available.")
            exit()

        if not self.ana_args.doData and not self.ana_args.MCflavour:
            print(f"----> ERROR: Requested MC run but did not specify --MCflavour. Please pick one..")
            exit()

        if not self.ana_args.doData and not self.ana_args.MCflavour in outnames_dict[self.ana_args.MCtype]:
            print(f"----> ERROR: Requested unknown --MCflavour for --MCtype {self.ana_args.MCtype}. Check the dictionary.")
            exit()

        #set the input/output directories:
        if self.ana_args.doData:
            self.input_dir = "/eos/experiment/fcc/ee/analyses/case-studies/aleph/LEP1_DATA/"
            self.output_dir = f"/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedData/{self.ana_args.year}/stage1_chargedpf/{self.ana_args.tag}"
            self.process_list = {
                "1994" : {"fraction" : self.ana_args.fraction},
            }

        else:
            self.input_dir = f"/eos/experiment/aleph/EDM4HEP/MC/{self.ana_args.year}/"
            output_name = outnames_dict[self.ana_args.MCtype][self.ana_args.MCflavour]
            self.output_dir = f"/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedMC/{self.ana_args.year}/{self.ana_args.MCtype}/stage1_chargedpf/{self.ana_args.tag}"
            self.process_list = {
                "QQB" : {"fraction" : self.ana_args.fraction, "output":output_name},
            }

        if self.ana_args.inputDir:
            self.input_dir = self.ana_args.inputDir
        if self.ana_args.outdir:
            self.output_dir = self.ana_args.outdir

        #set run options:
        self.n_threads = self.ana_args.nthreads
        self.include_paths = ["analyzer_chargedpf.h"]

    def analyzers(self, df):

        if self.ana_args.doData:
            df = df.Filter("AlephSelection::sel_class_filter(16)(ClassBitset) ")
            df = df.Define("jetPID", "AlephChargedPF::kNoTruthFloat")
        else:
            # Using Classbit to filter out QQbar samples and then get a specific flavor of jets
            # d-quark: 1, u-quark:2, s-quark:3, c-quark:4, b-quark: 5
            df = df.Define("jetPID", "AlephSelection::getJetPID(ClassBitset, MCParticles)")
            df = df.Filter(f"jetPID == {self.ana_args.MCflavour}")

        df = df.Define("event_number", "EventHeader.eventNumber")
        df = df.Define("run_number", "EventHeader.runNumber")

        ############################################# Candidate selection #######################################################
        # charged ReconstructedParticles with their own track, as a single group so that the
        # per-constituent helpers of analyzer.h can be reused unchanged
        df = df.Define("cpf_idx", "AlephChargedPF::select_indices(RecoParticles, _RecoParticles_tracks.index)")
        df = df.Define("cpf_group", "JetConstituentsUtils::build_constituents_cluster(RecoParticles, cpf_idx)")
        df = df.Define("cpf_rp", "AlephChargedPF::flatten(cpf_group)")
        df = df.Define("n_cpf", "int(cpf_rp.size())")

        ############################################# Kinematics #######################################################
        df = df.Define("cpf_px",     "ReconstructedParticle::get_px(cpf_rp)")
        df = df.Define("cpf_py",     "ReconstructedParticle::get_py(cpf_rp)")
        df = df.Define("cpf_pz",     "ReconstructedParticle::get_pz(cpf_rp)")
        df = df.Define("cpf_e",      "ReconstructedParticle::get_e(cpf_rp)")
        df = df.Define("cpf_p",      "ReconstructedParticle::get_p(cpf_rp)")
        df = df.Define("cpf_pt",     "ReconstructedParticle::get_pt(cpf_rp)")
        df = df.Define("cpf_theta",  "ReconstructedParticle::get_theta(cpf_rp)")
        df = df.Define("cpf_phi",    "ReconstructedParticle::get_phi(cpf_rp)")
        df = df.Define("cpf_charge", "ReconstructedParticle::get_charge(cpf_rp)")
        df = df.Define("cpf_mass",   "ReconstructedParticle::get_mass(cpf_rp)")
        df = df.Define("cpf_type",   "ReconstructedParticle::get_type(cpf_rp)")

        ############################################# PID flags #######################################################
        df = df.Define("cpf_types_group",  "AlephSelection::build_constituents_Types()(ParticleID, cpf_idx)")
        df = df.Define("cpf_pidType",      "AlephChargedPF::flatten(cpf_types_group)")
        df = df.Define("cpf_isMu",         "AlephChargedPF::flatten(AlephSelection::get_isType(cpf_types_group, 2))")
        df = df.Define("cpf_isEl",         "AlephChargedPF::flatten(AlephSelection::get_isType(cpf_types_group, 1))")
        df = df.Define("cpf_isChargedHad", "AlephChargedPF::flatten(AlephSelection::get_isType(cpf_types_group, 0))")

        ############################################# Track parameters #######################################################
        # tracks and track states re-indexed through the RecoParticle->Track relation (see analyzer.h)
        df = df.Define("TracksByRP", "AlephSelection::reindexByRPLink(Tracks, _RecoParticles_tracks.index)")
        df = df.Define("TrackStateFlipped", "AlephSelection::flipD0_copy(_Tracks_trackStates)")
        df = df.Define("TrackStateByRP", "AlephSelection::trackStatesByRPLink(TrackStateFlipped, Tracks, _RecoParticles_tracks.index)")

        df = df.Define("cpf_trackIdx", "AlephChargedPF::get_trackIndex(cpf_rp, _RecoParticles_tracks.index)")
        df = df.Define("cpf_d0",        "AlephChargedPF::flatten(AlephSelection::get_constituent_D0(cpf_group, TrackStateByRP))")
        df = df.Define("cpf_z0",        "AlephChargedPF::flatten(AlephSelection::get_constituent_Z0(cpf_group, TrackStateByRP))")
        df = df.Define("cpf_phi0",      "AlephChargedPF::get_trackStateParam(cpf_rp, TrackStateByRP, 0)")
        df = df.Define("cpf_omega",     "AlephChargedPF::get_trackStateParam(cpf_rp, TrackStateByRP, 1)")
        df = df.Define("cpf_tanLambda", "AlephChargedPF::get_trackStateParam(cpf_rp, TrackStateByRP, 2)")

        # covariance, lower-triangular in (d0, phi0, omega, z0, tanLambda): cov(a,b) at a*(a+1)/2 + b
        df = df.Define("cpf_d0d0",       "AlephChargedPF::flatten(AlephSelection::get_constituent_trackCov(cpf_group, TrackStateByRP, 0))")
        df = df.Define("cpf_phi0phi0",   "AlephChargedPF::flatten(AlephSelection::get_constituent_trackCov(cpf_group, TrackStateByRP, 2))")
        df = df.Define("cpf_omegaomega", "AlephChargedPF::flatten(AlephSelection::get_constituent_trackCov(cpf_group, TrackStateByRP, 5))")
        df = df.Define("cpf_z0z0",       "AlephChargedPF::flatten(AlephSelection::get_constituent_trackCov(cpf_group, TrackStateByRP, 9))")
        df = df.Define("cpf_tanLtanL",   "AlephChargedPF::flatten(AlephSelection::get_constituent_trackCov(cpf_group, TrackStateByRP, 14))")

        df = df.Define("cpf_trackChi2", "AlephChargedPF::flatten(AlephSelection::get_constituent_trackChi2(cpf_group, TracksByRP))")
        df = df.Define("cpf_trackNdof", "AlephChargedPF::flatten(AlephSelection::get_constituent_trackNdof(cpf_group, TracksByRP))")

        # subdetector hit counts (inside-out: VDET, ITC, TPC)
        df = df.Define("cpf_nTrackHits_VDET", "AlephChargedPF::flatten(AlephSelection::get_constituent_nTrackHits_VDET(cpf_group, TracksByRP, _Tracks_subdetectorHitNumbers))")
        df = df.Define("cpf_nTrackHits_ITC",  "AlephChargedPF::flatten(AlephSelection::get_constituent_nTrackHits_ITC(cpf_group, TracksByRP, _Tracks_subdetectorHitNumbers))")
        df = df.Define("cpf_nTrackHits_TPC",  "AlephChargedPF::flatten(AlephSelection::get_constituent_nTrackHits_TPC(cpf_group, TracksByRP, _Tracks_subdetectorHitNumbers))")

        ############################################# dE/dx and PID p-values #######################################################
        ## Pads
        df = df.Define("cpf_dEdx_pads_result", "AlephSelection::build_constituents_dEdx_PIDhypo()(RecoParticles, _RecoParticles_tracks.index, dEdxPads, _dEdxPads_track.index, cpf_idx, Tracks, _Tracks_trackStates, false, true)")
        df = df.Define("cpf_dEdx_pads_objs",  "cpf_dEdx_pads_result.dedx_constituents")
        df = df.Define("cpf_dEdx_pads_value", "AlephChargedPF::flatten(AlephSelection::get_dEdx_value(cpf_dEdx_pads_objs))")
        df = df.Define("cpf_dEdx_pads_error", "AlephChargedPF::flatten(AlephSelection::get_dEdx_error(cpf_dEdx_pads_objs))")
        df = df.Define("cpf_dEdx_pads_type",  "AlephChargedPF::flatten(AlephSelection::get_dEdx_type(cpf_dEdx_pads_objs))")

        ## p-values from Bethe-Bloch dE/dx vs p fits - order: ["e", "mu", "pi", "K", "p"]
        df = df.Define("cpf_pads_pvals", "cpf_dEdx_pads_result.pid_array_constituents")
        df = df.Define("cpf_PID_pval_pads_ele",    "AlephChargedPF::flatten(AlephSelection::get_PID_pvalue(cpf_pads_pvals, 0))")
        df = df.Define("cpf_PID_pval_pads_mu",     "AlephChargedPF::flatten(AlephSelection::get_PID_pvalue(cpf_pads_pvals, 1))")
        df = df.Define("cpf_PID_pval_pads_pi",     "AlephChargedPF::flatten(AlephSelection::get_PID_pvalue(cpf_pads_pvals, 2))")
        df = df.Define("cpf_PID_pval_pads_kaon",   "AlephChargedPF::flatten(AlephSelection::get_PID_pvalue(cpf_pads_pvals, 3))")
        df = df.Define("cpf_PID_pval_pads_proton", "AlephChargedPF::flatten(AlephSelection::get_PID_pvalue(cpf_pads_pvals, 4))")

        ## Wires
        df = df.Define("cpf_dEdx_wires_result", "AlephSelection::build_constituents_dEdx_PIDhypo()(RecoParticles, _RecoParticles_tracks.index, dEdxWires, _dEdxWires_track.index, cpf_idx, Tracks, _Tracks_trackStates, true, true)")
        df = df.Define("cpf_dEdx_wires_objs",  "cpf_dEdx_wires_result.dedx_constituents")
        df = df.Define("cpf_dEdx_wires_value", "AlephChargedPF::flatten(AlephSelection::get_dEdx_value(cpf_dEdx_wires_objs))")
        df = df.Define("cpf_dEdx_wires_error", "AlephChargedPF::flatten(AlephSelection::get_dEdx_error(cpf_dEdx_wires_objs))")
        df = df.Define("cpf_dEdx_wires_type",  "AlephChargedPF::flatten(AlephSelection::get_dEdx_type(cpf_dEdx_wires_objs))")

        df = df.Define("cpf_wires_pvals", "cpf_dEdx_wires_result.pid_array_constituents")
        df = df.Define("cpf_PID_pval_wires_ele",    "AlephChargedPF::flatten(AlephSelection::get_PID_pvalue(cpf_wires_pvals, 0))")
        df = df.Define("cpf_PID_pval_wires_mu",     "AlephChargedPF::flatten(AlephSelection::get_PID_pvalue(cpf_wires_pvals, 1))")
        df = df.Define("cpf_PID_pval_wires_pi",     "AlephChargedPF::flatten(AlephSelection::get_PID_pvalue(cpf_wires_pvals, 2))")
        df = df.Define("cpf_PID_pval_wires_kaon",   "AlephChargedPF::flatten(AlephSelection::get_PID_pvalue(cpf_wires_pvals, 3))")
        df = df.Define("cpf_PID_pval_wires_proton", "AlephChargedPF::flatten(AlephSelection::get_PID_pvalue(cpf_wires_pvals, 4))")

        ############################################# MC truth #######################################################
        if self.ana_args.doData:
            df = df.Define("cpf_mc_pdg",        "AlephChargedPF::constant_int(cpf_rp.size(), AlephChargedPF::kNoTruthInt)")
            df = df.Define("cpf_mc_parent_pdg", "AlephChargedPF::constant_int(cpf_rp.size(), AlephChargedPF::kNoTruthInt)")
            df = df.Define("cpf_mc_p",          "AlephChargedPF::constant(cpf_rp.size(), AlephChargedPF::kNoTruthFloat)")
            df = df.Define("cpf_mc_vertex_x",   "AlephChargedPF::constant(cpf_rp.size(), AlephChargedPF::kNoTruthFloat)")
            df = df.Define("cpf_mc_vertex_y",   "AlephChargedPF::constant(cpf_rp.size(), AlephChargedPF::kNoTruthFloat)")
            df = df.Define("cpf_mc_vertex_z",   "AlephChargedPF::constant(cpf_rp.size(), AlephChargedPF::kNoTruthFloat)")
            df = df.Define("cpf_mc_n_links",    "AlephChargedPF::constant_int(cpf_rp.size(), 0)")
        else:
            df = df.Define("cpf_mc", "AlephChargedPF::get_mcTruth(cpf_trackIdx, _trackMCLink_from.index, _trackMCLink_to.index, MCParticles, _MCParticles_parents.index)")
            df = df.Define("cpf_mc_pdg",        "cpf_mc.pdg")
            df = df.Define("cpf_mc_parent_pdg", "cpf_mc.parent_pdg")
            df = df.Define("cpf_mc_p",          "cpf_mc.p")
            df = df.Define("cpf_mc_vertex_x",   "cpf_mc.vx")
            df = df.Define("cpf_mc_vertex_y",   "cpf_mc.vy")
            df = df.Define("cpf_mc_vertex_z",   "cpf_mc.vz")
            df = df.Define("cpf_mc_n_links",    "cpf_mc.n_links")

        return df

    def output(self):

        return [
            # Event variables
            "run_number",
            "event_number",
            "n_cpf",
            "jetPID",

            # Kinematics
            "cpf_px",
            "cpf_py",
            "cpf_pz",
            "cpf_e",
            "cpf_p",
            "cpf_pt",
            "cpf_theta",
            "cpf_phi",
            "cpf_charge",
            "cpf_mass",
            "cpf_type",

            # PID flags
            "cpf_isMu",
            "cpf_isEl",
            "cpf_isChargedHad",
            "cpf_pidType",

            # dE/dx
            "cpf_dEdx_pads_value",
            "cpf_dEdx_pads_error",
            "cpf_dEdx_pads_type",
            "cpf_PID_pval_pads_ele",
            "cpf_PID_pval_pads_mu",
            "cpf_PID_pval_pads_pi",
            "cpf_PID_pval_pads_kaon",
            "cpf_PID_pval_pads_proton",
            "cpf_dEdx_wires_value",
            "cpf_dEdx_wires_error",
            "cpf_dEdx_wires_type",
            "cpf_PID_pval_wires_ele",
            "cpf_PID_pval_wires_mu",
            "cpf_PID_pval_wires_pi",
            "cpf_PID_pval_wires_kaon",
            "cpf_PID_pval_wires_proton",

            # Track parameters
            "cpf_d0",
            "cpf_z0",
            "cpf_phi0",
            "cpf_omega",
            "cpf_tanLambda",
            "cpf_d0d0",
            "cpf_z0z0",
            "cpf_phi0phi0",
            "cpf_omegaomega",
            "cpf_tanLtanL",
            "cpf_trackChi2",
            "cpf_trackNdof",
            "cpf_nTrackHits_VDET",
            "cpf_nTrackHits_ITC",
            "cpf_nTrackHits_TPC",
            "cpf_trackIdx",

            # MC truth
            "cpf_mc_pdg",
            "cpf_mc_parent_pdg",
            "cpf_mc_p",
            "cpf_mc_vertex_x",
            "cpf_mc_vertex_y",
            "cpf_mc_vertex_z",
            "cpf_mc_n_links",
            ]
