from argparse import ArgumentParser
import os 

class Analysis():
    # ------------------------------------------------------------------
    # Argument handling / I-O configuration
    # ------------------------------------------------------------------
    def __init__(self, cmdline_args):
        parser = ArgumentParser(
            description='Additional analysis arguments',
            usage='Provide additional arguments after analysis script path')
        parser.add_argument('--doData', action='store_true',
                            help='Run on data, instead of MC (which is the default behaviour).')
        parser.add_argument('--fraction', default=1.0, type=float,
                            help='Fraction of events to run, default is 1.0 = 100%')
        self.ana_args, _ = parser.parse_known_args(cmdline_args['remaining'])

        """
        unlike the Z->qq template, there is no truth-flavour split here
        (Z->tautau MC is not sub-divided by quark flavour), so the
        --MCtype/--MCflavour/outnames_dict machinery from the qq analysis
        has been dropped. If you do want to label 1-prong vs 3-prong output
        separately, that is done downstream from jet_nCharged, not here.
        """

        base_user_eos = "/eos/user/s/sjamoldi/stage1_Ztautau"
        if self.ana_args.doData:
            self.input_dir = "/eos/experiment/fcc/ee/analyses/case-studies/aleph/LEP1_DATA/"
            self.output_dir = os.path.join(base_user_eos, "data2_try") 
            self.process_list = {"1994": {"fraction": self.ana_args.fraction}}
            self.n_threads = 1
        else:
            self.input_dir = "/eos/user/s/sjamoldi/Ztautau/processed"
            self.output_dir = os.path.join(base_user_eos, "E_1p3p")
            self.process_list = {"MC": {"fraction": self.ana_args.fraction}}
            self.n_threads = 1
 
        self.include_paths = ["analyzer.h"]

    #------------------------------------------------------------------
    # Event loop
    # ------------------------------------------------------------------
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
        E_cm = 100
 
        df = df.Define("event_number", "EventHeader.eventNumber")
        df = df.Define("run_number", "EventHeader.runNumber")

        """
        MC truth info (kept for truth-matching / purity checks only -
        the kinematic and dE/dx variables below are all reco-level and
        identical in structure between data and MC).
        """

        #if self.ana_args.doData:
            #df = df.Filter("AlephSelection::sel_class_filter(16)(ClassBitset)   || AlephSelection::sel_class_filter(17)(ClassBitset) ")
            #df = df.Filter("AlephSelection::sel_class_filter(16)(ClassBitset) ") #basically each events are classified into distinct categories; this one says that we have hadrons in the final state so Z->qq
            #df = df.Define("jetPID", "-999")
            
        if not self.ana_args.doData:
            df = df.Define("MC_pdg", "MCParticle::get_pdg(MCParticles)")
            df = df.Define("MC_status", "MCParticle::get_genStatus(MCParticles)")
            df = df.Define("MC_momentum_x", "MCParticle::get_px(MCParticles)")
            df = df.Define("MC_momentum_y", "MCParticle::get_py(MCParticles)")
            df = df.Define("MC_momentum_z", "MCParticle::get_pz(MCParticles)")
            df = df.Define("RPtoMC", "FCCAnalyses::MCTruthMatching::build_RP2MC_via_track(RecoParticles, "
                                       "_RecoParticles_tracks.index, _trackMCLink_from.index, "
                                        "_trackMCLink_to.index, trackMCLink.weight)")
        # Flat RP kinematics (all PF candidates, unclustered)
        df = df.Define("RP_px", "ReconstructedParticle::get_px(RecoParticles)")
        df = df.Define("RP_py", "ReconstructedParticle::get_py(RecoParticles)")
        df = df.Define("RP_pz", "ReconstructedParticle::get_pz(RecoParticles)")
        df = df.Define("RP_e",  "ReconstructedParticle::get_e(RecoParticles)")
        df = df.Define("RP_m",  "ReconstructedParticle::get_mass(RecoParticles)")
        df = df.Define("RP_charge", "ReconstructedParticle::get_charge(RecoParticles)")


        # ------------------------------------------------------------
        # Two-jet clustering: each jet = one tau hemisphere
        # ------------------------------------------------------------
        # Anti-kt clustering and jet constituents ->bc we want to start with the soft ones

        df = df.Define("pjetc", "JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")
        df = df.Define("_jet", "JetClustering::clustering_ee_kt(2, 2, 1, 0)(pjetc)")
        df = df.Define("jets", "JetClusteringUtils::get_pseudoJets(_jet)")
        df = df.Define("_jetc", "JetClusteringUtils::get_constituents(_jet)")
        df = df.Define("jetc", "JetConstituentsUtils::build_constituents_cluster(RecoParticles, _jetc)")
        df = df.Define("jetConstituentsTypes", "AlephSelection::build_constituents_Types()(ParticleID, _jetc)")
        df = df.Define("event_njet", "JetConstituentsUtils::count_jets(jetc)")
        df = df.Filter("event_njet == 2", "2jets") #exactly two tau hemispheres


        if not self.ana_args.doData:
            df = df.Define("pfcand_mc_pdg", "JetConstituentsUtils::get_PIDs_cluster(RPtoMC.recind, RPtoMC.mcind, RecoParticles, MCParticles, _jetc)")
            df = df.Define("pfcand_mc_index", "AlephSelection::get_MCindex(RPtoMC.recind, RPtoMC.mcind, RecoParticles, _jetc)")
            df = df.Define("pfcand_mc_momentum", "AlephSelection::get_true_p_from_index(pfcand_mc_index, MCParticles)")
 
        # ------------------------------------------------------------
        # Particle-flow type flags, per jet per constituent
        # ------------------------------------------------------------
        df = df.Define("pfcand_isMu",         "AlephSelection::get_isType(jetConstituentsTypes,2)")
        df = df.Define("pfcand_isEl",         "AlephSelection::get_isType(jetConstituentsTypes,1)")
        df = df.Define("pfcand_isGamma",      "AlephSelection::get_isType(jetConstituentsTypes,4)")
        df = df.Define("pfcand_isChargedHad", "AlephSelection::get_isType(jetConstituentsTypes,0)")
        df = df.Define("pfcand_isNeutralHad", "AlephSelection::get_isType(jetConstituentsTypes,5)")
 

        # Charged multiplicity per jet == tau prong count.
        # AlephSelection::count_charged is NEW - see analyzer_additions.h
        df = df.Define("jet_nCharged", "AlephSelection::count_charged(pfcand_isChargedHad, pfcand_isMu, pfcand_isEl)")
 
        # Select one 1-prong hemisphere and one 3-prong hemisphere (either jet index)
        df = df.Filter("(jet_nCharged[0]==1 && jet_nCharged[1]==3) || (jet_nCharged[0]==3 && jet_nCharged[1]==1)", "1p3p_topology")

        """
        Prong count alone does not mean "pion": 1-prong is dominated by
        tau->e/mu nu-nu, and both prongs pick up extra pi0s from rho/a1
        decays. Neither cut below touches dE/dx, so the sample stays
        independent of what it is later used to calibrate.
        (a) lepton veto: every charged constituent in the jet must be a
        charged hadron, not e or mu - removes tau->e/mu nu-nu directly.
        """
        df = df.Define("jet_isHadronic", "AlephSelection::all_charged_are_hadronic(pfcand_isChargedHad, pfcand_isMu, pfcand_isEl)")
        # (b) no extra neutrals: total constituent count == charged count,
        # i.e. no accompanying pi0/photon - isolates pi-nu and 3pi-nu
        # specifically, excluding rho->pi pi0 and a1->3pi pi0.
        df = df.Define("j3", "jet_nCharged[0]==3 ? 0 : 1")
        df = df.Filter("jet_isHadronic[j3]", "3prong_hadronic")

        df = df.Define("j1", "jet_nCharged[0]==1 ? 0 : 1")

        # ------------------------------------------------------------
        # Vertex-proximity / TPC-crossing quality cuts
        # ------------------------------------------------------------
        # Raw (non-covariance) impact parameters. get_d0/get_z0, get_d0_cov/
        # get_z0_cov (same (jetc, TrackState) signature) - verify against your
        # build; if they are missing, AlephSelection::get_constituent_d0/z0 in
        # analyzer_additions.h are drop-in replacements with the same output shape.

        
        df = df.Define("TrackStateFlipped",f"AlephSelection::flipD0_copy( {coll['TrackState']} )")
        # re-index through the RecoParticle->Track relation so that .at(tracks_begin) picks the particle's own track
        df = df.Define("TrackStateByRP", "AlephSelection::reindexByRPLink(TrackStateFlipped, _RecoParticles_tracks.index)")
        df = df.Define("pfcand_d0",         "AlephSelection::get_constituent_trackParam(jetc, TrackStateByRP, 0)")
        df = df.Define("pfcand_z0",         "AlephSelection::get_constituent_trackParam(jetc, TrackStateByRP, 1)")
 
        # Per-jet worst-case |d0|, |z0| over its constituents - NEW helper.
        # Thresholds below (cm) are placeholders: tau flight cτ ≈ 87 μm, so
        # genuine 1-/3-prong decay products should sit very close to the
        # primary vertex; tune against your beamspot/IP resolution.
        df = df.Define("jet_max_absd0", "AlephSelection::max_abs_per_jet_charged(pfcand_d0, pfcand_isChargedHad, pfcand_isMu, pfcand_isEl)")
        df = df.Define("jet_max_absz0", "AlephSelection::max_abs_per_jet_charged(pfcand_z0, pfcand_isChargedHad, pfcand_isMu, pfcand_isEl)")

        
        df = df.Filter("jet_max_absd0[0]<2 && jet_max_absd0[1]<2 && "
               "jet_max_absz0[0]<10 && jet_max_absz0[1]<10", "IP_cuts")
        # Require every charged constituent to actually cross the TPC
        # (reuses the existing hit-counting helper from the qq analysis).
        df = df.Define("pfcand_nTrackHits_TPC", f"AlephSelection::get_constituent_nTrackHits_TPC(jetc, Tracks, _Tracks_subdetectorHitNumbers)")
        df = df.Define("jet_min_TPChits_charged", "AlephSelection::min_TPChits_charged(pfcand_nTrackHits_TPC, pfcand_isChargedHad, pfcand_isMu, pfcand_isEl)")
        df = df.Filter("jet_min_TPChits_charged[0] > 0 && jet_min_TPChits_charged[1] > 0", "TPC_hits")

        # ------------------------------------------------------------
        # Pion kinematics (= the tau daughters)
        # ------------------------------------------------------------
        df = df.Define("pfcand_e",  "JetConstituentsUtils::get_e(jetc)")

        df = df.Filter(f"All(pfcand_e[0] > 0.01*{E_cm}) && All(pfcand_e[1] > 0.01*{E_cm})", "per_track_energy")

        df = df.Filter(f"(Sum(pfcand_e[0] + Sum(pfcand_e[1])) > 0.1*{E_cm})", "total_energy")
        df = df.Define("pfcand_p",  "JetConstituentsUtils::get_p(jetc)")
        df = df.Define("pfcand_px", "JetConstituentsUtils::get_px(jetc)")
        df = df.Define("pfcand_py", "JetConstituentsUtils::get_py(jetc)")
        df = df.Define("pfcand_pz", "JetConstituentsUtils::get_pz(jetc)")
        df = df.Define("pfcand_charge", "JetConstituentsUtils::get_charge(jetc)")

         # ------------------------------------------------------------
        # dE/dx (unchanged from the Z->qq treemaker). These branches,
        # including the PID p-values, are the calibration TARGET - do not
        # filter on them above, or the sample purity becomes circular.
        # ------------------------------------------------------------
        df = df.Define("jet_constituents_dEdx_PIDhypo_pads_result", "AlephSelection::build_constituents_dEdx_PIDhypo()(RecoParticles, _RecoParticles_tracks.index, dEdxPads, _dEdxPads_track.index, _jetc, _Tracks_trackStates, false)")
        df = df.Define("jet_constituents_dEdx_pads_objs", "jet_constituents_dEdx_PIDhypo_pads_result.dedx_constituents")
        df = df.Define("pfcand_dEdx_pads_type",  "AlephSelection::get_dEdx_type(jet_constituents_dEdx_pads_objs)")
        df = df.Define("pfcand_dEdx_pads_value", "AlephSelection::get_dEdx_value(jet_constituents_dEdx_pads_objs)")
        df = df.Define("pfcand_dEdx_pads_error", "AlephSelection::get_dEdx_error(jet_constituents_dEdx_pads_objs)")
 
        df = df.Define("jet_constituents_PID_pvals_pads", "jet_constituents_dEdx_PIDhypo_pads_result.pid_array_constituents")
        df = df.Define("pfcand_PID_pval_pads_ele",    "AlephSelection::get_PID_pvalue(jet_constituents_PID_pvals_pads, 0)")
        df = df.Define("pfcand_PID_pval_pads_mu",     "AlephSelection::get_PID_pvalue(jet_constituents_PID_pvals_pads, 1)")
        df = df.Define("pfcand_PID_pval_pads_pi",     "AlephSelection::get_PID_pvalue(jet_constituents_PID_pvals_pads, 2)")
        df = df.Define("pfcand_PID_pval_pads_kaon",   "AlephSelection::get_PID_pvalue(jet_constituents_PID_pvals_pads, 3)")
        df = df.Define("pfcand_PID_pval_pads_proton", "AlephSelection::get_PID_pvalue(jet_constituents_PID_pvals_pads, 4)")
 
        df = df.Define("jet_constituents_dEdx_PIDhypo_wires_result", "AlephSelection::build_constituents_dEdx_PIDhypo()(RecoParticles, _RecoParticles_tracks.index, dEdxWires, _dEdxWires_track.index, _jetc, _Tracks_trackStates, true)")
        df = df.Define("jet_constituents_dEdx_wires_objs", "jet_constituents_dEdx_PIDhypo_wires_result.dedx_constituents")
        df = df.Define("pfcand_dEdx_wires_type",  "AlephSelection::get_dEdx_type(jet_constituents_dEdx_wires_objs)")
        df = df.Define("pfcand_dEdx_wires_value", "AlephSelection::get_dEdx_value(jet_constituents_dEdx_wires_objs)")
        df = df.Define("pfcand_dEdx_wires_error", "AlephSelection::get_dEdx_error(jet_constituents_dEdx_wires_objs)")
 
        df = df.Define("jet_constituents_PID_pvals_wires", "jet_constituents_dEdx_PIDhypo_wires_result.pid_array_constituents")
        df = df.Define("pfcand_PID_pval_wires_ele",    "AlephSelection::get_PID_pvalue(jet_constituents_PID_pvals_wires, 0)")
        df = df.Define("pfcand_PID_pval_wires_mu",     "AlephSelection::get_PID_pvalue(jet_constituents_PID_pvals_wires, 1)")
        df = df.Define("pfcand_PID_pval_wires_pi",     "AlephSelection::get_PID_pvalue(jet_constituents_PID_pvals_wires, 2)")
        df = df.Define("pfcand_PID_pval_wires_kaon",   "AlephSelection::get_PID_pvalue(jet_constituents_PID_pvals_wires, 3)")
        df = df.Define("pfcand_PID_pval_wires_proton", "AlephSelection::get_PID_pvalue(jet_constituents_PID_pvals_wires, 4)")

        # ------------------------------------------------------------
        # Invariant masses: tau visible mass (jet) and charged-only
        # ("daughters") mass. Both are NEW helpers.
        # ------------------------------------------------------------
        df = df.Define("jet_p4", "JetConstituentsUtils::compute_tlv_jets(jets)")
        df = df.Define("tau_visible_mass", "AlephSelection::get_masses(jet_p4)")
        df = df.Define("tau_charged_mass", "AlephSelection::get_charged_invmass(pfcand_px, pfcand_py, pfcand_pz, pfcand_e, "
                                             "pfcand_isChargedHad, pfcand_isMu, pfcand_isEl)")
        df = df.Define("tau3_visible_mass", "tau_visible_mass[j3]")
        df = df.Define("tau3_charged_mass", "tau_charged_mass[j3]")
 
        # Track-quality bookkeeping (kept from the qq treemaker)
        df = df.Define("pfcand_trackChi2",     "AlephSelection::get_constituent_trackChi2(jetc, Tracks)")
        df = df.Define("pfcand_trackNdof",     "AlephSelection::get_constituent_trackNdof(jetc, Tracks)")
        df = df.Define("pfcand_trackChi2Norm", "AlephSelection::get_constituent_trackChi2Norm(jetc, Tracks)")
        df = df.Define("pfcand_nTrackHits_VDET", "AlephSelection::get_constituent_nTrackHits_VDET(jetc, Tracks, _Tracks_subdetectorHitNumbers)")
        df = df.Define("pfcand_nTrackHits_ITC",  "AlephSelection::get_constituent_nTrackHits_ITC(jetc, Tracks, _Tracks_subdetectorHitNumbers)")
        # slice every per-jet / per-constituent branch down to just that hemisphere
        df = df.Define("pfcand3_isMu",         "pfcand_isMu[j3]")
        df = df.Define("pfcand3_isEl",         "pfcand_isEl[j3]")
        df = df.Define("pfcand3_isGamma",      "pfcand_isGamma[j3]")
        df = df.Define("pfcand3_isChargedHad", "pfcand_isChargedHad[j3]")
        df = df.Define("pfcand3_isNeutralHad", "pfcand_isNeutralHad[j3]")
        df = df.Define("pfcand1_isChargedHad", "pfcand_isChargedHad[j1]")

        df = df.Define("pfcand3_e",  "pfcand_e[j3]")
        df = df.Define("pfcand3_p",  "pfcand_p[j3]")
        df = df.Define("pfcand1_p",  "pfcand_p[j1]")
        df = df.Define("pfcand3_px", "pfcand_px[j3]")
        df = df.Define("pfcand3_py", "pfcand_py[j3]")
        df = df.Define("pfcand3_pz", "pfcand_pz[j3]")
        df = df.Define("pfcand3_charge", "pfcand_charge[j3]")

        df = df.Define("pfcand3_d0", "pfcand_d0[j3]")
        df = df.Define("pfcand3_z0", "pfcand_z0[j3]")

        df = df.Define("pfcand3_dEdx_wires_value",  "pfcand_dEdx_wires_value[j3]")
        df = df.Define("pfcand3_dEdx_wires_type",   "pfcand_dEdx_wires_type[j3]")
        df = df.Define("pfcand3_dEdx_wires_error",   "pfcand_dEdx_wires_error[j3]")

        df = df.Define("pfcand1_dEdx_wires_value",  "pfcand_dEdx_wires_value[j1]")
        df = df.Define("pfcand1_dEdx_wires_type",   "pfcand_dEdx_wires_type[j1]")
        df = df.Define("pfcand1_dEdx_wires_error",   "pfcand_dEdx_wires_error[j1]")
                
        

        
        if not self.ana_args.doData:
            df = df.Define("pfcand3_mc_pdg", "pfcand_mc_pdg[j3]")
            df = df.Define("pfcand1_mc_pdg", "pfcand_mc_pdg[j1]")
            """
            df = df.Define("jet_truth_decay_mode",
                "AlephSelection::get_jet_truth_decay_mode(pfcand_mc_index, MCParticles, "
                "_MCParticles_parents.index, _MCParticles_daughters.index)")

            df = df.Define("jet_truth_is_single_tau",
                "AlephSelection::jet_is_single_tau_matched(pfcand_mc_index, MCParticles, "
                "_MCParticles_parents.index)")

            df = df.Define("jet3_truth_decay_mode",    "jet_truth_decay_mode[j3]")
            df = df.Define("jet3_truth_is_single_tau", "jet_truth_is_single_tau[j3]")
            """

            df = df.Define("jet3_mc_index", "pfcand_mc_index[j3]")
            df = df.Define("MC_vertex_x", "MCParticle::get_vertex_x(MCParticles)")
            df = df.Define("MC_vertex_y", "MCParticle::get_vertex_y(MCParticles)")
            df = df.Define("MC_vertex_z", "MCParticle::get_vertex_z(MCParticles)")
            df = df.Define("jet3_same_vertex",
            "AlephSelection::jet_same_vertex(jet3_mc_index, MC_vertex_x, MC_vertex_y, MC_vertex_z, 1e-6f)")
            df = df.Define("truth_tau3pi_mass",
            "AlephSelection::truth_tau_mass_3pi_nu(jet3_mc_index, MCParticles, "
            "MC_vertex_x, MC_vertex_y, MC_vertex_z, 1e-3f)")

        report = df.Report()
        report.Print() 

        return df
        

    def output(self):
        branchList = [
            "event_number",
            "run_number",

            "jet_nCharged",
            "jet_max_absd0",
            "jet_max_absz0",
            "jet_min_TPChits_charged",
            "tau_visible_mass",
            "tau_charged_mass",

            "pfcand_isMu", "pfcand_isEl", "pfcand_isChargedHad", "pfcand_isGamma", "pfcand_isNeutralHad",
            "pfcand_e", "pfcand_p", "pfcand_px", "pfcand_py", "pfcand_pz", "pfcand_charge",
            "pfcand_d0", "pfcand_z0",
            "pfcand_nTrackHits_VDET", "pfcand_nTrackHits_ITC", "pfcand_nTrackHits_TPC",
            "pfcand_trackChi2", "pfcand_trackNdof", "pfcand_trackChi2Norm",

            "pfcand_dEdx_pads_type", "pfcand_dEdx_pads_value", "pfcand_dEdx_pads_error",
            "pfcand_PID_pval_pads_ele", "pfcand_PID_pval_pads_mu", "pfcand_PID_pval_pads_pi",
            "pfcand_PID_pval_pads_kaon", "pfcand_PID_pval_pads_proton",

            "pfcand_dEdx_wires_type", "pfcand_dEdx_wires_value", "pfcand_dEdx_wires_error",
            "pfcand_PID_pval_wires_ele", "pfcand_PID_pval_wires_mu", "pfcand_PID_pval_wires_pi",
            "pfcand_PID_pval_wires_kaon", "pfcand_PID_pval_wires_proton",

            "pfcand3_dEdx_wires_type", "pfcand3_dEdx_wires_value", "pfcand3_dEdx_wires_error",
            "pfcand3_p", "pfcand3_charge", "tau3_visible_mass", "tau3_charged_mass",
            "pfcand3_isMu", "pfcand3_isEl", "pfcand3_isGamma", "pfcand3_isChargedHad", "pfcand3_isNeutralHad",
            "pfcand1_isChargedHad", "pfcand1_dEdx_wires_value", "pfcand1_dEdx_wires_type",
            "pfcand1_dEdx_wires_error", "pfcand1_p",
        ]

        if not self.ana_args.doData:
            branchList += [
                "MC_pdg", "MC_status",
                "MC_momentum_x", "MC_momentum_y", "MC_momentum_z",

                "pfcand_mc_pdg", "pfcand_mc_index", "pfcand_mc_momentum",
                "pfcand3_mc_pdg", "pfcand1_mc_pdg",

                "jet3_mc_index",
                "MC_vertex_x", "MC_vertex_y", "MC_vertex_z",
                "jet3_same_vertex", "truth_tau3pi_mass",
            ]

        return branchList