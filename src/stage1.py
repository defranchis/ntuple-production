
import os
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
        parser.add_argument('--batch', action='store_true', 
                            help='Submit to HTCondor batch')
        parser.add_argument('--valid', action='store_true', 
                            help='Run tester file only for validation against Lukas ntuples.')
        parser.add_argument('--chunks', default=None, type=int,
                            help='Number of chunks per process/file')
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
        
        self.include_paths = ["analyzer.h"]

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


        # ==== Track selection (to harmonize with Luka's code)
        # Note: The selection strategy here only works if there is one trackstate stored pre track.
        # The code includes an assertion for that, if it is somehow not the case it will fail. 
        # df = df.Define("n_tracks_all", f"AlephSelection::select_tracks( {coll['PFTracks']} )")
        df = df.Define("n_tracks_all", "Tracks.size()")
        df = df.Define("chi2_tracks_all","AlephSelection::get_track_chi2( Tracks )") #TODO: use collection here
        df = df.Define("ndf_tracks_all","AlephSelection::get_track_ndf( Tracks )") #TODO: use collection here
        df = df.Define("chi2_o_ndf_tracks_all","AlephSelection::get_track_ndf( Tracks )") #TODO: use collection here
        
        # baseline track selection: positive definite cov matrix & chi2 < 10 
        df = df.Define("tracks_selected_baseline_result","AlephSelection::select_tracks_baseline( Tracks, _Tracks_trackStates )") #TODO: use collection here  0.75, 2.0
        df = df.Define("tracks_selected_baseline","tracks_selected_baseline_result.tracks") 
        df = df.Define("trackstates_selected_baseline","tracks_selected_baseline_result.trackStates") 

        # impose upper bounds on impact parameters to pre-select compatible tracks for the primary vertex fit 
        df = df.Define("tracks_selected_for_vertexfit_result","AlephSelection::select_tracks_impactparameters( tracks_selected_baseline_result, 0.75, 2.0 )") 
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

        # Luka's loose BS constraints from looking at data
        res_x_loose = 200. # in um
        res_y_loose = 100. # in um
        res_z_loose = 2. # in cm

        chi2max = 5. # the maximum chi2 under which tracks are compatible with vertex fit

        # Beamspot POSITION (the widths above are its size; this is its centre).
        # In simulation the beamspot is at the origin by construction. In data it is offset by
        # ~0.6 mm in x and ~0.2 mm in y, i.e. 2-3x the transverse widths used as the constraint,
        # so leaving it at 0 would bias the fit. Values are per-run, in the same 10um units as
        # the widths (see AlephSelection::get_beamspot in analyzer.h).
        # The json path is passed explicitly and lives on EOS: resolving it relative to the header
        # would break on condor, where analyzer.h is copied to the worker node and AFS may not be
        # readable. A copy is kept in the repo at Aleph/data/ as the version-controlled reference -
        # keep the two in sync. Override at runtime with $ALEPH_BEAMSPOT_JSON if needed.
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

        # Guard: with fewer than 2 IP-preselected tracks there is no meaningful primary vertex,
        # so return NO primary tracks (the PV fit then falls back to the dummy beamspot vertex).
        # FCCAnalyses' get_PrimaryTracks instead returns `seltracks` unchanged, i.e. the single
        # track - that is what the reference wrapper (getPrimaryTracks in analyzer_pvtools.cxx,
        # `if(tracksToUse.size() < 2){ return primaryTracks; }`) guards against. Without this we
        # get nPrim=1 where the reference has nPrim=0 (~1400 events / 1.05M in the full sweep).
        # note: the {{}} is an escaped literal {} for str.format - it is the empty RVec, not a placeholder
        df = df.Define("RecoedPrimaryTracks_looseBS", "trackstates_selected_for_vertexfit_flipped.size() < 2 ? ROOT::VecOps::RVec<edm4hep::TrackState>{{}} : VertexFitterSimple::get_PrimaryTracks(trackstates_selected_for_vertexfit_flipped, true, {},{},{}, Beamspot_x, Beamspot_y, Beamspot_z, {})".format(res_x_loose/10., res_y_loose/10., res_z_loose*1E03, chi2max)) # 10um as unit (x,y), 1cm as unit (z)
        df = df.Define("VertexObject_looseBS", "VertexFitterSimple::VertexFitter_Tk(1, RecoedPrimaryTracks_looseBS, true, {},{},{}, Beamspot_x, Beamspot_y, Beamspot_z)".format(res_x_loose/10., res_y_loose/10., res_z_loose*1E03)) # 10um as unit (x,y), 1cm as unit (z)
        df = df.Define("Vertex_refit_looseBS", "VertexingUtils::get_VertexData(VertexObject_looseBS)")
        df = df.Define("Vertex_refit_tlv", "TLorentzVector(Vertex_refit_looseBS.position.x, Vertex_refit_looseBS.position.y, Vertex_refit_looseBS.position.z, 0.)")
        # for retrieving secondary tracks, use the full list of selected tracks 
        df = df.Define("SecondaryTracks_looseBS", "VertexFitterSimple::get_NonPrimaryTracks(trackstates_selected_baseline_flipped, RecoedPrimaryTracks_looseBS)")

        df = df.Define("Vertex_refit_x", "Vertex_refit_looseBS.position.x")
        df = df.Define("Vertex_refit_y", "Vertex_refit_looseBS.position.y")
        df = df.Define("Vertex_refit_z", "Vertex_refit_looseBS.position.z")

        df = df.Define("n_primary_tracks", "ReconstructedParticle2Track::getTK_n(RecoedPrimaryTracks_looseBS)")
        df = df.Define("n_secondary_tracks", "ReconstructedParticle2Track::getTK_n(SecondaryTracks_looseBS)")

        # for comparison test, fit vertex with tracks all tracks:
        # df = df.Define("RecoedPrimaryTracks_looseBS_all_tracks", "VertexFitterSimple::get_PrimaryTracks(_Tracks_trackStates, true, {},{},{},0.,0.,0., {})".format(res_x_loose/10., res_y_loose/10., res_z_loose*1E03, chi2max)) # 10um as unit (x,y), 1cm as unit (z)
        # df = df.Define("VertexObject_looseBS_all_tracks", "VertexFitterSimple::VertexFitter_Tk(1, RecoedPrimaryTracks_looseBS_all_tracks, true, {},{},{},0.,0.,0.)".format(res_x_loose/10., res_y_loose/10., res_z_loose*1E03)) # 10um as unit (x,y), 1cm as unit (z)
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
        df = df.Define("SVs_looseBS", "FCCAnalyses::AlephSelection::get_SV_event_ALEPH("
            "SecondaryTracks_looseBS, "               # non-primary tracks
            "trackstates_selected_baseline_flipped, " # all tracks
            "VertexObject_looseBS, "                  # primary vertex
            "0.8, "                                   # dR prefilter cut
            "false)"                                  # exclusive V0 rejection (skip+break), matching FCCAnalyses@3a4de97 isV0 - the code that produced ntuples-withks
        )

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
        df = df.Define("pfcand_type",     "JetConstituentsUtils::get_type(jetc)") 
        df = df.Define("pfcand_erel",     "JetConstituentsUtils::get_erel_cluster(jets, jetc)")
        df = df.Define("pfcand_erel_log", "JetConstituentsUtils::get_erel_log_cluster(jets, jetc)")
        df = df.Define("pfcand_thetarel", "JetConstituentsUtils::get_thetarel_cluster(jets, jetc)")
        df = df.Define("pfcand_phirel",   "JetConstituentsUtils::get_phirel_cluster(jets, jetc)")

        # transverse momentum: ptrel is the ratio pT_constituent / pT_jet (same convention as erel)
        df = df.Define("pfcand_pt",        "JetConstituentsUtils::get_pt(jetc)")
        df = df.Define("pfcand_ptrel",     "AlephSelection::get_ptrel_cluster(jets, jetc)")
        df = df.Define("pfcand_ptrel_log", "AlephSelection::get_ptrel_log_cluster(jets, jetc)")

        # track fit quality per constituent (-1 for neutrals, which have no track)
        df = df.Define("pfcand_trackChi2",     "AlephSelection::get_constituent_trackChi2(jetc, Tracks)")
        df = df.Define("pfcand_trackNdof",     "AlephSelection::get_constituent_trackNdof(jetc, Tracks)")
        df = df.Define("pfcand_trackChi2Norm", "AlephSelection::get_constituent_trackChi2Norm(jetc, Tracks)")

        # subdetector hit counts per constituent (inside-out: VDET, ITC, TPC)
        df = df.Define("pfcand_nTrackHits_VDET", "AlephSelection::get_constituent_nTrackHits_VDET(jetc, Tracks, _Tracks_subdetectorHitNumbers)")
        df = df.Define("pfcand_nTrackHits_ITC",  "AlephSelection::get_constituent_nTrackHits_ITC(jetc, Tracks, _Tracks_subdetectorHitNumbers)")
        df = df.Define("pfcand_nTrackHits_TPC",  "AlephSelection::get_constituent_nTrackHits_TPC(jetc, Tracks, _Tracks_subdetectorHitNumbers)")

        df = df.Define("Bz", '1.5') # luka reads this from the event ? 

        ############################################# Track Parameters and Covariance #######################################################

        df = df.Define("TrackStateFlipped",f"AlephSelection::flipD0_copy( {coll['TrackState']} )")
        # re-index through the RecoParticle->Track relation so that .at(tracks_begin) picks the particle's own track
        df = df.Define("TrackStateByRP", "AlephSelection::reindexByRPLink(TrackStateFlipped, _RecoParticles_tracks.index)")

        df = df.Define("pfcand_dxy",        f'JetConstituentsUtils::XPtoPar_dxy(jetc, TrackStateByRP, Vertex_refit_tlv, Bz)') 
        df = df.Define("pfcand_dz",         f'JetConstituentsUtils::XPtoPar_dz(jetc, TrackStateByRP, Vertex_refit_tlv, Bz)') 
        df = df.Define("pfcand_phi0",       f'JetConstituentsUtils::XPtoPar_phi(jetc, TrackStateByRP, Vertex_refit_tlv, Bz)') 
        df = df.Define("pfcand_C",          f'JetConstituentsUtils::XPtoPar_C(jetc, TrackStateByRP, Bz)') 
        df = df.Define("pfcand_ct",         f'JetConstituentsUtils::XPtoPar_ct(jetc, TrackStateByRP, Bz)') 
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
        df = df.Define("jet_constituents_dEdx_PIDhypo_pads_result", "AlephSelection::build_constituents_dEdx_PIDhypo()(RecoParticles, _RecoParticles_tracks.index, dEdxPads, _dEdxPads_track.index, _jetc, false)" )
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
        df = df.Define("jet_constituents_dEdx_PIDhypo_wires_result", "AlephSelection::build_constituents_dEdx_PIDhypo()(RecoParticles, _RecoParticles_tracks.index, dEdxWires, _dEdxWires_track.index, _jetc, true)" )
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

        return [
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
            "pfcand_type",
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
