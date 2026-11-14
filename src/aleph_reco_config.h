#ifndef ALEPH_RECO_CONFIG_H
#define ALEPH_RECO_CONFIG_H
// Primary-vertex reconstruction constants shared by stage1.py, the standalone
// fitter (analyzer_pvnew.h), the joint vertex fitter (jointvtx_tables.h) and
// the standalone examples. Dependency-free; lengths in cm unless the name says
// otherwise (solenoid field and track conventions: aleph_units.h).
namespace FCCAnalyses {
namespace AlephReco {
// impact-parameter preselection of the primary-vertex candidate tracks
constexpr double kPVTrackD0Max = 0.75;   // |d0| [cm]
constexpr double kPVTrackZ0Max = 2.0;    // |z0| [cm]
constexpr double kPVTrackZ0MaxBeamspot = 5.0;  // |z0| [cm] from the run beam spot

// beam-spot constraint of the primary-vertex fit: loose widths from data
constexpr double kBeamSigmaX_um = 200.;
constexpr double kBeamSigmaY_um = 100.;
constexpr double kBeamSigmaZ_cm = 2.;
// VertexFitter_Tk scales these widths from um to mm by 1e-3, but the track
// states here are in cm, so its argument unit is 10 um. get_PrimaryTracks
// (--oldPV) scales them by the same 1e-3 but also scales the track lengths by
// 1e-3, so its constraint is 1000x too loose, i.e. effectively absent.
constexpr double kBeamSigmaXFit = kBeamSigmaX_um / 10.;
constexpr double kBeamSigmaYFit = kBeamSigmaY_um / 10.;
constexpr double kBeamSigmaZFit = kBeamSigmaZ_cm * 1e4 / 10.;

constexpr double kPVChi2Max = 5.;        // a track stays in the primary-vertex fit below this chi2
}  // namespace AlephReco
}  // namespace FCCAnalyses
#endif
