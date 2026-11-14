#ifndef JOINTVTX_EVENT_H
#define JOINTVTX_EVENT_H

/*
  Plain data of the joint PV + SV vertexer with V0s: the per-event input
  (Event), the output (Result) and the resolved configuration (Settings).

  Track arrays run over the event's baseline tracks in Tracks order; a "local"
  track index is a position in these arrays. Units: cm, GeV.
*/

#include <array>
#include <vector>

#include <Eigen/Dense>

#include "jointvtx_common.h"

namespace FCCAnalyses {
namespace AlephJointVtx {

struct Event {
  // per baseline track
  std::vector<Vec5> par;                    // raw state (D0, phi, omega, Z0, tanLambda)
  std::vector<std::array<double, 15>> cov;  // raw covariance, lower triangle (covMatrix[0:15])
  std::vector<Vec5> fpar;                   // fitter convention (D, phi0, C, z0, ct)
  std::vector<Mat5> fcov;
  std::vector<int> q;                       // -1 if omega_raw > 0, else +1
  std::vector<double> p;                    // ptConst / max(|omega_raw|, 1e-9) * sqrt(1 + tanLambda^2)
  std::vector<int> jet;                     // 0 / 1 = jet of the track's constituent, -1 none
  std::vector<int> nvdet, nitc, ntpc;       // subdetectorHitNumbers 0 / 1 / 2 (0 if absent)
  std::vector<int> origIdx;                 // index into Tracks
  // per event
  Vec3 bs = Vec3::Zero();                   // beam-spot centre
  std::array<double, 2> jetTheta{{JV_NAN, JV_NAN}};
  std::array<double, 2> jetPhi{{JV_NAN, JV_NAN}};
  std::array<double, 2> jetP{{JV_NAN, JV_NAN}};
  Vec3 thrust = Vec3::Zero();               // event thrust axis (only for non-finite jet axes)
  double ptConst = JV_PT_CONST;             // momentum constant p was computed with, GeV cm
  int nTracksAll = 0;                       // size of Tracks
  // V0s of the existing V0 finder, local track indices (-1 = not a baseline track)
  std::vector<int> v0T1, v0T2, v0Tight;
  int size() const { return static_cast<int>(q.size()); }
};

struct SvOut {
  Vec3 x = Vec3::Zero();
  Mat3 C = Mat3::Zero();
  double mass = 0.;        // pion hypothesis, at x
  Vec3 P = Vec3::Zero();   // summed track momenta at x
  std::vector<int> trk;    // local indices of the tracks assigned to this vertex
};

struct V0Out {
  int t1 = -1, t2 = -1;    // local indices: positive, negative track
  int pdg = 0;             // 310, 3122, 22
  int hyp = 0;             // V0Hyp
  Vec3 x = Vec3::Zero();
  Mat3 C = Mat3::Zero();
  double m = 0.;           // booked hypothesis; conversions: m_ee at the tangent point
  double w = 0.;           // final V0 weight
  double p = 0.;
  Vec3 P = Vec3::Zero();
};

// status of a Result. From 2 on, a fit or solve failed where the reference raises (its run ends): the value names
// the step, and the Result holds no vertices (pv and pvC NaN, no SVs, no V0s, empty per-track vectors)
enum ResultStatus : int {
  kResultOk = 0,
  kResultNotRun = 1,       // not filled by runEvent
  kResultFailRefPv = 2,    // reference PV: teardown fit, track chi2 at the jet-projected PV
  kResultFailSignif = 3,   // impact-parameter significances at the reference PV
  kResultFailSupply = 4,   // candidate supply: vertex-density maxima, candidate fits, merging
  kResultFailV0 = 5,       // V0 candidates
  kResultFailJoint = 6,    // first joint fit
  kResultFailSecond = 7,   // seeded second pass
  kResultFailSplit = 8,    // cascade split
  kResultFailOutput = 9    // significances at the final PV
};

struct Result {
  int status = kResultNotRun;
  Vec3 pv = Vec3::Constant(JV_NAN);
  Mat3 pvC = Mat3::Constant(JV_NAN);
  std::vector<SvOut> svs;
  std::vector<V0Out> v0s;  // weight above the output threshold
  // per baseline track
  std::vector<double> wpv;     // PV weight of the last PV refit (zero below the refit floor)
  std::vector<int> assign;     // -1 none, 0 PV, k >= 1 SV k - 1
  std::vector<double> wass;    // weight of the assignment
  std::vector<double> s3;      // 3D impact-parameter significance to the final PV
  std::vector<int> tv0;        // index into v0s when the V0 weight passes the hard threshold, -1 none
};

enum class BeamSpotSetting : int { S1 = 0, LR = 1 };

// widths and fitted tables of one beam-spot setting (generated tables header)
struct BeamSpotBundle;

// Resolved configuration; the nominal values are filled per beam-spot setting by nominalSettings().
struct Settings {
  const BeamSpotBundle* bundle = nullptr;
  Vec3 bsSig = Vec3::Constant(JV_NAN);   // beam-spot widths x, y, z
  double ptConst = JV_PT_CONST;          // momentum constant, GeV cm

  // track gate and reference PV
  int gateMinTpc = 0;
  bool hitRadius = true;
  bool removeTightV0 = false;
  double pvWindowD0 = JV_NAN, pvWindowZ0 = JV_NAN;
  double pvChi2max = JV_NAN;
  double pvTmin = JV_NAN;

  // seeds and acceptance
  double seedSig = JV_NAN, seedChi2 = JV_NAN, seedDr = JV_NAN, seedFsig = JV_NAN;
  bool seedCrossings = true, seedWindow = true, seedSameJet = false;
  bool accept = true, acceptFinal = true;
  double accCos = JV_NAN, accDisHi = JV_NAN, accDisLo = JV_NAN, accFsig2 = JV_NAN, accSigl = JV_NAN;
  int accMaxTrk = 0;
  double massMax = JV_NAN, vtxChi2ndf = JV_NAN;

  // growth and merging
  double growChi2 = JV_NAN;
  double mergeShare = JV_NAN, mergeSig = JV_NAN;
  bool mergeRequireShare = false;
  int minSvTracks = 0;

  // candidate scores
  double slChi2 = JV_NAN, slFsig = JV_NAN, slTrimChi2 = JV_NAN, slScoreMin = JV_NAN, slFinalMin = JV_NAN;
  bool slForward = true, slSymmetric = true, slScoreSeeded = true;

  // vertex-density supply
  double vdR0 = JV_NAN, vdFwdMargin = JV_NAN, vdGridChi2 = JV_NAN, vdGridCos = JV_NAN, vdGridPvChi2 = JV_NAN;
  double vdK0 = JV_NAN, vdMemberChi2 = JV_NAN, vdTube = JV_NAN, vdVmin = JV_NAN;
  int vdMaxStarts = 0, vdNClimb = 0;

  // joint soft-assignment fit
  std::vector<double> temps;
  double chi2Cut = JV_NAN, chi2CutNeg = JV_NAN, chi2CutSv = JV_NAN, assocChi2 = JV_NAN;
  double studentNu = JV_NAN, densK = JV_NAN, looEps = JV_NAN;
  double priorS = JV_NAN, svPriorS = JV_NAN, svPen = JV_NAN;
  bool lifetimePrior = true, priorT = false;
  bool freezeLast = true;
  int freezeN = 0;
  double svForward = JV_NAN, svMinP = JV_NAN, svRefitWmin = JV_NAN;
  int svMinVdet = 0;
  double hardMin = JV_NAN;

  // seeded second pass and cascade split
  int nIter = 0;
  double iterSig = JV_NAN;
  bool splitCasc = true, splitFinal = true;
  double splitDchi2 = JV_NAN, splitFsig = JV_NAN, splitOrder = JV_NAN, splitSig = JV_NAN;
  int splitNmax = 0, splitNmin = 0;

  // veto of V0-like pairs without the V0 class
  bool v0Veto = false;
  double ksWin = JV_NAN, lamWin = JV_NAN, eeWin = JV_NAN;

  // V0 class
  bool v0Class = false, v0Conv = false, v0Mass = true;
  double v0Admit = JV_NAN, v0Chi2Max = JV_NAN, v0ConvRmin = JV_NAN;
  std::array<double, 2> v0ConvHscale{{JV_NAN, JV_NAN}};
  double v0DzAbs = JV_NAN, v0DzSig = JV_NAN, v0HfChi2 = JV_NAN, v0HfSig = JV_NAN, v0HitScale = JV_NAN;
  int v0MaxComp = 0;
  double v0Nu = JV_NAN, v0OutMin = JV_NAN, v0PreKs = JV_NAN, v0PreLam = JV_NAN;
  std::array<double, 2> v0QtMin{{JV_NAN, JV_NAN}};
  double v0Rmax = JV_NAN, v0SmMin = JV_NAN, v0TwinExt = JV_NAN, v0TwinSig = JV_NAN;
  double v0WinKs = JV_NAN, v0WinLam = JV_NAN;
};

}  // namespace AlephJointVtx
}  // namespace FCCAnalyses

#endif
