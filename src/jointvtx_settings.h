#ifndef JOINTVTX_SETTINGS_H
#define JOINTVTX_SETTINGS_H

// GENERATED FILE, do not edit by hand (gen_tables.py, with jointvtx_tables.h).
// Nominal Settings of one beam-spot setting: every member from the generated constant of its configuration
// key, the beam-spot widths and tables from the bundle, the momentum constant of stage1.

#include "jointvtx_event.h"
#include "jointvtx_tables.h"

namespace FCCAnalyses {
namespace AlephJointVtx {

inline Settings nominalSettings(BeamSpotSetting b) {
  Settings s;
  s.bundle = (b == BeamSpotSetting::LR) ? &JV_BUNDLE_LR : &JV_BUNDLE_S1;
  for (int i = 0; i < 3; ++i) s.bsSig[i] = s.bundle->sigma[i];
  s.ptConst = JV_PT_CONST;
  s.gateMinTpc = JV_GATE_MIN_TPC;
  s.hitRadius = JV_HIT_RADIUS;
  s.removeTightV0 = JV_REMOVE_TIGHT_V0;
  s.pvWindowD0 = JV_PV_WINDOW_D0;
  s.pvWindowZ0 = JV_PV_WINDOW_Z0;
  s.pvChi2max = JV_PV_CHI2MAX;
  s.pvTmin = JV_PV_TMIN;
  s.seedSig = JV_SEED_SIG;
  s.seedChi2 = JV_SEED_CHI2;
  s.seedDr = JV_SEED_DR;
  s.seedFsig = JV_SEED_FSIG;
  s.seedCrossings = JV_SEED_CROSSINGS;
  s.seedWindow = JV_SEED_WINDOW;
  s.seedSameJet = JV_SEED_SAME_JET;
  s.accept = JV_ACCEPT;
  s.acceptFinal = JV_ACCEPT_FINAL;
  s.accCos = JV_ACC_COS;
  s.accDisHi = JV_ACC_DIS_HI;
  s.accDisLo = JV_ACC_DIS_LO;
  s.accFsig2 = JV_ACC_FSIG2;
  s.accSigl = JV_ACC_SIGL;
  s.accMaxTrk = JV_ACC_MAX_TRK;
  s.massMax = JV_MASS_MAX;
  s.vtxChi2ndf = JV_VTX_CHI2NDF;
  s.growChi2 = JV_GROW_CHI2;
  s.mergeShare = JV_MERGE_SHARE;
  s.mergeSig = JV_MERGE_SIG;
  s.mergeRequireShare = JV_MERGE_REQUIRE_SHARE;
  s.minSvTracks = JV_MIN_SV_TRACKS;
  s.slChi2 = JV_SL_CHI2;
  s.slFsig = JV_SL_FSIG;
  s.slTrimChi2 = JV_SL_TRIM_CHI2;
  s.slScoreMin = JV_SL_SCORE_MIN;
  s.slFinalMin = JV_SL_FINAL_MIN;
  s.slForward = JV_SL_FORWARD;
  s.slSymmetric = JV_SL_SYMMETRIC;
  s.slScoreSeeded = JV_SL_SCORE_SEEDED;
  s.vdR0 = JV_VD_R0;
  s.vdFwdMargin = JV_VD_FWD_MARGIN;
  s.vdGridChi2 = JV_VD_GRID_CHI2;
  s.vdGridCos = JV_VD_GRID_COS;
  s.vdGridPvChi2 = JV_VD_GRID_PVCHI2;
  s.vdK0 = JV_VD_K0;
  s.vdMemberChi2 = JV_VD_MEMBER_CHI2;
  s.vdTube = JV_VD_TUBE;
  s.vdVmin = JV_VD_VMIN;
  s.vdMaxStarts = JV_VD_MAX_STARTS;
  s.vdNClimb = JV_VD_N_CLIMB;
  s.temps.assign(JV_TEMPS.begin(), JV_TEMPS.end());
  s.chi2Cut = JV_CHI2_CUT;
  s.chi2CutNeg = JV_CHI2_CUT_NEG;
  s.chi2CutSv = JV_CHI2_CUT_SV;
  s.assocChi2 = JV_ASSOC_CHI2;
  s.studentNu = JV_STUDENT_NU;
  s.densK = JV_DENS_K;
  s.looEps = JV_LOO_EPS;
  s.priorS = JV_PRIOR_S;
  s.svPriorS = JV_SV_PRIOR_S;
  s.svPen = JV_SV_PEN;
  s.lifetimePrior = JV_LIFETIME_PRIOR;
  s.priorT = JV_PRIOR_T;
  s.freezeLast = JV_FREEZE_LAST;
  s.freezeN = JV_FREEZE_N;
  s.svForward = JV_SV_FORWARD;
  s.svMinP = JV_SV_MIN_P;
  s.svRefitWmin = JV_SV_REFIT_WMIN;
  s.svMinVdet = JV_SV_MIN_VDET;
  s.hardMin = JV_HARD_MIN;
  s.nIter = JV_N_ITER;
  s.iterSig = JV_ITER_SIG;
  s.splitCasc = JV_SPLIT_CASC;
  s.splitFinal = JV_SPLIT_FINAL;
  s.splitDchi2 = JV_SPLIT_DCHI2;
  s.splitFsig = JV_SPLIT_FSIG;
  s.splitOrder = JV_SPLIT_ORDER;
  s.splitSig = JV_SPLIT_SIG;
  s.splitNmax = JV_SPLIT_NMAX;
  s.splitNmin = JV_SPLIT_NMIN;
  s.v0Veto = JV_V0_VETO;
  s.ksWin = JV_KS_WIN;
  s.lamWin = JV_LAM_WIN;
  s.eeWin = JV_EE_WIN;
  s.v0Class = JV_V0_CLASS;
  s.v0Conv = JV_V0_CONV;
  s.v0Mass = JV_V0_MASS;
  s.v0Admit = JV_V0_ADMIT;
  s.v0Chi2Max = JV_V0_CHI2_MAX;
  s.v0ConvRmin = JV_V0_CONV_RMIN;
  s.v0ConvHscale = JV_V0_CONV_HSCALE;
  s.v0DzAbs = JV_V0_DZ_ABS;
  s.v0DzSig = JV_V0_DZ_SIG;
  s.v0HfChi2 = JV_V0_HF_CHI2;
  s.v0HfSig = JV_V0_HF_SIG;
  s.v0HitScale = JV_V0_HIT_SCALE;
  s.v0MaxComp = JV_V0_MAX_COMP;
  s.v0Nu = JV_V0_NU;
  s.v0OutMin = JV_V0_OUT_MIN;
  s.v0PreKs = JV_V0_PRE_KS;
  s.v0PreLam = JV_V0_PRE_LAM;
  s.v0QtMin = JV_V0_QT_MIN;
  s.v0Rmax = JV_V0_RMAX;
  s.v0SmMin = JV_V0_SM_MIN;
  s.v0TwinExt = JV_V0_TWIN_EXT;
  s.v0TwinSig = JV_V0_TWIN_SIG;
  s.v0WinKs = JV_V0_WIN_KS;
  s.v0WinLam = JV_V0_WIN_LAM;
  return s;
}

}  // namespace AlephJointVtx
}  // namespace FCCAnalyses

#endif
