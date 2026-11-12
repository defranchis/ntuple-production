"""Joint PV+SV vertexer for ALEPH, pluggable stages behind a config dict (DEFAULT = the first end-to-end version, row
hybrid_v0; other algorithms and the frozen final are sets of overrides).

S0 inputs      baseline tracks of stage1; optional quality gate; tight-V0 daughters removed; hit-radius rule: a track's
               maximal vertex radius from its VDET/ITC hits (tracks.max_radius_from_hits)
S1 reference   PV that claims no tracks: 'teardown' (beam-constrained fit, worst track removed while chi2 >= 5, = stage1),
               'jetproj_hard' (QFNDIP-like jet-projected constraint), 'jetproj_soft' (its soft correlated form),
               'nominal' (stage1's stored PV)
S2 seeds       pairs of tracks with lifetime-signed 3D IP significance > seed_sig w.r.t. the reference PV ('signed'),
               or of all non-PV tracks ('nonpv', LCFIPlus build-up), or stage1's svn vertices as candidates ('svn');
               pair fits started at the perigee and at the circle crossings, hit-radius veto, V0/conversion mass veto,
               flight significance and forward pointing
S3 growth      a track joins a candidate if chi2(candidate) < cap and (rule 'symmetric') chi2(candidate) < chi2(PV)
               [IVF] or (rule 'nonpv') the track is not in the reference PV [LCFIPlus]; hit-radius and mass veto;
               compatible candidates merged; optional exclusive arbitration (shared track -> lowest chi2)
S4 joint fit   annealed multi-vertex fit of PV + all candidates: soft weights w_kv = exp(-chi2_kv/2T) /
               (exp(-chi2_cut/2T) + sum_v' exp(-chi2_kv'/2T)); beam spot on the PV; lifetime-sign prior on the PV term
               (PV numerator x 1 / (1 + (s+/prior_s)^2), s+ = positive lifetime-signed IP significance); SV refits
               frozen at the last temperature; vertices within merge_sig merged, SVs within merge_sig of the PV dissolved;
               with v0_class, V0 candidates (K0S, Lambda, photon conversions; jv/v0.py) as further two-track
               hypotheses in the same weights
S5 iterate     unassigned tracks (or PV tracks with large positive lifetime sign) with s > seed_sig seed new pairs ->
               S3 -> S4 (n_iter rounds)
Alternatives and additions behind switches: seedless candidate supply instead of S2 + S3 ('supply': 'vd' vertex density,
'da' deterministic annealing, seedless.py); logistic scores of candidates and of final SVs (sl_score_min, sl_final_min);
S4 weight forms and outlier constants; B -> D cascade split of final SVs (split_casc); S6 acceptance window (accept).
Output: PV (x, C), SVs (x, C, mass, per-track weights), per-track PV weight and hard assignment (w > 0.5); with
v0_class also the V0 list (run_event).
Units cm, GeV.
"""
import time
import numpy as np

from . import fit as F
from . import tracks as TR
from . import io as IO
from . import engine as EN
from . import seedless as SL
from . import v0 as V0

M_KS, M_LAM, M_P = V0.M_KS, V0.M_LAM, V0.M_P
M_E = 0.000511

DEFAULT = dict(
    # S0
    gate_min_tpc=0, remove_tight_v0=True, hit_radius=True,
    # S1
    pv_ref='jetproj_soft', pv_window_d0=0.75, pv_window_z0=5.0, pv_chi2max=5.0,
    # S2
    seed_mode='signed', seed_sig=2.5, seed_chi2=9.0, seed_fsig=3.0, seed_crossings=True, v0_veto=True,
    # S2 additions (defaults reproduce hybrid_v0)
    seed_sig2=None, seed_abs=False, seed_min_vdet=0, seed_min_p=0., seed_max_sigd=None, seed_window=True, seed_same_jet=False, seed_sigl=None,
    veto_stage1_v0=False, ks_win=0.015, lam_win=0.006, ee_win=0.03, seed_sort='chi2',
    # S3 additions (defaults reproduce hybrid_v0)
    grow_pre_chi2=None, grow_sym_refit=False, grow_cos=None, grow_min_vdet=0, grow_sig_min=None, s3_merge_sig=None, merge_require_share=False,
    merge_chi2ndf=None, s4_cos=None,
    # kinematic acceptance of few-track vertices (n = 2, 3, >= 4): opening angle, momentum sum, softest track, pointing, sigma along p
    acc2_ang=None, acc2_psum=None, acc2_minp=None, acc2_cos=None, acc2_sigl=None, acc3_cos=None, acc3_sigl=None, acc3_psum=None, seed_kin=False,
    acc2_dis_hi=None, acc3_dis_lo=None, acc4_cos=None, acc4_dis_lo=None,
    # per-jet inclusive candidates from weak tracks
    jc=False, jc_sig=1.0, jc_chi2=5.0, jc_min=3, jc_cos=0.7, jc_rounds=1,
    # SV membership gate: minimum momentum
    sv_min_p=0.,
    # candidate supply: 'seeded' (S2 + S3) | 'vd' (vertex density) | 'da' (annealing with splitting) | unions 'seeded+vd' ... | 'truth'
    supply='seeded',
    vd_k0=1.0, vd_tube=1.0, vd_grid_chi2=16.0, vd_vmin=0.3, vd_R0=2.0, vd_n_climb=8, vd_max_starts=15, vd_fwd_margin=0.005,
    vd_member_chi2=9.0, vd_grid_pvchi2=4.0, vd_grid_cos=0.3,
    da_T0=128.0, da_cool=0.6, da_Tmin=1.0, da_chi2_cut=9.0, da_n_em=3, da_tol=2e-4, da_max_v=8, da_use_rho=False,
    da_prior=False, da_split_margin=1.0, da_merge_margin=0.8, da_min_side=2, da_fwd=True, da_Tpurge=2.0, da_wside=0.1,
    da_wmin=0.3, da_noise_T=4.0,
    sl_chi2=9.0, sl_symmetric=True, sl_forward=True, sl_trim_chi2=9.0, sl_fsig=3.0,
    # logistic scores (MC-fitted, jv/score_coefs.json candidates, jv/score_coefs_final.json final SVs); below the cut: dropped
    sl_score_min=-1e9, sl_score_seeded=False, sl_score_s5=None, sl_final_min=-1e9,
    # coefficient files of the two scores (in jv/), so refitted scores are config switches
    sl_coefs_cand='score_coefs.json', sl_coefs_final='score_coefs_final.json',
    # S4 forward rule cos >= 0 (mapped onto sv_forward = 0); S5 source 'seeded' | 'none'
    s4_forward=False, s5_mode='seeded',
    # beam constraint widths (cm) on every PV fit of S1 and S4 (default = stage1's; bs_sig_ref / s4_bs_sig override per stage);
    # chi2_cut_neg: PV outlier constant for tracks with negative lifetime sign (None: chi2_cut)
    bs_sig=(200e-4, 100e-4, 2.0), chi2_cut_neg=None,
    # S3
    grow_rule='symmetric', grow_chi2=6.0, vtx_chi2ndf=5.0, mass_max=6.0, merge_share=0.5, arbitration=False,
    # S4
    joint=True, temps=(64., 16., 4., 2., 1.), chi2_cut=5.0, assoc_chi2=16.0, lifetime_prior=True, prior_s=3.0, freeze_last=True,
    merge_sig=3.0, min_sv_tracks=2,
    # S4 engine switches (the defaults reproduce the baseline engine)
    chi2_mode_pv='point', chi2_mode_sv='point', loo_eps=1.0, dens_k=0.0, sv_pen=0.0, weight_form='gauss', student_nu=8.0,
    prior_form='lorentz', prior_T=False, prior_update=False, sv_prior_s=0.0, s4_bs_sig=None, freeze_n=1,
    merge_sig_pv=None, merge_sig_sv=None, merge_when='each', sv_forward=None, sv_min_vdet=0, cut_novdet=None, sv_refit_wmin=1e-4,
    pv_final='soft', sv_final='frozen', assoc_rule='any', step_guard=None, iter_update_s3=False, iter_pv_sig=None,
    sv_fwd_w=None, drop_veto=None, anneal_cand=False, sv_Tmax=None, pv_excl_sv=None, iter_sv_chi2=None, chi2_cut_sv=None, pv_Tmin=None, hard_min=0.5, iter_seed_sig=None, iter_fsig=None, accept_final=True, bs_sig_ref=None,
    # S5
    n_iter=1, iter_sig=4.0,
    # S6 candidate acceptance = the svn window (analyzer_svnew.h): displacement 0.03-3 cm, cos pointing >= 0.7,
    # sigma along the momentum <= 0.1 cm, <= 8 tracks, 2-track flight significance > 3, seed pair Delta R <= 0.8
    accept=True, acc_dis_lo=0.03, acc_dis_hi=3.0, acc_cos=0.7, acc_sigl=0.10, acc_max_trk=8, acc_fsig2=3.0, seed_dr=0.8,
    # SV source when no finder runs ('svn_asread' = stage1 svn as stored)
    sv_mode='finder',
    # B->D cascade split of final SVs (two-vertex hypothesis on the hard tracks, exhaustive bipartitions)
    split_casc=False, split_nmin=4, split_nmax=7, split_dchi2=9.0, split_sig=3.0, split_fsig=3.0, split_final=True,
    split_order=0.0, split_pick='chi2',     # split_pick 'valid': best chi2 among bipartitions passing the geometric cuts
    # flight-significance acceptance for every track count (None = off; acc_fsig2 keeps its 2-track role)
    acc_fsig_n=None,
    # stricter final-score cut for short flight: SVs with |x_SV - x_PV| < sl_short_fd (cm) need score >= sl_short_min
    sl_short_fd=None, sl_short_min=None,
    # momentum constant (jv/tracks.py PT_CONSTS): 'legacy' = 0.0045 GeV cm (frozen final) | 'stage1' = 0.0029979 x 1.5 T
    # (as stage1; the V0 configurations); the stored momenta t_p and the momenta of _mass / _kin_ok are rescaled to it
    pt_const='legacy',
    # V0 class (jv/v0.py, model in its docstring; off = frozen final): K0S, Lambda and (v0_conv) photon conversions as
    # two-track vertex hypotheses inside S4. Needs remove_tight_v0=False and v0_veto=False, else stage1's tight V0
    # daughters are removed and V0-like pairs vetoed before the V0 class sees them. Tables: file v0_params (in jv/).
    # Candidates: crossing radius < v0_rmax cm, |dz| < v0_dz_sig sigma or v0_dz_abs cm, mass at the crossing within
    # v0_pre_ks / v0_pre_lam GeV, two-track fit chi2 < v0_chi2_max, mass at the fit within v0_win_ks / v0_win_lam GeV,
    # Armenteros qT > v0_qt_min GeV (K0S, Lambda). Terms: Student-t track term (v0_nu dof); K0S / Lambda daughters with
    # inner hits inside the decay radius: rate eps(p) (table 'hit_eps') and chi2 / v0_hit_scale^2; mass term v0_mass
    # (error floor v0_sm_min GeV); PV pointing v0_point; prior odds 'prior'[v0_prior]. Admission log R >= v0_admit at the
    # reference PV; at most v0_max_comp candidates per connected component; output: V0s with weight > v0_out_min (their
    # tracks are V0-assigned above hard_min)
    v0_class=False, v0_params='v0_params.json', v0_prior=None,
    v0_rmax=180., v0_dz_sig=10., v0_dz_abs=1.0, v0_pre_ks=0.15, v0_pre_lam=0.08, v0_win_ks=0.06, v0_win_lam=0.03,
    v0_chi2_max=50., v0_nu=4.0, v0_hit_scale=3.0, v0_mass=True, v0_point=True,
    v0_admit=-4.6, v0_max_comp=12, v0_sm_min=0.0005, v0_qt_min=(0.015, 0.035), v0_out_min=0.05,
    # twin rule: an SV within v0_twin_sig standard deviations of a V0 candidate's vertex is the pair's own hypothesis
    # (None = off); v0_twin_ext: a number = such an SV holding more than this summed weight on other tracks is no twin
    # (None = no limit), 'frac' = the pair's share of the SV decides; v0_prior_nm = prior table of the weight without
    # the mass term (output v0_w0; None = not computed)
    v0_twin_sig=None, v0_twin_ext=None, v0_prior_nm=None,
    # photon conversions (v0_conv; jv/conv.py loose selection): pointing 'conv_point'[v0_conv_params]; hit term
    # 'hit_table'[v0_hit_table]['conversion'] (None = no hit term); track-term chi2 scales for one / both daughters with
    # inner hits inside the conversion radius v0_conv_hscale; conversion vertex radius from the axis >= v0_conv_rmin cm
    # (None = no bound)
    v0_conv=False, v0_conv_params=None, v0_hit_table=None, v0_conv_hscale=(1.0, 1.0), v0_conv_rmin=None,
    # V0-against-heavy-flavour term 'hfterm'[v0_hf] (None = off; the table also holds the PV-pointing factor
    # pv_wscale it was fitted with): further displaced tracks (IP significance > v0_hf_sig) at the V0 vertex (chi2 <
    # v0_hf_chi2) and in the daughters' jets, decay-length density, |cos theta*| when the table has it
    v0_hf=None, v0_hf_chi2=9.0, v0_hf_sig=3.0,
    # SV-pointing odds: the K0S / Lambda odds use the pointing to the best origin, PV or an SV, instead of the PV alone
    # (table 'pt'[v0_pt_params])
    v0_pt_point=False, v0_pt_params=None,
)

ROWS = {
    'hybrid_v0': {},
    'lcfiplus_like': dict(pv_ref='teardown', seed_mode='nonpv', grow_rule='nonpv', arbitration=True, joint=False, n_iter=0,
                          lifetime_prior=False),
    'mvf_like': dict(pv_ref='nominal', seed_mode='svn', joint=True, lifetime_prior=False, n_iter=0),
    'ivf_like': dict(pv_ref='nominal', seed_mode='signed', grow_rule='symmetric', arbitration=True, joint=False, n_iter=0,
                     lifetime_prior=False),
    'jetproj_pv_svn': dict(pv_ref='jetproj_hard', sv_mode='svn_asread', joint=False, n_iter=0),
    'jetproj_soft_svn': dict(pv_ref='jetproj_soft', sv_mode='svn_asread', joint=False, n_iter=0),
    'nominal': dict(pv_ref='nominal', sv_mode='svn_asread', joint=False, n_iter=0),
}


def config(name_or_overrides=None, **kw):
    c = dict(DEFAULT)
    if isinstance(name_or_overrides, str):
        c.update(ROWS[name_or_overrides])
    elif isinstance(name_or_overrides, dict):
        c.update(name_or_overrides)
    c.update(kw)
    bad = set(c) - set(DEFAULT)
    if bad:
        raise KeyError(f'unknown config keys {sorted(bad)}')
    # one forward rule, sv_forward = minimum cos(track direction at the SV, flight direction); the older switches s4_cos and
    # s4_forward (= cos >= 0) map onto it; the explicit sv_forward wins
    if c['sv_forward'] is None and c['s4_cos'] is not None:
        c['sv_forward'] = c['s4_cos']
    if c['sv_forward'] is None and c['s4_forward']:
        c['sv_forward'] = 0.0
    if c['v0_class']:
        # combinations that cannot run: the V0 class lives in the joint fit S4, and every table it selects must exist
        if c['sv_mode'] != 'finder' or not c['joint']:
            raise ValueError('v0_class needs the joint fit (sv_mode finder, joint true)')
        if (c['v0_prior'] is None or (c['v0_conv'] and c['v0_conv_params'] is None)
                or (c['v0_pt_point'] and c['v0_pt_params'] is None)):
            raise ValueError('v0_class needs v0_prior, v0_conv needs v0_conv_params, v0_pt_point needs v0_pt_params')
        V0.check_tables(c)
    return c


# ---------------------------------------------------------------------------------------------------------------- helpers
def _pre_window(E, idx, d0max=0.75, z0max=5.0):
    """stage1 PV pre-selection window (beam-spot-referenced |d0| <= 0.75 cm, |z0| <= 5 cm) on raw states"""
    p = E['t']['par'][idx]
    bx, by, bz = E['bs']
    D0r, phr, omr, Z0r, tlr = p[:, 0], p[:, 1], p[:, 2], p[:, 3], p[:, 4]
    s_ = bx * np.cos(phr) + by * np.sin(phr)
    d0p = D0r - bx * np.sin(phr) + by * np.cos(phr) - 0.5 * omr * s_ * s_
    z0p = Z0r - bz + tlr * s_
    return (np.abs(d0p) <= d0max) & (np.abs(z0p) <= z0max)


_PSCALE = 1.0      # momentum scale of the selected pt_const (tracks.pt_scale), set per event by run_event


def _mass(fpar, idx, x, masses=None):
    E_ = 0.; P = np.zeros(3)
    for k, t in enumerate(idx):
        m = TR.M_PI if masses is None else masses[k]
        p = TR.momentum_at(fpar[t], F._arc_near(fpar[t], x))
        if _PSCALE != 1.0:
            p = p * _PSCALE
        E_ += np.sqrt(p @ p + m * m); P += p
    return float(np.sqrt(max(E_ * E_ - P @ P, 0.))), P


def _v0_like(fpar, q, i, j, x, c=None):
    if q[i] * q[j] >= 0:
        return False
    ks, lam, ee = (0.015, 0.006, 0.03) if c is None else (c['ks_win'], c['lam_win'], c['ee_win'])
    if abs(_mass(fpar, [i, j], x)[0] - M_KS) < ks:
        return True
    for mm in ((M_P, TR.M_PI), (TR.M_PI, M_P)):
        if abs(_mass(fpar, [i, j], x, mm)[0] - M_LAM) < lam:
            return True
    if _mass(fpar, [i, j], x, (M_E, M_E))[0] < ee:
        return True
    return False


def _dsig2(x1, C1, x2, C2):
    d = x1 - x2
    try:
        return float(d @ np.linalg.solve(C1 + C2, d))
    except np.linalg.LinAlgError:
        return np.inf


def _accept(c, fpar, trk, x, C, pvx, pvC):
    """S6: svn candidate window on a vertex with tracks trk (indices into fpar)"""
    if not c['accept']:
        return True
    d = x - pvx
    dis = float(np.linalg.norm(d))
    if not (c['acc_dis_lo'] <= dis <= c['acc_dis_hi']) or len(trk) > c['acc_max_trk'] or len(trk) < 2:
        return False
    _, P = _mass(fpar, trk, x)
    pm = np.linalg.norm(P)
    if pm <= 0 or d @ P / (dis * pm) < c['acc_cos']:
        return False
    u = P / pm
    var = u @ C @ u
    if not (var > 0) or np.sqrt(var) > c['acc_sigl']:
        return False
    if len(trk) == 2 and _dsig2(x, C, pvx, pvC) <= c['acc_fsig2'] ** 2:
        return False
    if c['acc_fsig_n'] is not None and _dsig2(x, C, pvx, pvC) <= c['acc_fsig_n'] ** 2:
        return False
    return _kin_ok(c, fpar, trk, x, C, d, dis, P, pm, u)


def _kin_ok(c, fpar, trk, x, C, d, dis, P, pm, u):
    n = len(trk)
    cos_pf = d @ P / (dis * pm)
    if n > 3:
        if c['acc4_cos'] is not None and cos_pf < c['acc4_cos']:
            return False
        if c['acc4_dis_lo'] is not None and dis < c['acc4_dis_lo']:
            return False
        return True
    if n == 2:
        if c['acc2_dis_hi'] is not None and dis > c['acc2_dis_hi']:
            return False
        if c['acc2_cos'] is not None and cos_pf < c['acc2_cos']:
            return False
        if c['acc2_sigl'] is not None and u @ C @ u > c['acc2_sigl'] ** 2:
            return False
        if c['acc2_psum'] is not None and pm < c['acc2_psum']:
            return False
        if c['acc2_ang'] is not None or c['acc2_minp'] is not None:
            p1 = TR.momentum_at(fpar[trk[0]], F._arc_near(fpar[trk[0]], x))
            p2 = TR.momentum_at(fpar[trk[1]], F._arc_near(fpar[trk[1]], x))
            if _PSCALE != 1.0:
                p1, p2 = p1 * _PSCALE, p2 * _PSCALE
            n1, n2 = np.linalg.norm(p1), np.linalg.norm(p2)
            if c['acc2_minp'] is not None and min(n1, n2) < c['acc2_minp']:
                return False
            if c['acc2_ang'] is not None and p1 @ p2 < np.cos(c['acc2_ang']) * n1 * n2:
                return False
        return True
    if c['acc3_cos'] is not None and cos_pf < c['acc3_cos']:
        return False
    if c['acc3_dis_lo'] is not None and dis < c['acc3_dis_lo']:
        return False
    if c['acc3_sigl'] is not None and u @ C @ u > c['acc3_sigl'] ** 2:
        return False
    if c['acc3_psum'] is not None and pm < c['acc3_psum']:
        return False
    return True


def _dr(fpar, i, j):
    def eta_phi(p):
        th = np.arctan2(1., p[4])
        return -np.log(np.tan(0.5 * th)), p[1]
    e1, f1 = eta_phi(fpar[i]); e2, f2 = eta_phi(fpar[j])
    df = (f1 - f2 + np.pi) % (2 * np.pi) - np.pi
    return np.hypot(e1 - e2, df)


class Vtx:
    __slots__ = ('x', 'C', 'trk', 'chi2')

    def __init__(self, x, C, trk, chi2=0.):
        self.x, self.C, self.trk, self.chi2 = x, C, list(trk), chi2


# ---------------------------------------------------------------------------------------------------------------- stages
def stage0(E, c):
    t = E['t']
    n = len(t['q'])
    ok = np.isfinite(t['fpar']).all(1) & np.isfinite(t['fcov'].reshape(n, -1)).all(1)
    if c['gate_min_tpc'] > 0:
        ok &= t['ntpc'] >= c['gate_min_tpc']
    tightv0 = np.zeros(n, bool)
    v = E['v']
    for a, b, tg in zip(v['t1'], v['t2'], v['tight']):
        if tg and a >= 0 and b >= 0:
            tightv0[a] = tightv0[b] = True
    if c['remove_tight_v0']:
        ok &= ~tightv0
    rmax = TR.max_radius_from_hits(t['nvdet'], t['nitc']) if c['hit_radius'] else np.full(n, np.inf)
    return ok, rmax, tightv0


def stage1(E, c, ok):
    """reference PV: x, C, list of PV tracks (for the rows without a joint fit; the reference claims nothing)"""
    t = E['t']
    idx = np.nonzero(ok)[0]
    pre = idx[_pre_window(E, idx, c['pv_window_d0'], c['pv_window_z0'])] if len(idx) else idx
    bsr = np.asarray(c['bs_sig'], float) if c['bs_sig_ref'] is None else np.asarray(c['bs_sig_ref'], float)
    if c['pv_ref'] == 'nominal':
        return E['pv'].copy(), E['pvcov'].copy(), np.nonzero(t['mem'] & 1)[0]
    if c['pv_ref'] == 'teardown' or len(pre) < 2:
        x, C, kept, _ = F.teardown(t['fpar'], t['fcov'], pre, E['bs'], bs=E['bs'], bs_sig=bsr, chi2max=c['pv_chi2max'])
        return x, C, kept
    jet_of = t['jet'][pre]
    # teardown = fallback and the start of the soft form: its Gaussian weights vanish for all tracks when started at the
    # beam spot, up to ~1 cm away in z, and the fit then stays there (first 4000 events of mc_ZM4212.41: generated PV z
    # spread 0.74 cm; |dz| > 500 um in 5.9 % of events from the beam-spot start, 1.2 % from the teardown start)
    xt, Ct, kept, _ = F.teardown(t['fpar'], t['fcov'], pre, E['bs'], bs=E['bs'], bs_sig=bsr, chi2max=c['pv_chi2max'])
    try:
        x, C, nused = F.jetproj_pv(t['par'][pre], t['cov'][pre], jet_of, E['jet_theta'], E['jet_phi'], E['jet_p'], E['bs'],
                                   x0=xt, form='hard' if c['pv_ref'] == 'jetproj_hard' else 'soft', bs_sig=bsr)
    except Exception:
        return xt, Ct, kept
    if not np.all(np.isfinite(x)):
        return xt, Ct, kept
    chi2, _ = TR.chi2_to_point(t['fpar'][pre], t['fcov'][pre], x)
    return x, C, pre[chi2 < c['pv_chi2max']]


def stage2(E, c, ok, rmax, pvx, pvC, pvtrk, s3):
    """seed pairs -> list of (chi2, Vtx)"""
    t = E['t']
    fpar, fcov, q = t['fpar'], t['fcov'], t['q']
    if c['seed_mode'] == 'svn':
        out = []
        for k in range(len(E['s']['mass'])):
            tr = E['st']['t'][E['st']['sv'] == k]
            tr = [int(i) for i in tr if i >= 0 and ok[i]]
            if len(tr) >= 2:
                x, C, chi2, _ = F.fit_vertex(fpar, fcov, tr, E['s']['pos'][k], niter=15, tol=1e-6)
                out.append((float(chi2.sum()), Vtx(x, C, tr, float(chi2.sum()))))
        return out
    if c['seed_mode'] == 'truth':
        # diagnostic only (MC): perfect candidates = baseline tracks of each true weak-HF vertex
        out = []
        for h in range(len(E['h']['pos'])):
            tr = np.nonzero((t['hfv'] == h) & ok)[0]
            if len(tr) >= 2:
                x, C, chi2, _ = F.fit_vertex(fpar, fcov, tr, E['h']['pos'][h], niter=15, tol=1e-6)
                if np.all(np.isfinite(x)):
                    out.append((float(chi2.sum()), Vtx(x, C, [int(i) for i in tr], float(chi2.sum()))))
        return out
    if c['seed_mode'] == 'nonpv':
        inpv = np.zeros(len(q), bool); inpv[pvtrk] = True
        cand = np.nonzero(ok & ~inpv)[0]
        return _pairs(E, c, cand, rmax, pvx, pvC)
    ss = np.abs(s3) if c['seed_abs'] else s3
    q = ok.copy()
    if c['sv_min_p'] > 0:
        q &= t['p'] >= c['sv_min_p']
    if c['seed_min_vdet'] > 0:
        q &= t['nvdet'] >= c['seed_min_vdet']
    if c['seed_min_p'] > 0:
        q &= t['p'] >= c['seed_min_p']
    if c['seed_max_sigd'] is not None:
        q &= fcov[:, 0, 0] < c['seed_max_sigd'] ** 2
    lo = c['seed_sig'] if c['seed_sig2'] is None else c['seed_sig2']
    cand = np.nonzero(q & (ss > lo))[0]
    return _pairs(E, c, cand, rmax, pvx, pvC, strong=ss > c['seed_sig'], s=ss)


def _pairs(E, c, cand, rmax, pvx, pvC, strong=None, s=None):
    t = E['t']
    fpar, fcov, q = t['fpar'], t['fcov'], t['q']
    out = []
    vetoset = set()
    if c['veto_stage1_v0']:
        vetoset = {(min(int(a), int(b)), max(int(a), int(b))) for a, b in zip(E['v']['t1'], E['v']['t2']) if a >= 0 and b >= 0}
    for a in range(len(cand)):
        for b in range(a + 1, len(cand)):
            i, j = int(cand[a]), int(cand[b])
            if strong is not None and not (strong[i] or strong[j]):
                continue
            if c['accept'] and _dr(fpar, i, j) > c['seed_dr']:
                continue
            if c['seed_same_jet'] and t['jet'][i] != t['jet'][j]:
                continue
            if (min(i, j), max(i, j)) in vetoset:
                continue
            r = F.pair_fit(fpar, fcov, i, j, rmax[i], rmax[j], use_crossings=c['seed_crossings'])
            if r is None:
                continue
            x, C, chi2 = r
            if not chi2 < c['seed_chi2']:
                continue
            if np.hypot(x[0], x[1]) >= min(rmax[i], rmax[j]):
                continue
            if _dsig2(x, C, pvx, pvC) < c['seed_fsig'] ** 2:
                continue
            m, P = _mass(fpar, [i, j], x)
            if (x - pvx) @ P <= 0:
                continue
            if c['v0_veto'] and _v0_like(fpar, q, i, j, x, c):
                continue
            if c['seed_window'] and not _accept(c, fpar, [i, j], x, C, pvx, pvC):
                continue
            if c['seed_kin']:
                d_ = x - pvx; dis_ = np.linalg.norm(d_); pm_ = np.linalg.norm(P)
                if not _kin_ok(c, fpar, [i, j], x, C, d_, dis_, P, pm_, P / pm_):
                    continue
            if c['seed_sigl'] is not None:
                pm = np.linalg.norm(P)
                u = P / pm
                if u @ C @ u > c['seed_sigl'] ** 2:
                    continue
            key = chi2
            if c['seed_sort'] == 'sig' and s is not None:
                key = -min(s[i], s[j])
            out.append((key, Vtx(x, C, [i, j], chi2)))
    out.sort(key=lambda z: z[0])
    return out


def stage3(E, c, ok, rmax, seeds, pvx, pvtrk, pool=None, pvC=None, s3=None):
    """growth + merging (+ optional exclusive arbitration) -> list of Vtx"""
    t = E['t']
    fpar, fcov = t['fpar'], t['fcov']
    n = len(t['q'])
    if c['seed_mode'] in ('svn', 'truth'):
        cands = [v for _, v in seeds]
    else:
        inpv = np.zeros(n, bool); inpv[pvtrk] = True
        if pool is None:
            pool = np.nonzero(ok & ~inpv)[0] if c['grow_rule'] == 'nonpv' else np.nonzero(ok)[0]
        pool = np.asarray(pool, np.int64)
        if c['sv_min_p'] > 0 and len(pool):
            pool = pool[t['p'][pool] >= c['sv_min_p']]
        if c['grow_min_vdet'] > 0 and len(pool):
            pool = pool[t['nvdet'][pool] >= c['grow_min_vdet']]
        if c['grow_sig_min'] is not None and s3 is not None and len(pool):
            pool = pool[s3[pool] > c['grow_sig_min']]
        chi2pv, _ = TR.chi2_to_point(fpar[pool], fcov[pool], pvx) if len(pool) else (np.zeros(0), None)
        chi2pv_all = np.full(n, np.inf); chi2pv_all[pool] = chi2pv
        gate = c['grow_chi2'] if c['grow_pre_chi2'] is None else c['grow_pre_chi2']
        used = np.zeros(n, bool)
        cands = []
        for chi2s, v in seeds:
            if used[v.trk[0]] and used[v.trk[1]]:
                continue
            # growth
            ch, Lp = TR.chi2_to_point(fpar[pool], fcov[pool], v.x)
            rad = np.hypot(v.x[0], v.x[1])
            okj = (ch < gate) & (rad < rmax[pool]) & ~np.isin(pool, v.trk)
            if c['grow_rule'] == 'symmetric' and not c['grow_sym_refit']:
                okj &= ch < chi2pv
            if c['grow_cos'] is not None:
                okj &= _cos_flight(fpar[pool], Lp, v.x, pvx) > c['grow_cos']
            order = pool[okj][np.argsort(ch[okj])]
            for k in order:
                if c['accept'] and len(v.trk) >= c['acc_max_trk']:
                    break
                trial = v.trk + [int(k)]
                x, C, chi2, _ = F.fit_vertex(fpar, fcov, trial, v.x, niter=15, tol=1e-6)
                ndf = 2 * len(trial) - 3
                if c['grow_sym_refit'] and c['grow_rule'] == 'symmetric' and not chi2[-1] < chi2pv_all[k]:
                    continue
                if chi2[-1] < c['grow_chi2'] and chi2.sum() / ndf < c['vtx_chi2ndf'] and \
                        np.hypot(x[0], x[1]) < rmax[trial].min() and _mass(fpar, trial, x)[0] < c['mass_max']:
                    v = Vtx(x, C, trial, float(chi2.sum()))
            used[v.trk] = True
            cands.append(v)
        if c['jc'] and s3 is not None:
            cands += _jetclust(E, c, ok, rmax, used, pvx, pvC, s3)
    cands = _merge(E, c, cands)
    if c['arbitration']:
        cands = _arbitrate(E, c, cands)
    if pvC is not None:
        cands = [v for v in cands if _accept(c, fpar, v.trk, v.x, v.C, pvx, pvC)]
    return cands


def _jetclust(E, c, ok, rmax, used, pvx, pvC, s3):
    """per stage1 jet: all S0 tracks not used by a grown candidate with lifetime-signed significance > jc_sig pointing into the jet -> one inclusive vertex fit, worst track removed while
    chi2 >= jc_chi2 (tear-down of the jet's displaced tracks, ATLAS SV1 / LCFIPlus vertex refiner style); up to jc_rounds vertices per jet"""
    t = E['t']
    fpar, fcov = t['fpar'], t['fcov']
    axes = IO.jet_axes(E)
    out = []
    u = np.stack([np.cos(fpar[:, 1]), np.sin(fpar[:, 1]), fpar[:, 4]], 1)
    u = u / np.linalg.norm(u, axis=1)[:, None]
    for j in range(len(axes)):
        sel = ok & (t['jet'] == j) & (s3 > c['jc_sig']) & ~used & (u @ axes[j] > c['jc_cos'])
        for _ in range(c['jc_rounds']):
            idx = np.nonzero(sel)[0]
            if len(idx) < c['jc_min']:
                break
            x0 = pvx + 0.3 * axes[j]
            x, C, kept, chi2 = F.teardown(fpar, fcov, idx, x0, chi2max=c['jc_chi2'], nmin=2)
            if not np.all(np.isfinite(x)) or len(kept) < 2 or chi2.max() >= c['jc_chi2']:
                break
            if np.hypot(x[0], x[1]) >= rmax[kept].min():
                break
            m, P = _mass(fpar, kept, x)
            if (x - pvx) @ P <= 0 or m >= c['mass_max'] or _dsig2(x, C, pvx, pvC) < c['seed_fsig'] ** 2:
                break
            out.append(Vtx(x, C, [int(k) for k in kept], float(chi2.sum())))
            sel[kept] = False
    return out


def _merge(E, c, cands):
    t = E['t']
    fpar, fcov = t['fpar'], t['fcov']
    changed = True
    while changed and len(cands) > 1:
        changed = False
        for a in range(len(cands)):
            for b in range(a + 1, len(cands)):
                A, B = cands[a], cands[b]
                sh = len(set(A.trk) & set(B.trk))
                msig = c['merge_sig'] if c['s3_merge_sig'] is None else c['s3_merge_sig']
                close = _dsig2(A.x, A.C, B.x, B.C) < msig ** 2
                shared_big = sh >= c['merge_share'] * min(len(A.trk), len(B.trk))
                if (shared_big and close) if c['merge_require_share'] else (shared_big or close):
                    tr = sorted(set(A.trk) | set(B.trk))
                    x, C, chi2, _ = F.fit_vertex(fpar, fcov, tr, 0.5 * (A.x + B.x), niter=15, tol=1e-6)
                    mchi = c['vtx_chi2ndf'] if c['merge_chi2ndf'] is None else c['merge_chi2ndf']
                    if chi2.sum() / max(2 * len(tr) - 3, 1) < mchi and np.all(np.isfinite(x)):
                        cands[a] = Vtx(x, C, tr, float(chi2.sum()))
                        del cands[b]
                        changed = True
                        break
                    elif sh > 0:
                        # incompatible: the shared tracks stay with the larger candidate
                        small, big = (A, B) if len(A.trk) < len(B.trk) else (B, A)
                        small.trk = [k for k in small.trk if k not in big.trk]
                        if len(small.trk) < 2:
                            cands.remove(small)
                        else:
                            small.x, small.C, ch, _ = F.fit_vertex(fpar, fcov, small.trk, small.x, niter=15, tol=1e-6)
                        changed = True
                        break
            if changed:
                break
    return cands


def _arbitrate(E, c, cands):
    t = E['t']
    fpar, fcov = t['fpar'], t['fcov']
    owner = {}
    for vi, v in enumerate(cands):
        ch, _ = TR.chi2_to_point(fpar[v.trk], fcov[v.trk], v.x)
        for k, cc in zip(v.trk, ch):
            if k not in owner or cc < owner[k][1]:
                owner[k] = (vi, cc)
    out = []
    for vi, v in enumerate(cands):
        tr = [k for k in v.trk if owner[k][0] == vi]
        if len(tr) < 2:
            continue
        if tr != v.trk:
            x, C, chi2, _ = F.fit_vertex(fpar, fcov, tr, v.x, niter=15, tol=1e-6)
            v = Vtx(x, C, tr, float(chi2.sum()))
        out.append(v)
    return out


_SCORE = None


def cand_feats(E, trk, x, C, chi2, pvx, pvC, s3):
    """features of a candidate (tracks trk at x): candidate score and candidate table"""
    t = E['t']
    tr = np.asarray(trk)
    fpar = t['fpar']
    m, P = _mass(fpar, trk, x)
    d = x - pvx; dis = max(float(np.linalg.norm(d)), 1e-9); pm = max(float(np.linalg.norm(P)), 1e-9)
    u = P / pm
    ss = np.sort(s3[tr])
    drm = max(_dr(fpar, int(a), int(b)) for i, a in enumerate(tr) for b in tr[i + 1:])
    n = len(tr)
    ptrel = float(np.linalg.norm(np.cross(P, d / dis)))
    cosf = min(float(d @ P) / (dis * pm), 1.)
    axes = IO.jet_axes(E)
    jc = float(np.max(axes @ (d / dis)))
    jets = t['jet'][tr]
    return dict(n=n, mass=m, dis=dis, cos=cosf, sigl=float(np.sqrt(max(u @ C @ u, 0))), fsig2=_dsig2(x, C, pvx, pvC), chi2=chi2, n2v=int((t['nvdet'][tr] >= 2).sum()),
                n0v=int((t['nvdet'][tr] == 0).sum()), s0=ss[0], s1=ss[1], smax=ss[-1], drmax=drm, p=pm, r=float(np.hypot(x[0], x[1])), pmin=float(t['p'][tr].min()),
                nsig2=int((ss > 2).sum()), nneg=int((ss < -1).sum()), sumpos=float(np.minimum(np.maximum(ss, 0), 10).sum()), ptrel=ptrel, jetcos=jc, njetsame=int(max((jets == 0).sum(), (jets == 1).sum())))


def score_inputs(f):
    """engineered inputs of the logistic score from cand_feats (scalars or arrays)"""
    mc = np.sqrt(f['mass'] ** 2 + f['ptrel'] ** 2) + f['ptrel']
    return dict(n=np.minimum(f['n'], 6), lp=np.log(f['p'] + .1), lc=np.log(1 - np.minimum(f['cos'], 1) + 1e-4), dr=f['drmax'], s0=np.clip(f['s0'], -5, 10), s1=np.clip(f['s1'], -5, 15),
                sm=np.clip(f['smax'], -3, 20), n2v=np.minimum(f['n2v'], 4), n0v=np.minimum(f['n0v'], 3), lf=np.log(f['fsig2'] + 1), ls=np.log(f['sigl'] + 1e-3),
                ndf=f['chi2'] / (2 * f['n'] - 3), lm=np.log(f['mass'] + .05), ld=np.log(f['dis'] + 1e-3), lpm=np.log(f['pmin'] + .05), fsig2n=f['nsig2'] / f['n'], fneg=f['nneg'] / f['n'],
                sumpos=f['sumpos'] / f['n'], lmc=np.log(mc + .05), lptr=np.log(f['ptrel'] + .02), jc=np.log(1 - np.minimum(f['jetcos'], 1) + 1e-4), jsame=f['njetsame'] / f['n'],
                s0s=np.clip(f['s0'], -5, 3) ** 2, drs=f['drmax'] ** 2, lfs=np.log(f['fsig2'] + 1) ** 2, lps=np.log(f['p'] + .1) ** 2, ldlf=np.log(f['dis'] + 1e-3) * np.log(f['fsig2'] + 1),
                n2=(f['n'] == 2) * 1., n2lf=(f['n'] == 2) * np.log(f['fsig2'] + 1), n2s1=(f['n'] == 2) * np.clip(f['s1'], -5, 15))


_SCORES = {}
FINAL_DUMP = None
SPLIT_DUMP = None      # list -> _split_cascades records the best bipartition of each tried SV


def _score(E, trk, x, C, chi2, pvx, pvC, s3, which='score_coefs.json'):
    """logistic candidate score (coefficients: score_coefs.json for candidates, score_coefs_final.json for final SVs)"""
    S = _SCORES.get(which)
    if which.endswith('.pkl'):
        # gradient-boosted trees on the same inputs; returns the log-odds, so cuts read as for the logistic
        if S is None:
            import pickle, os
            S = pickle.load(open(os.path.join(os.path.dirname(__file__), which), 'rb'))
            _SCORES[which] = S
        F_ = score_inputs(cand_feats(E, trk, x, C, chi2, pvx, pvC, s3))
        p = float(S['model'].predict_proba(np.array([[F_[k] for k in S['names']]], float))[0, 1])
        p = min(max(p, 1e-9), 1 - 1e-9)
        return float(np.log(p / (1 - p)))
    if S is None:
        import json, os
        z = json.load(open(os.path.join(os.path.dirname(__file__), which)))
        S = {k: np.array(z[k]) if k != 'names' else z[k] for k in ('names', 'mu', 'sd', 'coef', 'b0')}
        _SCORES[which] = S
    F_ = score_inputs(cand_feats(E, trk, x, C, chi2, pvx, pvC, s3))
    Xv = np.array([F_[k] for k in S['names']], float)
    return float(((Xv - S['mu']) / S['sd']) @ S['coef'] + S['b0'])


def _flies_forward(fpar, L, d):
    ang = fpar[:, 1] + 2 * fpar[:, 2] * L
    return np.cos(ang) * d[0] + np.sin(ang) * d[1] + fpar[:, 4] * d[2] >= 0


def _candidate(E, c, trk, x0, rmax, pvx, pvC, s3):
    """candidate from a seedless maximum: fit the tracks, worst track out until all are within sl_trim_chi2; window as S3"""
    t = E['t']
    fpar, fcov, q = t['fpar'], t['fcov'], t['q']
    trk = [int(k) for k in trk]
    x = x0
    while len(trk) >= 2:
        x, C, chi2, _ = F.fit_vertex(fpar, fcov, trk, x, niter=15, tol=1e-6)
        if not np.all(np.isfinite(x)):
            return None
        j = int(np.argmax(chi2))
        rad = np.hypot(x[0], x[1])
        if len(trk) > 2 and (chi2[j] > c['sl_trim_chi2'] or rad >= rmax[trk].min() or len(trk) > c['acc_max_trk']):
            trk.pop(j)
            continue
        break
    if len(trk) < 2:
        return None
    if chi2.sum() / (2 * len(trk) - 3) >= c['vtx_chi2ndf'] or np.hypot(x[0], x[1]) >= rmax[trk].min():
        return None
    while len(trk) > 2 and _mass(fpar, trk, x)[0] >= c['mass_max']:
        trk.pop(int(np.argmax(chi2)))
        x, C, chi2, _ = F.fit_vertex(fpar, fcov, trk, x, niter=15, tol=1e-6)
    if _dsig2(x, C, pvx, pvC) < c['sl_fsig'] ** 2:
        return None
    m, P = _mass(fpar, trk, x)
    if (x - pvx) @ P <= 0 or m >= c['mass_max']:
        return None
    if len(trk) == 2 and c['v0_veto'] and _v0_like(fpar, q, trk[0], trk[1], x):
        return None
    if not _accept(c, fpar, trk, x, C, pvx, pvC):
        return None
    if c['sl_score_min'] > -1e8 and _score(E, trk, x, C, float(chi2.sum()), pvx, pvC, s3, c['sl_coefs_cand']) < c['sl_score_min']:
        return None
    return Vtx(x, C, trk, float(chi2.sum()))


def _tracks_for(c, fp, fc, x, sel_chi2, pvchi2, rm):
    """tracks a seedless vertex at x may claim (local indices)"""
    ch, L = TR.chi2_to_point(fp, fc, x)
    okt = (ch < sel_chi2) & (np.hypot(x[0], x[1]) < rm)
    if c['sl_symmetric']:
        okt &= ch < pvchi2
    return okt, L


def supply_vd(E, c, ok, rmax, pvx, pvC, s3):
    t = E['t']
    idx = np.nonzero(ok)[0]
    if len(idx) < 2:
        return []
    fp, fc = t['fpar'][idx], t['fcov'][idx]
    Ci = np.linalg.inv(pvC)
    mx, mv, dip, mc, ncl = SL.vd_maxima(fp, fc, rmax[idx], pvx, Ci, c['vd_k0'], c['vd_tube'], c['vd_grid_chi2'], c['vd_vmin'],
                                        c['vd_n_climb'], c['vd_max_starts'], c['vd_fwd_margin'], c['vd_member_chi2'], c['vd_grid_pvchi2'], c['vd_grid_cos'])
    pvchi2, _ = TR.chi2_to_point(fp, fc, pvx)
    out = []
    for m in np.argsort(-mv):
        if dip[m] >= c['vd_R0']:
            continue
        x = mx[m]
        okt, L = _tracks_for(c, fp, fc, x, c['sl_chi2'], pvchi2, rmax[idx])
        if c['sl_forward']:
            okt &= _flies_forward(fp, L, x - pvx)
        if okt.sum() < 2:
            continue
        v = _candidate(E, c, idx[okt], x, rmax, pvx, pvC, s3)
        if v is not None:
            out.append(v)
    return out


def supply_da(E, c, ok, rmax, pvx, pvC, s3, out_state=None):
    t = E['t']
    idx = np.nonzero(ok)[0]
    if len(idx) < 2:
        return []
    fp, fc = t['fpar'][idx], t['fcov'][idx]
    prior = np.ones(len(idx))
    if c['da_prior']:
        sp = np.maximum(s3[idx], 0.)
        prior = 1. / (1. + (sp / c['prior_s']) ** 2)
    Binv = np.diag(1. / np.asarray(c['bs_sig']) ** 2)
    x, C, w, ch, nsplit, nlin = SL.da_run(fp, fc, rmax[idx], pvx, Binv, E['bs'].astype(np.float64), c['da_T0'], c['da_cool'], c['da_Tmin'],
                                          c['da_chi2_cut'], c['da_n_em'], c['da_tol'], c['da_max_v'], c['da_use_rho'], prior,
                                          c['da_split_margin'], c['da_merge_margin'], c['da_min_side'], c['da_fwd'], c['da_Tpurge'],
                                          c['da_wside'], c['da_noise_T'])
    if out_state is not None:
        out_state.update(x=x, C=C, w=w, ch=ch, idx=idx, nsplit=nsplit, nlin=nlin)
    out = []
    for v in range(1, len(x)):
        okt = (w[v] > c['da_wmin']) & (ch[v] < c['sl_chi2'])
        if okt.sum() < 2:
            continue
        cand = _candidate(E, c, idx[okt], x[v], rmax, x[0], C[0], s3)
        if cand is not None:
            out.append(cand)
    return out


def supply_truth(E, c, ok, rmax, pvx, pvC):
    """oracle supply (MC only, ceiling study): the S0 tracks of each true weak-HF vertex as one candidate, fitted, no selection"""
    t = E['t']
    out = []
    for h in range(len(E['h']['isb'])):
        trk = np.nonzero(ok & (t['hfv'] == h))[0]
        if len(trk) < 2:
            continue
        x, C, chi2, _ = F.fit_vertex(t['fpar'], t['fcov'], trk, E['h']['pos'][h], niter=15, tol=1e-6)
        if np.all(np.isfinite(x)):
            out.append(Vtx(x, C, [int(k) for k in trk], float(chi2.sum())))
    return out


def supply(E, c, ok, rmax, pvx, pvC, pvtrk, s3):
    """S2 + S3 or their seedless replacements; a union goes through the S3 merging"""
    parts = c['supply'].split('+') if c['supply'] != 'all' else ['seeded', 'vd', 'da']
    cands = []
    for p in parts:
        if p == 'seeded':
            sd = stage3(E, c, ok, rmax, stage2(E, c, ok, rmax, pvx, pvC, pvtrk, s3), pvx, pvtrk, pvC=pvC, s3=s3)
            if c['sl_score_seeded'] and c['sl_score_min'] > -1e8:
                sd = [v for v in sd if _score(E, v.trk, v.x, v.C, v.chi2, pvx, pvC, s3, c['sl_coefs_cand']) >= c['sl_score_min']]
            cands += sd
        elif p == 'vd':
            cands += supply_vd(E, c, ok, rmax, pvx, pvC, s3)
        elif p == 'da':
            cands += supply_da(E, c, ok, rmax, pvx, pvC, s3)
        elif p == 'truth':
            cands += supply_truth(E, c, ok, rmax, pvx, pvC)
        else:
            raise KeyError(p)
    if parts[0] == 'truth':
        return cands
    if len(parts) > 1 or parts[0] != 'seeded':
        cands = _merge(E, c, cands)
        cands = [v for v in cands if _accept(c, t_fpar(E), v.trk, v.x, v.C, pvx, pvC)]
    return cands


def t_fpar(E):
    return E['t']['fpar']


MODE = {'point': 0, 'vcov': 1, 'loo': 2}
# in-process memo of S0-S3 per event, for tuning several S4 variants on one unit (None = off; CPU rows must run with it off)
UPCACHE = None
UP_KEYS = ('gate_min_tpc', 'remove_tight_v0', 'hit_radius', 'pv_ref', 'pv_window_d0', 'pv_window_z0', 'pv_chi2max', 'seed_mode', 'seed_sig',
           'seed_chi2', 'seed_fsig', 'seed_crossings', 'v0_veto', 'grow_rule', 'grow_chi2', 'vtx_chi2ndf', 'mass_max', 'merge_share',
           'arbitration', 'merge_sig', 'bs_sig_ref', 'accept', 'acc_dis_lo', 'acc_dis_hi', 'acc_cos', 'acc_sigl', 'acc_max_trk', 'acc_fsig2', 'seed_dr')


S4_ONLY = frozenset(('merge_sig_sv', 'merge_sig_pv', 'hard_min', 'n_iter', 'iter_sv_chi2', 'iter_sig', 'chi2_mode_pv', 'chi2_mode_sv',
                     'student_nu', 'chi2_cut', 'chi2_cut_sv', 'temps', 'freeze_n', 'freeze_last', 'pv_excl_sv', 'pv_final', 'sv_final',
                     'sv_forward', 'prior_T', 'weight_form', 'min_sv_tracks', 'assoc_chi2', 'sv_refit_wmin', 'loo_eps', 'pv_Tmin',
                     'chi2_cut_neg', 'merge_when', 'sv_Tmax', 'dens_k', 'sv_pen', 'sv_fwd_w', 'sv_min_vdet', 'sl_final_min', 'accept_final', 'drop_veto',
                     'split_casc', 'split_nmin', 'split_nmax', 'split_dchi2', 'split_sig', 'split_fsig', 'split_final', 'split_order',
                     'split_pick', 'sl_short_fd', 'sl_short_min'))


def _prior(c, s3f, T):
    if not c['lifetime_prior']:
        return np.ones(len(s3f))
    sp = np.maximum(s3f, 0.) / c['prior_s']
    f = c['prior_form']
    pr = 1. / (1. + sp ** 2) if f == 'lorentz' else (1. / (1. + sp ** 4) if f == 'lorentz4' else np.exp(-0.5 * sp ** 2))
    return pr ** (1. / T) if c['prior_T'] else pr


def _f(c, chi2, T):
    """soft-assignment function of the compatibility at temperature T"""
    if c['weight_form'] == 'gauss':
        return np.exp(-0.5 * np.minimum(chi2, 1e4) / T)
    nu = c['student_nu']
    return (1. + np.minimum(chi2, 1e6) / (nu * T)) ** (-0.5 * (nu + 2.))


def _weight_matrix(c, fpar, fcov, V, Wv, s3f, T, rmaxf, novdet, nvdet, extra=None, v0=None):
    """soft weights: V list of Vtx (0 = PV), Wv weights of the last fit per vertex (for the leave-one-out); v0 = V0Set
    (V0 class, jv/v0.py) or None. Returns W (V0 class on: non-V0 weights times 1 - V0 share; marginals in v0.P)"""
    nt = len(fpar)
    num = np.zeros((len(V), nt))
    modes = (MODE[c['chi2_mode_pv']], MODE[c['chi2_mode_sv']])
    pv = V[0].x
    for vi, v in enumerate(V):
        d = v.x - pv
        dn = np.linalg.norm(d)
        fdir = d / dn if dn > 0 else np.array([1., 0., 0.])
        mode = modes[0 if vi == 0 else 1]
        H = np.linalg.inv(v.C) if mode else np.zeros((3, 3))
        ch, dens, cosf, _ = EN.chi2_vertex(fpar, fcov, v.x, H, Wv[vi], mode, c['loo_eps'], fdir)
        eff = ch + c['dens_k'] * dens if mode else ch
        if vi == 0:
            num[0] = _f(c, eff, T)
        else:
            Tsv = T if c['sv_Tmax'] is None else min(T, c['sv_Tmax'])
            nn = _f(c, eff + c['sv_pen'], Tsv)
            allow = ch < c['assoc_chi2']
            allow &= np.hypot(v.x[0], v.x[1]) < rmaxf
            if c['sv_forward'] is not None and c['sv_fwd_w'] is None:
                allow &= cosf >= c['sv_forward']
            if c['sv_min_vdet'] > 0:
                allow &= nvdet >= c['sv_min_vdet']
            if extra is not None and extra[1] is not None:
                allow &= extra[1]       # sv_min_p
            if c['assoc_rule'] == 'cand' or (c['anneal_cand'] and T > c['temps'][-1]):
                am = np.zeros(nt, bool); am[v.trk] = True
                allow &= am
            nn[~allow] = 0.
            if c['sv_forward'] is not None and c['sv_fwd_w'] is not None:
                # soft forward rule, logistic in cos(track direction at the SV, flight direction)
                nn /= 1. + np.exp(-np.clip((cosf - c['sv_forward']) / c['sv_fwd_w'], -50., 50.))
            if c['sv_prior_s'] > 0:
                sp = np.maximum(s3f, 0.)
                nn *= sp ** 2 / (sp ** 2 + c['sv_prior_s'] ** 2)
            num[vi] = nn
    num[0] *= _prior(c, s3f, T)
    if extra is not None and len(extra) > 2 and extra[2] is not None:
        num[0][extra[2]] = 0.       # drop_veto: tracks of SVs rejected at the final selection do not return to the PV
    cutv = _f(c, np.full(nt, c['chi2_cut'], dtype=np.float64), T)
    if c['chi2_cut_neg'] is not None and extra is not None:
        # tracks with negative lifetime sign cannot come from a forward decay: their own PV outlier constant
        cutv[extra[0]] = _f(c, np.full(int(extra[0].sum()), c['chi2_cut_neg'], dtype=np.float64), T)
    if c['cut_novdet'] is not None:
        cutv[novdet] = _f(c, np.full(novdet.sum(), c['cut_novdet'], dtype=np.float64), T)
    sumn = num.sum(0)
    W = num / (cutv + sumn)
    if c['chi2_cut_sv'] is not None:
        W[1:] = num[1:] / (_f(c, np.full(nt, c['chi2_cut_sv'], dtype=np.float64), T) + sumn)
    if v0 is not None and v0.n:
        W = v0.apply(W, cutv + sumn, T, num, V)
    return W


def stage4(E, c, ok, rmax, pvx, pvC, cands, s3, fit_idx, axes=None, jet_of=None):
    """annealed joint fit. Returns pv (x, C), list of Vtx, weight matrix W (1 + nsv, len(fit_idx))"""
    t = E['t']
    fpar, fcov = t['fpar'][fit_idx], t['fcov'][fit_idx]
    nvdet = t['nvdet'][fit_idx]
    novdet = nvdet == 0
    nt = len(fit_idx)
    rmaxf = rmax[fit_idx]
    s3f = s3[fit_idx].copy()
    extra = (s3f < 0, (t['p'][fit_idx] >= c['sv_min_p']) if c['sv_min_p'] > 0 else None,
             E.get('_veto') if c['drop_veto'] is not None else None)
    V = [Vtx(pvx.copy(), pvC.copy(), [])] + [Vtx(v.x.copy(), v.C.copy(), v.trk) for v in cands]
    v0s = E.get('_v0') if c['v0_class'] else None
    bs_sig = np.asarray(c['bs_sig'], float) if c['s4_bs_sig'] is None else np.asarray(c['s4_bs_sig'], float)
    Binv = np.diag(1. / bs_sig ** 2)
    msig_pv = c['merge_sig'] if c['merge_sig_pv'] is None else c['merge_sig_pv']
    msig_sv = c['merge_sig'] if c['merge_sig_sv'] is None else c['merge_sig_sv']
    Wv = [np.zeros(nt)]
    for v in V[1:]:
        w0 = np.zeros(nt); w0[v.trk] = 1.
        Wv.append(w0)
    W = None
    temps = list(c['temps'])
    nT = len(temps)
    PVW = [None]      # weights of the last PV refit

    def guard(xn, xo):
        return np.all(np.isfinite(xn)) and (c['step_guard'] is None or np.linalg.norm(xn - xo) < c['step_guard'])

    def refit_pv(T=None):
        Wp = W
        if c['pv_Tmin'] is not None and T is not None and T < c['pv_Tmin']:
            Wp = _weight_matrix(c, fpar, fcov, V, Wv, s3f, c['pv_Tmin'], rmaxf, novdet, nvdet, extra, v0=v0s)
        if c['pv_excl_sv'] is not None and len(V) > 1:
            # tracks an SV holds at the current temperature do not enter the softened PV term
            Wp = Wp.copy()
            Wp[0][W[1:].max(0) > c['pv_excl_sv']] = 0.
        sel = np.nonzero(Wp[0] > 1e-4)[0]
        x, C, _, _ = F.fit_vertex_nb(fpar, fcov, sel.astype(np.int64), Wp[0][sel], V[0].x, Binv, E['bs'], np.zeros(nt), 15, 1e-6)
        if guard(x, V[0].x):
            V[0].x, V[0].C = x, C
            Wv[0] = Wp[0].copy(); Wv[0][Wp[0] <= 1e-4] = 0.
        PVW[0] = Wv[0]

    for it, T in enumerate(temps + [temps[-1]]):
        frozen = it >= nT - (c['freeze_n'] if c['freeze_last'] else 0)
        if c['prior_update'] and axes is not None and it > 0:
            _, s3n, _, _, _ = TR.ip_info(fpar, fcov, V[0].x, V[0].C, axes, jet_of[fit_idx])
            s3f = s3n
        W = _weight_matrix(c, fpar, fcov, V, Wv, s3f, T, rmaxf, novdet, nvdet, extra, v0=v0s)
        if it == nT:
            break
        refit_pv(T)
        if not frozen:
            for vi in range(1, len(V)):
                sel = np.nonzero(W[vi] > c['sv_refit_wmin'])[0]
                if W[vi][sel].sum() < 1.5 or len(sel) < 2:
                    continue
                x, C, _, _ = F.fit_vertex_nb(fpar, fcov, sel.astype(np.int64), W[vi][sel], V[vi].x, np.zeros((3, 3)),
                                             np.zeros(3), np.zeros(nt), 15, 1e-6)
                if guard(x, V[vi].x):
                    V[vi].x, V[vi].C = x, C
                    w = np.zeros(nt); w[sel] = W[vi][sel]
                    Wv[vi] = w
        # merging: SV-SV and SV-PV
        if c['merge_when'] == 'each' or it >= nT - 2:
            keep = [True] * len(V)
            for a in range(1, len(V)):
                if not keep[a]:
                    continue
                if _dsig2(V[a].x, V[a].C, V[0].x, V[0].C) < msig_pv ** 2:
                    keep[a] = False
                    continue
                for b in range(a + 1, len(V)):
                    if keep[b] and _dsig2(V[a].x, V[a].C, V[b].x, V[b].C) < msig_sv ** 2:
                        wa, wb = W[a].sum(), W[b].sum()
                        V[a].x = (wa * V[a].x + wb * V[b].x) / max(wa + wb, 1e-9)
                        if c['assoc_rule'] == 'cand' or c['anneal_cand']:
                            V[a].trk = sorted(set(V[a].trk) | set(V[b].trk))
                        keep[b] = False
            V = [v for v, k in zip(V, keep) if k]
            Wv = [w for w, k in zip(Wv, keep) if k]
    # drop SVs with fewer than min_sv_tracks hard tracks, then final weights
    for _ in range(3):
        hard = np.argmax(W, 0)
        hard[W.max(0) <= c['hard_min']] = -1
        good = [0] + [vi for vi in range(1, len(V)) if (hard == vi).sum() >= c['min_sv_tracks'] and
                      (not c['accept_final'] or _accept(c, fpar, np.nonzero(hard == vi)[0], V[vi].x, V[vi].C, V[0].x, V[0].C)) and
                      _final_ok(E, c, fit_idx[hard == vi], V[vi], V[0], s3)]
        if len(good) == len(V):
            break
        if c['drop_veto'] is not None:
            dm = (hard >= 1) & ~np.isin(hard, good) & (s3f > c['drop_veto'])
            extra = (extra[0], extra[1], dm if extra[2] is None else (extra[2] | dm))
        V = [V[g] for g in good]
        Wv = [Wv[g] for g in good]
        W = _weight_matrix(c, fpar, fcov, V, Wv, s3f, temps[-1], rmaxf, novdet, nvdet, extra, v0=v0s)
        refit_pv(temps[-1])
        W = _weight_matrix(c, fpar, fcov, V, Wv, s3f, temps[-1], rmaxf, novdet, nvdet, extra, v0=v0s)
    if v0s is not None:
        # V0 marginals of the final weights (the last call may have been the PV refit at pv_Tmin)
        W = _weight_matrix(c, fpar, fcov, V, Wv, s3f, temps[-1], rmaxf, novdet, nvdet, extra, v0=v0s)
        if c['v0_prior_nm'] is not None and v0s.n:
            v0s.weights_nomass()
    hard = np.argmax(W, 0)
    hard[W.max(0) <= c['hard_min']] = -1
    if c['sv_final'] == 'hard':
        # positions and covariances from the hard-assigned tracks (assignment unchanged)
        for vi in range(1, len(V)):
            sel = np.nonzero(hard == vi)[0]
            if len(sel) >= 2:
                x, C, _, _ = F.fit_vertex_nb(fpar, fcov, sel.astype(np.int64), np.ones(len(sel)), V[vi].x, np.zeros((3, 3)),
                                             np.zeros(3), np.zeros(nt), 15, 1e-6)
                if guard(x, V[vi].x):
                    V[vi].x, V[vi].C = x, C
    if c['pv_final'] == 'hard':
        sel = np.nonzero(hard == 0)[0]
        if len(sel) >= 1:
            x, C, _, _ = F.fit_vertex_nb(fpar, fcov, sel.astype(np.int64), np.ones(len(sel)), V[0].x, Binv, E['bs'], np.zeros(nt), 15, 1e-6)
            if guard(x, V[0].x):
                V[0].x, V[0].C = x, C
    for vi in range(1, len(V)):
        V[vi].trk = [int(k) for k in np.nonzero(hard == vi)[0]]
    E['_pvw'] = PVW[0]
    E['_veto'] = extra[2]       # carried into the S5 call of stage4 (same fit_idx)
    return V[0].x, V[0].C, V[1:], W


def _final_ok(E, c, trk, v, pv, s3):
    if c['sl_final_min'] < -1e8:
        return True
    tr = [int(k) for k in trk]
    ch, _ = TR.chi2_to_point(E['t']['fpar'][tr], E['t']['fcov'][tr], v.x)
    if FINAL_DUMP is not None:
        if len(tr) >= 2:
            FINAL_DUMP.append((tr, cand_feats(E, tr, v.x, v.C, float(ch.sum()), pv.x, pv.C, s3)))
        return True
    thr = c['sl_final_min']
    if c['sl_short_fd'] is not None and np.linalg.norm(v.x - pv.x) < c['sl_short_fd']:
        thr = max(thr, c['sl_short_min'])
    return len(tr) >= 2 and _score(E, tr, v.x, v.C, float(ch.sum()), pv.x, pv.C, s3, c['sl_coefs_final']) >= thr


def _split_cascades(E, c, fit_idx, V, hard, px, pC, s3):
    """try to split each final SV with split_nmin..split_nmax hard tracks into two vertices (B -> D cascade).
    Best bipartition (both halves >= 2 tracks) by chi2_A + chi2_B; accepted if it lowers the one-vertex chi2 by > split_dchi2,
    both halves fit with chi2/ndf < vtx_chi2ndf, the halves are > split_sig apart (combined covariance), both are > split_fsig
    from the PV, the far half lies downstream of the near one ((x_far - x_near) . (x_near - PV) >= split_order * norms), the
    hit-radius rule holds for each half and (split_final) both halves pass the final score. Returns new V and hard (local)."""
    from itertools import combinations
    t = E['t']
    fpar, fcov = t['fpar'], t['fcov']
    rmax = E.get('_rmax')
    V = list(V)
    hard = hard.copy()
    nv0 = len(V)
    for k in range(nv0):
        loc = np.nonzero(hard == k + 1)[0]
        n = len(loc)
        if n < c['split_nmin'] or n > c['split_nmax']:
            continue
        tr = fit_idx[loc]
        xa, Ca, cha, _ = F.fit_vertex(fpar, fcov, tr, V[k].x, niter=15, tol=1e-6)
        if not np.all(np.isfinite(xa)):
            continue
        chall = float(cha.sum())
        best = None
        rest = list(range(1, n))
        for na in range(1, n - 1):        # subset A always holds local track 0 -> each bipartition once
            if na + 1 < 2 or n - na - 1 < 2:
                continue
            for comb in combinations(rest, na):
                ia = [0] + list(comb)
                ib = [i for i in range(n) if i not in ia]
                x1, C1, ch1, _ = F.fit_vertex(fpar, fcov, tr[ia], V[k].x, niter=15, tol=1e-6)
                x2, C2, ch2, _ = F.fit_vertex(fpar, fcov, tr[ib], V[k].x, niter=15, tol=1e-6)
                if not (np.all(np.isfinite(x1)) and np.all(np.isfinite(x2))):
                    continue
                tot = float(ch1.sum() + ch2.sum())
                if c['split_pick'] == 'valid':
                    s1_, s2_ = float(ch1.sum()), float(ch2.sum())
                    if s1_ / max(2 * len(ia) - 3, 1) >= c['vtx_chi2ndf'] or s2_ / max(2 * len(ib) - 3, 1) >= c['vtx_chi2ndf']:
                        continue
                    if _dsig2(x1, C1, x2, C2) <= c['split_sig'] ** 2 or _dsig2(x1, C1, px, pC) <= c['split_fsig'] ** 2 or \
                            _dsig2(x2, C2, px, pC) <= c['split_fsig'] ** 2:
                        continue
                    e1_, e2_ = np.linalg.norm(x1 - px), np.linalg.norm(x2 - px)
                    xn_, xf_ = (x1, x2) if e1_ <= e2_ else (x2, x1)
                    if (xf_ - xn_) @ (xn_ - px) < c['split_order'] * np.linalg.norm(xf_ - xn_) * np.linalg.norm(xn_ - px):
                        continue
                if best is None or tot < best[0]:
                    best = (tot, ia, ib, x1, C1, float(ch1.sum()), x2, C2, float(ch2.sum()))
        if best is None:
            continue
        tot, ia, ib, x1, C1, c1, x2, C2, c2 = best
        if SPLIT_DUMP is not None:
            d1_, d2_ = np.linalg.norm(x1 - px), np.linalg.norm(x2 - px)
            xn_, xf_ = (x1, x2) if d1_ <= d2_ else (x2, x1)
            dn_, dd_ = xn_ - px, xf_ - xn_
            SPLIT_DUMP.append(dict(trk=tr.tolist(), ia=list(ia), ib=list(ib), dchi2=chall - tot, chall=chall, n=n,
                                   ndf1=c1 / max(2 * len(ia) - 3, 1), ndf2=c2 / max(2 * len(ib) - 3, 1), sep=_dsig2(x1, C1, x2, C2),
                                   fs1=_dsig2(x1, C1, px, pC), fs2=_dsig2(x2, C2, px, pC),
                                   order=float(dd_ @ dn_ / max(np.linalg.norm(dd_) * np.linalg.norm(dn_), 1e-12)),
                                   m1=_mass(fpar, tr[ia], x1)[0], m2=_mass(fpar, tr[ib], x2)[0], d1=d1_, d2=d2_, accepted=False))
        if chall - tot <= c['split_dchi2']:
            continue
        if c1 / max(2 * len(ia) - 3, 1) >= c['vtx_chi2ndf'] or c2 / max(2 * len(ib) - 3, 1) >= c['vtx_chi2ndf']:
            continue
        if _dsig2(x1, C1, x2, C2) <= c['split_sig'] ** 2:
            continue
        if _dsig2(x1, C1, px, pC) <= c['split_fsig'] ** 2 or _dsig2(x2, C2, px, pC) <= c['split_fsig'] ** 2:
            continue
        d1, d2 = np.linalg.norm(x1 - px), np.linalg.norm(x2 - px)
        (xn, xf) = (x1, x2) if d1 <= d2 else (x2, x1)
        dn, dd = xn - px, xf - xn
        if dd @ dn < c['split_order'] * np.linalg.norm(dd) * np.linalg.norm(dn):
            continue
        if rmax is not None and (np.hypot(x1[0], x1[1]) >= rmax[tr[ia]].min() or np.hypot(x2[0], x2[1]) >= rmax[tr[ib]].min()):
            continue
        ok_ = True
        for (xx, CC, ii) in ((x1, C1, ia), (x2, C2, ib)):
            m, P = _mass(fpar, tr[ii], xx)
            if (xx - px) @ P <= 0:
                ok_ = False
            if ok_ and c['split_final'] and c['sl_final_min'] > -1e8:
                ok_ = _final_ok(E, c, tr[ii], Vtx(xx, CC, []), Vtx(px, pC, []), s3)
        if not ok_:
            continue
        V[k] = Vtx(x1, C1, [int(i) for i in loc[ia]], c1)
        V.append(Vtx(x2, C2, [int(i) for i in loc[ib]], c2))
        if SPLIT_DUMP is not None:
            SPLIT_DUMP[-1]['accepted'] = True
        hard[loc[ib]] = len(V)
    return V, hard


# ---------------------------------------------------------------------------------------------------------------- driver
def run_event(E, c):
    """returns dict: pv (3), pvC (3x3), svs: list of (x, C, mass, ntrk_hard), wpv (n tracks), assign (n: -1 none,
    0 PV, k>=1 SV k), wass (weight of the assigned vertex), cpu (s, process time), s3 (lifetime-signed IP significance
    w.r.t. the output PV); with v0_class also v0s: list of (track a (positive), track b, booked PDG code, hypothesis,
    x, C, mass, final weight, momentum, weight without the mass term or -1) per V0 with weight > v0_out_min, tv0 (n:
    index in v0s of the track's V0 when its weight > hard_min, -1), and with v0_pt_point v0org = (most probable origin
    per output V0: 0 PV, k = SV k of the last S4, -1 none; its posterior; per output SV the summed V0 weight times
    origin posterior on it), at the last temperature of the last S4 (SVs added by the cascade split get 0)"""
    global _PSCALE
    t0 = time.process_time()
    E['_veto'] = None
    E['_v0'] = None
    t = E['t']
    n = len(t['q'])
    _PSCALE = TR.pt_scale(c['pt_const'])
    p_in = E.setdefault('_p_input', t['p'])      # stored momenta (harvest t_p, legacy constant), kept for a repeated call
    t['p'] = p_in * _PSCALE if _PSCALE != 1.0 else p_in
    axes = IO.jet_axes(E)
    jet_of = t['jet'].astype(np.int64)
    upk = None
    if UPCACHE is not None and c['sv_mode'] == 'finder':
        # the key holds every config key (UP_KEYS, the S0-S3 keys of the first version, misses later supply knobs) except
        # the keys read only by stage4 / the S5 loop (S4_ONLY), so S4-only variants reuse S0-S3
        upk = (E['row'], E.get('stem', ''), tuple(str(c[k]) for k in sorted(c) if k not in S4_ONLY))
    if upk is not None and upk in UPCACHE:
        ok, rmax, pvx, pvC, pvtrk, s3, cc = UPCACHE[upk]
        cands = [Vtx(v.x.copy(), v.C.copy(), list(v.trk), v.chi2) for v in cc]
    else:
        ok, rmax, tightv0 = stage0(E, c)
        pvx, pvC, pvtrk = stage1(E, c, ok)
        _, s3, _, _, _ = TR.ip_info(t['fpar'], t['fcov'], pvx, pvC, axes, jet_of)
        cands = None
    wpv = np.zeros(n); assign = np.full(n, -1, np.int64); wass = np.zeros(n)
    svs = []
    if c['sv_mode'] == 'svn_asread':
        wpv[pvtrk] = 1.; assign[pvtrk] = 0; wass[pvtrk] = 1.
        for k in range(len(E['s']['mass'])):
            tr = E['st']['t'][E['st']['sv'] == k]
            tr = tr[tr >= 0]
            assign[tr] = k + 1; wass[tr] = 1.
            svs.append((E['s']['pos'][k], IO.cov6_to_3(E['s']['cov6'][k]), float(E['s']['mass'][k]), len(tr)))
        return dict(pv=pvx, pvC=pvC, svs=svs, wpv=wpv, assign=assign, wass=wass, cpu=time.process_time() - t0, s3=s3)
    if cands is None:
        cands = supply(E, c, ok, rmax, pvx, pvC, pvtrk, s3)
        if upk is not None:
            UPCACHE[upk] = (ok, rmax, pvx, pvC, pvtrk, s3, [Vtx(v.x.copy(), v.C.copy(), list(v.trk), v.chi2) for v in cands])
    if not c['joint']:
        wpv[pvtrk] = 1.; assign[pvtrk] = 0; wass[pvtrk] = 1.
        for k, v in enumerate(cands):
            assign[v.trk] = k + 1; wass[v.trk] = 1.
            svs.append((v.x, v.C, _mass(t['fpar'], v.trk, v.x)[0], len(v.trk)))
        return dict(pv=pvx, pvC=pvC, svs=svs, wpv=wpv, assign=assign, wass=wass, cpu=time.process_time() - t0, s3=s3)
    # tracks entering the joint fit: S0 tracks in the PV window or in a candidate
    inc = np.zeros(n, bool)
    idx0 = np.nonzero(ok)[0]
    if len(idx0):
        inc[idx0[_pre_window(E, idx0, c['pv_window_d0'], c['pv_window_z0'])]] = True
    for v in cands:
        inc[v.trk] = True
    v0set = None
    if c['v0_class']:
        v0c = V0.candidates(E, c, ok, rmax, pvx, pvC, _v0_ref_denominators(c, t, ok, pvx, s3),
                            0.5 * TR.PT_CONSTS[c['pt_const']], s3=s3)
        for k in v0c:
            inc[k['a']] = inc[k['b']] = True
        v0set = V0.V0Set(c, v0c)
    fit_idx = np.nonzero(inc)[0]
    loc = np.full(n, -1); loc[fit_idx] = np.arange(len(fit_idx))
    for v in cands:
        v.trk = [int(loc[k]) for k in v.trk]
    if v0set is not None:
        v0set.localize(loc)
        E['_v0'] = v0set
    px, pC, V, W = stage4(E, c, ok, rmax, pvx, pvC, cands, s3, fit_idx, axes, jet_of)
    for itr in range(c['n_iter'] if c['s5_mode'] != 'none' else 0):
        hard = np.argmax(W, 0); hard[W.max(0) <= c['hard_min']] = -1
        if c['iter_update_s3']:
            _, s3, _, _, _ = TR.ip_info(t['fpar'], t['fcov'], px, pC, axes, jet_of)
        iss = c['seed_sig'] if c['iter_seed_sig'] is None else c['iter_seed_sig']
        lm = (hard < 0) | ((hard == 0) & (s3[fit_idx] > c['iter_sig']))
        if v0set is not None and v0set.n:
            lm &= v0set.share <= c['hard_min']        # tracks held by a V0 do not seed SVs
        if c['iter_sv_chi2'] is not None:
            for k, v in enumerate(V):
                m = hard == k + 1
                if m.sum() >= 4:
                    ch, _ = TR.chi2_to_point(t['fpar'][fit_idx], t['fcov'][fit_idx], v.x)
                    lm |= m & (ch > c['iter_sv_chi2'])
        left = fit_idx[lm & (s3[fit_idx] > iss)]
        if len(left) < 2:
            break
        new = _pairs(E, c if c['iter_fsig'] is None else dict(c, seed_fsig=c['iter_fsig'], acc_fsig2=c['iter_fsig']), left, rmax, px, pC)
        if not new:
            break
        newc = stage3(E, dict(c, arbitration=False), ok, rmax, new, px, np.zeros(0, np.int64), pool=left, pvC=pC)
        th5 = c['sl_score_s5'] if c['sl_score_s5'] is not None else (c['sl_score_min'] if c['sl_score_seeded'] else -1e9)
        if th5 > -1e8:
            newc = [v for v in newc if _score(E, v.trk, v.x, v.C, v.chi2, px, pC, s3, c['sl_coefs_cand']) >= th5]
        if not newc:
            break
        allc = [Vtx(v.x, v.C, v.trk) for v in V] + [Vtx(v.x, v.C, [int(loc[k]) for k in v.trk if loc[k] >= 0]) for v in newc]
        px, pC, V, W = stage4(E, c, ok, rmax, px, pC, allc, s3, fit_idx, axes, jet_of)
    wpv[fit_idx] = W[0] if E.get('_pvw') is None else E['_pvw']     # the weights the PV was actually fitted with
    hard = np.argmax(W, 0); hard[W.max(0) <= c['hard_min']] = -1
    if c['split_casc'] and len(V):
        E['_rmax'] = rmax
        V, hard = _split_cascades(E, c, fit_idx, V, hard, px, pC, s3)
    assign[fit_idx] = hard
    wass[fit_idx] = W.max(0)
    for k, v in enumerate(V):
        tr = fit_idx[hard == k + 1]
        svs.append((v.x, v.C, _mass(t['fpar'], tr, v.x)[0] if len(tr) else 0., len(tr)))
    _, s3o, _, _, _ = TR.ip_info(t['fpar'], t['fcov'], px, pC, axes, jet_of)
    out = dict(pv=px, pvC=pC, svs=svs, wpv=wpv, assign=assign, wass=wass, cpu=0., s3=s3o)
    if c['v0_class']:
        tv0 = np.full(n, -1, np.int64)
        v0s = []
        for k in range(v0set.n):
            if v0set.P[k] > c['v0_out_min']:
                K = v0set.k[k]
                w0 = float(v0set.P0[k]) if v0set.P0 is not None else -1.
                v0s.append((K['a'], K['b'], V0.SPEC_PDG[K['hyp']], K['hyp'], K['x'], K['C'], K['m'], float(v0set.P[k]), K['p'], w0))
                if v0set.P[k] > c['hard_min']:
                    tv0[K['a']] = tv0[K['b']] = len(v0s) - 1
        out['v0s'] = v0s
        out['tv0'] = tv0
        if c['v0_pt_point']:
            # origins of the output V0s and the V0 weight per SV (state of the last apply at the last temperature)
            st = v0set.org.get(c['temps'][-1])
            org, orgw = [], []
            ptw = np.zeros(len(V))
            for k in range(v0set.n):
                if v0set.P[k] > c['v0_out_min']:
                    if st is None:
                        org.append(-1); orgw.append(-1.)
                        continue
                    j = int(np.argmax(st[0][k]))
                    org.append(j); orgw.append(float(st[0][k, j]))
            if st is not None and st[1] is not None:
                for vi in range(min(len(V), st[0].shape[1] - 1)):
                    ptw[vi] = float((st[1] * st[0][:, vi + 1] * ~v0set.isconv).sum())
            out['v0org'] = (org, orgw, ptw)
    out['cpu'] = time.process_time() - t0
    return out


def _v0_ref_denominators(c, t, ok, pvx, s3):
    """S4 denominators (outlier constant + PV numerator, T = 1) of the S0 tracks at the reference PV: admission and
    calibration of the V0 class (inf for other tracks)"""
    A = np.full(len(t['q']), np.inf)
    oi = np.nonzero(ok)[0]
    if len(oi):
        ch, _ = TR.chi2_to_point(t['fpar'][oi], t['fcov'][oi], pvx)
        cut = _f(c, np.full(len(oi), c['chi2_cut'], dtype=np.float64), 1.)
        if c['chi2_cut_neg'] is not None:
            neg = s3[oi] < 0
            cut[neg] = _f(c, np.full(int(neg.sum()), c['chi2_cut_neg'], dtype=np.float64), 1.)
        A[oi] = cut + _f(c, ch, 1.) * _prior(c, s3[oi], 1.)
    return A


def oracle_event(E, prune=True, bs_sig=F.BS_SIG):
    """truth-primary oracle PV: tracks in the stage1 pre-selection window with truth class primary only (every track
    with b/c, V0 or other secondary ancestry and unmatched tracks dropped), beam constraint; prune=True: stage1's
    pruning (worst track removed while chi2 >= 5), prune=False: all truth-primary tracks kept"""
    t0 = time.process_time()
    t = E['t']
    n = len(t['q'])
    idx = np.arange(n)
    sel = idx[_pre_window(E, idx) & (t['cls'] == 0)]
    if len(sel) >= 1 and prune:
        x, C, sel, _ = F.teardown(t['fpar'], t['fcov'], sel, E['bs'], bs=E['bs'], bs_sig=np.asarray(bs_sig), chi2max=5.0)
    elif len(sel) >= 1:
        x, C, _, _ = F.fit_vertex(t['fpar'], t['fcov'], sel, E['bs'], bs=E['bs'], bs_sig=np.asarray(bs_sig), niter=15, tol=1e-6)
    else:
        x, C = E['bs'].copy(), np.diag(np.asarray(bs_sig) ** 2)
    wpv = np.zeros(n); wpv[sel] = 1.
    assign = np.full(n, -1, np.int64); assign[sel] = 0
    axes = IO.jet_axes(E)
    _, s3, _, _, _ = TR.ip_info(t['fpar'], t['fcov'], x, C, axes, t['jet'].astype(np.int64))
    return dict(pv=x, pvC=C, svs=[], wpv=wpv, assign=assign, wass=wpv.copy(), cpu=time.process_time() - t0, s3=s3)
