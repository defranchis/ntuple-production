"""Evaluator of the joint vertexer: the same metrics for every algorithm output on the same events.

Usage: python3 -m jv.evaluate <nproc> <out_tag> <sample> <rows.json> [--data <data sample>] [--ref nominal] [--nboot 100]
                              [--only <algo>,<algo>,...]
  reads  $JV_RUNS/<out_tag>/<sample>/<algo>/*.npz (jv.run) + the harvest inputs in $JV_INPUTS (truth)
  writes <rows.json> (metrics per algorithm, paired bootstrap differences vs --ref)
Definitions (all lengths um unless stated):
  PV residual r = PV - generated PV, per flavour (b / c / uds = event flavour from the primary quark, stage1 jetPID) and
    component: mean, median, 5 % trimmed mean (tmean), core width = IQR/1.349 (core), tail fraction = |r| > 3 x core of the
    reference algorithm (same threshold for every algorithm), pull core width/median (r / sigma); r_proj = r . u, u = unit
    vector from the generated PV to the least-displaced true weak-HF decay vertex of the event (HF bias direction)
  migration: truth class (primary / b / c / V0 / other secondary / unmatched) x hard assignment (PV only, SV, none) over
    all baseline tracks; mean PV weight per class, histograms of the PV weight and of the assigned-vertex weight per class
  SV truth: weak b/c decay vertices with >= 2 charged long-lived daughters and flight distance 30 um - 3 cm from the
    generated PV (denominator). Found = some reco SV holds >= 2 hard-assigned tracks from it; completeness = tracks of the
    vertex in its best SV / its baseline tracks (tracks of a vertex include its descendants, e.g. decays in flight);
    residual = |x_SV - x_true| (3D) of the best SV; efficiency, completeness and residual also vs true flight distance,
    hadron momentum, |cos theta|. Reco SV classes: HF (>= 2 tracks of truth class b or c, and b + c at least as many as V0 +
    other secondary and as primary + unmatched), V0/secondary (else, V0 + other secondary >= primary + unmatched), spurious
    (the rest). purity = HF fraction, also vs reconstructed flight distance, SV momentum (vector sum of the hard tracks at
    their perigee), |cos theta| of that momentum. uds fake rate = reco SVs per uds event, any class.
  V0 content of the SVs (MC inputs with the true-V0 arrays, t_v0t): a V0 daughter track = any baseline track whose
    largest-weight truth link is a daughter of a true V0 or a GALEPH continuation of one (v0t table: K0S -> pi pi,
    Lambda -> p pi, conversion -> e e; such tracks have t_cls = 3: jv/harvest.py); per 1000 events: HF-class SVs without
    any hard-assigned V0 daughter track (hf_clean), HF-class SVs with >= 1 (hf_v0dau), V0/secondary-class and spurious
    SVs (per algorithm: sv_per_1000_events[<class>]); clean-HF fraction = hf_clean / all SVs (per algorithm:
    clean_hf_fraction; purity counts hf_v0dau as pure).
  flight significance = sqrt(d^T (C_SV + C_PV)^-1 d), d = SV - PV.
  paired bootstrap: events resampled with replacement (same resample for both algorithms), std of the difference. Keys of
    diff_vs_ref: pv_<flavour>_<x|y|z>_<mean|median|tmean|core|tail|pull_core>, pv_<flavour>_proj_<mean|median|tmean>,
    sv_eff, sv_eff_b_hadron, sv_eff_c_hadron, sv_purity, sv_completeness, uds_fake_rate, sv_eff_vs_<var>_<bin>,
    sv_purity_vs_<var>_<bin>, cpu_ms_per_event, sv_per1000ev_<hf_clean|hf_v0dau|v0_or_secondary|spurious> (SVs per
    1000 events), sv_clean_hf_frac (the V0 content keys need inputs with t_v0t); each [difference, error].
"""
import glob
import json
import os
import sys
import numpy as np
from multiprocessing import Pool

from . import io as IO
from . import tracks as TR

FLAV = {'b': 5, 'c': 4, 'uds': 0}
CLS = ['primary', 'b', 'c', 'V0', 'other_sec', 'unmatched']
FD_BINS = [0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0]          # cm
P_BINS = [0, 5, 10, 15, 20, 30, 50]                        # GeV
CT_BINS = [0, 0.3, 0.5, 0.7, 0.85, 0.95, 1.0001]
RES_BINS = {'x': np.linspace(-300, 300, 81), 'y': np.linspace(-150, 150, 81), 'z': np.linspace(-300, 300, 81)}


def clean(x):
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, np.ndarray):
        return clean(x.tolist())
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        return float(x) if np.isfinite(x) else None
    return x


def load_algo(tag, sample, algo):
    fs = sorted(glob.glob(f'{IO.env_path("JV_RUNS")}/{tag}/{sample}/{algo}/*.npz'))
    parts = []
    for fn in fs:
        d = np.load(fn)
        parts.append({k: d[k] for k in d.files})
    parts.sort(key=lambda p: (str(p['stem']), int(p['a'])))
    return parts


def load_truth(parts):
    """per-event and per-track truth/input arrays aligned with the unit order of `parts`"""
    T = {k: [] for k in ('flav', 'gen', 'pv_nom', 'ntrk', 'thrust', 't_cls', 't_hfv', 'h_ev', 'h_pos', 'h_nch', 'h_fd',
                         'h_p3', 'h_isb', 'h_nmatched', 'near', 'unit', 't_ev', 'isMC', 't_p3', 'h_parent', 't_v0d')}
    files = {}
    ev0 = 0
    for p in parts:
        stem, a, b = str(p['stem']), int(p['a']), int(p['b'])
        if stem not in files:
            files = {stem: IO.File(stem)}
        f = files[stem]
        d, o = f.d, f.off
        sl = slice(a, b)
        T['flav'].append(d['ev_flav'][sl]); T['ntrk'].append(d['ev_ntrk'][sl]); T['thrust'].append(d['ev_thrust'][sl])
        T['pv_nom'].append(d['ev_pv'][sl])
        T['gen'].append(d['ev_gen'][sl])
        ta, tb = o['t'][a], o['t'][b]
        T['t_ev'].append(np.repeat(np.arange(b - a) + ev0, d['ev_ntrk'][sl]))
        par = d['t_par'][ta:tb]
        pt = TR.PT_CONST / np.maximum(np.abs(par[:, 2]), 1e-9)     # legacy constant: tables comparable across rows
        T['t_p3'].append(np.stack([pt * np.cos(par[:, 1]), pt * np.sin(par[:, 1]), pt * par[:, 4]], 1))
        if f.isMC:
            T['t_cls'].append(d['t_cls'][ta:tb])
            if 't_v0t' in d:                                  # true-V0 daughter track
                T['t_v0d'].append(d['t_v0t'][ta:tb] >= 0)
            # hfv local -> global index in the concatenated HF table
            ha, hb = o['h'][a], o['h'][b]
            hoff = np.repeat(o['h'][a:b] - ha, d['ev_ntrk'][sl])
            hv = d['t_hfv'][ta:tb].astype(np.int64)
            nh_before = sum(len(x) for x in T['h_nch'])
            T['t_hfv'].append(np.where(hv >= 0, hv + hoff + nh_before, -1))
            T['h_ev'].append(np.repeat(np.arange(b - a) + ev0, d['ev_nhf'][sl]))
            for k in ('pos', 'nch', 'fd', 'p3', 'isb', 'nmatched'):
                T['h_' + k].append(d['h_' + k][ha:hb])
            pa_ = d['h_parent'][ha:hb].astype(np.int64)
            T['h_parent'].append(np.where(pa_ >= 0, pa_ + np.repeat(o['h'][a:b] - ha, d['ev_nhf'][sl]) + nh_before, -1))
            # least-displaced true weak-HF decay vertex per event
            near = np.full((b - a, 3), np.nan)
            hev = np.repeat(np.arange(b - a), d['ev_nhf'][sl])
            fd = d['h_fd'][ha:hb]
            order = np.lexsort((fd, hev))
            first = np.r_[True, hev[order][1:] != hev[order][:-1]] if len(order) else np.zeros(0, bool)
            near[hev[order][first]] = d['h_pos'][ha:hb][order][first]
            T['near'].append(near)
        T['isMC'].append(np.full(b - a, f.isMC))
        ev0 += b - a
    out = {}
    for k, v in T.items():
        if v:
            out[k] = np.concatenate(v)
    if 'near' in out:
        u = out['near'] - out['gen']
        nu = np.linalg.norm(u, axis=1)
        out['unit'] = u / np.where(nu > 0, nu, np.nan)[:, None]
    return out


def concat(parts):
    o = {}
    ev0 = 0
    sv_ev = []
    for p in parts:
        sv_ev.append(p['sv_ev'] + ev0)
        ev0 += len(p['ev_row'])
    for k in parts[0]:
        if k in ('stem', 'a', 'b', 'config'):
            continue
        if k == 'sv_ev':
            o[k] = np.concatenate(sv_ev).astype(np.int64)
        else:
            o[k] = np.concatenate([p[k] for p in parts])
    o['config'] = str(parts[0]['config'])
    return o


# ------------------------------------------------------------------------------------------------------------- statistics
def core(v):
    v = v[np.isfinite(v)]
    if len(v) < 5:
        return np.nan
    q1, q3 = np.quantile(v, [0.25, 0.75])
    return (q3 - q1) / 1.349


def tmean(v, f=0.05):
    v = np.sort(v[np.isfinite(v)])
    k = int(f * len(v))
    return float(np.mean(v[k:len(v) - k])) if len(v) > 2 * k else np.nan


def stats_set(r):
    r = r[np.isfinite(r)]
    return dict(n=int(len(r)), mean=float(np.mean(r)) if len(r) else np.nan, median=float(np.median(r)) if len(r) else np.nan,
                tmean=tmean(r), core=core(r))


def boot_diff(fn, a, b, nboot, rng):
    """paired bootstrap: std of fn(a[idx]) - fn(b[idx]) over event resamples; returns (diff, err)"""
    n = len(a)
    d0 = fn(a) - fn(b)
    ds = []
    for _ in range(nboot):
        idx = rng.integers(0, n, n)
        ds.append(fn(a[idx]) - fn(b[idx]))
    return float(d0), float(np.nanstd(ds))


# ------------------------------------------------------------------------------------------------------------- metrics
def pv_metrics(A, T, ref_core=None):
    res = (A['ev_pv'] - T['gen']) * 1e4
    sig = np.sqrt(np.maximum(A['ev_pvc6'][:, :3], 1e-30)) * 1e4
    pull = res / sig
    proj = np.einsum('ij,ij->i', res, T['unit']) if 'unit' in T else None
    out = {}
    for fn, fv in FLAV.items():
        m = T['flav'] == fv
        o = {}
        for k, c in enumerate('xyz'):
            s = stats_set(res[m, k])
            thr = 3 * (ref_core[fn][c] if ref_core else s['core'])
            s['tail_thr_um'] = float(thr)
            s['tail'] = float(np.mean(np.abs(res[m, k]) > thr))
            s['pull_core'] = core(pull[m, k]); s['pull_median'] = float(np.median(pull[m, k]))
            hb = RES_BINS[c]
            s['hist'] = np.histogram(np.clip(res[m, k], hb[0], hb[-1] - 1e-9), hb)[0].tolist()
            s['hist_bins'] = hb.tolist()
            o[c] = s
        if proj is not None:
            ok = m & np.isfinite(proj)
            o['proj'] = stats_set(proj[ok])
        out[fn] = o
    return out, res, proj


def migration(A, T):
    cls = T['t_cls']
    asg = A['t_assign']
    col = np.where(asg == 0, 0, np.where(asg > 0, 1, 2))
    M = np.zeros((6, 3), np.int64)
    np.add.at(M, (cls, col), 1)
    wpv = {CLS[k]: float(np.mean(A['t_wpv'][cls == k])) if np.any(cls == k) else np.nan for k in range(6)}
    wb = np.linspace(0, 1, 21)
    hw = {CLS[k]: np.histogram(A['t_wpv'][cls == k], wb)[0].tolist() for k in range(6)}
    hs = {CLS[k]: np.histogram(A['t_wass'][(cls == k) & (asg > 0)], wb)[0].tolist() for k in range(6)}
    return dict(counts=M.tolist(), cols=['PV', 'SV', 'none'], rows=CLS,
                frac_row=(M / np.maximum(M.sum(1, keepdims=True), 1)).tolist(), mean_wpv=wpv, weight_bins=wb.tolist(),
                hist_wpv=hw, hist_wass_sv=hs)


def sv_tables(A, T):
    """per reco SV: class, best truth vertex; per truth vertex: found, completeness, residual"""
    nsv = len(A['sv_ev'])
    sv_off = np.zeros(len(A['ev_nsv']) + 1, np.int64)
    np.cumsum(A['ev_nsv'], out=sv_off[1:])
    asg = A['t_assign'].astype(np.int64)
    tev = T['t_ev']
    has = asg > 0
    gsv = np.where(has, sv_off[tev] + asg - 1, -1)        # global SV index of each track
    out = {'sv_class': np.full(nsv, 2, np.int8), 'sv_best_h': np.full(nsv, -1, np.int64), 'sv_nbest': np.zeros(nsv, np.int64)}
    P = np.zeros((nsv, 3))
    if nsv:
        np.add.at(P, gsv[has], T['t_p3'][has])
    pm = np.linalg.norm(P, axis=1)
    out['sv_p'] = pm
    out['sv_abscos'] = np.abs(P[:, 2]) / np.maximum(pm, 1e-9)
    out['sv_fd'] = np.linalg.norm(A['sv_pos'] - A['ev_pv'][A['sv_ev']], axis=1) if nsv else np.zeros(0)
    if 't_cls' not in T:
        return out, None
    cls = T['t_cls']; hfv = T['t_hfv']
    tk = np.nonzero(has)[0]
    # majority class per SV
    cnt = np.zeros((nsv, 6), np.int64)
    np.add.at(cnt, (gsv[tk], cls[tk]), 1)
    hfc = cnt[:, 1] + cnt[:, 2]; v0c = cnt[:, 3] + cnt[:, 4]; prc = cnt[:, 0] + cnt[:, 5]
    out['sv_class'] = np.where((hfc >= 2) & (hfc >= v0c) & (hfc >= prc), 0, np.where(v0c >= prc, 1, 2)).astype(np.int8)
    # (SV, truth vertex) pair counts
    th = tk[hfv[tk] >= 0]
    pairs = {}
    for s_, h_ in zip(gsv[th], hfv[th]):
        pairs[(s_, h_)] = pairs.get((s_, h_), 0) + 1
    nh = len(T['h_nch'])
    found = np.zeros(nh, bool); best_n = np.zeros(nh, np.int64); best_s = np.full(nh, -1, np.int64)
    for (s_, h_), n_ in pairs.items():
        if n_ > best_n[h_]:
            best_n[h_] = n_; best_s[h_] = s_
        if n_ > out['sv_nbest'][s_]:
            out['sv_nbest'][s_] = n_; out['sv_best_h'][s_] = h_
    found = best_n >= 2
    nbase = np.bincount(hfv[hfv >= 0], minlength=nh)
    resid = np.full(nh, np.nan)
    ok = best_s >= 0
    resid[ok] = np.linalg.norm(A['sv_pos'][best_s[ok]] - T['h_pos'][ok], axis=1) * 1e4
    one = np.zeros(nh, bool)
    one[out['sv_best_h'][out['sv_nbest'] >= 2]] = True
    H = dict(found=found, best_n=best_n, nbase=nbase, resid=resid, best_s=best_s, one=one)
    # track-level content of the SVs
    out['sv_cnt'] = cnt
    if 't_v0d' in T and len(T['t_v0d']) == len(cls):
        v0k = tk[T['t_v0d'][tk]]
        out['sv_nv0dau'] = np.bincount(gsv[v0k], minlength=nsv)
    if nsv:
        d = A['sv_pos'] - A['ev_pv'][A['sv_ev']]
        back = np.einsum('ij,ij->i', d[gsv[tk]], T['t_p3'][tk]) < 0
        out['sv_nback'] = np.bincount(gsv[tk][back], minlength=nsv)
    else:
        out['sv_nback'] = np.zeros(0, np.int64)
    return out, H


def sv_v0_split(S):
    """per reco SV: booleans hf_clean (HF class, no V0 daughter track), hf_v0dau (HF class, >= 1), v0_or_secondary, spurious"""
    sc, nv = S['sv_class'], S['sv_nv0dau']
    return {'hf_clean': (sc == 0) & (nv == 0), 'hf_v0dau': (sc == 0) & (nv > 0), 'v0_or_secondary': sc == 1, 'spurious': sc == 2}


def cascade_flags(H, T, den):
    """for c hadrons with a b-hadron parent, both in the denominator: BOTH (found in different SVs), MERGED (found, same best SV),
    B_ONLY, D_ONLY, NONE"""
    d = np.nonzero(den & (T['h_parent'] >= 0))[0]
    d = d[den[T['h_parent'][d]] & T['h_isb'][T['h_parent'][d]] & ~T['h_isb'][d]]
    b = T['h_parent'][d]
    fb, fd_ = H['found'][b], H['found'][d]
    same = H['best_s'][b] == H['best_s'][d]
    return d, dict(BOTH=fb & fd_ & ~same, MERGED=fb & fd_ & same, B_ONLY=fb & ~fd_, D_ONLY=~fb & fd_, NONE=~fb & ~fd_)


def sv_metrics(A, T, S, H):
    o = {}
    flav = T['flav']
    nsv_ev = A['ev_nsv']
    o['n_sv_per_event'] = {fn: float(np.mean(nsv_ev[flav == fv])) for fn, fv in FLAV.items()} if 'flav' in T else {}
    o['uds_fake_rate'] = float(np.mean(nsv_ev[flav == 0]))
    o['uds_frac_events_with_sv'] = float(np.mean(nsv_ev[flav == 0] > 0))
    if H is None:
        return o
    sc = S['sv_class']
    o['purity'] = float(np.mean(sc == 0)) if len(sc) else np.nan
    o['sv_class_frac'] = {k: float(np.mean(sc == i)) if len(sc) else np.nan for i, k in enumerate(['HF', 'V0_or_secondary', 'spurious'])}
    den = (T['h_nch'] >= 2) & (T['h_fd'] >= 0.003) & (T['h_fd'] <= 3.0)
    hfl = T['flav'][T['h_ev']]
    for nm, sel in (('all', den), ('b_hadron', den & T['h_isb']), ('c_hadron', den & ~T['h_isb'])):
        o[f'eff_{nm}'] = float(np.mean(H['found'][sel])) if sel.any() else np.nan
        o[f'n_den_{nm}'] = int(sel.sum())
    o['eff_in_b_events'] = float(np.mean(H['found'][den & (hfl == 5)]))
    o['eff_in_c_events'] = float(np.mean(H['found'][den & (hfl == 4)]))
    f = den & H['found']
    o['completeness_mean'] = float(np.mean(H['best_n'][f] / np.maximum(H['nbase'][f], 1)))
    for nm, sel in (('all', den), ('b_hadron', den & T['h_isb']), ('c_hadron', den & ~T['h_isb'])):
        o[f'eff_one2one_{nm}'] = float(np.mean(H['one'][sel])) if sel.any() else np.nan
    cnt = S['sv_cnt']
    nt = max(int(cnt.sum()), 1)
    o['trk_frac'] = {'hf': float((cnt[:, 1] + cnt[:, 2]).sum() / nt), 'primary': float(cnt[:, 0].sum() / nt), 'v0': float(cnt[:, 3].sum() / nt),
                     'other': float(cnt[:, 4].sum() / nt), 'unmatched': float(cnt[:, 5].sum() / nt), 'back': float(S['sv_nback'].sum() / nt)}
    o['tracks_per_sv'] = float(cnt.sum() / max(len(cnt), 1))
    if 'sv_nv0dau' in S:
        nev = len(nsv_ev)
        o['sv_per_1000_events'] = {k: float(1e3 * v.sum() / nev) for k, v in sv_v0_split(S).items()}
        o['clean_hf_fraction'] = float(np.mean(sv_v0_split(S)['hf_clean'])) if len(sc) else np.nan
    o['mass_ge_5p8'] = float(np.mean(A['sv_mass'] >= 5.8)) if len(A['sv_mass']) else np.nan
    d_, cf = cascade_flags(H, T, den)
    o['cascade'] = {k: float(np.mean(v)) if len(d_) else np.nan for k, v in cf.items()}
    o['cascade_n'] = int(len(d_))
    o['resid_median_um'] = float(np.nanmedian(H['resid'][f])); o['resid_core_um'] = core(H['resid'][f])
    p = np.linalg.norm(T['h_p3'], axis=1)
    ct = np.abs(T['h_p3'][:, 2]) / np.maximum(p, 1e-9)
    for nm, var, bins in (('fd_cm', T['h_fd'], FD_BINS), ('p_GeV', p, P_BINS), ('abscos', ct, CT_BINS)):
        effs, ns, comp, rs = [], [], [], []
        for lo, hi in zip(bins[:-1], bins[1:]):
            s = den & (var >= lo) & (var < hi)
            effs.append(float(np.mean(H['found'][s])) if s.any() else np.nan); ns.append(int(s.sum()))
            ff = s & H['found']
            comp.append(float(np.mean(H['best_n'][ff] / np.maximum(H['nbase'][ff], 1))) if ff.any() else np.nan)
            rs.append(float(np.nanmedian(H['resid'][ff])) if ff.any() else np.nan)
        o[f'eff_vs_{nm}'] = dict(bins=bins, eff=effs, n=ns, completeness=comp, resid_median_um=rs)
    for nm, var, bins in (('fd_cm', S['sv_fd'], FD_BINS), ('p_GeV', S['sv_p'], P_BINS), ('abscos', S['sv_abscos'], CT_BINS)):
        pur, ns = [], []
        for lo, hi in zip(bins[:-1], bins[1:]):
            s = (var >= lo) & (var < hi)
            pur.append(float(np.mean(sc[s] == 0)) if s.any() else np.nan); ns.append(int(s.sum()))
        o[f'purity_vs_{nm}'] = dict(bins=bins, purity=pur, n=ns)
    return o


def hemi_corr(A, T):
    """SV-mass hemisphere tag (an SV with mass > 1.8 GeV on the hemisphere side of the thrust axis) in b events:
    C = f_double / f_single^2 (sanity only)"""
    b = np.nonzero(T['flav'] == 5)[0]
    if len(b) == 0 or len(A['sv_ev']) == 0:
        return np.nan
    side = np.sign(np.einsum('ij,ij->i', A['sv_pos'] - A['ev_pv'][A['sv_ev']], T['thrust'][A['sv_ev']]))
    tag = np.zeros((len(A['ev_nsv']), 2), bool)
    m = A['sv_mass'] > 1.8
    tag[A['sv_ev'][m & (side > 0)], 0] = True
    tag[A['sv_ev'][m & (side <= 0)], 1] = True
    tb = tag[b]
    f1 = tb.mean()
    f2 = np.mean(tb[:, 0] & tb[:, 1])
    return float(f2 / f1 ** 2) if f1 > 0 else np.nan


def reco_dists(A, T):
    """reconstructed quantities for data/MC (arrays)"""
    npv = np.bincount(T['t_ev'][A['t_assign'] == 0], minlength=len(A['ev_nsv']))
    w = A['t_wpv'][A['t_wpv'] > 0.01]
    fd = np.linalg.norm(A['sv_pos'] - A['ev_pv'][A['sv_ev']], axis=1) if len(A['sv_ev']) else np.zeros(0)
    fsig = np.zeros(0)
    if len(A['sv_ev']):
        d = A['sv_pos'] - A['ev_pv'][A['sv_ev']]
        c6 = A['sv_c6'] + A['ev_pvc6'][A['sv_ev']]
        C = np.stack([c6[:, 0], c6[:, 3], c6[:, 4], c6[:, 3], c6[:, 1], c6[:, 5], c6[:, 4], c6[:, 5], c6[:, 2]], 1).reshape(-1, 3, 3)
        ok = np.linalg.det(C) > 0          # an SV without a valid covariance has no significance
        fsig = np.full(len(d), np.inf)
        fsig[ok] = np.sqrt(np.maximum(np.einsum('ki,ki->k', d[ok], np.linalg.solve(C[ok], d[ok][:, :, None])[:, :, 0]), 0.))
    return dict(n_pv_tracks=npv, pv_weight=w, n_sv=A['ev_nsv'], sv_mass=A['sv_mass'], sv_fd_cm=fd, sv_fsig=fsig,
                trk_ipsig=A['t_s3'])


DIST_BINS = dict(n_pv_tracks=np.arange(-0.5, 40.5, 1), pv_weight=np.linspace(0, 1, 26), n_sv=np.arange(-0.5, 8.5, 1),
                 sv_mass=np.linspace(0, 6, 31), sv_fd_cm=np.logspace(-3, 1, 33), sv_fsig=np.logspace(0, 3, 31),
                 trk_ipsig=np.linspace(-20, 40, 61))


def hist_dists(D):
    o = {}
    for k, b in DIST_BINS.items():
        v = np.asarray(D[k], float)
        h, _ = np.histogram(np.clip(v, b[0], b[-1] - 1e-9), b)
        o[k] = dict(bins=b.tolist(), counts=h.tolist(), mean=float(np.mean(v)) if len(v) else np.nan, n=int(len(v)))
    return o


# ------------------------------------------------------------------------------------------------------------- driver
def eval_algo(args):
    tag, sample, algo, ref_core = args
    parts = load_algo(tag, sample, algo)
    if not parts:
        return algo, None, None
    A = concat(parts)
    T = load_truth(parts)
    out = {'n_events': int(len(A['ev_nsv'])), 'cpu_ms_per_event': float(1e3 * np.mean(A['ev_cpu'])), 'config': A['config']}
    S, H = sv_tables(A, T)
    if T['isMC'].all():
        out['pv'], res, proj = pv_metrics(A, T, ref_core)
        out['migration'] = migration(A, T)
        out['hemi_corr_b'] = hemi_corr(A, T)
    out['sv'] = sv_metrics(A, T, S, H)
    out['dists'] = hist_dists(reco_dists(A, T))
    return algo, out, None


def ratio_diffs(items, nev, nboot, rng):
    """items: key -> (num_a, den_a, num_b, den_b), per-event arrays. Difference of the ratios of sums, paired bootstrap"""
    def rat(n, d, w=None):
        ds = d.sum() if w is None else d @ w
        return (n.sum() if w is None else n @ w) / ds if ds > 0 else np.nan
    d0 = {k: rat(v[0], v[1]) - rat(v[2], v[3]) for k, v in items.items()}
    bs = {k: [] for k in items}
    for _ in range(nboot):
        w = np.bincount(rng.integers(0, nev, nev), minlength=nev).astype(float)
        for k, v in items.items():
            bs[k].append(rat(v[0], v[1], w) - rat(v[2], v[3], w))
    return {k: (float(d0[k]), float(np.nanstd(bs[k]))) for k in items}


def paired(args):
    """paired bootstrap differences algo - ref on the same events"""
    tag, sample, algo, ref, nboot = args
    rng = np.random.default_rng(12345)
    PA, PB = load_algo(tag, sample, algo), load_algo(tag, sample, ref)
    A, B = concat(PA), concat(PB)
    T = load_truth(PA)
    nev = len(A['ev_nsv'])
    assert nev == len(B['ev_nsv']) and np.array_equal(A['ev_row'], B['ev_row'])
    o = {}
    items = {}
    if T['isMC'].all():
        ra = (A['ev_pv'] - T['gen']) * 1e4; rb = (B['ev_pv'] - T['gen']) * 1e4
        sa = np.sqrt(np.maximum(A['ev_pvc6'][:, :3], 1e-30)) * 1e4; sb = np.sqrt(np.maximum(B['ev_pvc6'][:, :3], 1e-30)) * 1e4
        pa = np.einsum('ij,ij->i', ra, T['unit']); pb = np.einsum('ij,ij->i', rb, T['unit'])
        for fn, fv in FLAV.items():
            m = T['flav'] == fv
            for k, c in enumerate('xyz'):
                for st, f in (('mean', np.mean), ('median', np.median), ('tmean', tmean), ('core', core)):
                    o[f'pv_{fn}_{c}_{st}'] = boot_diff(f, ra[m, k], rb[m, k], nboot, rng)
                o[f'pv_{fn}_{c}_pull_core'] = boot_diff(core, ra[m, k] / sa[m, k], rb[m, k] / sb[m, k], nboot, rng)
                thr = 3 * core(rb[m, k])        # tail threshold = 3 x core of the reference, as in pv_metrics
                items[f'pv_{fn}_{c}_tail'] = ((np.abs(ra[:, k]) > thr) & m, m, (np.abs(rb[:, k]) > thr) & m, m)
            ok = m & np.isfinite(pa) & np.isfinite(pb)
            for st, f in (('mean', np.mean), ('median', np.median), ('tmean', tmean)):
                o[f'pv_{fn}_proj_{st}'] = boot_diff(f, pa[ok], pb[ok], nboot, rng)
        SA, HA = sv_tables(A, T); SB, HB = sv_tables(B, T)
        den = (T['h_nch'] >= 2) & (T['h_fd'] >= 0.003) & (T['h_fd'] <= 3.0)
        hev = T['h_ev']

        def per_ev(ev, w):
            return np.bincount(ev, weights=np.asarray(w, float), minlength=nev)

        def eff_item(sel):
            nd = per_ev(hev[sel], np.ones(sel.sum()))
            return (per_ev(hev[sel], HA['found'][sel]), nd, per_ev(hev[sel], HB['found'][sel]), nd)
        items['sv_eff'] = eff_item(den)

        def one_item(sel):
            nd = per_ev(hev[sel], np.ones(sel.sum()))
            return (per_ev(hev[sel], HA['one'][sel]), nd, per_ev(hev[sel], HB['one'][sel]), nd)
        items['sv_eff_one2one'] = one_item(den)
        items['sv_eff_one2one_b_hadron'] = one_item(den & T['h_isb'])
        items['sv_eff_one2one_c_hadron'] = one_item(den & ~T['h_isb'])
        dA_, cA = cascade_flags(HA, T, den); dB_, cB = cascade_flags(HB, T, den)
        for k_ in cA:
            nd = per_ev(hev[dA_], np.ones(len(dA_)))
            items['sv_cascade_' + k_] = (per_ev(hev[dA_], cA[k_]), nd, per_ev(hev[dB_], cB[k_]), nd)

        def trkfrac(S_, X, idx):
            return per_ev(X['sv_ev'], S_['sv_cnt'][:, idx].sum(1)), per_ev(X['sv_ev'], S_['sv_cnt'].sum(1))
        items['sv_trk_hf_frac'] = trkfrac(SA, A, [1, 2]) + trkfrac(SB, B, [1, 2])
        items['sv_trk_primary_frac'] = trkfrac(SA, A, [0]) + trkfrac(SB, B, [0])
        items['sv_trk_back_frac'] = (per_ev(A['sv_ev'], SA['sv_nback']), per_ev(A['sv_ev'], SA['sv_cnt'].sum(1)),
                                     per_ev(B['sv_ev'], SB['sv_nback']), per_ev(B['sv_ev'], SB['sv_cnt'].sum(1)))
        items['sv_eff_b_hadron'] = eff_item(den & T['h_isb'])
        items['sv_eff_c_hadron'] = eff_item(den & ~T['h_isb'])
        p = np.linalg.norm(T['h_p3'], axis=1)
        ct = np.abs(T['h_p3'][:, 2]) / np.maximum(p, 1e-9)
        for nm, var, bins in (('fd_cm', T['h_fd'], FD_BINS), ('p_GeV', p, P_BINS), ('abscos', ct, CT_BINS)):
            for i, (lo, hi) in enumerate(zip(bins[:-1], bins[1:])):
                items[f'sv_eff_vs_{nm}_{i}'] = eff_item(den & (var >= lo) & (var < hi))

        def comp(H):
            f = den & H['found']
            return per_ev(hev[f], H['best_n'][f] / np.maximum(H['nbase'][f], 1)), per_ev(hev[f], np.ones(f.sum()))
        items['sv_completeness'] = comp(HA) + comp(HB)

        def pur(S, X, sel):
            return per_ev(X['sv_ev'][sel], S['sv_class'][sel] == 0), per_ev(X['sv_ev'][sel], np.ones(sel.sum()))
        items['sv_purity'] = pur(SA, A, np.ones(len(SA['sv_class']), bool)) + pur(SB, B, np.ones(len(SB['sv_class']), bool))
        for nm, key, bins in (('fd_cm', 'sv_fd', FD_BINS), ('p_GeV', 'sv_p', P_BINS), ('abscos', 'sv_abscos', CT_BINS)):
            for i, (lo, hi) in enumerate(zip(bins[:-1], bins[1:])):
                items[f'sv_purity_vs_{nm}_{i}'] = pur(SA, A, (SA[key] >= lo) & (SA[key] < hi)) + pur(SB, B, (SB[key] >= lo) & (SB[key] < hi))
        # track migration (truth class x PV / SV) and track purity (b/c fraction of SV tracks)
        tev = T['t_ev']; cls = T['t_cls']

        def mig(asg, k, sel):
            m = cls == k
            return per_ev(tev[m], sel(asg[m])), per_ev(tev[m], np.ones(m.sum()))
        for k, nm in ((0, 'primary'), (1, 'b'), (2, 'c'), (3, 'V0'), (4, 'other'), (5, 'unmatched')):
            items[f'mig_{nm}_in_pv'] = mig(A['t_assign'], k, lambda a: a == 0) + mig(B['t_assign'], k, lambda a: a == 0)
            items[f'mig_{nm}_in_sv'] = mig(A['t_assign'], k, lambda a: a > 0) + mig(B['t_assign'], k, lambda a: a > 0)

        def trkpur(X):
            s = X['t_assign'] > 0
            hf = (cls == 1) | (cls == 2)
            return per_ev(tev[s], hf[s]), per_ev(tev[s], np.ones(s.sum()))
        items['sv_trkpur_bc'] = trkpur(A) + trkpur(B)
        uds = (T['flav'] == 0).astype(float)
        items['uds_fake_rate'] = (A['ev_nsv'] * uds, uds, B['ev_nsv'] * uds, uds)
        if 'sv_nv0dau' in SA and 'sv_nv0dau' in SB:           # V0 content of the SVs (inputs with t_v0t)
            VA, VB = sv_v0_split(SA), sv_v0_split(SB)
            one = np.ones(nev)
            for k_ in VA:
                items['sv_per1000ev_' + k_] = (per_ev(A['sv_ev'], 1e3 * VA[k_]), one, per_ev(B['sv_ev'], 1e3 * VB[k_]), one)
            items['sv_clean_hf_frac'] = (per_ev(A['sv_ev'], VA['hf_clean']), A['ev_nsv'].astype(float),
                                         per_ev(B['sv_ev'], VB['hf_clean']), B['ev_nsv'].astype(float))
    o.update(ratio_diffs(items, nev, nboot, rng))
    ca, cb = 1e3 * A['ev_cpu'], 1e3 * B['ev_cpu']
    o['cpu_ms_per_event'] = boot_diff(np.mean, ca, cb, nboot, rng)
    return algo, o


def main():
    a = sys.argv[1:]
    nproc, tag, sample, out = int(a[0]), a[1], a[2], a[3]
    opts = dict(zip(a[4::2], a[5::2]))
    ref = opts.get('--ref', 'nominal')
    nboot = int(opts.get('--nboot', 100))
    data = opts.get('--data')
    algos = sorted(os.path.basename(p) for p in glob.glob(f'{IO.env_path("JV_RUNS")}/{tag}/{sample}/*'))
    if '--only' in opts:
        algos = [al for al in algos if al in opts['--only'].split(',') or al == ref]
    # reference core widths (tail threshold) first
    _, refo, _ = eval_algo((tag, sample, ref, None))
    ref_core = {fn: {c: refo['pv'][fn][c]['core'] for c in 'xyz'} for fn in FLAV} if 'pv' in refo else None
    with Pool(nproc) as p:
        res = dict((k, v) for k, v, _ in p.map(eval_algo, [(tag, sample, al, ref_core) for al in algos]))
        diffs = dict(p.map(paired, [(tag, sample, al, ref, nboot) for al in algos if al != ref]))
        dres = {}
        if data:
            dalgos = sorted(os.path.basename(x) for x in glob.glob(f'{IO.env_path("JV_RUNS")}/{tag}/{data}/*'))
            dres = dict((k, v) for k, v, _ in p.map(eval_algo, [(tag, data, al, None) for al in dalgos]))
    rows = {'sample': sample, 'tag': tag, 'ref': ref, 'nboot': nboot, 'definitions': __doc__, 'algos': {}}
    for al in algos:
        r = res[al]
        r['diff_vs_ref'] = diffs.get(al, {})
        rows['algos'][al] = r
    if data:
        rows['data_sample'] = data
        rows['data'] = dres
    json.dump(clean(rows), open(out, 'w'), indent=1)
    print('wrote', out)


if __name__ == '__main__':
    main()
