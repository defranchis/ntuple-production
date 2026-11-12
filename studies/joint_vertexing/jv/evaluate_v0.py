"""V0 evaluator of the joint vertexer: efficiency and purity of reconstructed V0 candidates (neutral two-prong decays:
K0S -> pi+ pi-, Lambda / anti-Lambda -> p pi, photon conversions gamma -> e+ e-) against the true-V0 table of the
harvested files (v0t_* and the per-track t_v0t / t_v0leg, jv/harvest.py), MC only.

Usage: python3 -m jv.evaluate_v0 <nproc> <sample> <out.json> [--rows stage1_loose,stage1_tight] [--ref stage1_loose]
                                 [--nboot 100]
  reads the harvested files in $JV_INPUTS (truth + stage1's V0 list v_*); rows: 'stage1_loose' = every stored stage1
  V0 candidate (loose tier, tight included), 'stage1_tight' = v_tight == 1, 'run:<tag>:<algo>[:<wmin>]' = the V0
  list of a jv.run output with the V0 class ($JV_RUNS/<tag>/<sample>/<algo>/<stem>_*.npz, v0_* arrays) with final V0
  weight v0_w > wmin (default 0.5; <tag> may be '<tag>/<run sample>', e.g. to evaluate DEV9 on a DEV run);
  'union:<row>+<row>' = both lists, a pair booked twice with the same species counted once; 'conv_loose' /
  'conv_tight' = the conversion finder jv.conv.find_conversions (sel 'loose' / 'tight') run here on every event (booked
  as conversions, vertex = CF vertex, mass = m_ee at the tangent point).
Truth: a track is a daughter track of a true V0 when its largest-weight link is the daughter or a GALEPH continuation
  of it (jv/harvest.py, continuation_parents).
Candidates (per file, flat over events): ev (event row), t1, t2 (local track indices, -1 = not in the track list), pdg
  (booked species: 310 K0S, 3122 Lambda or anti-Lambda, 22 conversion), pos (fitted vertex, cm), mass (GeV).
Definitions:
  species of a true V0: K0S (310), Lambda (|pdg| 3122, both charges), conversion (22).
  leg in the track list: the daughter is the largest-weight truth link of a baseline track (v0t_l1 / v0t_l2 >= 0).
  daughter track: EVERY baseline track whose largest-weight truth link is a daughter of a true V0 (t_v0t >= 0; also
    when several tracks share one daughter, not only the highest-weight one kept in v0t_l1 / v0t_l2).
  matched candidate: its two tracks are daughter tracks of the two different daughters of one true V0 of the booked
    species.
  efficiency: true V0s with >= 1 matched candidate / denominator; two denominators: 'all' = every true charged-mode decay
    of the species in the selected events (any radius), 'tracked' = those with both legs in the track list; versus the
    TRUE V0 momentum p (GeV), transverse decay radius rxy (cm, from the generated PV), |cos theta| of the V0 momentum,
    and split by origin (primary / b / c / other secondary, v0t_orig).
  purity: matched candidates / candidates booked as the species; versus the FITTED momentum (sum of the two track
    momenta at their point of closest approach to the fitted vertex, legacy momentum constant), fitted rxy (from the
    stage1 PV, ev_pv) and |cos theta|.
  candidate classes (every candidate, by the truth of its two tracks): matched; wrong_species (both tracks are the
    two daughters of one true V0 of another species); two_v0s (each track a daughter of a different true V0);
    one_leg_linked (one track a true daughter, the other with a truth link; also both tracks linked to the SAME
    daughter); one_leg_unlinked (one true daughter, the other track without truth link, truth class 5: e.g. soft K0S
    daughters whose truth link is lost); hf_pair (both tracks from b/c,
    truth classes 1-2); primary_pair (both truth class 0); other_secondary (both from other secondaries, class 4, or one
    class 4 + one class 1-2); other (the rest: mixed primary + secondary, unlinked pairs). Groups: matched; V0-related
    (wrong_species, two_v0s, one_leg_*); unmatched but physical (hf_pair, other_secondary); spurious (primary_pair, other).
  errors: binomial sqrt(e (1 - e) / n) for single values; differences of a row to --ref by a paired event bootstrap
    (events resampled with replacement, the same resample for both rows, std of the difference of ratios of sums).
"""
import json
import sys
import numpy as np
from multiprocessing import Pool

from . import io as IO
from . import tracks as TR

SPECIES = {'K0S': 310, 'Lambda': 3122, 'conversion': 22}
ORIG = ['primary', 'b', 'c', 'other_sec']
CLASSES = ['matched', 'wrong_species', 'two_v0s', 'one_leg_linked', 'one_leg_unlinked', 'hf_pair', 'primary_pair',
           'other_secondary', 'other']
GROUPS = {'matched': ['matched'], 'v0_related': ['wrong_species', 'two_v0s', 'one_leg_linked', 'one_leg_unlinked'],
          'physical_unmatched': ['hf_pair', 'other_secondary'], 'spurious': ['primary_pair', 'other']}
P_BINS = [0.2, 0.5, 1., 1.5, 2., 3., 4., 6., 9., 13., 20., 46.]                 # GeV
R_BINS = [0., 0.5, 1., 2., 3., 5., 8., 12., 20., 30., 45., 70., 100., 180.]     # cm
C_BINS = [0., 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0001]
VARS = (('p_GeV', P_BINS), ('rxy_cm', R_BINS), ('abscos', C_BINS))


def species_of(pdg):
    a = np.abs(np.asarray(pdg))
    return np.where(a == 310, 310, np.where(a == 3122, 3122, np.where(a == 22, 22, 0)))


def track_momenta(fpar, x):
    """momentum (n, 3) of tracks fpar (fitter convention, jv/tracks.py) at the xy point nearest to x (n, 3)"""
    D, p0, C = fpar[:, 0], fpar[:, 1], fpar[:, 2]
    Cs = np.where(C == 0, 1e-12, C)
    x0, y0 = -D * np.sin(p0), D * np.cos(p0)
    xc, yc = -(D + 0.5 / Cs) * np.sin(p0), (D + 0.5 / Cs) * np.cos(p0)
    dphi = (np.arctan2(x[:, 1] - yc, x[:, 0] - xc) - np.arctan2(y0 - yc, x0 - xc) + np.pi) % (2 * np.pi) - np.pi
    L = np.where(C == 0, 0., dphi / (2 * Cs))
    pt = 0.5 * TR.PT_CONST / np.maximum(np.abs(C), 1e-12)     # legacy constant: same binning for every row
    ang = p0 + 2 * C * L
    return np.stack([pt * np.cos(ang), pt * np.sin(ang), pt * fpar[:, 4]], 1)


def stage1_cands(f, tight_only):
    """stage1's stored V0 list (v_*) of one io.File as a candidate dict"""
    d = f.d
    ev = np.repeat(np.arange(len(d['ev_nv0'])), d['ev_nv0'])
    sel = d['v_tight'].astype(bool) if tight_only else np.ones(len(ev), bool)
    return dict(ev=ev[sel], t1=d['v_t1'][sel].astype(np.int64), t2=d['v_t2'][sel].astype(np.int64),
                pdg=species_of(d['v_pdg'][sel]), pos=d['v_pos'][sel], mass=d['v_mass'][sel])


def run_cands(f, sample, tag, algo, wmin):
    """V0 list of a jv.run output (V0 class on) for one io.File as a candidate dict"""
    import glob
    out = {k: [] for k in ('ev', 't1', 't2', 'pdg', 'pos', 'mass')}
    rd = tag if '/' in tag else f'{tag}/{sample}'          # '<tag>/<run sample>' reads another sample's run (DEV9 in DEV)
    fs = sorted(glob.glob(f'{IO.env_path("JV_RUNS")}/{rd}/{algo}/{f.stem}_*.npz'))
    if not fs:
        raise FileNotFoundError(f'{tag}/{sample}/{algo}: no output for {f.stem}')
    for fn in fs:
        d = np.load(fn)
        m = d['v0_w'] > wmin
        out['ev'].append(d['v0_ev'][m].astype(np.int64) + int(d['a'])); out['t1'].append(d['v0_t1'][m].astype(np.int64))
        out['t2'].append(d['v0_t2'][m].astype(np.int64)); out['pdg'].append(d['v0_pdg'][m].astype(np.int64))
        out['pos'].append(d['v0_pos'][m].reshape(-1, 3)); out['mass'].append(d['v0_mass'][m])
    return {k: np.concatenate(v) for k, v in out.items()}


def conv_cands(f, sel):
    """conversion candidates of jv.conv.find_conversions (sel 'loose' or 'tight') for one io.File, computed here"""
    from . import conv as CV
    d, o = f.d, f.off['t']
    out = {k: [] for k in ('ev', 't1', 't2', 'pdg', 'pos', 'mass')}
    for e in range(f.nev):
        a, b = o[e], o[e + 1]
        C = CV.find_conversions(f.par_fit[a:b], f.cov_fit[a:b], d['t_q'][a:b], d['t_nvdet'][a:b], d['t_nitc'][a:b],
                                d['ev_pv'][e], sel=sel)
        n = len(C['i'])
        out['ev'].append(np.full(n, e, np.int64)); out['t1'].append(C['i'].astype(np.int64)); out['t2'].append(C['j'].astype(np.int64))
        out['pdg'].append(np.full(n, 22, np.int64)); out['pos'].append(C['x'].reshape(-1, 3))
        out['mass'].append(C['Q'][:, CV.QI['tp_mee']])
    return {k: np.concatenate(v) for k, v in out.items()}


def row_cands(f, r, sample):
    if r in ('stage1_loose', 'stage1_tight'):
        return stage1_cands(f, r == 'stage1_tight')
    if r in ('conv_loose', 'conv_tight'):
        return conv_cands(f, r[5:])
    if r.startswith('union:'):
        A, B = (row_cands(f, x, sample) for x in r[6:].split('+'))
        C = {k: np.concatenate([A[k], B[k]]) for k in A}
        key = np.stack([C['ev'], np.minimum(C['t1'], C['t2']), np.maximum(C['t1'], C['t2']), C['pdg']], 1)
        _, first = np.unique(key, axis=0, return_index=True)
        return {k: v[np.sort(first)] for k, v in C.items()}
    if r.startswith('run:'):
        x = r.split(':')
        return run_cands(f, sample, x[1], x[2], float(x[3]) if len(x) > 3 else 0.5)
    raise KeyError(r)


def truth_tables(f):
    """per true V0 (v0t) and per track: species, legs as global track indices, origin, kinematics; leg / legside = per
    global track the true-V0 index it is a daughter track of and which daughter (1 positive, 2 negative), -1 if none"""
    d, o = f.d, f.off
    if 't_v0t' not in d:
        raise KeyError(f'{f.stem}: no per-track t_v0t (inputs without the true-V0 arrays, jv/harvest.py)')
    tev = np.repeat(np.arange(len(d['ev_nv0t'])), d['ev_nv0t'])
    g1 = np.where(d['v0t_l1'] >= 0, d['v0t_l1'] + o['t'][tev], -1)
    g2 = np.where(d['v0t_l2'] >= 0, d['v0t_l2'] + o['t'][tev], -1)
    trev = np.repeat(np.arange(len(d['ev_ntrk'])), d['ev_ntrk'])
    lt = d['t_v0t'].astype(np.int64)
    leg = np.where(lt >= 0, lt + o['v0t'][trev], -1)          # global track -> true V0 index (daughter track of)
    legside = np.where(lt >= 0, d['t_v0leg'].astype(np.int64), -1)
    p = np.linalg.norm(d['v0t_p3'], axis=1)
    return dict(ev=tev, g1=g1, g2=g2, spec=species_of(d['v0t_pdg']), orig=d['v0t_orig'].astype(np.int64), p=p,
                rxy=d['v0t_rxy'], abscos=np.abs(d['v0t_p3'][:, 2]) / np.maximum(p, 1e-9), leg=leg, legside=legside,
                cls=d['t_cls'], pos=d['v0t_pos'], p3=d['v0t_p3'])


def classify(f, T, C):
    """per candidate: class index (CLASSES), matched true-V0 index (-1), fitted p, rxy, |cos theta|"""
    o = f.off['t']
    n = len(C['ev'])
    cls = np.full(n, CLASSES.index('other'), np.int64)
    mt = np.full(n, -1, np.int64)
    ok = (C['t1'] >= 0) & (C['t2'] >= 0)
    g1 = np.where(ok, C['t1'] + o[C['ev']], 0)
    g2 = np.where(ok, C['t2'] + o[C['ev']], 0)
    l1 = np.where(ok, T['leg'][g1], -1); l2 = np.where(ok, T['leg'][g2], -1)
    s1 = np.where(ok, T['legside'][g1], -1); s2 = np.where(ok, T['legside'][g2], -1)
    c1 = np.where(ok, T['cls'][g1], 5); c2 = np.where(ok, T['cls'][g2], 5)
    samedau = (l1 >= 0) & (l1 == l2) & (s1 == s2)           # both tracks linked to the same daughter: not a V0 pair
    same = (l1 >= 0) & (l1 == l2) & ~samedau
    sp = np.where(l1 >= 0, T['spec'][np.maximum(l1, 0)], 0)
    hf1, hf2 = np.isin(c1, (1, 2)), np.isin(c2, (1, 2))
    one = (l1 >= 0) != (l2 >= 0)
    other_c = np.where(l1 >= 0, c2, c1)
    rules = [  # first match wins, from the most specific
        ('matched', same & (sp == C['pdg'])),
        ('wrong_species', same & (sp != C['pdg'])),
        ('two_v0s', (l1 >= 0) & (l2 >= 0) & (l1 != l2)),
        ('one_leg_unlinked', one & (other_c == 5)),
        ('one_leg_linked', one | samedau),
        ('hf_pair', hf1 & hf2),
        ('primary_pair', (c1 == 0) & (c2 == 0)),
        ('other_secondary', ((c1 == 4) & (c2 == 4)) | ((c1 == 4) & hf2) | (hf1 & (c2 == 4))),
    ]
    done = np.zeros(n, bool)
    for name, m in rules:
        m = m & ok & ~done
        cls[m] = CLASSES.index(name)
        done |= m
    mt[cls == 0] = l1[cls == 0]
    # fitted kinematics
    P = np.zeros((n, 3))
    if n and ok.any():
        P[ok] = track_momenta(f.par_fit[g1[ok]], C['pos'][ok]) + track_momenta(f.par_fit[g2[ok]], C['pos'][ok])
    pm = np.linalg.norm(P, axis=1)
    pv = f.d['ev_pv'][C['ev']] if n else np.zeros((0, 3))
    rxy = np.hypot(C['pos'][:, 0] - pv[:, 0], C['pos'][:, 1] - pv[:, 1]) if n else np.zeros(0)
    return dict(cls=cls, mt=mt, p=pm, rxy=rxy, abscos=np.abs(P[:, 2]) / np.maximum(pm, 1e-9), pdg=C['pdg'], ev=C['ev'],
                samedau=samedau & ok)


def file_job(args):
    stem, lim, rownames, sample = args
    f = IO.File(stem)
    T = truth_tables(f)
    if lim is not None:                    # FAST-like subsets: first lim events
        f.nev = min(f.nev, lim)
    keep = T['ev'] < f.nev
    out = {'stem': stem, 'nev': f.nev, 'truth': {k: v[keep] for k, v in T.items() if k not in ('leg', 'legside', 'cls')}}
    out['truth']['n1'] = f.d['v0t_n1'][keep]; out['truth']['n2'] = f.d['v0t_n2'][keep]
    for r in rownames:
        C = row_cands(f, r, sample)
        C = {k: v[C['ev'] < f.nev] for k, v in C.items()}
        K = classify(f, T, C)
        found = np.zeros(len(T['ev']), np.int64)
        np.add.at(found, K['mt'][K['mt'] >= 0], 1)
        out[r] = dict(cand=K, found=found[keep])
    return out


def merge(parts):
    ev0 = 0
    M = {'truth': {}, 'rows': {}}
    for p in parts:
        for k, v in p['truth'].items():
            M['truth'].setdefault(k, []).append(v + ev0 if k == 'ev' else v)
        for r in p:
            if r in ('stem', 'nev', 'truth'):
                continue
            R = M['rows'].setdefault(r, {'found': [], 'cand': {}})
            R['found'].append(p[r]['found'])
            for k, v in p[r]['cand'].items():
                R['cand'].setdefault(k, []).append(v + ev0 if k == 'ev' else v)
        ev0 += p['nev']
    M['nev'] = ev0
    M['truth'] = {k: np.concatenate(v) for k, v in M['truth'].items()}
    for R in M['rows'].values():
        R['found'] = np.concatenate(R['found'])
        R['cand'] = {k: np.concatenate(v) for k, v in R['cand'].items()}
    return M


def binom(k, n):
    e = k / n if n > 0 else np.nan
    return [float(e), float(np.sqrt(e * (1 - e) / n)) if n > 0 else None, int(n)]


def eff_tables(T, found):
    out = {}
    tracked = (T['g1'] >= 0) & (T['g2'] >= 0)
    hit = found > 0
    for sn, sv in SPECIES.items():
        s = T['spec'] == sv
        o = {}
        for dn, den in (('all', s), ('tracked', s & tracked)):
            od = {'total': binom(int(hit[den].sum()), int(den.sum()))}
            for vn, bins in VARS:
                v = T['p' if vn == 'p_GeV' else ('rxy' if vn == 'rxy_cm' else 'abscos')]
                od['vs_' + vn] = {'bins': bins, 'eff': [binom(int(hit[den & (v >= a) & (v < b)].sum()),
                                                              int((den & (v >= a) & (v < b)).sum()))
                                                        for a, b in zip(bins[:-1], bins[1:])]}
            od['by_origin'] = {on: binom(int(hit[den & (T['orig'] == k)].sum()), int((den & (T['orig'] == k)).sum()))
                               for k, on in enumerate(ORIG)}
            od['rxy_by_origin'] = {on: [binom(int(hit[den & (T['orig'] == k) & (T['rxy'] >= a) & (T['rxy'] < b)].sum()),
                                              int((den & (T['orig'] == k) & (T['rxy'] >= a) & (T['rxy'] < b)).sum()))
                                        for a, b in zip(R_BINS[:-1], R_BINS[1:])] for k, on in enumerate(ORIG)}
            o[dn] = od
        o['duplicates'] = int((found[s] > 1).sum())
        out[sn] = o
    return out


def pur_tables(K):
    out = {}
    for sn, sv in SPECIES.items():
        s = K['pdg'] == sv
        n = int(s.sum())
        o = {'n_cand': n, 'purity': binom(int((s & (K['cls'] == 0)).sum()), n),
             'classes': {c: binom(int((s & (K['cls'] == i)).sum()), n) for i, c in enumerate(CLASSES)},
             'groups': {g: binom(int((s & np.isin(K['cls'], [CLASSES.index(c) for c in cs])).sum()), n)
                        for g, cs in GROUPS.items()}}
        for vn, bins in VARS:
            v = K['p' if vn == 'p_GeV' else ('rxy' if vn == 'rxy_cm' else 'abscos')]
            o['vs_' + vn] = {'bins': bins, 'purity': [binom(int((s & (K['cls'] == 0) & (v >= a) & (v < b)).sum()),
                                                            int((s & (v >= a) & (v < b)).sum()))
                                                      for a, b in zip(bins[:-1], bins[1:])],
                             'classes': {c: [int((s & (K['cls'] == i) & (v >= a) & (v < b)).sum())
                                             for a, b in zip(bins[:-1], bins[1:])] for i, c in enumerate(CLASSES)}}
        out[sn] = o
    return out


def paired(M, row, ref, nboot, seed=12345):
    """paired event bootstrap of row - ref: efficiencies (both denominators) and purity per species"""
    rng = np.random.default_rng(seed)
    T, nev = M['truth'], M['nev']
    tracked = (T['g1'] >= 0) & (T['g2'] >= 0)
    A, B = M['rows'][row], M['rows'][ref]
    items = {}

    def pe(ev, w):
        return np.bincount(ev, weights=np.asarray(w, float), minlength=nev)
    for sn, sv in SPECIES.items():
        s = T['spec'] == sv
        for dn, den in (('all', s), ('tracked', s & tracked)):
            nd = pe(T['ev'][den], np.ones(den.sum()))
            items[f'eff_{sn}_{dn}'] = (pe(T['ev'][den], A['found'][den] > 0), nd, pe(T['ev'][den], B['found'][den] > 0), nd)
        for X in (A, B):
            k = X['cand']
            m = k['pdg'] == sv
            items.setdefault(f'pur_{sn}', []).extend([pe(k['ev'][m], k['cls'][m] == 0), pe(k['ev'][m], np.ones(m.sum()))])
    items = {k: tuple(v) for k, v in items.items()}

    def rat(n, d, w=None):
        ds = d.sum() if w is None else d @ w
        return (n.sum() if w is None else n @ w) / ds if ds > 0 else np.nan
    d0 = {k: rat(v[0], v[1]) - rat(v[2], v[3]) for k, v in items.items()}
    bs = {k: [] for k in items}
    for _ in range(nboot):
        w = np.bincount(rng.integers(0, nev, nev), minlength=nev).astype(float)
        for k, v in items.items():
            bs[k].append(rat(v[0], v[1], w) - rat(v[2], v[3], w))
    return {k: [float(d0[k]), float(np.nanstd(bs[k]))] for k in items}


def evaluate(sample, rownames, nproc, nboot=100, ref=None):
    units = []
    for s in IO.SAMPLES[sample]:
        stem, lim = (s.split(':')[0], int(s.split(':')[1])) if ':' in s else (s, None)
        units.append((stem, lim, rownames, sample.split(':')[0]))
    with Pool(min(nproc, len(units))) as p:
        parts = p.map(file_job, units)
    M = merge(parts)
    T = M['truth']
    res = {'sample': sample, 'n_events': M['nev'], 'definitions': __doc__, 'rows': {}}
    tracked = (T['g1'] >= 0) & (T['g2'] >= 0)
    res['truth'] = {sn: dict(n_all=int((T['spec'] == sv).sum()), n_tracked=int(((T['spec'] == sv) & tracked).sum()),
                             per_event=float((T['spec'] == sv).sum() / M['nev']),
                             by_origin={on: int(((T['spec'] == sv) & (T['orig'] == k)).sum()) for k, on in enumerate(ORIG)})
                    for sn, sv in SPECIES.items()}
    for r in rownames:
        R = M['rows'][r]
        res['rows'][r] = dict(n_cand=int(len(R['cand']['ev'])), cand_per_event=float(len(R['cand']['ev']) / M['nev']),
                              n_cand_same_daughter=int(R['cand']['samedau'].sum()),
                              eff=eff_tables(T, R['found']), purity=pur_tables(R['cand']))
        if ref is not None and r != ref:
            res['rows'][r]['diff_vs_ref'] = paired(M, r, ref, nboot)
    res['ref'] = ref
    return res, M


def clean(x):
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, np.ndarray):
        return clean(x.tolist())
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, (float, np.floating)):
        return float(x) if np.isfinite(x) else None
    return x


def main():
    a = sys.argv[1:]
    nproc, sample, out = int(a[0]), a[1], a[2]
    opts = dict(zip(a[3::2], a[4::2]))
    rows = opts.get('--rows', 'stage1_loose,stage1_tight').split(',')
    res, _ = evaluate(sample, rows, nproc, int(opts.get('--nboot', 100)), opts.get('--ref', rows[0]))
    json.dump(clean(res), open(out, 'w'), indent=1)
    print('wrote', out)


if __name__ == '__main__':
    main()
