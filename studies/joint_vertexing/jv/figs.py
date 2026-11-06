"""Standard figures of the joint-vertexing evaluator from a rows JSON (plot style: matplotlib defaults,
DejaVu Sans, ~4.5 x 4.05 in panels, tab10 colours, ratio panel when two things are compared, one plot per PNG + PDF).
The legend sits in head room above the data. Usage: python3 -m jv.figs <rows.json> <outdir> [<closure.json>]
(closure.json: output of a fitter-closure check, not included).
The command line draws the rows named in COL / LAB / ORDER and in the default `algos` of the functions (the comparison rows
nominal, hybrid_v0, oracle, lcfiplus_like, mvf_like, ivf_like, jetproj_*); other rows (e.g. final_LR) need entries in
COL / LAB and the functions called with `algos=`."""
import json
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({k: v for k, v in matplotlib.rcParamsDefault.items() if k != 'backend'})
SIZE = (4.5, 4.05)
DPI = 980 / 4.5
COL = {'nominal': 'k', 'hybrid_v0': 'tab:red', 'oracle': 'tab:green', 'oracle_noprune': 'tab:olive', 'lcfiplus_like': 'tab:blue',
       'mvf_like': 'tab:orange', 'ivf_like': 'tab:purple', 'jetproj_pv_svn': 'tab:brown', 'jetproj_soft_svn': 'tab:pink'}
LAB = {'nominal': 'nominal (stage1 PV + svn)', 'hybrid_v0': 'new algorithm v0', 'oracle': 'truth-primary oracle PV',
       'oracle_noprune': 'oracle PV, no pruning', 'lcfiplus_like': 'LCFIPlus-like', 'mvf_like': 'MVF-like',
       'ivf_like': 'IVF-like', 'jetproj_pv_svn': 'jet-proj. PV hard + svn', 'jetproj_soft_svn': 'jet-proj. PV soft + svn'}
ORDER = ['nominal', 'hybrid_v0', 'oracle', 'lcfiplus_like', 'mvf_like', 'ivf_like', 'jetproj_pv_svn', 'jetproj_soft_svn', 'oracle_noprune']


def headroom(ax, f=1.45, log=False):
    """extend the y axis so that a legend at the top stays clear of the data"""
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi * f if log else lo + (hi - lo) * f)


def save(fig, out, name):
    fig.savefig(f'{out}/{name}.png', dpi=DPI)
    fig.savefig(f'{out}/{name}.pdf')
    plt.close(fig)


def ratio_fig():
    fig, (ax, axr) = plt.subplots(2, 1, sharex=True, figsize=SIZE, gridspec_kw=dict(height_ratios=[3, 1]))
    fig.subplots_adjust(left=0.18, right=0.96, top=0.90, bottom=0.12, hspace=0.06)
    return fig, ax, axr


def pv_resid(R, out, fl, c, algos=('nominal', 'hybrid_v0', 'oracle')):
    fig, ax, axr = ratio_fig()
    ref = None
    for al in algos:
        o = R['algos'].get(al, {}).get('pv', {}).get(fl, {}).get(c)
        if not o:
            continue
        b = np.array(o['hist_bins']); h = np.array(o['hist'], float)
        n = h / h.sum() / np.diff(b)
        ax.stairs(n, b, color=COL[al], label=f"{LAB[al]} (core {o['core']:.0f} um)")
        if al == 'nominal':
            ref = (b, n, h)
        elif ref is not None:
            with np.errstate(divide='ignore', invalid='ignore'):
                axr.stairs(n / ref[1], b, color=COL[al])
    ax.set_title(f'PV residual {c}, {fl} events ({R["sample"]})', fontsize=10)
    ax.set_ylabel('normalised density [1/um]')
    ax.set_yscale('log')
    headroom(ax, 300, True)
    ax.legend(fontsize=7, loc='upper center')
    axr.axhline(1, color='r', ls='--', lw=0.8)
    axr.set_ylim(0, 2); axr.set_ylabel('ratio to\nnominal', fontsize=8)
    axr.set_xlabel(f'{c}(PV) - {c}(generated) [um]')
    save(fig, out, f'pv_resid_{fl}_{c}')


def diff_points(R, out, keyfmt, title, ylabel, name, flav=('b', 'c', 'uds')):
    fig, ax = plt.subplots(figsize=SIZE)
    fig.subplots_adjust(left=0.18, right=0.96, top=0.90, bottom=0.30)
    algs = [a for a in ORDER if a in R['algos'] and a != R['ref']]
    for k, fl in enumerate(flav):
        xs, ys, es = [], [], []
        for i, al in enumerate(algs):
            d = R['algos'][al].get('diff_vs_ref', {}).get(keyfmt.format(fl=fl))
            if d is None or d[0] is None:
                continue
            xs.append(i + 0.2 * (k - 1)); ys.append(d[0]); es.append(d[1] or 0)
        ax.errorbar(xs, ys, yerr=es, fmt='o', ms=4, color=f'C{k}', label=f'{fl} events')
    ax.axhline(0, color='k', lw=0.8)
    ax.set_xticks(range(len(algs))); ax.set_xticklabels([LAB[a] for a in algs], rotation=45, ha='right', fontsize=7)
    ax.set_ylabel(ylabel); ax.set_title(title, fontsize=10)
    headroom(ax, 1.3)
    ax.legend(fontsize=7, loc='upper center', ncol=3)
    save(fig, out, name)


def eff_vs(R, out, var, xlabel, algos=('nominal', 'hybrid_v0', 'lcfiplus_like', 'mvf_like', 'ivf_like'), logx=False):
    fig, ax, axr = ratio_fig()
    ref = None
    for al in algos:
        o = R['algos'].get(al, {}).get('sv', {}).get(f'eff_vs_{var}')
        if not o:
            continue
        b = np.array(o['bins'], float); e = np.array([np.nan if v is None else v for v in o['eff']]); n = np.array(o['n'], float)
        x = np.sqrt(b[:-1] * b[1:]) if logx else 0.5 * (b[:-1] + b[1:])
        err = np.sqrt(np.maximum(e * (1 - e), 0) / np.maximum(n, 1))
        ax.errorbar(x, e, yerr=err, fmt='o-', ms=3, lw=1, color=COL[al], label=LAB[al])
        if al == 'nominal':
            ref = e
        elif ref is not None:
            axr.plot(x, e / ref, 'o-', ms=3, lw=1, color=COL[al])
    if logx:
        ax.set_xscale('log')
    ax.set_ylim(0, 1); ax.set_ylabel('SV efficiency')
    ymax = np.nanmax([np.nanmax(np.array([np.nan if v is None else v for v in R['algos'][a]['sv'][f'eff_vs_{var}']['eff']], float))
                      for a in algos if a in R['algos']])
    ax.set_ylim(0, min(1.0, ymax + 0.35))
    ax.set_title(f'SV efficiency vs {xlabel} ({R["sample"]})', fontsize=10)
    ax.legend(fontsize=7, loc='upper left', ncol=2)
    axr.axhline(1, color='r', ls='--', lw=0.8); axr.set_ylim(0, 2); axr.set_ylabel('ratio to\nnominal', fontsize=8)
    axr.set_xlabel(xlabel)
    save(fig, out, f'sv_eff_vs_{var}')


def migration(R, out):
    fig, ax = plt.subplots(figsize=SIZE)
    fig.subplots_adjust(left=0.18, right=0.96, top=0.90, bottom=0.30)
    algs = [a for a in ORDER if a in R['algos'] and 'migration' in R['algos'][a]]
    rows = R['algos'][algs[0]]['migration']['rows']
    for k, cl in enumerate(('primary', 'b', 'c', 'V0')):
        r = rows.index(cl)
        ax.plot(range(len(algs)), [R['algos'][a]['migration']['frac_row'][r][0] for a in algs], 'o', color=f'C{k}',
                label=f'truth {cl} tracks')
    ax.set_xticks(range(len(algs))); ax.set_xticklabels([LAB[a] for a in algs], rotation=45, ha='right', fontsize=7)
    ax.set_ylabel('fraction assigned to the PV (w > 0.5)'); ax.set_ylim(0, 1.3); ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.set_title(f'Track migration into the PV ({R["sample"]})', fontsize=10)
    ax.legend(fontsize=7, loc='upper center', ncol=2)
    save(fig, out, 'migration_pv')


def datamc(R, out, key, xlabel, logx=False, algo='hybrid_v0'):
    if 'data' not in R or algo not in R['data']:
        return
    fig, ax, axr = ratio_fig()
    m = R['algos'][algo]['dists'][key]; d = R['data'][algo]['dists'][key]
    b = np.array(m['bins']); hm = np.array(m['counts'], float); hd = np.array(d['counts'], float)
    nm, nd = hm / max(hm.sum(), 1), hd / max(hd.sum(), 1)
    ax.stairs(nm, b, color='tab:blue', label=f'MC ({R["sample"]}), mean {m["mean"]:.3g}')
    x = 0.5 * (b[:-1] + b[1:]) if not logx else np.sqrt(b[:-1] * b[1:])
    ax.errorbar(x, nd, yerr=np.sqrt(hd) / max(hd.sum(), 1), fmt='o', ms=3, color='k', label=f'data, mean {d["mean"]:.3g}')
    if algo != 'nominal' and 'nominal' in R['data']:
        mn = np.array(R['algos']['nominal']['dists'][key]['counts'], float); dn = np.array(R['data']['nominal']['dists'][key]['counts'], float)
        ax.stairs(mn / max(mn.sum(), 1), b, color='tab:gray', ls='--', label='nominal, MC')
    with np.errstate(divide='ignore', invalid='ignore'):
        axr.errorbar(x, nd / nm, yerr=np.sqrt(hd) / max(hd.sum(), 1) / nm, fmt='o', ms=3, color='k')
    if logx:
        ax.set_xscale('log')
    ax.set_ylabel('fraction'); ax.set_title(f'{xlabel.split(" [")[0]}: data vs MC, {LAB[algo]}', fontsize=9)
    headroom(ax, 1.45)
    ax.legend(fontsize=7, loc='upper right')
    axr.axhline(1, color='r', ls='--', lw=0.8); axr.set_ylim(0.5, 1.5); axr.set_ylabel('data/MC', fontsize=8)
    axr.set_xlabel(xlabel)
    save(fig, out, f'datamc_{key}_{algo}')


def binned(R, out, key, field, ylabel, name, xlabel, algos=('nominal', 'hybrid_v0', 'lcfiplus_like', 'mvf_like', 'ivf_like'), logx=False,
           ratio=True, ymax=None):
    """a binned SV quantity (purity, completeness, residual) vs a variable, with a ratio to nominal"""
    fig, ax, axr = ratio_fig()
    ref = None
    top = 0.
    for al in algos:
        o = R['algos'].get(al, {}).get('sv', {}).get(key)
        if not o:
            continue
        b = np.array(o['bins'], float); v = np.array([np.nan if x is None else x for x in o[field]], float)
        x = np.sqrt(b[:-1] * b[1:]) if logx else 0.5 * (b[:-1] + b[1:])
        ax.plot(x, v, 'o-', ms=3, lw=1, color=COL[al], label=LAB[al])
        top = max(top, np.nanmax(v))
        if al == 'nominal':
            ref = v
        elif ref is not None:
            axr.plot(x, v / ref, 'o-', ms=3, lw=1, color=COL[al])
    if logx:
        ax.set_xscale('log')
    ax.set_ylim(0, ymax if ymax else top * 1.5)
    ax.set_ylabel(ylabel); ax.set_title(f'{ylabel.split(" (")[0]} vs {xlabel.split(" [")[0]} ({R["sample"]})', fontsize=10)
    if ymax:
        ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.legend(fontsize=7, loc='upper left', ncol=2)
    axr.axhline(1, color='r', ls='--', lw=0.8); axr.set_ylim(0, 2); axr.set_ylabel('ratio to\nnominal', fontsize=8)
    axr.set_xlabel(xlabel)
    save(fig, out, name)


def weights(R, out, algo='hybrid_v0'):
    m = R['algos'].get(algo, {}).get('migration')
    if not m or 'hist_wpv' not in m:
        return
    fig, ax = plt.subplots(figsize=SIZE)
    fig.subplots_adjust(left=0.18, right=0.96, top=0.90, bottom=0.12)
    b = np.array(m['weight_bins'])
    for k, cl in enumerate(('primary', 'b', 'c', 'V0')):
        h = np.array(m['hist_wpv'][cl], float)
        ax.stairs(h / max(h.sum(), 1), b, color=f'C{k}', label=f'truth {cl} tracks')
    ax.set_yscale('log'); headroom(ax, 30, True)
    ax.set_xlabel('PV weight per track'); ax.set_ylabel('fraction of the tracks of the class')
    ax.set_title(f'PV track weight, {LAB[algo]} ({R["sample"]})', fontsize=10)
    ax.legend(fontsize=7, loc='upper center', ncol=2)
    save(fig, out, f'weights_pv_{algo}')


def closure(path, out):
    c = json.load(open(path))
    b = np.array(c['hist_bins_um'])
    for key, nm, tt in (('hist_pv_abs_um', 'closure_pv', 'stage1 PV'), ('hist_sv_abs_um', 'closure_sv', 'svn')):
        fig, ax = plt.subplots(figsize=SIZE)
        fig.subplots_adjust(left=0.18, right=0.96, top=0.90, bottom=0.12)
        for k, cc in enumerate('xyz'):
            ax.stairs(np.array(c[key][cc], float), b, color=f'C{k}', label=cc)
        ax.set_xscale('log'); ax.set_xlim(1e-4, 1e5); ax.set_yscale('symlog', linthresh=1); headroom(ax, 30, True)
        ax.set_xlabel('|refit - stored position| [um]'); ax.set_ylabel('vertices')
        ax.set_title(f'Fitter closure on {tt} track lists ({c["sample"]})', fontsize=10)
        ax.legend(fontsize=7, loc='upper right')
        save(fig, out, nm)


def main():
    R = json.load(open(sys.argv[1]))
    out = sys.argv[2]
    os.makedirs(out, exist_ok=True)
    if len(sys.argv) > 3:
        closure(sys.argv[3], out)
    weights(R, out)
    for var, xl, lg in (('fd_cm', 'flight distance [cm]', True), ('p_GeV', 'momentum [GeV]', False), ('abscos', '|cos theta|', False)):
        binned(R, out, f'purity_vs_{var}', 'purity', 'SV purity (HF fraction)', f'sv_purity_vs_{var}', 'reconstructed ' + xl, logx=lg, ymax=1.5)
        binned(R, out, f'eff_vs_{var}', 'completeness', 'track completeness', f'sv_completeness_vs_{var}', 'true ' + xl, logx=lg, ymax=1.5)
        binned(R, out, f'eff_vs_{var}', 'resid_median_um', 'median |SV - true| [um]', f'sv_resid_vs_{var}', 'true ' + xl, logx=lg)
    for fl in ('b', 'c', 'uds'):
        for c in 'xyz':
            pv_resid(R, out, fl, c)
    diff_points(R, out, 'pv_{fl}_proj_median', 'PV bias along the HF direction, vs nominal', 'change of median r.u [um]', 'd_pv_proj_median')
    diff_points(R, out, 'pv_{fl}_proj_tmean', 'PV bias along the HF direction (5% trimmed mean), vs nominal', 'change of trimmed mean r.u [um]', 'd_pv_proj_tmean')
    for c in 'xyz':
        diff_points(R, out, 'pv_{fl}_%s_core' % c, f'PV core width {c}, vs nominal', f'change of core width {c} [um]', f'd_pv_core_{c}')
    eff_vs(R, out, 'fd_cm', 'true flight distance [cm]', logx=True)
    eff_vs(R, out, 'p_GeV', 'true hadron momentum [GeV]')
    eff_vs(R, out, 'abscos', 'true hadron |cos theta|')
    migration(R, out)
    for algo in ('hybrid_v0', 'nominal'):
        for key, xl, lg in (('n_pv_tracks', 'tracks in the PV (w > 0.5)', False), ('pv_weight', 'PV track weight', False),
                            ('n_sv', 'SVs per event', False), ('sv_mass', 'SV mass [GeV]', False),
                            ('sv_fd_cm', 'SV flight distance [cm]', True), ('sv_fsig', 'SV flight significance', True),
                            ('trk_ipsig', 'track signed 3D IP significance', False)):
            datamc(R, out, key, xl, lg, algo)
    print('figures in', out)


if __name__ == '__main__':
    main()
