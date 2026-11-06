"""Run vertexing algorithms (hybrid.ROWS names, 'oracle', or a JSON config file) over a sample, multiprocessing.

Usage: python3 -m jv.run <nproc> <sample> <out_tag> <algo> [<algo> ...]
  sample : FAST | DEV | HELDOUT | DATA (jv.io.SAMPLES), optionally ':<maxev per unit>' for a smoke test (e.g. FAST:50)
  algo   : a hybrid.ROWS key (hybrid_v0, lcfiplus_like, mvf_like, ivf_like, jetproj_pv_svn, nominal), 'oracle' (MC only),
           or name=path.json (config overrides of hybrid.DEFAULT; a name starting with 'oracle' runs the oracle PV with
           the config's bs_sig)
  Unit boundaries depend on nproc (3 x nproc units, rounded up to a multiple of the file count); results per event do not.
Environment: JV_INPUTS (harvest npz directory, read), JV_RUNS (output root directory); both required, no default.
Output : $JV_RUNS/<out_tag>/<sample>/<algo>/<stem>_<a>_<b>.npz
  ev_pv (3), ev_pvc6 (xx yy zz yx zx zy), ev_nsv, ev_cpu (s, process time), ev_row; sv_ev (local event index in the unit),
  sv_pos, sv_c6, sv_mass, sv_ntrk; t_wpv, t_assign (-1 none, 0 PV, k>=1 SV k), t_wass, t_s3 (per baseline track, input
  order of the unit's events)
"""
import json
import os
import sys
import time
import traceback
import numpy as np
from multiprocessing import Pool

from . import io as IO
from . import hybrid as HY


def _c6(C):
    return [C[0, 0], C[1, 1], C[2, 2], C[1, 0], C[2, 0], C[2, 1]]


def unit(args):
    stem, a, b, algos, tag, sample, maxev = args
    try:
        f = IO.File(stem)
        HY.UPCACHE = {} if os.environ.get('JV_UPCACHE') else None
        if maxev:
            b = min(b, a + maxev)
        res = {}
        for name, cfg in algos:      # warm-up (numba compilation outside the CPU measurement)
            E0 = f.event(a)
            HY.oracle_event(E0, bs_sig=cfg['bs_sig'] if cfg else HY.F.BS_SIG) if name.startswith('oracle') else HY.run_event(E0, cfg)
        for name, cfg in algos:
            o = {k: [] for k in ('ev_pv', 'ev_pvc6', 'ev_nsv', 'ev_cpu', 'ev_row', 'sv_ev', 'sv_pos', 'sv_c6', 'sv_mass',
                                 'sv_ntrk', 't_wpv', 't_assign', 't_wass', 't_s3')}
            for e in range(a, b):
                E = f.event(e)
                try:
                    r = HY.oracle_event(E, bs_sig=cfg['bs_sig'] if cfg else HY.F.BS_SIG) if name.startswith('oracle') and not name.startswith('oracle_noprune') else (HY.oracle_event(E, prune=False, bs_sig=cfg['bs_sig'] if cfg else HY.F.BS_SIG) if name.startswith('oracle_noprune')
                                                                     else HY.run_event(E, cfg))
                except Exception:
                    raise RuntimeError(f'{name} {stem} event {e}: ' + traceback.format_exc()[-1500:])
                o['ev_pv'].append(r['pv']); o['ev_pvc6'].append(_c6(r['pvC'])); o['ev_nsv'].append(len(r['svs']))
                o['ev_cpu'].append(r['cpu']); o['ev_row'].append(e)
                for (x, C, m, nt) in r['svs']:
                    o['sv_ev'].append(e - a); o['sv_pos'].append(x); o['sv_c6'].append(_c6(C)); o['sv_mass'].append(m)
                    o['sv_ntrk'].append(nt)
                o['t_wpv'].append(r['wpv']); o['t_assign'].append(r['assign']); o['t_wass'].append(r['wass'])
                o['t_s3'].append(r['s3'])
            arr = {}
            for k, v in o.items():
                if k.startswith('t_'):
                    arr[k] = np.concatenate(v) if v else np.zeros(0)
                elif k in ('sv_pos', 'sv_c6'):
                    arr[k] = np.asarray(v, np.float64).reshape(-1, 3 if k == 'sv_pos' else 6)
                else:
                    arr[k] = np.asarray(v)
            arr['t_assign'] = arr['t_assign'].astype(np.int16)
            for k in ('t_wpv', 't_wass', 't_s3'):
                arr[k] = arr[k].astype(np.float32)
            d = f'{IO.env_path("JV_RUNS")}/{tag}/{sample}/{name}'
            os.makedirs(d, exist_ok=True)
            fn = f'{d}/{stem}_{a}_{b}.npz'
            np.savez_compressed(fn + '.tmp.npz', stem=stem, a=a, b=b, config=json.dumps(cfg if cfg else {}), **arr)
            os.replace(fn + '.tmp.npz', fn)
            res[name] = (b - a, float(np.sum(arr['ev_cpu'])))
        return stem, a, b, 'OK', res
    except Exception:
        return stem, a, b, 'FAILED', traceback.format_exc()[-3000:]


def parse_algos(names):
    out = []
    for nm in names:
        if '=' in nm:
            k, p = nm.split('=', 1)
            out.append((k, HY.config(json.load(open(p)))))
        elif nm in ('oracle', 'oracle_noprune'):
            out.append((nm, None))
        else:
            out.append((nm, HY.config(nm)))
    return out


if __name__ == '__main__':
    nproc, sample, tag = int(sys.argv[1]), sys.argv[2], sys.argv[3]
    maxev = None
    if ':' in sample:
        sample, maxev = sample.split(':')[0], int(sample.split(':')[1])
    algos = parse_algos(sys.argv[4:])
    IO.env_path('JV_RUNS')
    units = IO.sample_units(sample, nchunk=max(1, -(-3 * nproc // len(IO.SAMPLES[sample]))))
    t0 = time.time()
    tot = {}
    nfail = 0
    with Pool(nproc) as p:
        for r in p.imap_unordered(unit, [(s, a, b, algos, tag, sample, maxev) for s, a, b in units]):
            if r[3] != 'OK':
                nfail += 1
                print(r, flush=True)
                continue
            for k, (n, c) in r[4].items():
                tot.setdefault(k, [0, 0.]); tot[k][0] += n; tot[k][1] += c
    for k, (n, c) in tot.items():
        print(f'CPU {k}: {n} events, {1e3 * c / max(n, 1):.2f} ms/event (process time)', flush=True)
    print(f'wall {time.time() - t0:.0f} s, units {len(units)}, failed {nfail}', flush=True)
    sys.exit(1 if nfail else 0)
