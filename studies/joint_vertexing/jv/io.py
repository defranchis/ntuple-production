"""Input access of the joint vertexer: samples, file lists, per-event views of the harvest npz files.

Directories are passed by environment variable (no default): JV_INPUTS = harvest npz files (<stem>.npz)."""
import os
import numpy as np

from . import tracks as TR


def env_path(name):
    """directory given by the environment variable name (required, no default)"""
    p = os.environ.get(name)
    if not p:
        raise RuntimeError(f'environment variable {name} is not set (see README.md)')
    return p


SAMPLES = {
    # tuning sample; FAST = first 20000 events of its first file
    'DEV': [f'mc_ZM4212.{i}' for i in range(41, 51)],
    'FAST': ['mc_ZM4212.41:20000'],
    'MID': [f'mc_ZM4212.{i}' for i in range(41, 45)],
    # training file of a refitted score and of the V0 tables (last DEV file) and DEV without it
    'TRAIN': ['mc_ZM4212.50'], 'DEV9': [f'mc_ZM4212.{i}' for i in range(41, 50)],
    # two full DEV files for regression runs against stored DEV outputs
    'REG2': ['mc_ZM4212.41', 'mc_ZM4212.42'],
    'DATA1': ['data_Y15224.10'],
    # untouched until the final integration (frozen config, evaluated once)
    'HELDOUT': [f'mc_ZM4212.{i}' for i in range(51, 61)],
    'DATA': [f'data_Y15224.{i}' for i in range(10, 20)],
}


def sample_units(name, nchunk=4):
    """work units (stem, first event, last event) of a sample: each file split into nchunk event ranges"""
    units = []
    for s in SAMPLES[name]:
        stem, lim = (s.split(':')[0], int(s.split(':')[1])) if ':' in s else (s, None)
        n = int(np.load(f'{env_path("JV_INPUTS")}/{stem}.npz')['ev_run'].shape[0]) if lim is None else lim
        edges = np.linspace(0, n, nchunk + 1).astype(int)
        units += [(stem, int(a), int(b)) for a, b in zip(edges[:-1], edges[1:]) if b > a]
    return units


class File:
    """lazy per-event views of one harvest npz"""

    def __init__(self, stem):
        d = np.load(f'{env_path("JV_INPUTS")}/{stem}.npz')
        self.d = {k: d[k] for k in d.files}
        self.stem = stem
        self.isMC = bool(self.d['isMC'])
        self.nev = len(self.d['ev_run'])
        self.off = {}
        for key, cnt in (('t', 'ev_ntrk'), ('v', 'ev_nv0'), ('s', 'ev_nsv'), ('st', 'ev_nst'), ('h', 'ev_nhf'), ('v0t', 'ev_nv0t')):
            if key == 'v0t' and cnt not in self.d:       # offsets of the true-V0 table (MC; jv.evaluate_v0)
                continue
            c = self.d[cnt].astype(np.int64)
            o = np.zeros(len(c) + 1, np.int64)
            np.cumsum(c, out=o[1:])
            self.off[key] = o
        self.par_fit, self.cov_fit = TR.to_fit(self.d['t_par'], self.d['t_cov'])

    def event(self, e):
        d, o = self.d, self.off
        E = {k[3:]: d[k][e] for k in d if k.startswith('ev_')}
        a, b = o['t'][e], o['t'][e + 1]
        E['t'] = {k[2:]: d[k][a:b] for k in d if k.startswith('t_')}
        E['t']['fpar'] = self.par_fit[a:b]
        E['t']['fcov'] = self.cov_fit[a:b]
        E['t_off'] = a
        for key in ('v', 's', 'st', 'h'):
            a2, b2 = o[key][e], o[key][e + 1]
            E[key] = {k[len(key) + 1:]: d[k][a2:b2] for k in d if k.startswith(key + '_') and not
                      (key == 's' and k.startswith('st_'))}
        E['isMC'] = self.isMC
        E['pvcov'] = cov6_to_3(E['pvcov6'])
        E['row'] = e
        return E


def cov6_to_3(c6):
    """(xx yy zz yx zx zy) -> 3x3"""
    xx, yy, zz, yx, zx, zy = c6
    return np.array([[xx, yx, zx], [yx, yy, zy], [zx, zy, zz]], np.float64)


def jet_axes(E):
    th, ph = E['jet_theta'], E['jet_phi']
    ax = np.array([[np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)] for t, p in zip(th, ph)])
    if not np.all(np.isfinite(ax)):
        t = E['thrust']
        ax = np.array([t, -t])
    return ax
