"""Joint vertexer: per-event harvest of one input file (stage1 output + the raw EDM4hep file it was produced from).

Inputs (read-only; directories by environment variable, required, no default):
  stage1 : $JV_STAGE1/<stem>.root   (stage1 output: PV, V0s, svn SVs, jets, track membership bits)
  raw    : $JV_RAW/<stem>.root      (EDM4hep input of stage1: track states + covariance, hits, MC; stage1 does not store
                                     the track covariance)
Joined by (run, event); the event list is stage1's (stage1 event selection applied). Output: one npz per stem in
  $JV_INPUTS/<stem>.npz
Layout (flat arrays + per-event counts; local indices refer to the event's baseline-track list, Tracks order):
  ev_*  per event: run, event, flav (5 b, 4 c, 0 uds/other: stage1 jetPID), pv_good, pv (stage1 Vertex_refit, cm),
        pvcov6 (xx yy zz yx zx zy), bs (beam spot, cm), gen (generated PV, MC), thrust (unit axis of the analysis
        particles, recomputed: stage1 does not store it), jet_theta/phi/p (2 exclusive Durham jets), ntrk, nv0, nsv, nhf
  t_*   per baseline track (trk_member bit 9): idx (Tracks index), par (raw ALEPH state D0 phi omega Z0 tanLambda),
        cov (15, lower triangle of the 5x5 of the raw state), mem (trk_member bits; bit0 = in stage1 PV), jet (0/1, -1),
        nvdet/nitc/ntpc (subdetectorHitNumbers 0/1/2), q (charge), p (GeV);
        MC: cls (0 primary, 1 from b, 2 from c not via b, 3 V0 (Ks/Lambda/conversion leg), 4 other secondary,
        5 unmatched), prod (true production vertex), hfv (local index of the weak-HF truth vertex the track comes from,
        -1), ptrue (MC momentum of the linked particle), v0t (local index in the event's v0t table of the true V0
        whose daughter is the track's largest-weight truth link or a continuation of it, see below; -1 if none; every
        such track, also when several tracks share the daughter), v0leg (1 positive, 2 negative daughter, -1 if none),
        lpdg (PDG code of the largest-weight link, 0 none), chain (continuation steps from the link to its chain head:
        0 = the link is no continuation, -1 = no link), cls_chain (every class judged on the chain head: e.g. a
        restarted primary electron is class 0, not 4)
  v_*   per stage1 V0 (v0n): t1, t2 (local), tight, pdg, pos, mass
  s_*   per stage1 SV (svn): pos (absolute = PV + svn_d*), cov6, mass, ntrk; st_sv/st_t: constituent (SV local, track local)
  h_*   per weak-HF truth vertex (every weakly decaying b/c hadron that flew, cascade-linked):
        pos (decay point = production vertex of its first daughter), isb, pdg, nch (charged long-lived particles produced
        at the decay point), nmatched (of them linked to a baseline track), parent (local HF index, -1), p3 (momentum),
        fd; MC: nmatched_chain (as nmatched, a daughter counted when a track links to it or to one of its continuations)
  v0t_* (MC, count ev_nv0t) per true V0 with a charged two-prong decay recorded in the MC (exactly two daughters:
        K0S -> pi+ pi-, Lambda / anti-Lambda -> p pi, photon -> e+ e- with the conversion point > 1 um from the
        generated PV): pdg (310, +-3122, 22), pos (decay point = production vertex of the first daughter), p3 (V0
        momentum), rxy (transverse distance of the decay point from the generated PV), prod (production vertex of the
        V0), l1 / l2 (local index of the baseline track linked to the positive / negative daughter as for t_v0t,
        highest link weight if several, -1 if none), n1 / n2 (number of such baseline tracks), lp1 / lp2 (MC momentum
        of the positive / negative daughter), orig (0 primary: produced at the generated PV, strong and
        electromagnetic decays there included; 1 from a weak b-hadron decay chain; 2 from a weak c-hadron chain not via
        b; 3 other secondary: produced away from the PV, e.g. hadronic interactions, decays of other long-lived
        particles), mc (MC index of the V0).
  Continuation chains (MC): GALEPH ends the MC record of a charged particle at some interactions and continues it as a
        NEW MC particle (continuation_parents: same PDG code, starts at the parent's end point, the parent's only
        daughter of that code there, every other daughter there a photon or none: bremsstrahlung or a bare restart);
        the chain head is the first particle of the chain (chain_heads, any length). A track whose largest-weight link
        is a continuation of a V0 daughter counts as that daughter (t_v0t, t_v0leg, v0t_l1 / l2, v0t_n1 / n2) and has
        t_cls = 3 and t_hfv = -1; without following the chains it would get the class of the link itself (1, 2, 4 or 0)
        and, for 1 / 2, its HF vertex. Every other t_cls / t_hfv entry is judged on the link itself.
Usage: python3 -m jv.harvest <nproc> <stem> [<stem> ...]   (<stem>:<n> = first n events only -> test_<stem>.npz)"""
import os
import sys
import time
import traceback
import numpy as np
from multiprocessing import Pool

from .io import env_path
from . import tracks as TR
TSN = '_Tracks_trackStates/_Tracks_trackStates.'
TOL = 1e-4   # 1 um: same-point tolerance of the truth vertices
WEAK_B = {511, 521, 531, 541, 5122, 5132, 5232, 5332, 5142, 5242, 5342, 5512, 5532, 5542, 5312, 5322}
WEAK_C = {411, 421, 431, 4122, 4132, 4232, 4332, 4412, 4422, 4432, 4312, 4322}
LONG_CH = np.array([11, 13, 211, 321, 2212, 3112, 3222, 3312, 3334])


def thrust_axis(px, py, pz):
    P = np.stack([px, py, pz], 1).astype(np.float64)
    if len(P) < 2:
        return np.array([0., 0., 1.])
    pm = np.linalg.norm(P, axis=1)
    best, bt = None, -1.
    starts = list(P[np.argsort(-pm)[:4]])
    for t in starts:
        t = t / max(np.linalg.norm(t), 1e-12)
        for _ in range(30):
            s = np.sign(P @ t)
            s[s == 0] = 1
            tn = (s[:, None] * P).sum(0)
            tn /= max(np.linalg.norm(tn), 1e-12)
            if np.allclose(tn, t, atol=1e-10):
                break
            t = tn
        T = np.abs(P @ t).sum()
        if T > bt:
            bt, best = T, t
    return best


def true_v0s(pdg, chg, vtx, p3, dbeg, dend, didx, par1, gen):
    """true V0s with a charged two-prong decay in the MC record (v0t_* in the docstring): list of
    (mc index, positive daughter, negative daughter, origin)"""
    apdg = np.abs(pdg)
    nmc = len(pdg)
    out = []
    for i in np.nonzero(np.isin(apdg, (310, 3122, 22)))[0]:
        if dend[i] - dbeg[i] != 2 or dbeg[i] + 1 >= len(didx):
            continue
        d0, d1 = didx[dbeg[i]], didx[dbeg[i] + 1]
        if not (0 <= d0 < nmc and 0 <= d1 < nmc) or chg[d0] * chg[d1] >= 0:
            continue
        a, a0, a1 = int(apdg[i]), int(apdg[d0]), int(apdg[d1])
        if a == 310 and not (a0 == 211 and a1 == 211):
            continue
        if a == 3122 and sorted((a0, a1)) != [211, 2212]:
            continue
        if a == 22 and not (a0 == 11 and a1 == 11 and np.sum((vtx[d0] - gen) ** 2) > TOL * TOL):
            continue
        dp, dn = (d0, d1) if chg[d0] > 0 else (d1, d0)
        org, j, steps = 0, par1[i], 0
        while j >= 0 and steps < 60:
            aj = int(apdg[j])
            if aj in WEAK_B:
                org = 1
                break
            if aj in WEAK_C and org == 0:
                org = 2
            j = par1[j]; steps += 1
        if org == 0 and np.sum((vtx[i] - gen) ** 2) > TOL * TOL:
            org = 3
        out.append((int(i), int(dp), int(dn), org))
    return out


def continuation_parents(pdg, chg, vtx, endp, dbeg, dend, didx, par1):
    """per MC particle the MC index of the particle it CONTINUES, -1 if none. GALEPH ends the record of a
    charged particle at some interactions (bremsstrahlung, also with no recorded product) and starts a NEW MC particle for
    the same physical particle there. Rule: particle i continues its first parent j if (1) same PDG code (same species
    and charge), charged; (2) i starts at j's end point (within TOL); (3) i is j's only daughter with that PDG code
    starting at that point; (4) every other daughter of j starting at that point is a photon (or there is none). Hadronic
    interactions (other hadrons, nuclear fragments produced there) and delta rays (a second e- of the same code) fail
    (3)-(4)."""
    nmc = len(pdg)
    cp = np.full(nmc, -1, np.int64)
    if nmc == 0:
        return cp
    j = par1
    ok = (j >= 0) & (chg != 0)
    jj = np.where(ok, j, 0)
    ok &= pdg == pdg[jj]
    ok &= np.sum((vtx - endp[jj]) ** 2, 1) < TOL * TOL
    for i in np.nonzero(ok)[0]:
        p_ = int(jj[i])
        nsame, good = 0, True
        for q in range(dbeg[p_], dend[p_]):
            if q >= len(didx):
                good = False
                break
            d = int(didx[q])
            if not (0 <= d < nmc) or np.sum((vtx[d] - endp[p_]) ** 2) >= TOL * TOL:
                continue
            if pdg[d] == pdg[p_]:
                nsame += 1
            elif pdg[d] != 22:
                good = False
        if good and nsame == 1:
            cp[i] = p_
    return cp


def chain_heads(cp):
    """per MC particle the head of its continuation chain (itself if it continues nothing) and the number
    of continuation steps to it (chains of any length)"""
    head = np.arange(len(cp), dtype=np.int64)
    nst = np.zeros(len(cp), np.int64)
    for i in np.nonzero(cp >= 0)[0]:
        h, s = i, 0
        while cp[h] >= 0 and s < 1000:
            h = cp[h]; s += 1
        head[i], nst[i] = h, s
    return head, nst


def job(stem):
    try:
        return _job(stem)
    except Exception:
        return stem, 'FAILED', traceback.format_exc()[-3000:]


def _job(stem):
    import uproot
    import awkward as ak
    t0 = time.time()
    maxev = None
    if ':' in stem:                      # test mode: <stem>:<maxev> -> test_<stem>.npz
        stem, maxev = stem.split(':')[0], int(stem.split(':')[1])
    isMC = stem.startswith('mc_')
    sbr = ['run_number', 'event_number', 'jetPID', 'pv_good', 'Vertex_refit_x', 'Vertex_refit_y', 'Vertex_refit_z',
           'Vertex_refit_cov_xx', 'Vertex_refit_cov_yy', 'Vertex_refit_cov_zz', 'Vertex_refit_cov_yx', 'Vertex_refit_cov_zx',
           'Vertex_refit_cov_zy', 'gen_vertex_x', 'gen_vertex_y', 'gen_vertex_z', 'trk_member', 'Beamspot_x_cm',
           'Beamspot_y_cm', 'Beamspot_z_cm', 'jet_p', 'jet_theta', 'jet_phi', 'pfcand_trackIdx', 'pfcand_px', 'pfcand_py',
           'pfcand_pz', 'v0n_trk1_origIdx', 'v0n_trk2_origIdx', 'v0n_tight', 'v0n_pdg', 'v0n_vx', 'v0n_vy', 'v0n_vz',
           'v0n_invM', 'v0n_trk1_q', 'svn_dx', 'svn_dy', 'svn_dz', 'svn_cov_xx', 'svn_cov_yy', 'svn_cov_zz', 'svn_cov_yx',
           'svn_cov_zx', 'svn_cov_zy', 'svn_mass', 'svn_ntracks', 'svn_trk_sv', 'svn_trk_origIdx']
    S = uproot.open(f'{env_path("JV_STAGE1")}/{stem}.root')['events'].arrays(sbr, library='ak')
    rbr = ['EventHeader/EventHeader.runNumber', 'EventHeader/EventHeader.eventNumber', 'Tracks/Tracks.trackStates_begin',
           TSN + 'D0', TSN + 'phi', TSN + 'omega', TSN + 'Z0', TSN + 'tanLambda', TSN + 'covMatrix.values[21]',
           'Tracks/Tracks.subdetectorHitNumbers_begin', 'Tracks/Tracks.subdetectorHitNumbers_end',
           '_Tracks_subdetectorHitNumbers']
    mcb = ['MCParticles/MCParticles.PDG', 'MCParticles/MCParticles.charge', 'MCParticles/MCParticles.vertex.x',
           'MCParticles/MCParticles.vertex.y', 'MCParticles/MCParticles.vertex.z', 'MCParticles/MCParticles.momentum.x',
           'MCParticles/MCParticles.momentum.y', 'MCParticles/MCParticles.momentum.z', 'MCParticles/MCParticles.parents_begin',
           'MCParticles/MCParticles.parents_end', 'MCParticles/MCParticles.daughters_begin',
           'MCParticles/MCParticles.daughters_end', '_MCParticles_parents/_MCParticles_parents.index',
           '_MCParticles_daughters/_MCParticles_daughters.index', 'trackMCLink/trackMCLink.weight',
           '_trackMCLink_from/_trackMCLink_from.index', '_trackMCLink_to/_trackMCLink_to.index',
           'MCParticles/MCParticles.endpoint.x', 'MCParticles/MCParticles.endpoint.y', 'MCParticles/MCParticles.endpoint.z']
    R = uproot.open(f'{env_path("JV_RAW")}/{stem}.root')['events'].arrays(rbr + (mcb if isMC else []), library='ak')
    sc = lambda a: ak.to_numpy(ak.firsts(a)) if a.ndim > 1 else ak.to_numpy(a)
    rkey = {k: i for i, k in enumerate(zip(ak.to_numpy(R['EventHeader/EventHeader.runNumber'][:, 0]).tolist(),
                                           ak.to_numpy(R['EventHeader/EventHeader.eventNumber'][:, 0]).tolist()))}
    srun, sevt = sc(S.run_number), sc(S.event_number)
    ri = np.array([rkey[k] for k in zip(srun.tolist(), sevt.tolist())])
    R = R[ri]
    if maxev is not None:
        S, R, srun, sevt, ri = S[:maxev], R[:maxev], srun[:maxev], sevt[:maxev], ri[:maxev]
    nev = len(S)
    L = {k: ak.to_list(S[k]) for k in sbr if k not in ('run_number', 'event_number')}
    tsb = R['Tracks/Tracks.trackStates_begin']
    st = {k: ak.to_list(R[TSN + k][tsb]) for k in ('D0', 'phi', 'omega', 'Z0', 'tanLambda')}
    cvl = ak.to_list(R[TSN + 'covMatrix.values[21]'][tsb])
    hbl = ak.to_list(R['Tracks/Tracks.subdetectorHitNumbers_begin'])
    hel = ak.to_list(R['Tracks/Tracks.subdetectorHitNumbers_end'])
    hnl = ak.to_list(R['_Tracks_subdetectorHitNumbers'])
    if isMC:
        M = {k.split('MCParticles.', 1)[1] if 'MCParticles/MCParticles.' in k else k: ak.to_list(R[k]) for k in mcb}
    E = {k: [] for k in ('run', 'event', 'flav', 'pv_good', 'pv', 'pvcov6', 'bs', 'gen', 'thrust', 'jet_theta', 'jet_phi',
                         'jet_p', 'ntrk', 'nv0', 'nsv', 'nst', 'nhf') + (('nv0t',) if isMC else ())}
    V0T = {k: [] for k in ('pdg', 'pos', 'p3', 'rxy', 'prod', 'l1', 'l2', 'n1', 'n2', 'lp1', 'lp2', 'orig', 'mc')}
    T = {k: [] for k in ('idx', 'par', 'cov', 'mem', 'jet', 'nvdet', 'nitc', 'ntpc', 'q', 'p', 'cls', 'prod', 'hfv', 'ptrue',
                         'v0t', 'v0leg') + (('lpdg', 'chain', 'cls_chain') if isMC else ())}
    V = {k: [] for k in ('t1', 't2', 'tight', 'pdg', 'pos', 'mass')}
    SV = {k: [] for k in ('pos', 'cov6', 'mass', 'ntrk')}
    ST = {k: [] for k in ('sv', 't')}
    H = {k: [] for k in ('pos', 'isb', 'pdg', 'nch', 'nmatched', 'parent', 'p3', 'fd') + (('nmatched_chain',) if isMC else ())}
    qchk = [0, 0]
    chain_stats = dict(hfv_chain_diff=0)    # tracks whose HF truth vertex would change if found from the chain head
    for e in range(nev):
        mem = np.array(L['trk_member'][e], np.int64)
        base = np.nonzero(mem & 512)[0]
        loc = np.full(len(mem), -1, np.int64)
        loc[base] = np.arange(len(base))
        pv = np.array([L['Vertex_refit_x'][e], L['Vertex_refit_y'][e], L['Vertex_refit_z'][e]], np.float64)
        E['run'].append(int(srun[e])); E['event'].append(int(sevt[e]))
        jp = L['jetPID'][e]
        E['flav'].append(int(jp) if jp in (4, 5) else 0)
        E['pv_good'].append(int(L['pv_good'][e]))
        E['pv'].append(pv)
        E['pvcov6'].append([L[f'Vertex_refit_cov_{c}'][e] for c in ('xx', 'yy', 'zz', 'yx', 'zx', 'zy')])
        E['bs'].append([L['Beamspot_x_cm'][e], L['Beamspot_y_cm'][e], L['Beamspot_z_cm'][e]])
        E['gen'].append([L['gen_vertex_x'][e], L['gen_vertex_y'][e], L['gen_vertex_z'][e]])
        px = np.concatenate([np.asarray(j, np.float64) for j in L['pfcand_px'][e]]) if len(L['pfcand_px'][e]) else np.zeros(0)
        py = np.concatenate([np.asarray(j, np.float64) for j in L['pfcand_py'][e]]) if len(L['pfcand_py'][e]) else np.zeros(0)
        pz = np.concatenate([np.asarray(j, np.float64) for j in L['pfcand_pz'][e]]) if len(L['pfcand_pz'][e]) else np.zeros(0)
        E['thrust'].append(thrust_axis(px, py, pz))
        th, ph, pj = list(L['jet_theta'][e]), list(L['jet_phi'][e]), list(L['jet_p'][e])
        while len(th) < 2:
            th.append(np.nan); ph.append(np.nan); pj.append(np.nan)
        E['jet_theta'].append(th[:2]); E['jet_phi'].append(ph[:2]); E['jet_p'].append(pj[:2])
        # tracks
        jet = np.full(len(mem), -1, np.int64)
        for j, lst in enumerate(L['pfcand_trackIdx'][e][:2]):
            for ti in lst:
                if 0 <= ti < len(mem):
                    jet[ti] = j
        hb, he, hn = hbl[e], hel[e], hnl[e]
        for o in base:
            T['idx'].append(o)
            T['par'].append([st['D0'][e][o], st['phi'][e][o], st['omega'][e][o], st['Z0'][e][o], st['tanLambda'][e][o]])
            T['cov'].append(cvl[e][o][:15])
            T['mem'].append(mem[o]); T['jet'].append(jet[o])
            nh = [hn[hb[o] + k] if hb[o] + k < he[o] else 0 for k in range(3)]
            T['nvdet'].append(nh[0]); T['nitc'].append(nh[1]); T['ntpc'].append(nh[2])
            om = st['omega'][e][o]
            T['q'].append(-1 if om > 0 else 1)   # = stage1 v0n_trk1_q (checked: q_match)
            T['p'].append(TR.PT_CONST / max(abs(om), 1e-9) * np.sqrt(1 + st['tanLambda'][e][o] ** 2))
        E['ntrk'].append(len(base))
        # stage1 V0s (charge convention check on leg 1)
        nv = 0
        for k in range(len(L['v0n_tight'][e])):
            a, b = L['v0n_trk1_origIdx'][e][k], L['v0n_trk2_origIdx'][e][k]
            la = loc[a] if 0 <= a < len(loc) else -1
            lb = loc[b] if 0 <= b < len(loc) else -1
            if la >= 0:
                qchk[0] += 1
                qchk[1] += int((-1 if st['omega'][e][a] > 0 else 1) == L['v0n_trk1_q'][e][k])
            V['t1'].append(la); V['t2'].append(lb); V['tight'].append(L['v0n_tight'][e][k]); V['pdg'].append(L['v0n_pdg'][e][k])
            V['pos'].append([L['v0n_vx'][e][k], L['v0n_vy'][e][k], L['v0n_vz'][e][k]]); V['mass'].append(L['v0n_invM'][e][k])
            nv += 1
        E['nv0'].append(nv)
        nsv = len(L['svn_mass'][e])
        for k in range(nsv):
            SV['pos'].append(pv + np.array([L['svn_dx'][e][k], L['svn_dy'][e][k], L['svn_dz'][e][k]]))
            SV['cov6'].append([L[f'svn_cov_{c}'][e][k] for c in ('xx', 'yy', 'zz', 'yx', 'zx', 'zy')])
            SV['mass'].append(L['svn_mass'][e][k]); SV['ntrk'].append(L['svn_ntracks'][e][k])
        nst = 0
        for s, o in zip(L['svn_trk_sv'][e], L['svn_trk_origIdx'][e]):
            ST['sv'].append(s); ST['t'].append(loc[o] if 0 <= o < len(loc) else -1); nst += 1
        E['nsv'].append(nsv); E['nst'].append(nst)
        # MC truth
        nhf = 0
        if isMC:
            pdg = np.array(M['PDG'][e], np.int64)
            apdg = np.abs(pdg)
            nmc = len(pdg)
            vtx = np.stack([M['vertex.' + c][e] for c in 'xyz'], 1).astype(np.float64) if nmc else np.zeros((0, 3))
            p3 = np.stack([M['momentum.' + c][e] for c in 'xyz'], 1).astype(np.float64) if nmc else np.zeros((0, 3))
            chg = np.asarray(M['charge'][e], np.float64)
            pbeg, pend, pidx = M['parents_begin'][e], M['parents_end'][e], M['_MCParticles_parents/_MCParticles_parents.index'][e]
            dbeg, dend, didx = M['daughters_begin'][e], M['daughters_end'][e], M['_MCParticles_daughters/_MCParticles_daughters.index'][e]
            par1 = np.array([pidx[pbeg[i]] if pend[i] > pbeg[i] and pbeg[i] < len(pidx) else -1 for i in range(nmc)], np.int64)
            gen = np.array(E['gen'][-1], np.float64)
            # GALEPH continuation chains (continuation_parents) -> chain head of every MC particle
            endp = np.stack([M['endpoint.' + c][e] for c in 'xyz'], 1).astype(np.float64) if nmc else np.zeros((0, 3))
            head, nstep = chain_heads(continuation_parents(pdg, chg, vtx, endp, dbeg, dend, didx, par1))
            # weak HF hadrons that flew
            hfl = []
            seen = []
            for i in range(nmc):
                a = int(apdg[i])
                if a not in WEAK_B and a not in WEAK_C:
                    continue
                if dend[i] <= dbeg[i] or dbeg[i] >= len(didx):
                    continue
                d0 = didx[dbeg[i]]
                if not (0 <= d0 < nmc):
                    continue
                dec = vtx[d0]
                if np.sum((dec - vtx[i]) ** 2) < TOL * TOL:
                    continue
                if any(sa == a and np.sum((sd - dec) ** 2) < TOL * TOL for sa, sd in seen):
                    continue
                seen.append((a, dec))
                hfl.append(i)
            hpos = np.array([vtx[par_d] for par_d in [didx[dbeg[i]] for i in hfl]]) if hfl else np.zeros((0, 3))
            # per-track LW link
            lw = np.full(len(mem), -1, np.int64)
            bw = np.full(len(mem), -1e30)
            for f_, t_, w_ in zip(M['_trackMCLink_from/_trackMCLink_from.index'][e], M['_trackMCLink_to/_trackMCLink_to.index'][e],
                                  M['trackMCLink/trackMCLink.weight'][e]):
                if 0 <= f_ < len(mem) and 0 <= t_ < nmc and w_ > bw[f_]:
                    bw[f_] = w_; lw[f_] = t_
            linked_mc = set(int(lw[o]) for o in base if lw[o] >= 0)
            linked_head = set(int(head[m]) for m in linked_mc)          # their chain heads (for t_cls_chain)
            # ancestry flags
            isb_anc = np.zeros(nmc, bool); isc_anc = np.zeros(nmc, bool)
            for i in linked_mc | linked_head:
                j, steps = par1[i], 0
                while j >= 0 and steps < 60:
                    aj = int(apdg[j])
                    if aj in WEAK_B:
                        isb_anc[i] = True
                    elif aj in WEAK_C:
                        isc_anc[i] = True
                    j = par1[j]; steps += 1
            hidx = {}
            for k, i in enumerate(hfl):
                dec = hpos[k]
                d2 = np.sum((vtx - dec) ** 2, 1)
                at = np.nonzero((d2 < TOL * TOL) & (np.abs(chg) > 0) & np.isin(apdg, LONG_CH))[0]
                H['pos'].append(dec); H['isb'].append(int(apdg[i]) in WEAK_B); H['pdg'].append(int(pdg[i]))
                H['nch'].append(len(at)); H['nmatched'].append(int(sum(int(m) in linked_mc for m in at)))
                H['nmatched_chain'].append(int(sum(int(m) in linked_head for m in at)))
                H['p3'].append(p3[i]); H['fd'].append(float(np.linalg.norm(dec - gen)))
                pr = vtx[i]
                pa = -1
                for k2 in range(len(hfl)):
                    if k2 != k and np.sum((hpos[k2] - pr) ** 2) < TOL * TOL:
                        pa = k2; break
                H['parent'].append(pa)
                hidx[k] = i
            nhf = len(hfl)

            def _isv0(m_):
                pa_ = par1[m_]
                disp_ = np.linalg.norm(vtx[m_] - gen) > TOL
                return pa_ >= 0 and (int(apdg[pa_]) in (310, 3122) or (int(apdg[pa_]) == 22 and int(apdg[m_]) == 11 and disp_))

            def _cls(m_, isv0_):
                if isv0_:
                    return 3
                if isb_anc[m_]:
                    return 1
                if isc_anc[m_]:
                    return 2
                return 4 if np.linalg.norm(vtx[m_] - gen) > TOL else 0

            def _hfv(m_):
                hv_, j, steps = -1, m_, 0
                while j >= 0 and steps < 60 and hv_ < 0:
                    d2 = np.sum((hpos - vtx[j]) ** 2, 1)
                    kk = int(np.argmin(d2))
                    if d2[kk] < TOL * TOL:
                        hv_ = kk
                    j = par1[j]; steps += 1
                return hv_

            for o in base:
                m = lw[o]
                if m < 0:
                    T['cls'].append(5); T['prod'].append([np.nan] * 3); T['hfv'].append(-1); T['ptrue'].append(np.nan)
                    T['lpdg'].append(0); T['chain'].append(-1); T['cls_chain'].append(5)
                    continue
                pr = vtx[m]
                mh = int(head[m])
                # class 3 (V0 daughter) judged on the chain head; classes 1, 2, 4, 0 on the link
                c = _cls(m, _isv0(mh))
                hv = _hfv(m) if c in (1, 2) and nhf else -1
                ch = _cls(mh, _isv0(mh))             # t_cls_chain: every class judged on the chain head
                if ch in (1, 2) and nhf and _hfv(mh) != hv:
                    chain_stats['hfv_chain_diff'] += 1
                T['cls'].append(c); T['prod'].append(pr); T['hfv'].append(hv); T['ptrue'].append(float(np.linalg.norm(p3[m])))
                T['lpdg'].append(int(pdg[m])); T['chain'].append(int(nstep[m])); T['cls_chain'].append(ch)
            # true V0 table: baseline track(s) linked to each daughter, highest link weight kept; a track counts for a
            # daughter when its largest-weight link is the daughter or one of its continuations
            leg = {}
            for o in base:
                m = int(head[lw[o]]) if lw[o] >= 0 else -1
                if m >= 0:
                    n_, w_, l_ = leg.get(m, (0, -1e30, -1))
                    leg[m] = (n_ + 1, bw[o], int(loc[o])) if bw[o] > w_ else (n_ + 1, w_, l_)
            tv = true_v0s(pdg, chg, vtx, p3, dbeg, dend, didx, par1, gen)
            for i, dp, dn, org in tv:
                V0T['pdg'].append(int(pdg[i])); V0T['pos'].append(vtx[dp]); V0T['p3'].append(p3[i])
                V0T['rxy'].append(float(np.hypot(*(vtx[dp][:2] - gen[:2])))); V0T['prod'].append(vtx[i])
                for s_, d_ in (('1', dp), ('2', dn)):
                    n_, _, l_ = leg.get(d_, (0, 0., -1))
                    V0T['l' + s_].append(l_); V0T['n' + s_].append(n_); V0T['lp' + s_].append(p3[d_])
                V0T['orig'].append(org); V0T['mc'].append(i)
            E['nv0t'].append(len(tv))
            # every baseline track whose largest-weight link is a true-V0 daughter or a continuation of one -> (V0 index,
            # which daughter)
            dau = {}
            for k_, (i, dp, dn, org) in enumerate(tv):
                dau[dp] = (k_, 1); dau[dn] = (k_, 2)
            for o in base:
                k_, s_ = dau.get(int(head[lw[o]]) if lw[o] >= 0 else -1, (-1, -1))
                T['v0t'].append(k_); T['v0leg'].append(s_)
        E['nhf'].append(nhf)
    out = {}
    for pre, D in (('ev_', E), ('t_', T), ('v_', V), ('s_', SV), ('st_', ST), ('h_', H)) + ((('v0t_', V0T),) if isMC else ()):
        for k, v in D.items():
            if len(v) == 0 and not isMC and pre in ('t_', 'h_') and k in ('cls', 'prod', 'hfv', 'ptrue', 'pos', 'isb', 'pdg',
                                                                       'nch', 'nmatched', 'parent', 'p3', 'fd', 'v0t', 'v0leg'):
                continue
            out[pre + k] = np.asarray(v)
    for k in ('ev_pv', 'ev_pvcov6', 'ev_bs', 'ev_gen', 'ev_thrust', 'ev_jet_theta', 'ev_jet_phi', 'ev_jet_p', 't_par', 't_cov',
              't_p', 't_prod', 't_ptrue', 'v_pos', 'v_mass', 's_pos', 's_cov6', 's_mass', 'h_pos', 'h_p3', 'h_fd',
              'v0t_pos', 'v0t_p3', 'v0t_rxy', 'v0t_prod', 'v0t_lp1', 'v0t_lp2'):
        if k in out:
            out[k] = out[k].astype(np.float64)
    out['isMC'] = np.array(isMC)
    out['stem'] = np.array(stem)
    OUT = env_path('JV_INPUTS')
    os.makedirs(OUT, exist_ok=True)
    nm = stem if maxev is None else f'test_{stem}'
    np.savez_compressed(f'{OUT}/{nm}.npz.tmp.npz', **out)
    os.replace(f'{OUT}/{nm}.npz.tmp.npz', f'{OUT}/{nm}.npz')
    msg = dict(nev=nev, ntrk=int(sum(E['ntrk'])), nv0=int(sum(E['nv0'])), nsv=int(sum(E['nsv'])),
               st_unmapped=int(np.sum(np.asarray(ST['t']) < 0)) if ST['t'] else 0,
               v0_unmapped=int(np.sum((np.asarray(V['t1']) < 0) | (np.asarray(V['t2']) < 0))) if V['t1'] else 0,
               q_match=f'{qchk[1]}/{qchk[0]}', sec=round(time.time() - t0, 1))
    if isMC:
        msg['cls_counts'] = np.bincount(np.asarray(T['cls']), minlength=6).tolist()
        msg['nhf'] = int(sum(E['nhf']))
        msg['nv0t'] = int(sum(E['nv0t']))
        msg['v0t_legs_tracked'] = int(sum(n > 0 for n in V0T['n1'] + V0T['n2']))
        msg['v0t_leg_tracks'] = int(sum(V0T['n1']) + sum(V0T['n2']))
        msg['t_v0t_tracks'] = int(sum(k_ >= 0 for k_ in T['v0t']))
        msg['t_chain_tracks'] = int(sum(s_ > 0 for s_ in T['chain']))
        msg['cls_chain_counts'] = np.bincount(np.asarray(T['cls_chain']), minlength=6).tolist()
        msg['hfv_chain_diff'] = chain_stats['hfv_chain_diff']
    return stem, 'OK', msg


if __name__ == '__main__':
    nproc = int(sys.argv[1])
    with Pool(nproc) as p:
        for r in p.imap_unordered(job, sys.argv[2:]):
            print(r, flush=True)
    print('DONE', flush=True)
