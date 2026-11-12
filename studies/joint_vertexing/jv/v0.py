"""V0 class of the joint vertexer: K0S -> pi+ pi-, Lambda / anti-Lambda -> p pi and photon conversions gamma -> e+ e-
as two-track vertex hypotheses inside the S4 annealing (jv.hybrid.stage4), competing with the PV, the SV candidates and
the outlier term for their two tracks. Switched on by the config key v0_class (default off: the frozen final is
unchanged). Parameter tables: jv/v0_params.json, selected by the config keys named below.

Candidates (raw_candidates, conv_raw): opposite-charge pairs of the S0 tracks.
  K0S / Lambda: pre-selection at the xy crossings of the two track circles: crossing radius < v0_rmax, z of the two
  tracks compatible (v0_dz_sig standard deviations or v0_dz_abs cm), mass at the crossing inside the pre-window of a
  species (v0_pre_ks / v0_pre_lam GeV). Two-track fit (fit.fit_vertex_nb) from each passing crossing, lowest chi2 kept;
  chi2 < v0_chi2_max; mass at the fitted vertex inside the wide window (v0_win_ks / v0_win_lam); Armenteros qT of the
  daughters relative to the V0 momentum > v0_qt_min (K0S, Lambda). K0S = pi pi; Lambda: both assignments (positive track
  = proton: Lambda, negative track = antiproton: anti-Lambda), the one with the smaller mass pull kept. No pointing cut.
  Conversions (v0_conv): the loose selection of jv/conv.py (constrained conversion fit CF), conversion vertex radius
  from the detector axis >= v0_conv_rmin.
Pair hypothesis k = (tracks a, b; species s) in the weight matrix (V0Set.apply):
  log R_k = b_s(r) + sum_{t = a, b} [log f(chi2_tk / h_t^2, T) - log A_t] + [log L_m + log L_p + log L_h + log L_x] / T
  f       = Student-t soft-assignment function of the track term (nu = v0_nu);
  chi2_tk = chi2 of track t to the V0 vertex (pair fit, track covariance); h_t = v0_hit_scale for a K0S / Lambda
            daughter with VDET or ITC hits inside the decay radius (wrongly attached inner hits: biased, over-confident
            parameters), else 1; conversions: v0_conv_hscale (conv_terms);
  A_t     = S4 denominator of track t (outlier constant + PV numerator + SV numerators): the pair competes with
            everything else its two tracks could be; twin rule (V0Set.denominators): an SV within v0_twin_sig standard
            deviations of the V0 vertex is the pair's own hypothesis, so its numerators leave the pair's denominators;
  L_m     = exp(-z^2 / 2) / (sigma_m / 1 MeV), z = (m - m_s) / sigma_m; sigma_m propagated from the two track
            covariances at fixed arcs (no vertex error), floored at v0_sm_min (K0S / Lambda only);
  L_p     = (1 + q_p / nu)^(-w_s): q_p = 2D significance^2 of the miss distance of the V0 flight line from the
            reference PV (V0 vertex + PV covariance + momentum-direction error x decay length); nu = tail of the true-V0
            q_p distribution (2D Student-t fit); w_s = strength per species fitted with the prior ('point' table) times
            the PV-pointing factor pv_wscale of the HF-term table (K0S / Lambda; conversions: conv_terms); with
            v0_pt_point the pointing to the best origin, PV or SV (V0Set.origins);
  L_h     = K0S / Lambda: product over daughters with inner hits inside the decay radius of eps(p) (rate of wrongly
            attached inner hits per track momentum, table 'hit_eps'); conversions: conv_terms;
  L_x     = V0-against-heavy-flavour term (v0_hf; hf_features, hf_term), 0 when off;
  b_s(r)  = log prior odds per species and transverse distance r from the reference PV (table 'prior'[v0_prior]),
            calibrated so that R_k / (1 + R_k) at the reference PV (T = 1, no SVs) is the probability that the
            candidate is a true V0 of the species.
  Admission: candidates with log R_k < v0_admit at the reference PV are dropped before the fit; at most v0_max_comp
  candidates per connected component (the lowest log R at the reference PV dropped).
  Coupling: the two daughters enter or leave together and a track belongs to at most one V0: P_k = exact marginal of
  candidate k over the matchings (sets of candidates without a shared track) weighted by prod R_k, enumerated per
  connected component; every non-V0 weight of track t is multiplied by 1 - sum_{k contains t} P_k. A V0 vertex is the
  fixed pair-fit point (its tracks never change, so it is not refitted). Output: V0s with P_k > v0_out_min.
The tables were fitted on simulation (file mc_ZM4212.50); the final weight P_k is not calibrated as a probability.
Units cm, GeV.
"""
import json
import os
import numpy as np
from numba import njit

from . import fit as F
from . import tracks as TR

M_PI = TR.M_PI
M_P, M_KS, M_LAM = 0.938272, 0.497611, 1.115683
SPEC_PDG = {1: 310, 2: 3122, 3: 3122, 4: 22}    # hypothesis code -> booked species (2: Lambda, 3: anti-Lambda, 4: conversion)
NF = 30                                          # float columns of the kernel output (layout in pair_kernel)


@njit(cache=True)
def _circle(p):
    D, phi0, C = p[0], p[1], p[2]
    xc = -(D + 0.5 / C) * np.sin(phi0)
    yc = (D + 0.5 / C) * np.cos(phi0)
    return xc, yc, 0.5 / abs(C)


@njit(cache=True)
def _crossings(p1, p2, out):
    """xy crossings of the two circles (as tracks.circle_crossings): number of points written to out (2, 2)"""
    if p1[2] == 0 or p2[2] == 0:
        return 0
    x0, y0, r0 = _circle(p1)
    x1, y1, r1 = _circle(p2)
    dx, dy = x1 - x0, y1 - y0
    d = np.sqrt(dx * dx + dy * dy)
    if not np.isfinite(d) or not d > 0:
        return 0
    ux, uy = dx / d, dy / d
    if d >= r0 + r1:
        t = r0 + 0.5 * (d - r0 - r1)
        out[0, 0] = x0 + t * ux; out[0, 1] = y0 + t * uy
        return 1
    if d <= abs(r0 - r1):
        s = 1. if r0 > r1 else -1.
        out[0, 0] = 0.5 * (x0 + x1 + s * (r0 + r1) * ux); out[0, 1] = 0.5 * (y0 + y1 + s * (r0 + r1) * uy)
        return 1
    a = (r0 * r0 - r1 * r1 + d * d) / (2. * d)
    h2 = r0 * r0 - a * a
    X, Y = x0 + a * ux, y0 + a * uy
    if h2 > 0:
        h = np.sqrt(h2)
        out[0, 0] = X - h * uy; out[0, 1] = Y + h * ux
        out[1, 0] = X + h * uy; out[1, 1] = Y - h * ux
        return 2
    out[0, 0] = X; out[0, 1] = Y
    return 1


@njit(cache=True)
def _arc(p, x, y):
    """transverse arc of the circle point nearest to (x, y) (= tracks.arc_to_xy)"""
    D, phi0, C = p[0], p[1], p[2]
    if C == 0:
        return 0.
    x0, y0 = -D * np.sin(phi0), D * np.cos(phi0)
    xc, yc = -(D + 0.5 / C) * np.sin(phi0), (D + 0.5 / C) * np.cos(phi0)
    a0 = np.arctan2(y0 - yc, x0 - xc)
    a1 = np.arctan2(y - yc, x - xc)
    dphi = (a1 - a0 + np.pi) % (2 * np.pi) - np.pi
    return dphi / (2 * C)


@njit(cache=True)
def _mom_jac(p, L, K, mom, J):
    """momentum at arc L (fixed arc) and its Jacobian d mom / d par (3 x 5); K = PT_CONST / 2 (pT = K / |C|)"""
    C, ct = p[2], p[4]
    pt = K / max(abs(C), 1e-12)
    ang = p[1] + 2 * C * L
    ca, sa = np.cos(ang), np.sin(ang)
    mom[0] = pt * ca; mom[1] = pt * sa; mom[2] = pt * ct
    for i in range(3):
        for j in range(5):
            J[i, j] = 0.
    J[0, 1] = -pt * sa; J[1, 1] = pt * ca
    Cs = C if C != 0 else 1e-12
    J[0, 2] = -pt / Cs * ca - 2 * L * pt * sa
    J[1, 2] = -pt / Cs * sa + 2 * L * pt * ca
    J[2, 2] = -pt / Cs * ct
    J[2, 4] = pt


@njit(cache=True)
def _mass_err(pa, Ja, ca, pb, Jb, cb, ma, mb):
    """invariant mass of momenta pa, pb with masses ma, mb and its error from the track covariances (Jacobians Ja, Jb)"""
    Ea = np.sqrt(pa @ pa + ma * ma); Eb = np.sqrt(pb @ pb + mb * mb)
    P = pa + pb
    m2 = (Ea + Eb) ** 2 - P @ P
    if not m2 > 0:
        return 0., 1.
    m = np.sqrt(m2)
    ga = ((Ea + Eb) * pa / Ea - P) / m
    gb = ((Ea + Eb) * pb / Eb - P) / m
    va = ga @ Ja
    vb = gb @ Jb
    s2 = va @ ca @ va + vb @ cb @ vb
    return m, np.sqrt(max(s2, 1e-18))


@njit(cache=True)
def pair_kernel(par, cov, q, idx, K, rmax_v, dz_sig, dz_abs, pre_ks, pre_lam, win_ks, win_lam, chi2max, qt_ks, qt_lam, oi, of):
    """V0 pair candidates among the tracks idx (local indices into par). Writes rows to oi (n, 3) = (a, b, hypothesis:
    1 K0S, 2 Lambda (a or b positive = proton), 3 anti-Lambda (negative = antiproton)) and of (n, NF) =
    x(0:3) C(3:12) chi2_a chi2_b(12:14) m sigma_m(14:16) P(16:19) CP(19:28) |p_a| |p_b|(28:30); a = positive track.
    Returns the number of rows."""
    n = len(idx)
    nout = 0
    cr = np.zeros((2, 2))
    ma = np.zeros(3); mb = np.zeros(3); Ja = np.zeros((3, 5)); Jb = np.zeros((3, 5))
    Binv = np.zeros((3, 3)); c0 = np.zeros(3)
    two = np.zeros(2, np.int64); w2 = np.ones(2)
    Ls = np.zeros(par.shape[0])
    for u in range(n):
        for v in range(u + 1, n):
            i0, j0 = idx[u], idx[v]
            if q[i0] * q[j0] >= 0:
                continue
            i, j = (i0, j0) if q[i0] > 0 else (j0, i0)       # i positive, j negative
            pi_, pj_ = par[i], par[j]
            nc = _crossings(pi_, pj_, cr)
            best_chi2 = 1e30
            bx = np.zeros(3); bC = np.zeros((3, 3)); bch = np.zeros(2)
            for k in range(nc):
                cx, cy = cr[k, 0], cr[k, 1]
                if cx * cx + cy * cy >= rmax_v * rmax_v:
                    continue
                Li, Lj = _arc(pi_, cx, cy), _arc(pj_, cx, cy)
                zi = pi_[3] + pi_[4] * Li; zj = pj_[3] + pj_[4] * Lj
                vi = cov[i][3, 3] + 2 * Li * cov[i][3, 4] + Li * Li * cov[i][4, 4]
                vj = cov[j][3, 3] + 2 * Lj * cov[j][3, 4] + Lj * Lj * cov[j][4, 4]
                dz = abs(zi - zj)
                if dz > dz_abs and dz > dz_sig * np.sqrt(max(vi + vj, 1e-18)):
                    continue
                _mom_jac(pi_, Li, K, ma, Ja)
                _mom_jac(pj_, Lj, K, mb, Jb)
                Ea, Eb = np.sqrt(ma @ ma + M_PI * M_PI), np.sqrt(mb @ mb + M_PI * M_PI)
                P = ma + mb
                mks = np.sqrt(max((Ea + Eb) ** 2 - P @ P, 0.))
                Eap, Ebp = np.sqrt(ma @ ma + M_P * M_P), np.sqrt(mb @ mb + M_P * M_P)
                ml1 = np.sqrt(max((Eap + Eb) ** 2 - P @ P, 0.))
                ml2 = np.sqrt(max((Ea + Ebp) ** 2 - P @ P, 0.))
                if abs(mks - M_KS) > pre_ks and abs(ml1 - M_LAM) > pre_lam and abs(ml2 - M_LAM) > pre_lam:
                    continue
                Ls[i] = Li; Ls[j] = Lj
                two[0] = i; two[1] = j
                x0 = np.array([cx, cy, 0.5 * (zi + zj)])
                x, C, ch, conv = F.fit_vertex_nb(par, cov, two, w2, x0, Binv, c0, Ls, 15, 1e-6)
                tot = ch[0] + ch[1]
                if not (np.isfinite(tot) and np.all(np.isfinite(x))):
                    continue
                if tot < best_chi2:
                    best_chi2 = tot
                    bx[:] = x; bC[:, :] = C; bch[0] = ch[0]; bch[1] = ch[1]
            if not best_chi2 < chi2max or bx[0] * bx[0] + bx[1] * bx[1] >= rmax_v * rmax_v:
                continue
            Li, Lj = _arc(pi_, bx[0], bx[1]), _arc(pj_, bx[0], bx[1])
            _mom_jac(pi_, Li, K, ma, Ja)
            _mom_jac(pj_, Lj, K, mb, Jb)
            P = ma + mb
            pP = np.sqrt(P @ P)
            qt2 = 0.
            if pP > 0:
                pl = (ma @ P) / pP
                qt2 = ma @ ma - pl * pl                      # Armenteros qT^2: photon conversions have qT ~ 0
            CP = Ja @ cov[i] @ Ja.T + Jb @ cov[j] @ Jb.T
            for hyp in range(1, 4):
                m1 = M_PI if hyp != 2 else M_P
                m2 = M_PI if hyp != 3 else M_P
                m, sm = _mass_err(ma, Ja, cov[i], mb, Jb, cov[j], m1, m2)
                if hyp == 1 and (abs(m - M_KS) > win_ks or qt2 < qt_ks * qt_ks):
                    continue
                if hyp > 1 and (abs(m - M_LAM) > win_lam or qt2 < qt_lam * qt_lam):
                    continue
                if nout >= oi.shape[0]:
                    return nout
                oi[nout, 0] = i; oi[nout, 1] = j; oi[nout, 2] = hyp
                of[nout, 0:3] = bx
                for r in range(3):
                    for s in range(3):
                        of[nout, 3 + 3 * r + s] = bC[r, s]
                        of[nout, 19 + 3 * r + s] = CP[r, s]
                of[nout, 12] = bch[0]; of[nout, 13] = bch[1]
                of[nout, 14] = m; of[nout, 15] = sm
                of[nout, 16:19] = P
                of[nout, 28] = np.sqrt(ma @ ma); of[nout, 29] = np.sqrt(mb @ mb)
                nout += 1
    return nout


@njit(cache=True)
def _component_marginals(logR, conf):
    """exact marginals of the candidates of one component: sets without conflicts (conf[k] = bitmask of the candidates
    sharing a track with k) weighted by exp(sum logR); returns P (n)"""
    n = len(logR)
    nset = 1 << n
    lw = np.full(nset, -np.inf)
    mx = -np.inf
    for S in range(nset):
        ok = True
        s = 0.
        for k in range(n):
            if S & (1 << k):
                if conf[k] & S:
                    ok = False
                    break
                s += logR[k]
        if ok:
            lw[S] = s
            if s > mx:
                mx = s
    Z = 0.
    P = np.zeros(n)
    for S in range(nset):
        if lw[S] > -np.inf:
            w = np.exp(lw[S] - mx)
            Z += w
            for k in range(n):
                if S & (1 << k):
                    P[k] += w
    return P / Z


def _bin(edges, v):
    return int(np.clip(np.searchsorted(edges, v, side='right') - 1, 0, len(edges) - 2))


_PAR = {}


def load_json(name):
    """parameter file name (in jv/), read once per process"""
    if name not in _PAR:
        _PAR[name] = json.load(open(os.path.join(os.path.dirname(__file__), name)))
    return _PAR[name]


def check_tables(c):
    """KeyError naming the first missing entry of c['v0_params'] that the V0 configuration c selects: checks that the
    selected tables exist (point strength, prior, hit_eps, conversion pointing and hit tables, the HF table with its
    K0S, Lambda and conversion entries and the K0S / Lambda pv_wscale) and, for v0_pt_point, the fields logrho, w, nu;
    the table contents are not checked"""
    P = load_json(c['v0_params'])
    need = [('point', 'w', c['v0_prior']), ('prior', c['v0_prior']), ('hit_eps', 'eps')]
    if c['v0_prior_nm'] is not None:
        need.append(('prior', c['v0_prior_nm']))
    if c['v0_conv']:
        need.append(('conv_point', c['v0_conv_params']))
        if c['v0_hit_table'] is not None:
            need.append(('hit_table', c['v0_hit_table'], 'conversion'))
    if c['v0_hf'] is not None:
        need += [('hfterm', c['v0_hf'], s) for s in ('K0S', 'Lambda', 'conversion')]
        need += [('hfterm', c['v0_hf'], s, 'pv_wscale') for s in ('K0S', 'Lambda')]
    if c['v0_pt_point']:
        need += [('pt', c['v0_pt_params'], s, f) for s in ('K0S', 'Lambda') for f in ('logrho', 'w', 'nu')]
    for path in need:
        x = P
        for key in path:
            if not isinstance(x, dict) or key not in x:
                raise KeyError(f"{c['v0_params']}: no table {'/'.join(str(p) for p in path)}")
            x = x[key]


def track_logf(c, chi2, T):
    """log of the V0 track-term soft-assignment function (Student-t, v0_nu dof)"""
    nu = c['v0_nu']
    return -0.5 * (nu + 2.) * np.log1p(np.minimum(chi2, 1e6) / (nu * T))


def hit_eps(c, p):
    H = load_json(c['v0_params'])['hit_eps']
    return float(H['eps'][_bin(H['p_bins'], p)])


def point_term(c, q, spec):
    """log L_p = -w log(1 + q / nu) (PV pointing, relative to q = 0); 0 when v0_point is off. nu = tail of the true-V0
    q_p distribution (2D Student-t; w = (nu + 2) / 2 would be its density alone); w = strength per species fitted
    together with the prior, because the other candidates also spread to large q_p ('point' 'w'[v0_prior]), times the
    PV-pointing factor the HF-term table was fitted with ('hfterm'[v0_hf][species]['pv_wscale'], when v0_hf is set)"""
    if not c['v0_point']:
        return 0. * q
    P = load_json(c['v0_params'])
    nus = P['point']['nu_signal']
    w = P['point']['w'][c['v0_prior']][spec]
    if c['v0_hf'] is not None:
        w = w * P['hfterm'][c['v0_hf']][spec]['pv_wscale']
    return -w * np.log1p(np.maximum(q, 0.) / nus)


def pointing_q(x, Cx, P, CP, pvx, pvC):
    """2D significance^2 of the miss distance of the V0 flight line (x, direction P) from the PV"""
    pm = np.linalg.norm(P)
    if not pm > 0:
        return 0.
    u = P / pm
    d = x - pvx
    L = d @ u
    m = d - L * u
    e1 = np.cross(u, [0., 0., 1.]) if abs(u[2]) < 0.9 else np.cross(u, [1., 0., 0.])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(u, e1)
    B = np.stack([e1, e2])
    Cm = B @ (Cx + pvC + (L * L / (pm * pm)) * CP) @ B.T
    mm = B @ m
    try:
        return max(float(mm @ np.linalg.solve(Cm, mm)), 0.)
    except np.linalg.LinAlgError:
        return 0.


def flight_sig(x, C, pvx, pvC):
    """3D flight significance of the V0 vertex from the reference PV (V0 + PV covariance)"""
    d = x - pvx
    try:
        return float(np.sqrt(max(d @ np.linalg.solve(C + pvC, d), 0.)))
    except np.linalg.LinAlgError:
        return 0.


def raw_candidates(E, c, ok, rmax, pvx, pvC, K):
    """kernel candidates of the event with their static quantities (dicts, global track indices); Lambda: both
    assignments kept here"""
    t = E['t']
    idx = np.nonzero(ok)[0].astype(np.int64)
    n = len(idx)
    if n < 2:
        return []
    nmax = 3 * n * (n - 1) // 2
    oi = np.zeros((nmax, 3), np.int64)
    of = np.zeros((nmax, NF))
    q = t['q'].astype(np.float64)
    nc = pair_kernel(t['fpar'], t['fcov'], q, idx, K, c['v0_rmax'], c['v0_dz_sig'], c['v0_dz_abs'], c['v0_pre_ks'],
                     c['v0_pre_lam'], c['v0_win_ks'], c['v0_win_lam'], c['v0_chi2_max'], c['v0_qt_min'][0], c['v0_qt_min'][1], oi, of)
    out = []
    for r in range(nc):
        a, b, hyp = int(oi[r, 0]), int(oi[r, 1]), int(oi[r, 2])
        f = of[r]
        x = f[0:3].copy(); C = f[3:12].reshape(3, 3).copy(); P = f[16:19].copy(); CP = f[19:28].reshape(3, 3).copy()
        rv = float(np.hypot(x[0], x[1]))
        pV = float(np.linalg.norm(P))
        spec = 'K0S' if hyp == 1 else 'Lambda'
        out.append(dict(a=a, b=b, hyp=hyp, spec=spec, x=x, C=C, chi2a=float(f[12]), chi2b=float(f[13]), m=float(f[14]),
                        sm=float(f[15]), P=P, CP=CP, pa=float(f[28]), pb=float(f[29]), p=pV,
                        viola=bool(rv >= rmax[a]), violb=bool(rv >= rmax[b]),
                        r=float(np.hypot(x[0] - pvx[0], x[1] - pvx[1])), qp=pointing_q(x, C, P, CP, pvx, pvC),
                        fsig=flight_sig(x, C, pvx, pvC), rv=rv))
    return out


def conv_raw(E, c, ok, rmax, pvx, pvC, K):
    """photon-conversion candidates of the S0 tracks: the loose selection of jv/conv.py (CF chi2 < 200, m_ee at the
    tangent point < 50 MeV, cos of the pointing angle to the PV > 0.9; momenta with K = PT_CONST / 2 of pt_const); per
    candidate the CF vertex and covariance, the CF chi2 split per track, the z-side chi2 of the tangent point (dz,
    dtanLambda with their correlation, 2 dof), cos of the pointing angle"""
    from . import conv as CV
    t = E['t']
    idx = np.nonzero(ok)[0]
    if len(idx) < 2:
        return []
    R = CV.find_conversions(t['fpar'][idx], t['fcov'][idx], t['q'][idx], t['nvdet'][idx], t['nitc'][idx], pvx, sel='loose',
                            k_pt=K)
    Q = R['Q']
    QI = CV.QI
    out = []
    for r in range(len(R['i'])):
        i, j = int(idx[R['i'][r]]), int(idx[R['j'][r]])
        a, b = (i, j) if t['q'][i] > 0 else (j, i)
        q_ = Q[r]
        x = R['x'][r].copy(); C = R['cov'][r].copy()
        chi2 = float(q_[QI['cf_chi2']])
        c1, c2 = float(q_[QI['cf_chi2_1']]), float(q_[QI['cf_chi2_2']])
        ca, cb = (c1, c2) if a == i else (c2, c1)
        sdz, sdc, rzc = q_[QI['tp_sdz']], q_[QI['tp_sdct']], q_[QI['tp_rzc']]
        u, v = q_[QI['tp_dz']] / max(sdz, 1e-12), q_[QI['tp_dct']] / max(sdc, 1e-12)
        zs = float((u * u - 2 * rzc * u * v + v * v) / max(1. - rzc * rzc, 1e-6))
        P = np.array([q_[QI['cf_px']], q_[QI['cf_py']], q_[QI['cf_pz']]])
        rv = float(np.hypot(x[0], x[1]))
        pa, pb = (float(q_[QI['p1']]), float(q_[QI['p2']])) if a == i else (float(q_[QI['p2']]), float(q_[QI['p1']]))
        out.append(dict(a=a, b=b, hyp=4, spec='conversion', x=x, C=C, chi2a=ca, chi2b=cb, chi2cf=chi2, zside=zs,
                        m=float(q_[QI['tp_mee']]), sm=0., P=P, CP=np.zeros((3, 3)), pa=pa, pb=pb, p=float(np.linalg.norm(P)),
                        viola=bool(rv >= rmax[a]), violb=bool(rv >= rmax[b]), cosp=float(q_[QI['cf_cosp']]),
                        r=float(np.hypot(x[0] - pvx[0], x[1] - pvx[1])), qp=0., fsig=flight_sig(x, C, pvx, pvC), rv=rv))
    return out


def conv_terms(c, k):
    """static terms of a conversion candidate: no mass term (m_ee is the collinearity the CF chi2 already holds);
    pointing log L_p = -w log(1 + max(1 - cos, 0) / tau) (one-sided, 'conv_point'[v0_conv_params]); hit term =
    measured log likelihood ratio of 0 / 1 / 2 daughters with inner hits inside the CF radius ('hit_table'
    [v0_hit_table]['conversion'], radius bins from the axis; 0 when v0_hit_table is None); chi2 of the track term: the
    CF chi2 split per track divided by h^2 for a daughter with such hits (v0_conv_hscale: one / both); with both, half
    the z-side chi2 of the tangent point per track replaces the CF chi2 (the transverse gap is biased then)"""
    Pp = load_json(c['v0_params'])
    cp = Pp['conv_point'][c['v0_conv_params']]
    lp = -cp['w'] * np.log1p(max(1. - k['cosp'], 0.) / cp['tau'])
    n = int(k['viola']) + int(k['violb'])
    lh = 0.
    if c['v0_hit_table'] is not None:
        H = Pp['hit_table'][c['v0_hit_table']]['conversion']
        lh = float(H['log_lr'][_bin(H['r_bins'], k['rv'])][n])
    h1, h2 = c['v0_conv_hscale']
    ha = hb = 1.
    if n == 2:
        k['chi2a'] = k['chi2b'] = 0.5 * k['zside']
        ha = hb = h2
    elif n == 1:
        ha, hb = (h1, 1.) if k['viola'] else (1., h1)
    return 0., 0., float(lp), lh, ha, hb


def static_terms(c, k):
    """T-independent pieces of candidate k: mass pull z, log L_m, log L_p, log L_h, hit scales h_a, h_b"""
    if k['hyp'] == 4:
        return conv_terms(c, k)
    m0 = M_KS if k['hyp'] == 1 else M_LAM
    sm = max(k['sm'], c['v0_sm_min'])
    z = (k['m'] - m0) / sm
    lm = (-0.5 * z * z - np.log(sm / 0.001)) if c['v0_mass'] else 0.
    lp = float(point_term(c, np.float64(k['qp']), k['spec']))
    lh, ha, hb = 0., 1., 1.
    if k['viola']:
        lh += np.log(hit_eps(c, k['pa'])); ha = c['v0_hit_scale']
    if k['violb']:
        lh += np.log(hit_eps(c, k['pb'])); hb = c['v0_hit_scale']
    return z, lm, lp, lh, ha, hb


CTAU = {'K0S': 2.6844, 'Lambda': 7.89}           # c tau [cm] (PDG)
M_SPEC = {'K0S': M_KS, 'Lambda': M_LAM}


def hf_features(E, c, ok, cands, s3):
    """static features of the V0-against-heavy-flavour (HF) term per candidate (written into the dicts):
    n_ext = other S0 tracks compatible with the V0 vertex (chi2 to the point < v0_hf_chi2, track covariance) and displaced
            from the reference PV (lifetime-signed 3D IP significance s3 > v0_hf_sig): further tracks of an HF vertex;
    n_jet = other S0 tracks of the jet(s) of the two daughters with s3 > v0_hf_sig (jet context);
    dl    = log density of the transverse decay length r for the species' lifetime: -log(lam) - r / lam,
            lam = (pT / m) c tau (0 for conversions);
    ct    = |cos theta*| of the decay (abs_cos_star; 0 for conversions)"""
    t = E['t']
    idx = np.nonzero(ok)[0]
    disp = idx[s3[idx] > c['v0_hf_sig']]
    jet = t['jet'].astype(np.int64)
    for k in cands:
        o = disp[(disp != k['a']) & (disp != k['b'])]
        k['n_ext'] = 0; k['n_jet'] = 0
        if len(o):
            ch, _ = TR.chi2_to_point(t['fpar'][o], t['fcov'][o], k['x'])
            k['n_ext'] = int((ch < c['v0_hf_chi2']).sum())
            ja, jb = jet[k['a']], jet[k['b']]
            k['n_jet'] = int((((jet[o] == ja) & (ja >= 0)) | ((jet[o] == jb) & (jb >= 0))).sum())
        if k['spec'] in CTAU:
            lam = max(float(np.hypot(k['P'][0], k['P'][1])) / M_SPEC[k['spec']] * CTAU[k['spec']], 1e-3)
            k['dl'] = float(-np.log(lam) - k['r'] / lam)
            k['ct'] = abs_cos_star(k)
        else:
            k['dl'] = 0.
            k['ct'] = 0.


def abs_cos_star(k):
    """|cos theta*| of the positive daughter in the V0 rest frame (hypothesis masses, fitted pair mass), from |p_a|,
    |p_b| and the pair momentum: flat for a true two-body decay, peaked at 1 for combinatorial pairs"""
    ma = M_P if k['hyp'] == 2 else M_PI
    mb = M_P if k['hyp'] == 3 else M_PI
    P = float(np.linalg.norm(k['P'])); M = max(k['m'], ma + mb + 1e-6)
    if not P > 0:
        return 0.
    pl = (P * P + k['pa'] ** 2 - k['pb'] ** 2) / (2 * P)
    E = np.sqrt(P * P + M * M)
    g, bt = E / M, P / E
    ps = np.sqrt(max((M * M - (ma + mb) ** 2) * (M * M - (ma - mb) ** 2), 0.)) / (2 * M)
    es = (M * M + ma * ma - mb * mb) / (2 * M)
    if not ps > 0:
        return 1.
    return float(min(abs((pl / g - bt * es) / ps), 1.))


def hf_term(c, k):
    """log likelihood-ratio term V0 against the other candidates (HF two-prongs above all) from the hf_features,
    'hfterm'[v0_hf][species]: b0 + ext[min(n_ext, 2) - 1] + jet[min(n_jet, 3) - 1] + wdl * dl (+ ct[|cos theta*| bin]
    when the table has 'ct'); logistic fit with the reference-stage log odds as offset"""
    H = load_json(c['v0_params'])['hfterm'][c['v0_hf']][k['spec']]
    lx = H['b0'] + H['wdl'] * k['dl']
    if k['n_ext'] > 0:
        lx += H['ext'][min(k['n_ext'], len(H['ext'])) - 1]
    if k['n_jet'] > 0:
        lx += H['jet'][min(k['n_jet'], len(H['jet'])) - 1]
    if 'ct' in H and k['spec'] in CTAU:
        lx += H['ct'][_bin(H['ct_bins'], k['ct'])]
    return float(lx)


def prior_b(c, k, key='v0_prior'):
    """log prior odds per species and transverse distance r from the reference PV; tables with 'fs_bins' are binned
    in (r, flight significance) below r = 'fs_rmax' (key '<species>_fs': [r bin][fsig bin])"""
    if c[key] is None:
        return 0.
    B = load_json(c['v0_params'])['prior'][c[key]]
    rb = _bin(B['r_bins'], k['r'])
    if 'fs_bins' in B and k['r'] < B['fs_rmax']:
        return float(B[k['spec'] + '_fs'][rb][_bin(B['fs_bins'], k['fsig'])])
    return float(B[k['spec']][rb])


class V0Set:
    """the admitted V0 candidates of one event inside S4 (track indices local to fit_idx after localize)"""

    def __init__(self, c, cands):
        self.c = c
        self.k = cands
        self.n = len(cands)
        self.a = np.array([k['a'] for k in cands], np.int64)
        self.b = np.array([k['b'] for k in cands], np.int64)
        self.chi2a = np.array([k['chi2a'] / k['ha'] ** 2 for k in cands])
        self.chi2b = np.array([k['chi2b'] / k['hb'] ** 2 for k in cands])
        self.stat = np.array([k['lm'] + k['lp'] + k['lh'] for k in cands])
        self.bprior = np.array([k['bp'] for k in cands])
        self.stat_nm = np.array([k['lp'] + k['lh'] for k in cands])
        self.bprior_nm = np.array([k.get('bp_nm', 0.) for k in cands])
        if any('lx' in k for k in cands):
            # V0-against-HF term (v0_hf): a static term as the mass, pointing and hit terms
            lx = np.array([k.get('lx', 0.) for k in cands])
            self.stat = self.stat + lx
            self.stat_nm = self.stat_nm + lx
        self.x = [k['x'] for k in cands]
        self.C = [k['C'] for k in cands]
        # SV-pointing odds (v0_pt_point): each K0S / Lambda candidate as a neutral pseudo-track (vertex x, C and momentum
        # P, CP of the pair fit); dstat = log of the pointing to the best origin minus the PV pointing term (None = PV
        # only); org[T] = (origin posteriors, V0 weights) of the last apply at temperature T (see origins)
        self.Pm = [k['P'] for k in cands]
        self.CPm = [k['CP'] for k in cands]
        self.lp_pv = np.array([k['lp'] for k in cands])
        self.isconv = np.array([k['hyp'] == 4 for k in cands], bool)
        self.dstat = None
        self.org = {}
        self.last = None
        self.P0 = None
        self.P = np.zeros(self.n)
        self.comps = []
        self._components()

    def _components(self):
        tr = {}
        parent = list(range(self.n))

        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i
        for k in range(self.n):
            for t in (self.a[k], self.b[k]):
                if t in tr:
                    parent[find(k)] = find(tr[t])
                else:
                    tr[t] = k
        groups = {}
        for k in range(self.n):
            groups.setdefault(find(k), []).append(k)
        for g in groups.values():
            g = np.array(sorted(g), np.int64)
            conf = np.zeros(len(g), np.int64)
            for u in range(len(g)):
                for v in range(len(g)):
                    if u != v and len({self.a[g[u]], self.b[g[u]]} & {self.a[g[v]], self.b[g[v]]}):
                        conf[u] |= 1 << v
            self.comps.append((g, conf))

    def localize(self, loc):
        self.a = loc[self.a]; self.b = loc[self.b]

    def logR(self, A, T, Aa=None, Ab=None, nomass=False):
        """log R per candidate; Aa, Ab = denominators of the two tracks per candidate (default A of the tracks);
        nomass: the weight without the mass term (log L_m neutral, own prior b_nm)"""
        c = self.c
        if Aa is None:
            Aa, Ab = A[self.a], A[self.b]
        if self.dstat is not None:
            # v0_pt_point: the PV pointing term replaced by the pointing to the best origin (PV or SV)
            if nomass:
                return (self.bprior_nm + track_logf(c, self.chi2a, T) + track_logf(c, self.chi2b, T) - np.log(Aa) -
                        np.log(Ab) + (self.stat_nm + self.dstat) / T)
            return (self.bprior + track_logf(c, self.chi2a, T) + track_logf(c, self.chi2b, T) - np.log(Aa) -
                    np.log(Ab) + (self.stat + self.dstat) / T)
        if nomass:
            return (self.bprior_nm + track_logf(c, self.chi2a, T) + track_logf(c, self.chi2b, T) - np.log(Aa) -
                    np.log(Ab) + self.stat_nm / T)
        return (self.bprior + track_logf(c, self.chi2a, T) + track_logf(c, self.chi2b, T) - np.log(Aa) -
                np.log(Ab) + self.stat / T)

    def denominators(self, A, num, V):
        """per-candidate denominators of its two tracks: A, minus (v0_twin_sig) the numerators of the twin SVs, i.e. the
        SVs whose current position lies within v0_twin_sig standard deviations (V0 + SV covariance) of the V0 vertex: a
        twin is the same vertex hypothesis for the pair, so the V0 odds (prior, mass, pointing) decide between 'V0' and
        'not a V0' (a two-prong SV included, as in the calibration at the reference PV). v0_twin_ext: a number = an SV
        whose current soft weights on tracks other than the pair sum to more than it is no twin (a further-track HF
        vertex); 'frac' = the twin's numerators leave with the factor phi = (w_a + w_b) / (sum of the SV's current soft
        weights), 1 for a pure two-track SV"""
        Aa, Ab = A[self.a].copy(), A[self.b].copy()
        ts = self.c['v0_twin_sig']
        te = self.c['v0_twin_ext']
        if ts is None or num is None or V is None or len(V) < 2:
            return Aa, Ab
        for k in range(self.n):
            xk, Ck = self.x[k], self.C[k]
            for vi in range(1, len(V)):
                d = V[vi].x - xk
                try:
                    q = float(d @ np.linalg.solve(V[vi].C + Ck, d))
                except np.linalg.LinAlgError:
                    continue
                if q < ts * ts and te == 'frac':
                    wo = num[vi] / np.maximum(A, 1e-300)
                    phi = float(wo[self.a[k]] + wo[self.b[k]]) / max(float(wo.sum()), 1e-300)
                    Aa[k] -= phi * num[vi, self.a[k]]
                    Ab[k] -= phi * num[vi, self.b[k]]
                    continue
                if q < ts * ts and te is not None:
                    wo = num[vi] / np.maximum(A, 1e-300)
                    if float(wo.sum() - wo[self.a[k]] - wo[self.b[k]]) > te:
                        continue
                if q < ts * ts:
                    Aa[k] -= num[vi, self.a[k]]
                    Ab[k] -= num[vi, self.b[k]]
        return np.maximum(Aa, 1e-300), np.maximum(Ab, 1e-300)

    def marginals(self, lr):
        P = np.zeros(self.n)
        for g, conf in self.comps:
            if len(g) == 1:
                P[g[0]] = 1. / (1. + np.exp(-np.clip(lr[g[0]], -700, 700)))
            else:
                P[g] = _component_marginals(lr[g], conf)
        return P

    def apply(self, W, A, T, num=None, V=None):
        """V0 marginals at temperature T against the S4 denominators A (per local track; num = numerators per vertex and
        track, V = vertices, for the twin rule); scales the non-V0 weights"""
        Aa, Ab = self.denominators(A, num, V)
        if self.c['v0_pt_point']:
            self.origins(V, T)
        self.last = (Aa, Ab, T)
        self.P = self.marginals(self.logR(A, T, Aa, Ab))
        if T in self.org:
            self.org[T] = (self.org[T][0], self.P.copy())
        s = np.zeros(W.shape[1])
        np.add.at(s, self.a, self.P)
        np.add.at(s, self.b, self.P)
        self.share = s
        return W * (1. - np.minimum(s, 1.))[None, :]

    def weights_nomass(self):
        """weights without the mass term at the state of the last apply (same denominators, temperature and
        exclusivity; log L_m neutral; prior b_nm calibrated for it, v0_prior_nm): stored in P0"""
        if self.last is None:
            self.P0 = np.zeros(self.n)
            return self.P0
        Aa, Ab, T = self.last
        self.P0 = self.marginals(self.logR(None, T, Aa, Ab, nomass=True))
        return self.P0

    def origins(self, V, T):
        """SV-pointing odds (v0_pt_point): origins of the neutral pseudo-track of each K0S / Lambda candidate among the
        PV and the current SVs V[1:]. log f_0 = the PV pointing term (log L_p at the reference PV, unchanged); log f_v =
        log rho_s - w_s log(1 + q_v / nu_s) for an SV v that the V0 vertex lies downstream of ((x_V0 - x_v) . P > 0) and
        is no twin of (more than v0_twin_sig standard deviations away), q_v = 2D significance^2 of the miss distance of
        the V0 flight line from v (pointing_q with the SV covariance); rho_s, w_s, nu_s = 'pt'[v0_pt_params][species].
        Conversions keep the PV only. Sets dstat = log sum_v f_v - log f_0 (>= 0) and org[T] = (omega, None): omega
        (n, len(V)) = origin posteriors softmax(log f / T); apply adds the V0 weights"""
        c = self.c
        nv = len(V)
        lf = np.full((self.n, nv), -np.inf)
        lf[:, 0] = self.lp_pv
        if nv > 1 and self.n:
            Pt = load_json(c['v0_params'])['pt'][c['v0_pt_params']]
            ts = c['v0_twin_sig'] if c['v0_twin_sig'] is not None else 0.
            for k in range(self.n):
                if self.isconv[k]:
                    continue
                prm = Pt[self.k[k]['spec']]
                for vi in range(1, nv):
                    d = self.x[k] - V[vi].x
                    if d @ self.Pm[k] <= 0:
                        continue
                    try:
                        q2 = float(d @ np.linalg.solve(V[vi].C + self.C[k], d))
                    except np.linalg.LinAlgError:
                        continue
                    if q2 < ts * ts:
                        continue
                    q = pointing_q(self.x[k], self.C[k], self.Pm[k], self.CPm[k], V[vi].x, V[vi].C)
                    lf[k, vi] = prm['logrho'] - prm['w'] * np.log1p(q / prm['nu'])
        mx = lf.max(1)
        lse = mx + np.log(np.exp(lf - mx[:, None]).sum(1))
        self.dstat = np.where(np.isfinite(lf[:, 1:]).any(1), lse - self.lp_pv, 0.) if nv > 1 else np.zeros(self.n)
        lt = lf / T
        mt = lt.max(1)
        om = np.exp(lt - mt[:, None])
        om /= om.sum(1)[:, None]
        self.org[T] = (om, None)


def candidates(E, c, ok, rmax, pvx, pvC, A_ref, K, s3=None):
    """admitted candidates of the event (global track indices): kernel and conversion supply (v0_conv, radius bound
    v0_conv_rmin) + one Lambda assignment per pair + static terms + prior + HF term (v0_hf; needs s3, the lifetime-signed
    IP significances at the reference PV) + admission (log R at the reference PV, T = 1, A_ref = reference denominators
    per track) + the component-size limit"""
    raw = raw_candidates(E, c, ok, rmax, pvx, pvC, K)
    if c['v0_conv']:
        cr = conv_raw(E, c, ok, rmax, pvx, pvC, K)
        if c['v0_conv_rmin'] is not None:
            cr = [k for k in cr if k['rv'] >= c['v0_conv_rmin']]
        raw = raw + cr
    if c['v0_hf'] is not None:
        hf_features(E, c, ok, raw, s3)
    best = {}
    for k in raw:
        z, lm, lp, lh, ha, hb = static_terms(c, k)
        k.update(z=z, lm=lm, lp=lp, lh=lh, ha=ha, hb=hb, bp=prior_b(c, k), bp_nm=prior_b(c, k, 'v0_prior_nm'))
        if c['v0_hf'] is not None:
            k['lx'] = hf_term(c, k)
        key = (k['a'], k['b'], k['spec'])
        if key not in best or abs(z) < abs(best[key]['z']):
            best[key] = k
    out = []
    for k in best.values():
        lr = (k['bp'] + track_logf(c, k['chi2a'] / k['ha'] ** 2, 1.) + track_logf(c, k['chi2b'] / k['hb'] ** 2, 1.) -
              np.log(A_ref[k['a']]) - np.log(A_ref[k['b']]) + k['lm'] + k['lp'] + k['lh'])
        if 'lx' in k:
            lr = lr + k['lx']
        k['lr_ref'] = float(lr)
        if lr >= c['v0_admit']:
            out.append(k)
    # at most v0_max_comp candidates per connected component (exact enumeration): the lowest log R_ref dropped
    out.sort(key=lambda k: -k['lr_ref'])
    keep, cnt, comp_of = [], {}, {}
    for k in out:
        ca, cb = comp_of.get(k['a']), comp_of.get(k['b'])
        cid = ca if ca is not None else (cb if cb is not None else len(cnt))
        if ca is not None and cb is not None and ca != cb:
            n = cnt[ca] + cnt[cb]
            if n + 1 > c['v0_max_comp']:
                continue
            for t, g in list(comp_of.items()):
                if g == cb:
                    comp_of[t] = ca
            cnt[ca] = n; cnt[cb] = 0; cid = ca
        if cnt.get(cid, 0) + 1 > c['v0_max_comp']:
            continue
        cnt[cid] = cnt.get(cid, 0) + 1
        comp_of[k['a']] = cid; comp_of[k['b']] = cid
        keep.append(k)
    return keep
