"""Vertex fits of the joint vertexer.

fit_vertex   weighted Gauss-Newton fit of one vertex to tracks (linearised helix model of tracks.py, relinearised every
             iteration), optional Gaussian constraint (beam spot) -- with all weights 1 and the beam constraint
             (200, 100, 2e4 um) this is the stage1 PV fit (analyzer_pvnew.h)
teardown     fit + removal of the worst track while its chi2 >= chi2max (= stage1 PV pruning, chi2max 5)
pair_fit     two-track fit started at the perigee and at the circle crossings; solutions beyond the maximal
             vertex radius of either track excluded (hit-radius rule); lowest chi2 kept
jetproj_pv   jet-projected PV (QFNDIP-like, line model at the origin-referenced perigee), hard or soft form (jetproj.py)
Units: cm.
"""
import numpy as np
from numba import njit

from . import tracks as TR
from . import jetproj as JP

BS_SIG = np.array([200e-4, 100e-4, 2.0])    # stage1 PV beam constraint (cm) = aleph_reco_config.h


@njit(cache=True)
def fit_vertex_nb(par, cov, idx, w, x0, Binv, c, Ls, niter, tol):
    """weighted fit of the tracks idx (weights w, same length) from x0; Binv = constraint information (zeros = none),
    centre c. Ls: per-track start arcs (all tracks, updated). Returns x, C (3x3), per-track chi2 (len idx), converged"""
    n = len(idx)
    x = x0.copy()
    X = np.zeros((n, 3)); D = np.zeros((n, 3, 3)); chi2 = np.zeros(n)
    conv = False
    for it in range(niter):
        TR.lin_all(par, cov, idx, x, Ls, X, D, chi2)
        H = Binv.copy()
        g = Binv @ (x - c)
        for k in range(n):
            H += w[k] * D[k]
            g += w[k] * (D[k] @ (x - X[k]))
        dx = -np.linalg.solve(H, g)
        if not np.all(np.isfinite(dx)):
            break
        x = x + dx
        if np.max(np.abs(dx)) < tol:
            conv = True
            break
    TR.lin_all(par, cov, idx, x, Ls, X, D, chi2)
    H = Binv.copy()
    for k in range(n):
        H += w[k] * D[k]
    C = np.linalg.inv(H)
    return x, C, chi2, conv


@njit(cache=True)
def _phi(z):
    return np.exp(-0.5 * z * z)


@njit(cache=True)
def fit_pv_eu_nb(par, cov, wU, wEb, x0, Binv, c, axes, jet_of, sj2, kf, kb, ke, niter, tol, useeu):
    """PV fit with each track's information split in E = t x J / |t x J| (perpendicular to track and jet) and U = E x t (in the
    track-jet plane, positive along the jet). Base weights wU (U component), wEb (E component) per track. Jet-angle variance
    min(l4^2, sj2 (l4 / nu)^2) added to E (l4 = offset of the track from the vertex along U, nu = track-jet angle; decay length x
    jet angle error); U gets the one-sided weight phi(z) / (phi(z) + phi(kf)) in front (l4 > 0), kb behind, z = l4 / sigma_U;
    E gets the outlier weight with ke. The track nearest to the jet axis in angle when jet_of < 0; tracks with t.J <= 0 keep
    the full information with weight wU, and so do tracks with useeu False. Weights recomputed at every iteration. Returns x, C"""
    n = par.shape[0]
    x = x0.copy()
    Ls = np.zeros(n)
    X = np.zeros(3); D = np.zeros((3, 3)); a = np.zeros(3); A = np.zeros((3, 5)); X2 = np.zeros(3)
    H = Binv.copy()
    for it in range(niter):
        H = Binv.copy()
        g = Binv @ (x - c)
        for t in range(n):
            if wU[t] < 1e-4 and wEb[t] < 1e-4:
                continue
            c2, L = TR.lin1(par[t], cov[t], x, Ls[t], X, D)
            Ls[t] = L
            TR.helix1(par[t], L, X2, a, A)
            tv = a / np.sqrt(a @ a)
            j = jet_of[t]
            if j < 0:
                best = -2.
                for k in range(axes.shape[0]):
                    d = axes[k] @ tv
                    if d > best:
                        best = d; j = k
            J = axes[j]
            cn = tv @ J
            s2_ = 1. - cn * cn
            r = X - x
            if cn <= 0 or s2_ < 1e-10 or not useeu[t]:
                H += wU[t] * D
                g += wU[t] * (D @ (x - X))
                continue
            s_ = np.sqrt(s2_)
            nu = np.arctan2(s_, cn)
            Ev = np.cross(tv, J) / s_
            Uv = (J - tv * cn) / s_
            dee = Ev @ D @ Ev; deu = Ev @ D @ Uv; duu = Uv @ D @ Uv
            det = dee * duu - deu * deu
            if not det > 0:
                H += wU[t] * D
                g += wU[t] * (D @ (x - X))
                continue
            cee = duu / det; cuu = dee / det; ceu = -deu / det
            l3 = r @ Ev; l4 = r @ Uv
            jet2 = min(l4 * l4, sj2[j] * (l4 / nu) ** 2)
            cee2 = cee + jet2
            zu = l4 / np.sqrt(cuu)
            zk = kf if l4 > 0 else kb
            wu = wU[t] * _phi(zu) / (_phi(zu) + _phi(zk))
            ze = l3 / np.sqrt(cee2)
            we = wEb[t] * _phi(ze) / (_phi(ze) + _phi(ke))
            d2 = cee2 * cuu - ceu * ceu
            if not d2 > 0:
                continue
            iee = cuu / d2; iuu = cee2 / d2; ieu = -ceu / d2
            sw = np.sqrt(we * wu)
            I = (we * iee) * np.outer(Ev, Ev) + (wu * iuu) * np.outer(Uv, Uv) + (sw * ieu) * (np.outer(Ev, Uv) + np.outer(Uv, Ev))
            H += I
            g += I @ (x - X)
        dx = -np.linalg.solve(H, g)
        if not np.all(np.isfinite(dx)):
            break
        x = x + dx
        if np.max(np.abs(dx)) < tol:
            break
    C = np.linalg.inv(H)
    return x, C


def fit_vertex(par, cov, idx, x0, w=None, bs=None, bs_sig=BS_SIG, Ls=None, niter=40, tol=1e-10):
    idx = np.asarray(idx, np.int64)
    w = np.ones(len(idx)) if w is None else np.asarray(w, np.float64)
    if bs is not None:
        Binv = np.diag(1. / np.asarray(bs_sig) ** 2); c = np.asarray(bs, np.float64)
    else:
        Binv = np.zeros((3, 3)); c = np.zeros(3)
    Ls = np.zeros(len(par)) if Ls is None else Ls
    return fit_vertex_nb(par, cov, idx, w, np.asarray(x0, np.float64), Binv, c, Ls, niter, tol)


def teardown(par, cov, idx, x0, bs=None, bs_sig=BS_SIG, chi2max=5.0, nmin=2):
    """fit, then drop the worst track while its chi2 >= chi2max and more than nmin tracks remain. Returns x, C, kept idx,
    chi2 of kept"""
    idx = list(np.asarray(idx, np.int64))
    Ls = np.zeros(len(par))
    if len(idx) == 0:
        if bs is None:
            return np.asarray(x0, float), np.eye(3), np.zeros(0, np.int64), np.zeros(0)
        return np.asarray(bs, float), np.diag(np.asarray(bs_sig) ** 2), np.zeros(0, np.int64), np.zeros(0)
    x = np.asarray(x0, float)
    while True:
        x, C, chi2, _ = fit_vertex(par, cov, idx, x, bs=bs, bs_sig=bs_sig, Ls=Ls)
        j = int(np.argmax(chi2))
        if not (chi2[j] >= chi2max) or len(idx) - 1 < nmin:
            return x, C, np.array(idx, np.int64), chi2
        idx.pop(j)


def pair_fit(par, cov, i, j, rmax_i=np.inf, rmax_j=np.inf, x0=None, use_crossings=True):
    """two-track vertex: fits from the perigee start x0 (default: midpoint of the perigee points) and from each circle
    crossing; returns (x, C, chi2) of the lowest-chi2 solution inside min(rmax_i, rmax_j) (transverse radius from the
    detector axis), else the perigee-start solution; None if nothing finite"""
    rmax = min(rmax_i, rmax_j)
    if x0 is None:
        x0 = 0.5 * (TR.point_at(par[i], 0.) + TR.point_at(par[j], 0.))
    starts = [(np.asarray(x0, float), 0., 0.)]
    if use_crossings:
        for (cx, cy) in TR.circle_crossings(par[i], par[j]):
            Li, Lj = TR.arc_to_xy(par[i], cx, cy), TR.arc_to_xy(par[j], cx, cy)
            zi, zj = par[i][3] + par[i][4] * Li, par[j][3] + par[j][4] * Lj
            starts.append((np.array([cx, cy, 0.5 * (zi + zj)]), Li, Lj))
    best, first = None, None
    idx = np.array([i, j], np.int64)
    for k, (xs, Li, Lj) in enumerate(starts):
        Ls = np.zeros(len(par)); Ls[i] = Li; Ls[j] = Lj
        x, C, chi2, conv = fit_vertex(par, cov, idx, xs, Ls=Ls, niter=15, tol=1e-6)
        tot = float(chi2.sum())
        if not np.isfinite(tot) or not np.all(np.isfinite(x)):
            continue
        sol = (x, C, tot)
        if k == 0:
            first = sol
        if np.hypot(x[0], x[1]) < rmax and (best is None or tot < best[2]):
            best = sol
    return best if best is not None else first


def jetproj_pv(par_raw, cov15, jet_of, jet_theta, jet_phi, jet_p, bs, x0=None, form='hard', bs_sig=BS_SIG, kf=1.5, kb=3.0, ke=3.0):
    """jet-projected PV of the tracks given (raw states), QFNDIP rule (findip.F; hard) or its soft correlated form
    (soft_corr, kf 1.5 kb 3 ke 3). Returns x, C (3x3), number of tracks used"""
    L = JP
    q = L.line_params(np.asarray(par_raw), np.asarray(cov15))
    trk = tuple(q[:, k] for k in range(7))
    axes = L.jet_axes(jet_theta, jet_phi)
    sg = [L.jet_sig2_aleph(p) for p in jet_p]
    bsW = np.diag(1. / np.asarray(bs_sig) ** 2)
    bsc = np.asarray(bs, float)
    Q = L.prep_qfn(trk, axes, np.asarray(jet_of, np.int64))
    if form == 'hard':
        return L.qfndip(Q, sg, bsW, bsc, 0.0)
    Ceu = L.eu_cov(trk, Q)
    xs = bsc if x0 is None else np.asarray(x0, float)
    return L.soft_corr(Q, Ceu, sg, bsW, bsc, xs, kf=kf, kb=kb, ke=ke)


def inv_mass(par, idx, x, Ls=None, masses=None):
    """invariant mass (GeV) of the tracks idx with momenta at their PCA to x (pion mass unless masses given)"""
    E = 0.; P = np.zeros(3)
    X = np.zeros(3); D = np.zeros((3, 3))
    for k, t in enumerate(idx):
        m = TR.M_PI if masses is None else masses[k]
        _, L = TR.lin1(par[t], np.eye(5) * 1e-6 if Ls is None else Ls, x, 0., X, D) if False else (0., None)
        L = _arc_near(par[t], x)
        p = TR.momentum_at(par[t], L)
        E += np.sqrt(p @ p + m * m); P += p
    return float(np.sqrt(max(E * E - P @ P, 0.))), P


def _arc_near(p, x):
    """transverse arc of the track's xy point nearest to x (geometric, no covariance)"""
    if p[2] == 0:
        return 0.
    return TR.arc_to_xy(p, x[0], x[1])
