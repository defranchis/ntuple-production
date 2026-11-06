"""Seedless candidate supply for the joint vertexer: drop-in replacements of the stages S2 + S3 of jv.hybrid.

Two engines, both on all S0 tracks (no IP-significance selection of tracks, no pair seed that must survive):

vd   vertex-density function (ZVTOP idea). Track tube f_i(r) = exp(-chi2_i(r) / (2 tube)), chi2_i = tube chi2 of track i at the
     point r (tracks.lin1), IP term f_0(r) = k0 exp(-(r - PV)^T C^-1 (r - PV) / 2). V(r) = sum f - sum f^2 / sum f. Pair crossings
     only serve as the search grid: every grid point climbs V by a fixed-point iteration, and the climbs give the maxima.
     A maximum is resolved from the PV if the lowest V on the straight segment to the PV over the lower of the two maxima is
     below vd_R0 (ZVTOP's resolution criterion).
da   3-D deterministic annealing (CMS DA generalised to anisotropic track metrics). The event starts as ONE vertex (beam prior).
     Free energy F = -T sum_k log(Z + sum_v rho_v exp(-chi2_kv / 2T)). A vertex is unstable against splitting along delta when
     delta^T (A - B / T) delta < 0, A = sum_k w_kv D_k (+ prior), B = sum_k w_kv g_k g_k^T, g_k = D_k (x_v - X_k): critical temperature
     T_c = largest eigenvalue of A^-1 B. On reaching T_c the vertex splits along the eigenvector; the two copies start at the
     1-D optima of the tracks pulling to either side. Copies merge again if the critical temperature of their union is below T.

Units cm, GeV.
"""
import numpy as np
from numba import njit

from . import tracks as TR

BIG = 1e6


@njit(cache=True)
def _vd_value(chi2, tube, f0):
    s1 = f0; s2 = f0 * f0
    for k in range(len(chi2)):
        f = np.exp(-0.5 * chi2[k] / tube)
        s1 += f; s2 += f * f
    return (s1 - s2 / s1) if s1 > 1e-100 else 0.


@njit(cache=True)
def vd_maxima(par, cov, rmax, pvx, pvCi, k0, tube, grid_chi2, vmin, n_climb, max_starts, fwd_margin, member_chi2, grid_pvchi2, grid_cos):
    """Returns maxima positions (m, 3), V values (m), dip ratio to the PV (m), chi2 (m, nt) of every track at the maximum,
    number of climbs"""
    nt = par.shape[0]
    Ls0 = np.zeros(nt); X0 = np.zeros((nt, 3)); D0 = np.zeros((nt, 3, 3)); c0 = np.zeros(nt); a0 = np.zeros((nt, 3))
    a_ = np.zeros(3); A_ = np.zeros((3, 5)); Xt = np.zeros(3)
    for k in range(nt):
        c0[k], Ls0[k] = TR.lin1(par[k], cov[k], pvx, 0., X0[k], D0[k])
        TR.helix1(par[k], Ls0[k], Xt, a_, A_)
        a0[k] = a_ / np.sqrt(a_ @ a_)
    # search grid: linearised pair crossings
    nmax = nt * (nt - 1) // 2
    si = np.zeros(nmax, np.int64); sj = np.zeros(nmax, np.int64); sx = np.zeros((nmax, 3)); sc = np.zeros(nmax)
    ns = 0
    eps = np.eye(3)     # (1 cm)^-2: keeps nearly parallel pairs from running away
    for i in range(nt):
        for j in range(i + 1, nt):
            if max(c0[i], c0[j]) < grid_pvchi2 or a0[i] @ a0[j] < grid_cos:
                continue        # two tracks inside the PV tube cannot climb away from it; back-to-back pairs give no vertex
            H = D0[i] + D0[j] + eps
            x = np.linalg.solve(H, D0[i] @ (X0[i] - pvx) + D0[j] @ (X0[j] - pvx)) + pvx
            ri = x - X0[i]; rj = x - X0[j]
            ch = ri @ D0[i] @ ri + rj @ D0[j] @ rj
            if ch > grid_chi2:
                continue
            d = x - pvx
            if d @ d > 3.5 * 3.5:
                continue
            if np.hypot(x[0], x[1]) >= min(rmax[i], rmax[j]):
                continue
            if a0[i] @ d < -fwd_margin or a0[j] @ d < -fwd_margin:
                continue
            # score: V with the tubes linearised at the PV
            s1 = 0.; s2 = 0.
            for k in range(nt):
                r = x - X0[k]
                f = np.exp(-0.5 * (r @ D0[k] @ r) / tube)
                s1 += f; s2 += f * f
            si[ns] = i; sj[ns] = j; sx[ns] = x; sc[ns] = (s1 - s2 / s1) if s1 > 1e-100 else 0.
            ns += 1
    order = np.argsort(-sc[:ns])
    covered = np.zeros((nt, nt), np.bool_)
    mx = np.zeros((max_starts, 3)); mv = np.zeros(max_starts); mc = np.zeros((max_starts, nt))
    nm = 0; nclimb = 0
    Ls = np.zeros(nt); X = np.zeros((nt, 3)); D = np.zeros((nt, 3, 3)); ch = np.zeros(nt); f = np.zeros(nt)
    for s_ in order:
        i, j = si[s_], sj[s_]
        if covered[i, j]:
            continue
        if nclimb >= max_starts:
            break
        nclimb += 1
        x = sx[s_].copy()
        Ls[:] = Ls0
        for it in range(n_climb):
            s1 = 0.; s2 = 0.
            d = x - pvx
            f0 = k0 * np.exp(-0.5 * (d @ pvCi @ d))
            for k in range(nt):
                ch[k], Ls[k] = TR.lin1(par[k], cov[k], x, Ls[k], X[k], D[k])
                f[k] = np.exp(-0.5 * ch[k] / tube)
                s1 += f[k]; s2 += f[k] * f[k]
            s1 += f0; s2 += f0 * f0
            if s1 <= 1e-100:
                break
            H = np.zeros((3, 3)); g = np.zeros(3)
            for k in range(nt):
                p = f[k] / s1
                a = f[k] * (1. - 2. * p + s2 / (s1 * s1)) / tube
                if a > 1e-8:
                    H += a * D[k]; g += a * (D[k] @ (X[k] - pvx))
            p = f0 / s1
            a = f0 * (1. - 2. * p + s2 / (s1 * s1))
            H += a * pvCi
            H += 1e-6 * eps
            xn = np.linalg.solve(H, g) + pvx
            if not np.all(np.isfinite(xn)):
                break
            dx = xn - x
            nd = np.sqrt(dx @ dx)
            if nd > 0.3:        # no step control in a fixed point: limit the step (3 mm)
                xn = x + dx * 0.3 / nd
            x = xn
            if nd < 2e-4:
                break
        d = x - pvx
        f0 = k0 * np.exp(-0.5 * (d @ pvCi @ d))
        for k in range(nt):
            ch[k], Ls[k] = TR.lin1(par[k], cov[k], x, Ls[k], X[k], D[k])
        v = _vd_value(ch, tube, f0)
        mem = np.nonzero(ch < member_chi2)[0]
        for a1 in mem:
            for a2 in mem:
                covered[a1, a2] = True
        covered[i, j] = True
        if v < vmin or len(mem) < 2:
            continue
        # duplicate of an earlier maximum?
        dup = False
        for m in range(nm):
            dm = x - mx[m]
            if np.sqrt(dm @ dm) < 5e-4:
                dup = True
                break
        if dup:
            continue
        mx[nm] = x; mv[nm] = v; mc[nm] = ch; nm += 1
        if nm == max_starts:
            break
    # PV value and dip along the straight segment to the PV
    vpv = _vd_value(c0, tube, k0)
    dip = np.zeros(nm)
    for m in range(nm):
        lo = min(mv[m], vpv)
        vm = mv[m]
        for s, frac in enumerate((0.15, 0.3, 0.5, 0.7, 0.85)):
            xs = mx[m] + frac * (pvx - mx[m])
            d = xs - pvx
            f0 = k0 * np.exp(-0.5 * (d @ pvCi @ d))
            for k in range(nt):
                ch[k], _ = TR.lin1(par[k], cov[k], xs, Ls0[k], X[k], D[k])
            vm = min(vm, _vd_value(ch, tube, f0))
        dip[m] = vm / lo if lo > 1e-100 else 1.
    return mx[:nm], mv[:nm], dip, mc[:nm], nclimb


@njit(cache=True)
def _tc(A, B):
    """largest eigenvalues of A^-1 B and the matching directions in space (columns, largest first); A symmetric, floored"""
    S = 0.5 * (A + A.T)
    if not (np.all(np.isfinite(S)) and np.all(np.isfinite(B))):
        return np.zeros(3), np.eye(3)
    ea, Ua = np.linalg.eigh(S)
    fl = max(1e-12 * ea[2], 1e-12)
    Li = np.zeros((3, 3))
    for i in range(3):
        Li[i] = Ua[:, i] / np.sqrt(max(ea[i], fl))
    M = Li @ B @ Li.T
    M = 0.5 * (M + M.T)
    e, U = np.linalg.eigh(M)
    dirs = Li.T @ U
    return e[::-1], dirs[:, ::-1]


@njit(cache=True)
def da_run(par, cov, rmax, x0, Binv, cbs, T0, cool, Tmin, chi2_cut, n_em, tol, max_v, use_rho, prior, split_margin,
           merge_margin, min_side, fwd, Tpurge, wside, noise_T):
    """Deterministic annealing from one vertex. Returns x (nv, 3), C (nv, 3, 3), w (nv, nt), chi2 (nv, nt), n splits, n lin"""
    nt = par.shape[0]
    x = np.zeros((max_v, 3)); rho = np.zeros(max_v)
    Ls = np.zeros((max_v, nt)); X = np.zeros((max_v, nt, 3)); D = np.zeros((max_v, nt, 3, 3)); ch = np.zeros((max_v, nt))
    w = np.zeros((max_v, nt)); A = np.zeros((max_v, 3, 3))
    x[0] = x0; rho[0] = 1.; nv = 1
    a_ = np.zeros(3); A_ = np.zeros((3, 5))
    eps = 1e-2 * np.eye(3)   # prior of (10 cm)^-2 on SVs: keeps one-track phases of a vertex finite
    T = T0
    nsplit = 0; nlin = 0
    while True:
        Z = np.exp(-0.5 * chi2_cut / T) if T <= noise_T else 1e-30   # as CMS DA: no outlier term during the splitting phase
        for em in range(n_em):
            # E step
            den = np.full(nt, Z)
            for v in range(nv):
                rad = np.hypot(x[v, 0], x[v, 1])
                dv = x[v] - x[0]
                for k in range(nt):
                    c2, L = TR.lin1(par[k], cov[k], x[v], Ls[v, k], X[v, k], D[v, k])
                    Ls[v, k] = L
                    nlin += 1
                    if v > 0 and rad >= rmax[k]:
                        c2 = BIG
                    elif v > 0 and fwd:
                        ang = par[k, 1] + 2 * par[k, 2] * L
                        if np.cos(ang) * dv[0] + np.sin(ang) * dv[1] + par[k, 4] * dv[2] < 0:
                            c2 = BIG
                    ch[v, k] = c2
                    e = np.exp(-0.5 * min(c2, 1e4) / T)
                    if use_rho:
                        e *= rho[v]
                    if v == 0:
                        e *= prior[k]
                    w[v, k] = e
                    den[k] += e
            # M step
            shift = 0.
            tot = 0.
            for v in range(nv):
                H = (Binv.copy() if v == 0 else eps.copy()); g = (Binv @ (cbs - x[0]) if v == 0 else np.zeros(3))
                sw = 0.
                for k in range(nt):
                    w[v, k] /= den[k]
                    wk = w[v, k]
                    if wk > 1e-5:
                        H += wk * D[v, k]; g += wk * (D[v, k] @ (X[v, k] - x[v]))
                        sw += wk
                dx = np.linalg.solve(H, g)
                A[v] = H
                if v == 0 and False:
                    pass
                nd = np.sqrt(dx @ dx)
                if nd > 0.3:
                    dx *= 0.3 / nd
                    nd = 0.3
                if np.all(np.isfinite(dx)):
                    x[v] += dx
                shift = max(shift, nd)
                rho[v] = sw
                tot += sw
            if tot > 0:
                rho[:nv] /= tot
            if shift < tol:
                break
        # merge: union of two vertices stable at this T
        changed = True
        while changed and nv > 1:
            changed = False
            for a in range(nv):
                for b in range(a + 1, nv):
                    ra, rb = rho[a], rho[b]
                    if ra + rb <= 0:
                        xm = 0.5 * (x[a] + x[b])
                    elif a == 0:
                        xm = x[0].copy()
                    else:
                        xm = (ra * x[a] + rb * x[b]) / (ra + rb)
                    d = x[a] - x[b]
                    if d @ d > 0.04:     # beyond 2 mm no merge by this test (S3/S4 merging decides)
                        continue
                    Hm = (Binv + eps) if a == 0 else eps.copy()
                    Bm = np.zeros((3, 3))
                    for k in range(nt):
                        for v_, wk in ((a, w[a, k]), (b, w[b, k])):
                            if wk > 1e-5:
                                r = xm - X[v_, k]
                                gk = D[v_, k] @ r
                                Hm += wk * D[v_, k]
                                Bm += wk * np.outer(gk, gk)
                    e, _d = _tc(Hm, Bm)
                    if e[0] < T * merge_margin:
                        x[a] = xm
                        rho[a] = ra + rb
                        w[a] = w[a] + w[b]
                        last = nv - 1
                        if b != last:
                            x[b] = x[last]; rho[b] = rho[last]; Ls[b] = Ls[last]; X[b] = X[last]; D[b] = D[last]
                            ch[b] = ch[last]; w[b] = w[last]; A[b] = A[last]
                        nv -= 1
                        changed = True
                        break
                if changed:
                    break
        # purge: below Tpurge a non-PV vertex needs 2 tracks
        if T <= Tpurge:
            v = 1
            while v < nv:
                n2 = 0
                for k in range(nt):
                    if w[v, k] > 0.5:
                        n2 += 1
                if n2 < 2:
                    last = nv - 1
                    if v != last:
                        x[v] = x[last]; rho[v] = rho[last]; Ls[v] = Ls[last]; X[v] = X[last]; D[v] = D[last]
                        ch[v] = ch[last]; w[v] = w[last]; A[v] = A[last]
                    nv -= 1
                else:
                    v += 1
        if T <= Tmin:
            break
        # split
        nv0 = nv
        for v in range(nv0):
            if nv >= max_v:
                break
            Bv = np.zeros((3, 3))
            gk = np.zeros((nt, 3))
            for k in range(nt):
                wk = w[v, k]
                if wk > 1e-5:
                    gk[k] = D[v, k] @ (X[v, k] - x[v])   # pull of track k on the vertex
                    Bv += wk * np.outer(gk[k], gk[k])
            e, dirs = _tc(A[v], Bv)
            nextT = max(T * cool, Tmin)
            for q in range(3):
                if e[q] < nextT * split_margin:
                    break
                de = dirs[:, q]
                # the two sides
                n1 = 0; n2 = 0; a1n = 0.; a1d = 0.; a2n = 0.; a2d = 0.; w1 = 0.; w2 = 0.
                for k in range(nt):
                    wk = w[v, k]
                    if wk > 1e-5:
                        gd = gk[k] @ de
                        dd = de @ D[v, k] @ de
                        if gd > 0:
                            a1n += wk * gd; a1d += wk * dd; w1 += wk
                            if wk > wside:
                                n1 += 1
                        else:
                            a2n += wk * gd; a2d += wk * dd; w2 += wk
                            if wk > wside:
                                n2 += 1
                if n1 < min_side or n2 < min_side or a1d <= 0 or a2d <= 0:
                    continue
                a1 = a1n / a1d; a2 = a2n / a2d
                nw = nv
                # the heavier side stays in place of v (and keeps the prior when v is the PV)
                xo = x[v].copy()
                if w1 >= w2:
                    x[v] = xo + a1 * de; x[nw] = xo + a2 * de
                    rho[nw] = rho[v] * w2 / (w1 + w2); rho[v] = rho[v] * w1 / (w1 + w2)
                else:
                    x[v] = xo + a2 * de; x[nw] = xo + a1 * de
                    rho[nw] = rho[v] * w1 / (w1 + w2); rho[v] = rho[v] * w2 / (w1 + w2)
                Ls[nw] = Ls[v]
                nv += 1
                nsplit += 1
                break
        T = max(T * cool, Tmin)
    # final state at the final T
    C = np.zeros((nv, 3, 3))
    for v in range(nv):
        C[v] = np.linalg.inv(A[v])
    return x[:nv].copy(), C, w[:nv].copy(), ch[:nv].copy(), nsplit, nlin
