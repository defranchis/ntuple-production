"""Photon-conversion finder of the joint vertexer (gamma -> e+ e- in the detector material), numba kernels.

A photon converts into two electrons that start at the same point with nearly the same direction (opening angle of
order m_e / E). In the transverse plane the two track circles (opposite curvature) are then tangent at the conversion
point. Two formulations of the conversion vertex:

TP (tangent point, closed form): the point on the line joining the two circle centres, in the middle of the gap (or of
   the overlap) between the circles. There the two track directions are parallel in xy by construction. Quality: the
   signed transverse gap between the circles (> 0 separated, < 0 overlapping), the difference of the z of the two tracks
   at that point and the difference of their dip angles (tanLambda); errors propagated from the track covariances
   (numerical Jacobian); vertex z = error-weighted mean of the two z.
CF (constrained conversion fit): 7 parameters (vertex x, y, z, common direction phi, common tanLambda, the two
   curvatures) fitted to the 10 measured helix parameters of the two tracks (Levenberg-Marquardt with a second-order
   check at the solution, numerical derivatives, started at TP; see pair_cf): a common vertex with collinear daughters,
   3 degrees of freedom, m_ee = 0 by construction. The chi2 splits exactly into the two tracks' contributions.

Track model and conventions: jv/tracks.py (fitter convention p = (D, phi0, C, z0, ct), helix point at transverse arc L,
direction phi0 + 2 C L, circle centre (-(D + 0.5 / C) sin phi0, (D + 0.5 / C) cos phi0), radius 0.5 / |C|); the helix
kernels are repeated here as numba functions of this module.
Momentum constant: pT = K / |C| = 2 K / |omega_raw|; K is an argument (k_pt; jv.v0 passes PT_CONST / 2 of the config key
pt_const), default K_PT = PT_CONSTS['stage1'] / 2 = 0.0029979 x 1.5 / 2 (jv/tracks.py).
Units: cm, GeV.

Interface (per event; par, cov = (n, 5), (n, 5, 5) in the fitter convention of jv.tracks.to_fit; q = charges):
  pair_tp(p1, c1, p2, c2, K)                 -> TP result vector (TP_NAMES; ok = 0 when undefined (same curvature
                                                sign: no tangent point), all other entries NaN then)
  pair_cf(p1, c1, p2, c2, tp, K)             -> CF result vector (CF_NAMES; ok 1 converged to a local minimum, 2 not
                                                converged, 0 undefined / failed with NaN entries; chi2 = chi2_1 +
                                                chi2_2, the two tracks' contributions), vertex covariance (3, 3; NaN
                                                for ok = 0)
  find_pairs(par, cov, q, ok=None, pre=PRE)  -> (i, j) index arrays of opposite-charge pairs passing the prefilter PRE
  pair_quality(par, cov, i, j, pv, rmax, k_pt=K_PT) -> quality matrix (n_pairs, Q_NAMES) + CF vertex covariances
                                                (n_pairs, 3, 3)
  select(Q, LOOSE | TIGHT | dict)            -> mask of the rows passing a selection
  arbitrate(ev, i, j, score)                 -> mask: greedy one-candidate-per-track (lowest score first)
  find_conversions(par, cov, q, nvdet, nitc, pv, sel='loose' | 'tight' | dict | None, k_pt=K_PT)
                                             -> dict i, j, x (CF vertex), cov (3, 3), Q of one event's candidates
Quality matrix columns (Q_NAMES): tp_* (TP vector), cf_* (CF vector), tp_r / cf_r (vertex radius from the detector axis),
  tp_chi2 (gap, dz, dct with their propagated covariance, 3 dof), tp_dgam / cf_dgam (signed transverse distance of the
  photon line from the PV), tp_cosp / cf_cosp (cos of the 3D angle between PV -> vertex and the photon momentum), cf_sl
  (CF vertex error along the photon), nhx1 / nhx2 (1 when the TP radius is above the hit-radius bound of track i / j,
  jv.tracks.max_radius_from_hits), p1 / p2 (track momenta at the TP). Every column is NaN where its TP / CF is undefined
  (tp_ok / cf_ok = 0).
Selections (chosen on the simulation file mc_ZM4212.50): LOOSE = supply to the joint fit (CF chi2 < 200, m_ee < 50 MeV,
  cos pointing > 0.9; no radius prior on the detector material); TIGHT = stand-alone list, followed by arbitrate() on the
  CF chi2.
"""
import numpy as np
from numba import njit

from . import tracks as TR

K_PT = 0.5 * TR.PT_CONSTS['stage1']    # GeV cm: default pT = K_PT / |C| (C = omega_raw / 2 in the fitter convention)
M_E = 0.000510999
TWO_PI = 2. * np.pi

TP_NAMES = ['ok', 'x', 'y', 'z', 'gap', 'dz', 'dct', 'L1', 'L2', 'phi', 'ct', 'sxx', 'syy', 'szz', 'sxy', 'sxz', 'syz',
            'sgap', 'sdz', 'sdct', 'rgz', 'rgc', 'rzc', 'mee', 'px', 'py', 'pz']
CF_NAMES = ['ok', 'x', 'y', 'z', 'phi', 'ct', 'C1', 'C2', 'chi2', 'niter', 'px', 'py', 'pz', 'chi2_1', 'chi2_2', 'dec']
N_TP, N_CF = len(TP_NAMES), len(CF_NAMES)


# ----------------------------------------------------------------------------------------------- helix geometry
@njit(cache=True)
def _wrap(a):
    return (a + np.pi) % TWO_PI - np.pi


@njit(cache=True)
def _centre(p):
    rs = 0.5 / p[2]
    return -(p[0] + rs) * np.sin(p[1]), (p[0] + rs) * np.cos(p[1]), rs


@njit(cache=True)
def _arc_to(p, x, y):
    """transverse arc L of the circle point of track p nearest to (x, y) (|2 C L| < pi)"""
    xc, yc, rs = _centre(p)
    x0, y0 = -p[0] * np.sin(p[1]), p[0] * np.cos(p[1])
    return _wrap(np.arctan2(y - yc, x - xc) - np.arctan2(y0 - yc, x0 - xc)) / (2. * p[2])


@njit(cache=True)
def _point(p, L, X):
    D, p0, C = p[0], p[1], p[2]
    u = C * L
    f = L * (np.sin(u) / u if abs(u) > 1e-6 else 1. - u * u / 6.)
    ang = p0 + u
    X[0] = -D * np.sin(p0) + f * np.cos(ang)
    X[1] = D * np.cos(p0) + f * np.sin(ang)
    X[2] = p[3] + p[4] * L


@njit(cache=True)
def _mom(p, L, K, P):
    pt = K / max(abs(p[2]), 1e-12)
    ang = p[1] + 2. * p[2] * L
    P[0] = pt * np.cos(ang); P[1] = pt * np.sin(ang); P[2] = pt * p[4]


@njit(cache=True)
def perigee_from_ref(xv, yv, zv, phi, ct, C, Lref, out):
    """fitter-convention parameters (D, phi0, C, z0, ct) of the helix through (xv, yv, zv) with direction phi there; the
    transverse arc L from the perigee to the vertex is taken on the turn of the helix nearest to the reference arc Lref:
    continuous in the vertex position near a half turn"""
    rs = 0.5 / C
    xc = xv - rs * np.sin(phi); yc = yv + rs * np.cos(phi)
    ra = np.hypot(xc, yc)
    rho = ra if C > 0 else -ra
    phi0 = np.arctan2(-xc / rho, yc / rho)
    D = rho - rs
    x0, y0 = -D * np.sin(phi0), D * np.cos(phi0)
    a = np.arctan2(yv - yc, xv - xc) - np.arctan2(y0 - yc, x0 - xc)
    a0 = 2. * C * Lref
    a = a0 + _wrap(a - a0)
    L = a / (2. * C)
    out[0] = D; out[1] = phi0; out[2] = C; out[3] = zv - ct * L; out[4] = ct


# ----------------------------------------------------------------------------------------------- TP: tangent point
@njit(cache=True)
def _tp_raw(p1, p2, o):
    """o = (x, y, z1, z2, gap, dct, L1, L2); returns False if undefined (zero curvature, same curvature sign, d = 0)"""
    if p1[2] == 0. or p2[2] == 0. or p1[2] * p2[2] > 0.:
        return False
    x1, y1, r1 = _centre(p1)
    x2, y2, r2 = _centre(p2)
    R1, R2 = abs(r1), abs(r2)
    dx, dy = x2 - x1, y2 - y1
    d = np.hypot(dx, dy)
    if not d > 0.:
        return False
    ux, uy = dx / d, dy / d
    gap = d - R1 - R2
    t = R1 + 0.5 * gap
    x, y = x1 + t * ux, y1 + t * uy
    L1, L2 = _arc_to(p1, x, y), _arc_to(p2, x, y)
    o[0] = x; o[1] = y
    o[2] = p1[3] + p1[4] * L1; o[3] = p2[3] + p2[4] * L2
    o[4] = gap; o[5] = p1[4] - p2[4]; o[6] = L1; o[7] = L2
    return True


@njit(cache=True)
def _step(p, k):
    if k == 2:
        return 1e-6 * max(abs(p[2]), 1e-5)
    if k == 1 or k == 4:
        return 1e-7
    return 1e-6


@njit(cache=True)
def pair_tp(p1, c1, p2, c2, K):
    """TP result vector (TP_NAMES): vertex, gap, dz, dct, arcs, photon direction, vertex covariance (6 numbers), errors
    of gap / dz / dct and their correlations, m_ee and photon momentum at the TP (track momenta at the TP arcs)"""
    r = np.full(N_TP, np.nan)
    r[0] = 0.
    o = np.zeros(8)
    if not _tp_raw(p1, p2, o):
        return r
    # Jacobian of (x, y, z1, z2, gap) w.r.t. the 10 parameters (dct is linear)
    J1 = np.zeros((5, 5)); J2 = np.zeros((5, 5))
    q = np.zeros(8); qm = np.zeros(8)
    for side in range(2):
        for k in range(5):
            pp = (p1 if side == 0 else p2).copy()
            pm = pp.copy()
            h = _step(pp, k)
            pp[k] += h; pm[k] -= h
            if side == 0:
                okp = _tp_raw(pp, p2, q); okm = _tp_raw(pm, p2, qm)
            else:
                okp = _tp_raw(p1, pp, q); okm = _tp_raw(p1, pm, qm)
            if not (okp and okm):
                return r
            for a in range(5):
                v = (q[a] - qm[a]) / (2. * h)
                if side == 0:
                    J1[a, k] = v
                else:
                    J2[a, k] = v
    V = J1 @ c1 @ J1.T + J2 @ c2 @ J2.T              # cov of (x, y, z1, z2, gap)
    # dct = ct1 - ct2: covariances with the other outputs
    g1 = np.zeros(5); g1[4] = 1.
    g2 = np.zeros(5); g2[4] = -1.
    Vc = J1 @ c1 @ g1 + J2 @ c2 @ g2                 # cov(out_a, dct)
    vcc = g1 @ c1 @ g1 + g2 @ c2 @ g2
    # vertex z = w z1 + (1 - w) z2 with minimal variance
    s11, s22, s12 = V[2, 2], V[3, 3], V[2, 3]
    den = s11 + s22 - 2. * s12
    w = (s22 - s12) / den if den > 0 else 0.5
    T = np.zeros((3, 5)); T[0, 0] = 1.; T[1, 1] = 1.; T[2, 2] = w; T[2, 3] = 1. - w
    Cx = T @ V @ T.T
    # dz = z1 - z2
    gz = np.array([0., 0., 1., -1., 0.])
    vzz = gz @ V @ gz
    vgg = V[4, 4]
    vgz = V[4] @ gz
    vgc = Vc[4]
    vzc = gz @ Vc
    L1, L2 = o[6], o[7]
    P1 = np.zeros(3); P2 = np.zeros(3)
    _mom(p1, L1, K, P1); _mom(p2, L2, K, P2)
    E = np.sqrt(P1 @ P1 + M_E * M_E) + np.sqrt(P2 @ P2 + M_E * M_E)
    P = P1 + P2
    r[0] = 1.
    r[1] = o[0]; r[2] = o[1]; r[3] = w * o[2] + (1. - w) * o[3]
    r[4] = o[4]; r[5] = o[2] - o[3]; r[6] = o[5]; r[7] = L1; r[8] = L2
    r[9] = np.arctan2(P[1], P[0]); r[10] = P[2] / max(np.hypot(P[0], P[1]), 1e-12)
    r[11] = Cx[0, 0]; r[12] = Cx[1, 1]; r[13] = Cx[2, 2]; r[14] = Cx[0, 1]; r[15] = Cx[0, 2]; r[16] = Cx[1, 2]
    r[17] = np.sqrt(max(vgg, 0.)); r[18] = np.sqrt(max(vzz, 0.)); r[19] = np.sqrt(max(vcc, 0.))
    r[20] = vgz / max(r[17] * r[18], 1e-300); r[21] = vgc / max(r[17] * r[19], 1e-300); r[22] = vzc / max(r[18] * r[19], 1e-300)
    r[23] = np.sqrt(max(E * E - P @ P, 0.))
    r[24] = P[0]; r[25] = P[1]; r[26] = P[2]
    return r


# ----------------------------------------------------------------------------------------------- CF: constrained fit
@njit(cache=True)
def _cf_pred(th, Lr, h):
    """predicted parameters of the two tracks (h[0:5], h[5:10]) from th = (x, y, z, phi, ct, C1, C2)"""
    a = np.zeros(5); b = np.zeros(5)
    perigee_from_ref(th[0], th[1], th[2], th[3], th[4], th[5], Lr[0], a)
    perigee_from_ref(th[0], th[1], th[2], th[3], th[4], th[6], Lr[1], b)
    for k in range(5):
        h[k] = a[k]; h[5 + k] = b[k]


@njit(cache=True)
def _cf_res(th, m, Lr, r):
    _cf_pred(th, Lr, r)
    for k in range(10):
        r[k] -= m[k]
    r[1] = _wrap(r[1]); r[6] = _wrap(r[6])


CF_MAXIT = 500          # maximum Levenberg-Marquardt iterations of the CF (large-residual pairs converge linearly)
CF_NESCAPE = 5          # maximum restarts of the CF from a saddle point (negative curvature of the chi2)
CF_NEG = 1e-3          # negative curvature: a Hessian eigenvalue of chi2 / 2 below -CF_NEG (parameters in units of their errors)
CF_TOL = 1e-6           # convergence: predicted chi2 decrease of a full Gauss-Newton step (g^T H^-1 g) below this


@njit(cache=True)
def _cf_jac(th, m, Lr, J, rp, rm):
    """numerical Jacobian (10, 7) of the CF residuals (central differences)"""
    for k in range(7):
        if k < 3:
            h = 1e-6
        elif k < 5:
            h = 1e-7
        else:
            h = 1e-6 * max(abs(th[k]), 1e-5)
        tq = th.copy(); tq[k] += h
        tm = th.copy(); tm[k] -= h
        _cf_res(tq, m, Lr, rp); _cf_res(tm, m, Lr, rm)
        for a in range(10):
            J[a, k] = _wrap(rp[a] - rm[a]) / (2. * h) if (a == 1 or a == 6) else (rp[a] - rm[a]) / (2. * h)


@njit(cache=True)
def _cf_decrement(H, g):
    """Newton decrement g^T H^-1 g (predicted chi2 decrease of the full Gauss-Newton step; inf if H is singular)"""
    try:
        d = np.linalg.solve(H, g)
    except Exception:
        return np.inf
    v = g @ d
    return v if np.isfinite(v) else np.inf


@njit(cache=True)
def _cf_lm(th, m, Lr, W, chi2):
    """Levenberg-Marquardt from th (chi2 = its chi2): returns (th, chi2, iterations)"""
    r = np.zeros(10); rp = np.zeros(10); rm = np.zeros(10)
    J = np.zeros((10, 7))
    _cf_res(th, m, Lr, r)
    lam = 1e-3
    it = 0
    for it in range(CF_MAXIT):
        _cf_jac(th, m, Lr, J, rp, rm)
        H = J.T @ W @ J
        g = J.T @ W @ r
        if _cf_decrement(H, g) < CF_TOL:
            return th, chi2, it + 1
        accepted = False
        for _ in range(40):
            A = H.copy()
            for k in range(7):
                A[k, k] *= 1. + lam
            try:
                d = -np.linalg.solve(A, g)
            except Exception:
                lam *= 10.
                continue
            tn = th + d
            tn[3] = _wrap(tn[3])
            _cf_res(tn, m, Lr, rp)
            cn = rp @ W @ rp
            if np.isfinite(cn) and cn < chi2:
                th = tn; chi2 = cn
                for a in range(10):
                    r[a] = rp[a]
                lam = max(lam * 0.2, 1e-9)
                accepted = True
                break
            lam *= 10.
        if not accepted:
            return th, chi2, it + 1
    return th, chi2, CF_MAXIT


@njit(cache=True)
def _cf_chi2(th, m, Lr, W, r):
    _cf_res(th, m, Lr, r)
    return r @ W @ r


@njit(cache=True)
def _cf_escape(th, m, Lr, W, chi2):
    """full numerical Hessian of the chi2 at th in units of the Gauss-Newton parameter errors s_k (step 0.01 s_k);
    returns (negative curvature found: lowest eigenvalue of chi2 / 2 < -CF_NEG, best point, its chi2) where the best
    point is the lowest chi2 of steps t s v (t = +-0.5, 1, 2, 4, 8, 16) along each negative-curvature eigenvector v"""
    r = np.zeros(10); rp = np.zeros(10); rm = np.zeros(10)
    J = np.zeros((10, 7))
    _cf_jac(th, m, Lr, J, rp, rm)
    H = J.T @ W @ J
    try:
        C = np.linalg.inv(H)
    except Exception:
        return False, th, chi2
    sg = np.zeros(7)
    for k in range(7):
        if not C[k, k] > 0:
            return False, th, chi2
        sg[k] = np.sqrt(C[k, k])
    h = 0.01
    F = np.zeros((7, 7))
    for a in range(7):
        for b in range(a, 7):
            if a == b:
                tq = th.copy(); tq[a] += h * sg[a]
                tm = th.copy(); tm[a] -= h * sg[a]
                v = (_cf_chi2(tq, m, Lr, W, r) - 2. * chi2 + _cf_chi2(tm, m, Lr, W, r)) / (h * h)
            else:
                t1 = th.copy(); t1[a] += h * sg[a]; t1[b] += h * sg[b]
                t2 = th.copy(); t2[a] += h * sg[a]; t2[b] -= h * sg[b]
                t3 = th.copy(); t3[a] -= h * sg[a]; t3[b] += h * sg[b]
                t4 = th.copy(); t4[a] -= h * sg[a]; t4[b] -= h * sg[b]
                v = (_cf_chi2(t1, m, Lr, W, r) - _cf_chi2(t2, m, Lr, W, r) - _cf_chi2(t3, m, Lr, W, r)
                     + _cf_chi2(t4, m, Lr, W, r)) / (4. * h * h)
            F[a, b] = 0.5 * v; F[b, a] = 0.5 * v          # Hessian of chi2 / 2 (the GN matrix in these units ~ 1)
    if not np.all(np.isfinite(F)):
        return False, th, chi2
    ev, V = np.linalg.eigh(F)
    if not ev[0] < -CF_NEG:
        return False, th, chi2
    best, cb = th, chi2
    for q in range(7):
        if not ev[q] < -CF_NEG:
            break
        for t in (0.5, -0.5, 1., -1., 2., -2., 4., -4., 8., -8., 16., -16.):
            tn = th + t * sg * V[:, q]
            tn[3] = _wrap(tn[3])
            c = _cf_chi2(tn, m, Lr, W, r)
            if np.isfinite(c) and c < cb:
                best, cb = tn, c
    return True, best, cb


@njit(cache=True)
def pair_cf(p1, c1, p2, c2, tp, K):
    """CF result vector (CF_NAMES) and the 3x3 vertex covariance; tp = TP result vector (start at the TP vertex,
    photon direction at the TP and the two measured curvatures; the helix turn of each track is the one of its TP arc).
    Minimiser: Levenberg-Marquardt (Gauss-Newton normal matrix H = J^T W J, numerical Jacobian, Marquardt damping
    H_kk (1 + lambda); a step is accepted only if it lowers the chi2), at most CF_MAXIT iterations; at the stationary
    point the full numerical Hessian of the chi2 is checked (_cf_escape): along a direction of negative curvature (a
    saddle point, possible for large residuals, invisible to H) the fit steps downhill and restarts (at most CF_NESCAPE
    times). ok = 1: converged = the predicted chi2 decrease of a full Gauss-Newton step g^T H^-1 g (g = J^T W r) is below
    CF_TOL and the Hessian has no negative curvature (a local minimum; several minima can exist for badly matching
    pairs); ok = 2: not converged (decrement above CF_TOL, or a saddle that could not be left); ok = 0: undefined (no
    tangent point: same curvature sign) or failed (singular matrices, non-finite numbers): every other entry and the
    covariance are NaN. niter = total LM iterations; chi2_1 / chi2_2 = contributions of track 1 / 2 (r_k^T C_k^-1 r_k
    of the track's 5 helix residuals; the track covariances are independent, so chi2 = chi2_1 + chi2_2 exactly);
    dec = final decrement."""
    res = np.full(N_CF, np.nan); Cx = np.full((3, 3), np.nan)
    res[0] = 0.
    if not tp[0] > 0.:
        return res, Cx
    W = np.zeros((10, 10))
    try:
        W[:5, :5] = np.linalg.inv(c1); W[5:, 5:] = np.linalg.inv(c2)
    except Exception:
        return res, Cx
    m = np.zeros(10); m[:5] = p1; m[5:] = p2
    th = np.array([tp[1], tp[2], tp[3], tp[9], tp[10], p1[2], p2[2]])
    Lr = np.array([tp[7], tp[8]])           # arcs of the two tracks to the TP: the helix turn of the model
    r = np.zeros(10); rp = np.zeros(10); rm = np.zeros(10)
    J = np.zeros((10, 7))
    _cf_res(th, m, Lr, r)
    chi2 = r @ W @ r
    if not np.isfinite(chi2):
        return res, Cx
    it = 0
    ntot = 0
    neg = False
    for _ in range(CF_NESCAPE + 1):
        th, chi2, it = _cf_lm(th, m, Lr, W, chi2)
        ntot += it
        # second-order check at the stationary point: the Gauss-Newton matrix cannot see negative curvature of the
        # chi2 (large residuals), so a saddle point would pass the decrement test; step along negative-curvature
        # directions and restart while that lowers the chi2
        neg, tn, cn = _cf_escape(th, m, Lr, W, chi2)
        if not (neg and cn < chi2):
            break
        th = tn; chi2 = cn
    it = ntot
    # decrement and covariance at the returned point
    _cf_res(th, m, Lr, r)
    chi2 = r @ W @ r
    _cf_jac(th, m, Lr, J, rp, rm)
    H = J.T @ W @ J
    g = J.T @ W @ r
    dec = _cf_decrement(H, g)
    try:
        Cth = np.linalg.inv(H)
    except Exception:
        return res, Cx
    if not (np.all(np.isfinite(Cth)) and np.all(np.isfinite(th)) and np.isfinite(chi2)):
        return res, Cx
    pt1, pt2 = K / max(abs(th[5]), 1e-12), K / max(abs(th[6]), 1e-12)
    pt = pt1 + pt2
    res[0] = 1. if (dec < CF_TOL and not neg) else 2.
    for k in range(7):
        res[1 + k] = th[k]
    res[8] = chi2; res[9] = it
    res[10] = pt * np.cos(th[3]); res[11] = pt * np.sin(th[3]); res[12] = pt * th[4]
    res[13] = r[:5] @ W[:5, :5] @ r[:5]; res[14] = r[5:] @ W[5:, 5:] @ r[5:]; res[15] = dec
    for a in range(3):
        for b in range(3):
            Cx[a, b] = Cth[a, b]
    return res, Cx


# ----------------------------------------------------------------------------------------------- pair loops
@njit(cache=True)
def _prefilter(p1, p2, pre, o):
    """cheap TP-based prefilter: pre = (max |gap|, max |dz|, max |dct|, min arc, max vertex radius)"""
    if not _tp_raw(p1, p2, o):
        return False
    if abs(o[4]) > pre[0] or abs(o[2] - o[3]) > pre[1] or abs(o[5]) > pre[2]:
        return False
    if o[6] < pre[3] or o[7] < pre[3]:
        return False
    if np.hypot(o[0], o[1]) > pre[4]:
        return False
    return True


@njit(cache=True)
def find_pairs_nb(par, q, ok, pre):
    n = par.shape[0]
    I = np.empty(n * n // 2 + 1, np.int64); Jx = np.empty(n * n // 2 + 1, np.int64)
    k = 0
    o = np.zeros(8)
    for i in range(n):
        if not ok[i]:
            continue
        for j in range(i + 1, n):
            if not ok[j] or q[i] * q[j] >= 0:
                continue
            if _prefilter(par[i], par[j], pre, o):
                I[k] = i; Jx[k] = j; k += 1
    return I[:k], Jx[:k]


Q_NAMES = (['tp_' + n for n in TP_NAMES] + ['cf_' + n for n in CF_NAMES] +
           ['tp_r', 'cf_r', 'tp_chi2', 'tp_dgam', 'tp_cosp', 'cf_dgam', 'cf_cosp', 'cf_sl', 'nhx1', 'nhx2', 'p1', 'p2'])
QI = {n: k for k, n in enumerate(Q_NAMES)}
N_Q = len(Q_NAMES)


@njit(cache=True)
def _quality(p1, c1, p2, c2, pv, rmax1, rmax2, K, out, Cx):
    ntp, ncf = N_TP, N_CF
    tp = pair_tp(p1, c1, p2, c2, K)
    cf, C = pair_cf(p1, c1, p2, c2, tp, K)
    for k in range(ntp):
        out[k] = tp[k]
    for k in range(ncf):
        out[ntp + k] = cf[k]
    for a in range(3):
        for b in range(3):
            Cx[a, b] = C[a, b]
    b0 = ntp + ncf
    out[b0] = np.hypot(tp[1], tp[2]) if tp[0] > 0 else np.nan
    out[b0 + 1] = np.hypot(cf[1], cf[2]) if cf[0] > 0 else np.nan
    # TP conversion chi2 of (gap, dz, dct) with their correlations (3 constraints)
    if tp[0] > 0:
        s = np.array([tp[17], tp[18], tp[19]])
        R = np.array([[1., tp[20], tp[21]], [tp[20], 1., tp[22]], [tp[21], tp[22], 1.]])
        v = np.array([tp[4], tp[5], tp[6]]) / np.maximum(s, 1e-12)
        try:
            out[b0 + 2] = v @ np.linalg.solve(R, v)
        except Exception:
            out[b0 + 2] = np.nan
    else:
        out[b0 + 2] = np.nan
    # pointing of the photon to the PV: signed transverse distance of the photon line from the PV, cos(3D angle between
    # PV -> vertex and the photon momentum)
    for which in range(2):
        if which == 0:
            ok_ = tp[0] > 0; x = tp[1]; y = tp[2]; z = tp[3]; px = tp[24]; py = tp[25]; pz = tp[26]
        else:
            ok_ = cf[0] > 0; x = cf[1]; y = cf[2]; z = cf[3]; px = cf[10]; py = cf[11]; pz = cf[12]
        if ok_:
            dx, dy, dz = x - pv[0], y - pv[1], z - pv[2]
            pt = max(np.hypot(px, py), 1e-12)
            dg = (dx * py - dy * px) / pt
            cp = (dx * px + dy * py + dz * pz) / max(np.sqrt(dx * dx + dy * dy + dz * dz) * np.sqrt(px * px + py * py + pz * pz), 1e-300)
        else:
            dg = np.nan; cp = np.nan
        out[b0 + 3 + 2 * which] = dg; out[b0 + 4 + 2 * which] = cp
    # CF position error along the photon momentum
    if cf[0] > 0:
        u = np.array([cf[10], cf[11], cf[12]]); u = u / max(np.sqrt(u @ u), 1e-300)
        out[b0 + 7] = np.sqrt(max(u @ C @ u, 0.))
    else:
        out[b0 + 7] = np.nan
    # inner hits impossible for the vertex radius (hit-radius bound of jv.tracks.max_radius_from_hits), TP radius
    if tp[0] > 0:
        r = out[b0]
        out[b0 + 8] = 1. if r > rmax1 else 0.
        out[b0 + 9] = 1. if r > rmax2 else 0.
    else:                                   # no tangent point: undefined
        out[b0 + 8] = np.nan
        out[b0 + 9] = np.nan
    P = np.zeros(3)
    _mom(p1, tp[7], K, P); out[b0 + 10] = np.sqrt(P @ P)
    _mom(p2, tp[8], K, P); out[b0 + 11] = np.sqrt(P @ P)


@njit(cache=True)
def pair_quality_nb(par, cov, I, J, pv_of, pv, rmax, K):
    n = len(I)
    Q = np.full((n, N_Q), np.nan)
    CX = np.zeros((n, 3, 3))
    for k in range(n):
        i, j = I[k], J[k]
        _quality(par[i], cov[i], par[j], cov[j], pv[pv_of[k]], rmax[i], rmax[j], K, Q[k], CX[k])
    return Q, CX


# ----------------------------------------------------------------------------------------------- selections
# prefilter (cheap, TP without errors); values in cm, GeV
PRE = np.array([3.0, 15.0, 0.2, -5.0, 200.])   # |gap| < 3 cm, |dz| < 15 cm, |dct| < 0.2, arcs > -5 cm, R < 200 cm


def find_pairs(par, cov, q, ok=None, pre=PRE):
    """opposite-charge pairs (i < j) of usable tracks passing the prefilter"""
    if ok is None:
        ok = np.isfinite(par).all(1) & np.isfinite(cov.reshape(len(par), -1)).all(1)
    return find_pairs_nb(np.ascontiguousarray(par), np.asarray(q, np.int64), np.asarray(ok, np.bool_), np.asarray(pre, float))


def pair_quality(par, cov, I, J, pv, rmax, k_pt=K_PT):
    """quality matrix (Q_NAMES) and CF vertex covariances of the pairs (I, J) of one event, PV pv (3,), hit-radius
    bounds rmax per track, momentum constant k_pt (pT = k_pt / |C|)"""
    pv = np.asarray(pv, float).reshape(1, 3)
    return pair_quality_nb(np.ascontiguousarray(par), np.ascontiguousarray(cov), np.asarray(I, np.int64),
                           np.asarray(J, np.int64), np.zeros(len(I), np.int64), pv, np.asarray(rmax, float), float(k_pt))


def select(Q, sel):
    """boolean mask of the candidate rows of Q passing the selection dict sel (keys: see LOOSE / TIGHT)"""
    m = (Q[:, QI['tp_ok']] > 0) & (Q[:, QI['cf_ok']] > 0)
    for k, v in sel.items():
        if v is None:
            continue
        lo, hi = v
        col = Q[:, QI[k]]
        if k in ('tp_dgam', 'cf_dgam', 'tp_gap', 'tp_dz', 'tp_dct'):   # symmetric variables: cut on |value|
            col = np.abs(col)
        if lo is not None:
            m &= col >= lo
        if hi is not None:
            m &= col < hi
    return m


def find_conversions(par, cov, q, nvdet, nitc, pv, sel='loose', pre=PRE, k_pt=K_PT):
    """conversion candidates of one event: dict with i, j (track indices), x (CF vertex), cov (3x3), Q (quality
    matrix, Q_NAMES); sel = 'loose' (supply to the joint fit), 'tight' (stand-alone list, arbitrated: each track in at
    most one candidate), a selection dict, or None (every prefiltered pair); k_pt = momentum constant (pT = k_pt / |C|)"""
    I, J = find_pairs(par, cov, q, pre=pre)
    rmax = TR.max_radius_from_hits(nvdet, nitc)
    Q, CX = pair_quality(par, cov, I, J, pv, rmax, k_pt)
    if sel is None:
        m = np.ones(len(I), bool)
    else:
        d = LOOSE if sel == 'loose' else (TIGHT if sel == 'tight' else sel)
        m = select(Q, d)
        if sel == 'tight':
            k = np.nonzero(m)[0]
            m[k] = arbitrate(np.zeros(len(k)), I[k], J[k], Q[k, QI['cf_chi2']])
    x = Q[:, [QI['cf_x'], QI['cf_y'], QI['cf_z']]]
    return dict(i=I[m], j=J[m], x=x[m], cov=CX[m], Q=Q[m])


def arbitrate(ev, I, J, score):
    """greedy one-candidate-per-track arbitration: candidates in increasing score (e.g. CF chi2) are kept unless one of
    their tracks is already used in the same event; returns the keep mask"""
    ev, I, J = np.asarray(ev, np.int64), np.asarray(I, np.int64), np.asarray(J, np.int64)
    keep = np.zeros(len(I), bool)
    used = set()
    for k in np.argsort(np.asarray(score, float), kind='stable'):
        a, b = (ev[k], I[k]), (ev[k], J[k])
        if a in used or b in used:
            continue
        keep[k] = True
        used.add(a); used.add(b)
    return keep


# the two selections (cuts on the quality matrix; tuple = (lower, upper), None = open)
LOOSE = dict(cf_chi2=(None, 200.), tp_mee=(None, 0.05), cf_cosp=(0.9, None))
TIGHT = dict(cf_chi2=(None, 15.), tp_mee=(None, 0.02), cf_cosp=(0.995, None), cf_r=(2., None))   # + arbitrate()
