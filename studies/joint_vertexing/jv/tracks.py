"""Track model of the joint vertexer (numba).

Helix: the model of the stage1 PV fit (analyzer_pvnew.h), re-implemented in Python:
parameters p = (D, phi0, C, z0, ct) in the vertex-fitter convention (D = -D0_raw, C = omega_raw / 2, ct = tanLambda;
conversion `to_fit` from the raw ALEPH state stored in the EDM4hep input), arc parameter L = transverse arc length from the perigee.
Linearised track model: at a point x the track is replaced by its point of closest approach X (weighted metric) and
the rank-2 weight matrix D = W - W a a^T W / (a^T W a), W = (A Cov A^T)^-1; chi2(track, x') = (x' - X)^T D (x' - X).
Units: cm, GeV; B = 1.5 T (pT = PT_CONST / |omega_raw| = PT_CONST / 2 / |C|).
Momentum constant: PT_CONSTS['legacy'] = 0.0045 (the frozen final, used by the functions of this module) or
PT_CONSTS['stage1'] = kPtPerTeslaCm x kBz = 0.0029979 x 1.5 = 0.00449685 (stage1, src/aleph_units.h), selected by the
config key pt_const of jv.hybrid: the stored momenta and the momenta of hybrid._mass / _kin_ok are rescaled by
pt_scale(name); jv.v0 and jv.conv take PT_CONSTS[pt_const] / 2 as an argument.
"""
import numpy as np
from numba import njit

KTRI = np.array([[0, 1, 3, 6, 10], [1, 2, 4, 7, 11], [3, 4, 5, 8, 12], [6, 7, 8, 9, 13], [10, 11, 12, 13, 14]])
KSC = np.array([1., 1., -0.5, 1., 1.])
FLIP = np.array([1, 4, 6, 8, 10, 12])
PT_CONST_LEGACY = 0.0045                    # GeV cm: pT = PT_CONST / |omega_raw|, frozen final (c = 0.3, B = 1.5 T)
PT_PER_TESLA_CM = 0.0029979                 # GeV / (T cm): stage1 src/aleph_units.h kPtPerTeslaCm
BZ_TESLA = 1.5                              # T: stage1 src/aleph_units.h kBz
PT_CONSTS = {'legacy': PT_CONST_LEGACY, 'stage1': PT_PER_TESLA_CM * BZ_TESLA}
PT_CONST = PT_CONST_LEGACY                  # the constant of the functions below (momentum_at); 0.5 * PT_CONST per |C|
_PT_HALF = 0.5 * PT_CONST_LEGACY            # = 0.00225 exactly (halving is exact in floating point)
M_PI = 0.13957039


def to_fit(par_raw, cov15):
    """raw ALEPH state (D0, phi, omega, Z0, tanLambda) + 15-element lower-triangle covariance -> fitter convention
    (n, 5), (n, 5, 5)"""
    p = np.asarray(par_raw, np.float64).copy()
    p[:, 0] *= -1.
    p[:, 2] *= -1.
    p *= KSC[None, :]
    c = np.asarray(cov15, np.float64).copy()
    if c.shape[1] < 15:
        raise ValueError('cov needs 15 elements')
    c = c[:, :15].copy()
    c[:, FLIP] *= -1.
    cc = c[:, KTRI] * KSC[None, :, None] * KSC[None, None, :]
    return p, cc


@njit(cache=True)
def _sinc(u):
    if abs(u) < 1e-4:
        return 1. - u * u / 6. + u ** 4 / 120.
    return np.sin(u) / u


@njit(cache=True)
def _dsinc(u):
    if abs(u) < 1e-4:
        return -u / 3. + u ** 3 / 30.
    return (u * np.cos(u) - np.sin(u)) / (u * u)


@njit(cache=True)
def helix1(p, L, X, a, A):
    D, p0, C, z0, ct = p[0], p[1], p[2], p[3], p[4]
    u = C * L
    f = L * _sinc(u)
    ang = p0 + u
    ca, sa, cp, sp = np.cos(ang), np.sin(ang), np.cos(p0), np.sin(p0)
    X[0] = -D * sp + f * ca
    X[1] = D * cp + f * sa
    X[2] = z0 + ct * L
    ang2 = p0 + 2 * C * L
    a[0] = np.cos(ang2); a[1] = np.sin(ang2); a[2] = ct
    dfdC = L * L * _dsinc(u)
    for i in range(3):
        for j in range(5):
            A[i, j] = 0.
    A[0, 0] = -sp; A[1, 0] = cp
    A[0, 1] = -D * cp - f * sa; A[1, 1] = -D * sp + f * ca
    A[0, 2] = dfdC * ca - f * L * sa; A[1, 2] = dfdC * sa + f * L * ca
    A[2, 3] = 1.; A[2, 4] = L


@njit(cache=True)
def reginv3(M):
    """inverse of a symmetric 3x3 (analytic cofactors); eigenvalue-floored inverse when ill-conditioned"""
    a, b, c = M[0, 0], 0.5 * (M[0, 1] + M[1, 0]), 0.5 * (M[0, 2] + M[2, 0])
    d, e, f = M[1, 1], 0.5 * (M[1, 2] + M[2, 1]), M[2, 2]
    A = d * f - e * e; B = c * e - b * f; Cc = b * e - c * d
    det = a * A + b * B + c * Cc
    tr = abs(a) + abs(d) + abs(f)
    out = np.zeros((3, 3))
    if det > 1e-9 * tr * tr * tr and det > 0:
        out[0, 0] = A / det; out[0, 1] = B / det; out[0, 2] = Cc / det
        out[1, 1] = (a * f - c * c) / det; out[1, 2] = (b * c - a * e) / det; out[2, 2] = (a * d - b * b) / det
        out[1, 0] = out[0, 1]; out[2, 0] = out[0, 2]; out[2, 1] = out[1, 2]
        return out
    S = 0.5 * (M + M.T)
    w, U = np.linalg.eigh(S)
    wmax = np.max(np.abs(w))
    fl = 1e-12 * wmax
    for k in range(3):
        wc = w[k] if w[k] > fl else fl
        if wc <= 0:
            wc = 1e-300
        for i in range(3):
            for j in range(3):
                out[i, j] += U[i, k] * U[j, k] / wc
    return out


@njit(cache=True)
def lin1(p, cov, x, L0, X, D):
    """linearise one track at point x (start arc L0): fills X (PCA), D (3x3); returns (chi2, L)"""
    a = np.zeros(3); A = np.zeros((3, 5)); L = L0
    W = np.zeros((3, 3)); aw = np.zeros(3)
    for it in range(12):
        helix1(p, L, X, a, A)
        W = reginv3(A @ cov @ A.T)
        aw = W @ a
        den = a @ aw
        if not (den > 0) or not np.isfinite(den):
            break
        dL = ((x - X) @ aw) / den
        if abs(dL) > 50.:
            dL = 50. * np.sign(dL)
        L = L + dL
        if abs(dL) < 1e-10:
            break
    helix1(p, L, X, a, A)
    W = reginv3(A @ cov @ A.T)
    aw = W @ a
    den = a @ aw
    for i in range(3):
        for j in range(3):
            D[i, j] = W[i, j] - (aw[i] * aw[j] / den if (den > 0 and np.isfinite(den)) else 0.)
    for i in range(3):
        for j in range(i + 1, 3):
            s = 0.5 * (D[i, j] + D[j, i]); D[i, j] = s; D[j, i] = s
    r = x - X
    return r @ D @ r, L


@njit(cache=True)
def lin_all(par, cov, idx, x, Ls, Xo, Do, chi2o):
    """linearise tracks idx at the point x; Ls (per track, all tracks) holds the start arc and is updated"""
    for k in range(len(idx)):
        t = idx[k]
        c2, L = lin1(par[t], cov[t], x, Ls[t], Xo[k], Do[k])
        chi2o[k] = c2
        Ls[t] = L


@njit(cache=True)
def chi2_to_point(par, cov, x):
    """chi2 of every track to the point x (track covariance only) and the arc L at the PCA"""
    n = par.shape[0]
    out = np.zeros(n); Lo = np.zeros(n)
    X = np.zeros(3); D = np.zeros((3, 3))
    for t in range(n):
        out[t], Lo[t] = lin1(par[t], cov[t], x, 0., X, D)
    return out, Lo


@njit(cache=True)
def momentum_at(p, L):
    """momentum vector (GeV) of the track at arc L"""
    pt = _PT_HALF / max(abs(p[2]), 1e-12)
    ang = p[1] + 2 * p[2] * L
    return np.array([pt * np.cos(ang), pt * np.sin(ang), pt * p[4]])


@njit(cache=True)
def point_at(p, L):
    X = np.zeros(3); a = np.zeros(3); A = np.zeros((3, 5))
    helix1(p, L, X, a, A)
    return X


@njit(cache=True)
def ip_info(par, cov, pv, pvcov, axis_dirs, jet_of):
    """per track w.r.t. the vertex pv (cov pvcov) and a lifetime-sign axis (axis_dirs[jet_of[t]], or the track's own
    nearest axis when jet_of < 0): returns ip3 (cm), sig3 (signed by lifetime sign), ip2 (transverse, signed by lifetime
    sign), sig2, lifetime sign (+1/-1)"""
    n = par.shape[0]
    ip3 = np.zeros(n); s3 = np.zeros(n); ip2 = np.zeros(n); s2 = np.zeros(n); sg = np.zeros(n)
    X = np.zeros(3); D = np.zeros((3, 3)); a = np.zeros(3); A = np.zeros((3, 5))
    for t in range(n):
        c2, L = lin1(par[t], cov[t], pv, 0., X, D)
        helix1(par[t], L, X, a, A)
        r = X - pv
        # track position covariance at the PCA plus the vertex covariance, projected on r
        Ct = A @ cov[t] @ A.T
        d3 = np.sqrt(r @ r)
        if d3 > 0:
            u = r / d3
            var3 = u @ (Ct + pvcov) @ u
        else:
            var3 = 1.
        mom = momentum_at(par[t], L)
        j = jet_of[t]
        if j < 0:
            best = -2.
            for k in range(axis_dirs.shape[0]):
                d = axis_dirs[k] @ mom
                if d > best:
                    best = d; j = k
        # lifetime sign: closest approach of the track line and the axis, in front of or behind the vertex. r comes from
        # the chi2-metric PCA and has a component along the track, which is removed first (with it: b tracks 62 % positive,
        # without: 79 %, primaries unchanged; MC truth classes, w.r.t. the generated PV, first 4000 events of mc_ZM4212.41)
        re_ = r - ((r @ a) / (a @ a)) * a
        sgn = 1. if (re_ @ axis_dirs[j]) >= 0 else -1.
        sg[t] = sgn
        ip3[t] = d3
        s3[t] = sgn * d3 / np.sqrt(max(var3, 1e-14))
        # transverse: PCA in xy (approx. same L), signed distance
        r2 = np.array([r[0], r[1], 0.])
        d2 = np.sqrt(r2 @ r2)
        if d2 > 0:
            u2 = r2 / d2
            var2 = u2 @ (Ct + pvcov) @ u2
        else:
            var2 = 1.
        ip2[t] = sgn * d2
        s2[t] = sgn * d2 / np.sqrt(max(var2, 1e-14))
    return ip3, s3, ip2, s2, sg


def pt_scale(name):
    """momenta of the selected constant (config pt_const: 'legacy' | 'stage1') over those of PT_CONST (all momenta
    scale with the constant: pT and pz = pT tanLambda)"""
    return PT_CONSTS[name] / PT_CONST


def max_radius_from_hits(nvdet, nitc):
    """upper bound [cm] on the transverse radius of the track's production vertex from its hits (layer radius +
    tolerance): 6.5 + 1.1 with >= 2 VDET hits, 11.3 + 1.8 with 1 VDET hit, at most 26 + 2.9 with ITC hits; inf otherwise"""
    nvdet = np.asarray(nvdet); nitc = np.asarray(nitc)
    r = np.full(len(nvdet), np.inf)
    r[nvdet >= 2] = 6.5 + 1.1
    r[nvdet == 1] = 11.3 + 1.8
    r = np.where(nitc >= 1, np.minimum(r, 26. + 2.9), r)
    return r


def circle_crossings(p1, p2):
    """xy crossing points of the two track circles (fitter convention); 0, 1 or 2 points
    (disjoint / nested circles: the midpoint of the gap)"""
    xc, yc, rc = np.zeros(2), np.zeros(2), np.zeros(2)
    for k, p in enumerate((p1, p2)):
        D, phi0, C = p[0], p[1], p[2]
        if C == 0:
            return []
        xc[k] = -(D + 0.5 / C) * np.sin(phi0)
        yc[k] = (D + 0.5 / C) * np.cos(phi0)
        rc[k] = 0.5 / abs(C)
    dx, dy = xc[1] - xc[0], yc[1] - yc[0]
    d = np.hypot(dx, dy)
    if not np.isfinite(d) or not d > 0:
        return []
    ux, uy = dx / d, dy / d
    if d >= rc[0] + rc[1]:
        t = rc[0] + 0.5 * (d - rc[0] - rc[1])
        return [(xc[0] + t * ux, yc[0] + t * uy)]
    if d <= abs(rc[0] - rc[1]):
        s = 1. if rc[0] > rc[1] else -1.
        return [(0.5 * (xc[0] + xc[1] + s * (rc[0] + rc[1]) * ux), 0.5 * (yc[0] + yc[1] + s * (rc[0] + rc[1]) * uy))]
    a = (rc[0] ** 2 - rc[1] ** 2 + d * d) / (2. * d)
    h2 = rc[0] ** 2 - a * a
    X, Y = xc[0] + a * ux, yc[0] + a * uy
    if h2 > 0:
        h = np.sqrt(h2)
        return [(X - h * uy, Y + h * ux), (X + h * uy, Y - h * ux)]
    return [(X, Y)]


def arc_to_xy(p, x, y):
    """transverse arc length L of the track (fitter convention) at the circle point nearest to (x, y)"""
    D, phi0, C = p[0], p[1], p[2]
    x0, y0 = -D * np.sin(phi0), D * np.cos(phi0)
    xc, yc = -(D + 0.5 / C) * np.sin(phi0), (D + 0.5 / C) * np.cos(phi0)
    a0 = np.arctan2(y0 - yc, x0 - xc)
    a1 = np.arctan2(y - yc, x - xc)
    dphi = (a1 - a0 + np.pi) % (2 * np.pi) - np.pi
    return dphi / (2 * C)
