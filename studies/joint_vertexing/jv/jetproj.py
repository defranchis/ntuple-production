"""Jet-projected primary vertex (QFNDIP-like) on a straight-line track model.

Each track is a straight line at its origin-referenced perigee with the full (d0, z0) covariance. Its offset from the
vertex is split along E (perpendicular to the track and its jet axis) and U (in the track-jet plane); the jet-angle error
enters the E variance. Two forms, both used by jv.fit.jetproj_pv:
  hard  qfndip: the ALEPH QFNDIP rule (findip.F): 4 sigma windows, component-wise outlier removal at dchi2 >= 9
  soft  soft_corr: one-sided lifetime weight on U and an outlier weight on E, iterated to convergence,
        with the correlated (E, U) covariance (eu_cov)
Units: cm (lengths), GeV.
"""
import numpy as np

MINERR2 = 6.25e-6                              # QFNDIP minimum component error (25 um)^2, trot.F:32
NSIG, DCHI2 = 4.0, 9.0                         # findip.F: 4 sigma windows, component-wise outlier dchi2 >= 9


def jet_sig2_aleph(p):
    """ALEPH JETERR (jeterr.F:38-39, 68): jet angular error^2 [rad^2], p capped at 45 GeV"""
    p = min(p, 45.0)
    return 4.0 * (0.2862 - 0.01205 * p + 0.000135 * p * p) ** 2


def line_params(par_raw, cov_raw):
    """raw ALEPH-convention track state -> line-model inputs in the LCIO convention of pfcand_*:
    d0 = -D0, z0, phi0, ct = tanLambda, cov(d0,d0), cov(z0,z0), cov(d0,z0) = -cov(D0,Z0)"""
    d0 = -par_raw[:, 0].astype(np.float64)
    z0 = par_raw[:, 3].astype(np.float64)
    phi = par_raw[:, 1].astype(np.float64)
    ct = par_raw[:, 4].astype(np.float64)
    c00 = cov_raw[:, 0].astype(np.float64)
    c11 = cov_raw[:, 9].astype(np.float64)
    c01 = -cov_raw[:, 6].astype(np.float64)
    return np.stack([d0, z0, phi, ct, c00, c11, c01], 1)


def prep_qfn(trk, axes, kj_force=None, minerr2=MINERR2, t3sign=1.0):
    """QFNDIP track preparation: per track the directions H0 (theta), H1 (transverse), E (perpendicular to track and jet)
    and U (in the track-jet plane), their error^2 and offsets, the track-jet angle nu and the jet index. Each track goes
    with the jet of largest t.j; kj_force (per track, -1 = not forced) overrides that association. minerr2 = component
    error floor^2 (ALEPH 25 um, trot.F:32). t3sign: sign of the (H1, H0) offset covariance T_3 = t3sign * cov(d0, z0) cos(l);
    +1 (default) is consistent with T0 = (-d0 sin phi, d0 cos phi, z0), whose offsets along H1 / H0 are +d0 / +z0 cos(l)."""
    d0, z0, phi, ct, c00, c11, c01 = trk
    cosl = 1.0 / np.sqrt(1.0 + ct * ct); sinl = ct * cosl
    T1 = np.stack([cosl * np.cos(phi), cosl * np.sin(phi), sinl], 1)
    T0 = np.stack([-d0 * np.sin(phi), d0 * np.cos(phi), z0], 1)
    T_1 = c00; T_2 = c11 * cosl ** 2; T_3 = t3sign * c01 * cosl
    N = len(d0)
    H = np.zeros((N, 4, 3)); E2 = np.zeros((N, 5))
    norm = np.sqrt(T1[:, 0] ** 2 + T1[:, 1] ** 2)
    H[:, 1] = np.stack([-T1[:, 1] / norm, T1[:, 0] / norm, np.zeros(N)], 1)
    H[:, 0] = np.stack([-T1[:, 0] * T1[:, 2] / norm, -T1[:, 1] * T1[:, 2] / norm, norm], 1)
    E2[:, 1] = np.maximum(T_1, minerr2); E2[:, 0] = np.maximum(T_2, minerr2)
    dots = T1 @ np.array(axes).T
    kj = np.argmax(dots, 1)
    if kj_force is not None:
        # forced (jet-constituent) association, except for a track pointing away from its own jet (t.j <= 0)
        kf_ok = (kj_force >= 0) & (dots[np.arange(N), np.maximum(kj_force, 0)] > 0)
        kj = np.where(kf_ok, kj_force, kj)
    md = dots[np.arange(N), kj]
    nu = np.arccos(np.minimum(md, 1.0))
    for i in range(N):
        J1 = axes[kj[i]]
        if nu[i] > 0:
            s = np.sin(nu[i]); t = T1[i]
            E = np.array([t[1] * J1[2] - t[2] * J1[1], t[2] * J1[0] - t[0] * J1[2], t[0] * J1[1] - t[1] * J1[0]]) / s
            U = np.array([-t[1] * E[2] + t[2] * E[1], -t[2] * E[0] + t[0] * E[2], -t[0] * E[1] + t[1] * E[0]])
            H[i, 2] = E; H[i, 3] = U
            cose = E @ H[i, 1]; sine = E @ H[i, 0]
            E2[i, 2] = max(T_1[i] * cose ** 2 + T_2[i] * sine ** 2 + 2 * T_3[i] * cose * sine, minerr2)
            E2[i, 3] = max(T_2[i] * cose ** 2 + T_1[i] * sine ** 2 - 2 * T_3[i] * cose * sine, minerr2)
        else:
            H[i, 2] = H[i, 0]; H[i, 3] = H[i, 1]; E2[i, 2] = E2[i, 0]; E2[i, 3] = E2[i, 1]
    E2[:, 4] = E2[:, 2]
    DT0 = np.einsum('kta,ka->kt', H, T0)
    return H, E2, DT0, nu, kj


def qfndip(Q, axes_sig2, bsW, bsc, window, allfull=False):
    """ALEPH findip.F logic (QFNDIP jet-projected PV); returns x, covariance (inverse information of the last solve), tracks used"""
    H, E2, DT0, nu, kj = Q
    E2 = E2.copy(); N = len(nu)
    flag = np.full(N, 4 if not allfull else 3)
    old = None
    bsb = bsW @ bsc
    for meta in range(6):
        ii, tt = np.nonzero(((flag[:, None] >> np.arange(3)[None, :]) & 1).astype(bool))
        while True:
            hs = H[ii, tt]; e2 = E2[ii, tt]; dt = DT0[ii, tt]
            A = bsW + (hs.T / e2) @ hs; b = bsb + hs.T @ (dt / e2)
            Ai = np.linalg.inv(A); x = Ai @ b
            if len(ii) == 0:
                break
            r = dt - hs @ x; den = e2 - np.einsum('ka,ab,kb->k', hs, Ai, hs)
            d = np.where(den > 0, r * r / np.where(den > 0, den, 1.0), 0.0)
            kb = int(np.argmax(d))
            if d[kb] >= DCHI2:
                ii = np.delete(ii, kb); tt = np.delete(tt, kb)
            else:
                break
        nused = len(np.unique(ii))
        if old is not None:
            dx = x - old; d2 = dx @ dx
            if d2 == 0 or not (d2 > (dx @ Ai @ dx) / d2):
                return x, Ai, nused
            if allfull:
                return x, Ai, nused
        if meta == 5:
            return x, Ai, nused
        old = x.copy()
        L = DT0 - np.einsum('kta,a->kt', H, x)
        S = np.sqrt(E2[:, :4])
        if allfull:
            flag = np.where((np.abs(L[:, 0]) < NSIG * S[:, 0]) & (np.abs(L[:, 1]) < NSIG * S[:, 1]), 3, 0)
            continue
        full = (L[:, 3] < window) & (np.abs(L[:, 0]) < NSIG * S[:, 0]) & (np.abs(L[:, 1]) < NSIG * S[:, 1])
        eonly = (~full) & (np.abs(L[:, 2]) < NSIG * np.sqrt(E2[:, 2]))
        flag = np.where(full, 3, np.where(eonly, 4, 0))
        s2 = np.array([axes_sig2[k] for k in kj])
        nus = np.where(nu > 0, nu, 1.0)
        jerr2 = np.where(nu > 0, np.minimum(L[:, 3] ** 2, s2 * (L[:, 3] / nus) ** 2), L[:, 3] ** 2)
        E2[:, 2] = np.where(eonly, E2[:, 4] + jerr2, E2[:, 2])
    return x, Ai, nused


def phi_(z):
    return np.exp(-0.5 * z * z)


TOL_CM = 1e-6       # convergence of the soft fits: position step below 0.01 um
MAXITER = 500


def _iters(niter, maxiter):
    return niter if niter is not None else maxiter


def _done(niter, xn, x, tol):
    return niter is None and np.max(np.abs(xn - x)) < tol


def jet_axes(theta, phi):
    return [np.array([np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)]) for t, p in zip(theta, phi)]


def eu_cov(trk, Q, minerr2=MINERR2, t3sign=1.0):
    """2x2 covariance of the (E, U) track offsets, rotated from the (phi, theta)-direction offsets with their correlation
    (prep_qfn keeps only the two diagonal terms). E = a1 H1 + a0 H0, U = b1 H1 + b0 H0 (H1 = transverse direction
    perpendicular to the track: variance T_1 = var(d0); H0 = theta direction: T_2 = var(z0) cos^2(l); covariance T_3)."""
    d0, z0, phi, ct, c00, c11, c01 = trk
    H, E2, DT0, nu, kj = Q
    cosl = 1.0 / np.sqrt(1.0 + ct * ct)
    T_1 = c00; T_2 = c11 * cosl ** 2; T_3 = t3sign * c01 * cosl      # t3sign: see prep_qfn
    a1 = np.einsum('ka,ka->k', H[:, 2], H[:, 1]); a0 = np.einsum('ka,ka->k', H[:, 2], H[:, 0])
    b1 = np.einsum('ka,ka->k', H[:, 3], H[:, 1]); b0 = np.einsum('ka,ka->k', H[:, 3], H[:, 0])
    C = np.zeros((len(d0), 2, 2))
    C[:, 0, 0] = np.maximum(a1 * a1 * T_1 + a0 * a0 * T_2 + 2 * a1 * a0 * T_3, minerr2)
    C[:, 1, 1] = np.maximum(b1 * b1 * T_1 + b0 * b0 * T_2 + 2 * b1 * b0 * T_3, minerr2)
    C[:, 0, 1] = C[:, 1, 0] = a1 * b1 * T_1 + a0 * b0 * T_2 + (a1 * b0 + a0 * b1) * T_3
    return C


def soft_corr(Q, Ceu, axes_sig2, bsW, bsc, x0, kf=1.5, kb=3.0, ke=3.0, niter=None, tol=TOL_CM, maxiter=MAXITER, info=None):
    """soft form of the jet-projected PV with the full (E, U) covariance: per track information H2^T W^1/2 C'^-1 W^1/2 H2,
    H2 = (E, U), W = diag(wE, wu), C' = Ceu with the jet-angle term added to the E variance. wu = one-sided lifetime weight
    on the in-plane offset u (phi(z) / (phi(z) + phi(k)), k = kf in front of the vertex, kb behind), wE = outlier weight on
    the E offset (ke). With all weights 1 and no jet term it is the full straight-line track information.
    niter = None: iterate until the position step is < tol (cm), at most maxiter iterations; otherwise niter iterations.
    info (dict, optional) receives 'niter' and 'converged'. Returns x, covariance (inverse weighted information at the last
    iteration), number of tracks with E weight > 0.5."""
    H, E2, DT0, nu, kj = Q
    x = x0.copy(); N = len(nu)
    s2j = np.array([axes_sig2[k] for k in kj])
    H2 = H[:, 2:4]                      # (N, 2, 3)
    D2 = DT0[:, 2:4]
    conv = False
    for it in range(_iters(niter, maxiter)):
        L_ = D2 - np.einsum('kta,a->kt', H2, x)
        l3, l4 = L_[:, 0], L_[:, 1]
        jet2 = np.where(nu > 0, np.minimum(l4 * l4, s2j * (l4 / np.where(nu > 0, nu, 1.0)) ** 2), l4 * l4)
        Cp = Ceu.copy(); Cp[:, 0, 0] += jet2
        sE = np.sqrt(Cp[:, 0, 0]); su = np.sqrt(Cp[:, 1, 1])
        z = l4 / su
        wu = np.where(l4 <= 0, phi_(z) / (phi_(z) + phi_(kb)), phi_(z) / (phi_(z) + phi_(kf)))
        wE = phi_(l3 / sE) / (phi_(l3 / sE) + phi_(ke))
        sw = np.stack([np.sqrt(wE), np.sqrt(wu)], 1)
        Wt = np.linalg.inv(Cp) * sw[:, :, None] * sw[:, None, :]
        A = bsW + np.einsum('kia,kij,kjb->ab', H2, Wt, H2)
        b = bsW @ bsc + np.einsum('kia,kij,kj->a', H2, Wt, D2)
        xn = np.linalg.solve(A, b)
        conv = _done(niter, xn, x, tol)
        x = xn
        if conv:
            break
    if info is not None:
        info['niter'] = it + 1; info['converged'] = conv
    return x, np.linalg.inv(A), int((wE > 0.5).sum())
