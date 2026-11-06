"""Joint-fit engine helpers: track-vertex compatibility with the vertex covariance and leave-one-out.

chi2_vertex(par, cov, x, H, w, mode, eps, fdir): per-track compatibility with the vertex at x with information matrix H (inverse covariance)
  mode 0 'point': chi2 = r^T D r, track covariance only (baseline)
  mode 1 'vcov' : chi2 = r^T (D - D (H + D)^-1 D) r, residual covariance = track + vertex (the track's own contribution stays in H)
  mode 2 'loo'  : leave-one-out: the track's contribution w D is removed from H and from the position first (H_i = H - w D + eps 1,
                  x_i = x + H_i^-1 w D (x - X)), then the same formula with H_i, r_i = x_i - X. This equals the chi2 increase of adding
                  the track to the vertex fitted without it.
  dens = ln det(H_i + D) - ln det(H_i) >= 0: minus twice the log of the volume factor of the marginal likelihood of the track given the
  vertex (matrix determinant lemma); eps is a weak prior on the vertex position (information, cm^-2) that makes H_i invertible for a
  vertex left with a single track. Also returned: cosf = cosine between the track direction at its closest approach and fdir, L = arc.
Units cm.
"""
import numpy as np
from numba import njit

from . import tracks as TR


@njit(cache=True)
def chi2_vertex(par, cov, x, H, w, mode, eps, fdir):
    n = par.shape[0]
    chi2 = np.zeros(n); dens = np.zeros(n); cosf = np.zeros(n); Lo = np.zeros(n)
    X = np.zeros(3); D = np.zeros((3, 3))
    for t in range(n):
        c2, L = TR.lin1(par[t], cov[t], x, 0., X, D)
        Lo[t] = L
        ang = par[t][1] + 2 * par[t][2] * L
        ct = par[t][4]
        nrm = np.sqrt(1. + ct * ct)
        cosf[t] = (np.cos(ang) * fdir[0] + np.sin(ang) * fdir[1] + ct * fdir[2]) / nrm
        if mode == 0:
            chi2[t] = c2
            continue
        Hi = H + eps * np.eye(3)
        xi = x
        if mode == 2 and w[t] > 0:
            Hi = Hi - w[t] * D
            Ci = TR.reginv3(Hi)
            xi = x + Ci @ (w[t] * (D @ (x - X)))
        r = xi - X
        M = Hi + D
        Mi = TR.reginv3(M)
        Dr = D @ r
        chi2[t] = max(r @ Dr - Dr @ (Mi @ Dr), 0.)
        dH = np.linalg.det(Hi)
        dM = np.linalg.det(M)
        dens[t] = np.log(dM / dH) if (dH > 0 and dM > 0) else 0.
    return chi2, dens, cosf, Lo
