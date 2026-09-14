"""Phenotype shapes, discretisation maps, and analytic estimands.

Generative model (latent, continuous):

    Z ~ F0                      base shape: normal / lognormal / bimodal
    Y*| g=k  =  k*a + (1 + k*s) * Z
    Y        =  D(Y*)           deterministic, monotone non-decreasing map

g in {0,1,2} with HWE frequencies at minor allele frequency p.

Three estimands are computed by quadrature, all for the RIF slope
    b_tau = -Cov(g, F_{Y|g}(q_tau)) / ( f(q_tau) * Var(g) )

  b_latent : uses the latent Y* (its own quantile, its own true density).
             This is what a user believes they are estimating when they round,
             censor or bin a continuous trait.
  b_obs    : uses the observed Y and its TRUE density.  Only defined when Y has
             no atom at q_tau.
  b_ship   : uses the observed Y and the density the SHIPPED estimator actually
             converges to, namely the population Gaussian-KDE convolution
             (f * K_h)(q_tau) at Silverman bandwidth h = 1.06 sd(Y) N^{-1/5}.
             This is the asymptotic target of Stage 1 at that N.
"""
from __future__ import annotations

import numpy as np
from scipy import stats as sps

# ----------------------------------------------------------------- base shapes

class Base:
    """Latent base distribution Z, standardised to mean 0 / sd 1."""

    def __init__(self, name: str):
        self.name = name
        if name == "normal":
            self._comp_w = np.array([1.0])
            self._comp_m = np.array([0.0])
            self._comp_s = np.array([1.0])
            self._kind = "mix"
        elif name == "bimodal":
            # two well separated modes -> deep trough at the median
            w = np.array([0.5, 0.5]); m = np.array([-1.5, 1.5]); s = np.array([0.5, 0.5])
            mu = float(w @ m)
            var = float(w @ (s ** 2 + (m - mu) ** 2))
            self._comp_w, self._comp_m, self._comp_s = w, (m - mu) / np.sqrt(var), s / np.sqrt(var)
            self._kind = "mix"
        elif name == "bimodal_uneven":
            w = np.array([0.75, 0.25]); m = np.array([-1.0, 2.0]); s = np.array([0.45, 0.45])
            mu = float(w @ m)
            var = float(w @ (s ** 2 + (m - mu) ** 2))
            self._comp_w, self._comp_m, self._comp_s = w, (m - mu) / np.sqrt(var), s / np.sqrt(var)
            self._kind = "mix"
        elif name == "lognormal":
            self._sdlog = 0.9
            self._kind = "lognormal"
            m = np.exp(self._sdlog ** 2 / 2.0)
            v = (np.exp(self._sdlog ** 2) - 1.0) * np.exp(self._sdlog ** 2)
            self._ln_m, self._ln_sd = m, np.sqrt(v)
        else:
            raise ValueError(name)

    # -- cdf / pdf / ppf of the standardised base
    def cdf(self, x):
        x = np.asarray(x, dtype=float)
        if self._kind == "mix":
            return sum(w * sps.norm.cdf(x, m, s)
                       for w, m, s in zip(self._comp_w, self._comp_m, self._comp_s))
        u = x * self._ln_sd + self._ln_m
        out = np.zeros_like(u)
        pos = u > 0
        out[pos] = sps.norm.cdf(np.log(u[pos]) / self._sdlog)
        return out

    def pdf(self, x):
        x = np.asarray(x, dtype=float)
        if self._kind == "mix":
            return sum(w * sps.norm.pdf(x, m, s)
                       for w, m, s in zip(self._comp_w, self._comp_m, self._comp_s))
        u = x * self._ln_sd + self._ln_m
        out = np.zeros_like(u)
        pos = u > 0
        out[pos] = sps.norm.pdf(np.log(u[pos]) / self._sdlog) / (u[pos] * self._sdlog) * self._ln_sd
        return out

    def rvs(self, n, rng):
        if self._kind == "mix":
            j = rng.choice(len(self._comp_w), size=n, p=self._comp_w)
            return rng.normal(self._comp_m[j], self._comp_s[j])
        return (np.exp(rng.normal(0.0, self._sdlog, size=n)) - self._ln_m) / self._ln_sd

    def support(self):
        if self._kind == "lognormal":
            return (-self._ln_m / self._ln_sd + 1e-12, 60.0)
        return (-12.0, 12.0)


# ------------------------------------------------------------ discretisation D

class Dmap:
    """Monotone non-decreasing map applied to the latent.

    Needs two things: apply() for simulation, and ulim(y) = sup{u : D(u) <= y}
    so that  P(Y <= y | g) = P(Y* <= ulim(y) | g).
    """

    def __init__(self, kind: str, **kw):
        self.kind = kind
        self.kw = kw

    def apply(self, y):
        k = self.kind
        if k == "identity":
            return y
        if k == "round":
            r = self.kw["step"]
            return np.round(y / r) * r
        if k == "floor":
            return np.maximum(y, self.kw["at"])
        if k == "ceiling":
            return np.minimum(y, self.kw["at"])
        if k == "clip_round":
            r = self.kw["step"]; lo = self.kw["lo"]; hi = self.kw["hi"]
            return np.clip(np.round(y / r) * r, lo, hi)
        raise ValueError(k)

    def ulim(self, y):
        """sup{u : D(u) <= y}; -inf if no such u."""
        k = self.kind
        y = np.asarray(y, dtype=float)
        if k == "identity":
            return y
        if k == "round":
            r = self.kw["step"]
            # D(u) <= y  <=>  round(u/r) <= floor(y/r + eps)  <=> u < (m+0.5) r
            m = np.floor(y / r + 1e-9)
            return (m + 0.5) * r
        if k == "floor":
            c = self.kw["at"]
            return np.where(y >= c - 1e-8, y, -np.inf)
        if k == "ceiling":
            c = self.kw["at"]
            return np.where(y >= c - 1e-8, np.inf, y)
        if k == "clip_round":
            r = self.kw["step"]; lo = self.kw["lo"]; hi = self.kw["hi"]
            m = np.floor(y / r + 1e-9)
            u = (m + 0.5) * r
            u = np.where(y < lo - 1e-8, -np.inf, u)
            u = np.where(y >= hi - 1e-8, np.inf, u)
            return u
        raise ValueError(k)


# --------------------------------------------------------------- genotype grid

def geno_grid(maf: float):
    p = maf
    probs = np.array([(1 - p) ** 2, 2 * p * (1 - p), p ** 2])
    g = np.array([0.0, 1.0, 2.0])
    gbar = float(probs @ g)
    var = float(probs @ (g - gbar) ** 2)
    return g, probs, gbar, var


# ------------------------------------------------------------------- quadrature

def _latent_grid(base: Base, a: float, s: float, n=200001):
    lo, hi = base.support()
    # widest possible observed range over the three genotypes
    los = [k * a + (1 + k * s) * lo for k in (0, 1, 2)]
    his = [k * a + (1 + k * s) * hi for k in (0, 1, 2)]
    return np.linspace(min(los), max(his), n)


class Truth:
    """All population quantities for one (base, D, maf, a, s, N) cell."""

    def __init__(self, base: Base, D: Dmap, maf: float, a: float, s: float,
                 taus: np.ndarray, N: int, ngrid: int = 400001,
                 h_override: float | None = None):
        self.base, self.D, self.maf, self.a, self.s = base, D, maf, a, s
        self.taus = np.asarray(taus, dtype=float)
        self.N = N
        g, probs, gbar, gvar = geno_grid(maf)
        self.g, self.probs, self.gbar, self.gvar = g, probs, gbar, gvar

        yg = _latent_grid(base, a, s, ngrid)
        self.yg = yg
        self.Fstar = np.stack([base.cdf((yg - k * a) / (1 + k * s)) for k in g])
        self.fstar = np.stack([base.pdf((yg - k * a) / (1 + k * s)) / (1 + k * s) for k in g])
        self.Fstar_mix = probs @ self.Fstar
        self.fstar_mix = probs @ self.fstar

        self.yobs = D.apply(yg)

        w = self.fstar_mix
        m1 = np.trapezoid(self.yobs * w, yg)
        m2 = np.trapezoid(self.yobs ** 2 * w, yg)
        self.sd_obs = float(np.sqrt(max(m2 - m1 ** 2, 1e-300)))
        self.h = 1.06 * self.sd_obs * N ** (-0.2) if h_override is None else float(h_override)

        # ---------------- latent estimand
        self.q_latent = np.array([self._inv(self.Fstar_mix, yg, t) for t in self.taus])
        self.f_latent = np.interp(self.q_latent, yg, self.fstar_mix)
        Fk_l = np.stack([np.interp(self.q_latent, yg, self.Fstar[k]) for k in range(3)])
        self.b_latent = -(probs @ ((g - gbar)[:, None] * Fk_l)) / (self.f_latent * gvar)

        # ---------------- observed distribution
        sup = np.unique(np.round(self.yobs, 10))
        self.sup = sup
        Fsup_k = np.stack([self._cdf_obs_k(k, sup) for k in range(3)])
        Fsup = probs @ Fsup_k
        self.Fsup = Fsup
        pmf = np.diff(np.concatenate([[0.0], Fsup]))          # mass at each support value
        self.pmf = np.clip(pmf, 0.0, None)

        if D.kind == "identity":
            self.q_obs = self.q_latent.copy()
        else:
            idx = np.searchsorted(Fsup, self.taus - 1e-12)
            idx = np.clip(idx, 0, len(sup) - 1)
            self.q_obs = sup[idx]

        Fk_o = np.stack([self._cdf_obs_k(k, self.q_obs) for k in range(3)])
        self.F_obs = probs @ Fk_o

        # Atom detection on the latent grid: a support value v carries an atom
        # iff {u : D(u) = v} has positive Lebesgue measure, i.e. iff many grid
        # points map to it.  Mass is integrated over that latent interval.
        inv = np.searchsorted(sup, np.round(self.yobs, 10))
        inv = np.clip(inv, 0, len(sup) - 1)
        counts = np.bincount(inv, minlength=len(sup))
        dy = yg[1] - yg[0]
        mass = np.bincount(inv, weights=self.fstar_mix * dy, minlength=len(sup))
        self.is_atom_sup = counts > 3
        ai = np.clip(np.searchsorted(sup, self.q_obs - 1e-9), 0, len(sup) - 1)
        self.atom = np.where(self.is_atom_sup[ai], mass[ai], 0.0)
        self.atom_all = np.where(self.is_atom_sup, mass, 0.0)

        ul = D.ulim(self.q_obs)
        with np.errstate(invalid="ignore"):
            f_true = np.where(np.isfinite(ul), np.interp(ul, yg, self.fstar_mix), np.nan)
        self.f_obs_true = np.where(self.atom > 0, np.nan, f_true)

        # ---------------- shipped estimand
        self.f_ship = np.array([self._kde_pop(q) for q in self.q_obs])
        self.b_ship = -(probs @ ((g - gbar)[:, None] * Fk_o)) / (self.f_ship * gvar)
        with np.errstate(invalid="ignore", divide="ignore"):
            self.b_obs = -(probs @ ((g - gbar)[:, None] * Fk_o)) / (self.f_obs_true * gvar)
        self.tie_frac = self.atom

    @staticmethod
    def _inv(F, y, t):
        i = int(np.searchsorted(F, t))
        i = min(max(i, 1), len(y) - 1)
        f0, f1 = F[i - 1], F[i]
        if f1 == f0:
            return float(y[i])
        return float(y[i - 1] + (t - f0) * (y[i] - y[i - 1]) / (f1 - f0))

    def _cdf_obs_k(self, k, y):
        """P(Y <= y | g = k), vectorised."""
        y = np.atleast_1d(np.asarray(y, dtype=float))
        ul = np.atleast_1d(self.D.ulim(y)).astype(float)
        out = np.empty_like(ul)
        fin = np.isfinite(ul)
        if np.any(fin):
            out[fin] = self.base.cdf((ul[fin] - self.g[k] * self.a) / (1 + self.g[k] * self.s))
        out[~fin & (ul < 0)] = 0.0
        out[~fin & (ul > 0)] = 1.0
        return out

    def _kde_pop(self, q):
        z = (q - self.yobs) / self.h
        kk = np.exp(-0.5 * z ** 2) / (self.h * np.sqrt(2 * np.pi))
        return float(np.trapezoid(kk * self.fstar_mix, self.yg))


def simulate(base: Base, D: Dmap, maf_vec, a, s, N, rng):
    """Return (G (N,M), Y (N,))."""
    M = len(maf_vec)
    G = np.empty((N, M))
    for j, p in enumerate(maf_vec):
        G[:, j] = rng.binomial(2, p, size=N).astype(float)
    return G


def make_pheno(base: Base, D: Dmap, g, a, s, N, rng):
    Z = base.rvs(N, rng)
    Ystar = g * a + (1.0 + g * s) * Z
    return D.apply(Ystar), Ystar
