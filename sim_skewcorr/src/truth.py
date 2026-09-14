"""Exact population quantities for the simulation cells.

Everything here is KNOWN TRUTH: the data-generating law is written down, so the
unconditional CDF, density, quantiles and the population RIF-regression estimand
b(tau) are all computable to machine precision by 3-point sums over g plus a
1-D root find.  No simulation is used to define truth.

Data-generating model, for genotype g in {0,1,2} at HWE frequency p:

    gc  = g - 2p
    s_g = exp(gamma * gc)
    Y   = delta * gc + s_g * e ,      e ~ F0 standardised to mean 0, variance 1

so that delta is a pure location effect in SD/allele and gamma is a pure
log-scale (dispersion) effect in log-SD/allele.  Covariates, when present, add
an independent normal term and are handled by `covar_sd` (a convolution that is
closed form for the normal base and done by quadrature for the skewed one).

The RIF-regression estimand at level tau is

    b(tau) = - E[ gc * F(q_tau | g) ] / ( f(q_tau) * Var(g) )

with q_tau the tau-quantile of the UNCONDITIONAL Y and f its density there.
"""

from __future__ import annotations

import numpy as np
from scipy import optimize, stats

TAUS = np.linspace(0.02, 0.98, 45)


# --------------------------------------------------------------------------
# base (standardised) error distributions
# --------------------------------------------------------------------------
class NormalBase:
    name = "normal"

    def cdf(self, x):
        return stats.norm.cdf(x)

    def pdf(self, x):
        return stats.norm.pdf(x)

    def rvs(self, rng, n):
        return rng.standard_normal(n)


class SkewBase:
    """Standardised chi-square with 4 df: skewness 1.414, excess kurtosis 3."""

    name = "skew"
    df = 4.0

    def __init__(self):
        self.m = self.df
        self.s = np.sqrt(2.0 * self.df)

    def cdf(self, x):
        return stats.chi2.cdf(x * self.s + self.m, self.df)

    def pdf(self, x):
        return stats.chi2.pdf(x * self.s + self.m, self.df) * self.s

    def rvs(self, rng, n):
        return (rng.chisquare(self.df, n) - self.m) / self.s



class HeavyBase:
    """Standardised Student t with 6 df: skewness 0, excess kurtosis 3."""

    name = "heavy"
    df = 6.0

    def __init__(self):
        self.s = np.sqrt(self.df / (self.df - 2.0))

    def cdf(self, x):
        return stats.t.cdf(x * self.s, self.df)

    def pdf(self, x):
        return stats.t.pdf(x * self.s, self.df) * self.s

    def rvs(self, rng, n):
        return rng.standard_t(self.df, n) / self.s


class VeryHeavyBase:
    """Standardised Student t with 4 df: excess kurtosis infinite (df=4 -> kurtosis undefined)."""

    name = "vheavy"
    df = 4.5

    def __init__(self):
        self.s = np.sqrt(self.df / (self.df - 2.0))

    def cdf(self, x):
        return stats.t.cdf(x * self.s, self.df)

    def pdf(self, x):
        return stats.t.pdf(x * self.s, self.df) * self.s

    def rvs(self, rng, n):
        return rng.standard_t(self.df, n) / self.s


class MildSkewBase:
    """Standardised chi-square with 20 df: skewness 0.632, excess kurtosis 0.6."""

    name = "mildskew"
    df = 20.0

    def __init__(self):
        self.m = self.df
        self.s = np.sqrt(2.0 * self.df)

    def cdf(self, x):
        return stats.chi2.cdf(x * self.s + self.m, self.df)

    def pdf(self, x):
        return stats.chi2.pdf(x * self.s + self.m, self.df) * self.s

    def rvs(self, rng, n):
        return (rng.chisquare(self.df, n) - self.m) / self.s


BASES = {"normal": NormalBase(), "skew": SkewBase(), "heavy": HeavyBase(),
         "vheavy": VeryHeavyBase(), "mildskew": MildSkewBase()}


# --------------------------------------------------------------------------
# genotype law
# --------------------------------------------------------------------------
def geno_law(p):
    gvals = np.array([0.0, 1.0, 2.0])
    pr = np.array([(1 - p) ** 2, 2 * p * (1 - p), p ** 2])
    gc = gvals - 2.0 * p
    var_g = 2.0 * p * (1 - p)
    return gc, pr, var_g


# --------------------------------------------------------------------------
# unconditional law of Y
# --------------------------------------------------------------------------
class Cell:
    """One simulation cell: MAF, base shape, location effect, dispersion effect."""

    def __init__(self, maf, shape, delta, gamma, covar_sd=0.0, n_quad=4001):
        self.maf = float(maf)
        self.shape = shape
        self.delta = float(delta)
        self.gamma = float(gamma)
        self.covar_sd = float(covar_sd)
        self.base = BASES[shape]
        self.gc, self.pr, self.var_g = geno_law(self.maf)
        self.n_quad = int(n_quad)
        if self.covar_sd > 0:
            # Gauss-Hermite nodes for the independent normal covariate term.
            x, w = np.polynomial.hermite_e.hermegauss(61)
            self._cx = x * self.covar_sd
            self._cw = w / np.sqrt(2.0 * np.pi)
        else:
            self._cx = np.array([0.0])
            self._cw = np.array([1.0])

    # conditional on g
    def _cdf_g(self, y, k):
        s = np.exp(self.gamma * self.gc[k])
        mu = self.delta * self.gc[k]
        y = np.atleast_1d(np.asarray(y, dtype=float))
        # convolve over the covariate term
        z = (y[:, None] - mu - self._cx[None, :]) / s
        return (self.base.cdf(z) * self._cw[None, :]).sum(axis=1)

    def _pdf_g(self, y, k):
        s = np.exp(self.gamma * self.gc[k])
        mu = self.delta * self.gc[k]
        y = np.atleast_1d(np.asarray(y, dtype=float))
        z = (y[:, None] - mu - self._cx[None, :]) / s
        return (self.base.pdf(z) * self._cw[None, :]).sum(axis=1) / s

    def cdf(self, y):
        out = np.zeros(np.atleast_1d(y).shape)
        for k in range(3):
            out = out + self.pr[k] * self._cdf_g(y, k)
        return out

    def pdf(self, y):
        out = np.zeros(np.atleast_1d(y).shape)
        for k in range(3):
            out = out + self.pr[k] * self._pdf_g(y, k)
        return out

    def quantile(self, taus):
        lo, hi = -60.0, 60.0
        out = np.empty(len(taus))
        for i, t in enumerate(taus):
            out[i] = optimize.brentq(lambda y: self.cdf(y)[0] - t, lo, hi, xtol=1e-12, rtol=1e-14)
        return out

    def moments(self):
        """Mean and variance of the unconditional Y (closed form)."""
        # E[Y] = sum_k pr_k * delta*gc_k  (base has mean 0, covariate mean 0)
        m = float((self.pr * self.delta * self.gc).sum())
        # E[Y^2 | g] = (delta gc)^2 + s^2 + covar_sd^2
        s2 = np.exp(2.0 * self.gamma * self.gc)
        ey2 = float((self.pr * ((self.delta * self.gc) ** 2 + s2 + self.covar_sd ** 2)).sum())
        return m, ey2 - m * m

    def b_tau(self, taus=TAUS):
        """Population RIF-regression estimand b(tau) -- exact."""
        q = self.quantile(taus)
        f = self.pdf(q)
        num = np.zeros(len(taus))
        for k in range(3):
            num += self.pr[k] * self.gc[k] * self._cdf_g(q, k)
        return -num / (f * self.var_g), q, f

    def sigma0_diag(self, taus=TAUS):
        """Population diagonal of Sigma_0 = Var(psi_tau) (used for |Z| scaling).

        psi_tau = q_tau + (tau - 1{Y<=q_tau})/f(q_tau); Var = tau(1-tau)/f^2.
        """
        q = self.quantile(taus)
        f = self.pdf(q)
        return taus * (1.0 - taus) / (f * f)


def design_W(cell, taus=TAUS):
    """Location / log-scale weight matrix, built from POPULATION quantities.

    Column 1 is the location column (all ones): a location shift of size delta
    gives b(tau) = delta for every tau.
    Column 2 is z_tau/(2 sd) for a normal baseline -- exactly
    make_weights_vqtl() -- and its distribution-free generalisation
    (q_tau - m)/(2 var) otherwise.  Under a log-scale perturbation of size gamma
    the first-order estimand is b(tau) = gamma (q_tau - m), so the second
    parameter has true value 2 * var * gamma.
    """
    q = cell.quantile(taus)
    m, v = cell.moments()
    W = np.column_stack([np.ones(len(taus)), (q - m) / (2.0 * v)])
    return W, q, m, v
