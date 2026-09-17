"""Simulation and known-truth machinery for the G-MoQ validation suite.

WHAT THIS IS FOR
----------------
Tiers 1-3 need simulated data whose Stage 2 answer is known in advance. Two
different "known answers" are needed and they are not the same thing:

  theta_pop     the population projection A b(tau) of Stage 1's OWN estimand.
                b(tau) is the population OLS slope of the RIF on the genotype,
                computed here EXACTLY (no Monte Carlo) from the conditional
                CDF. Coverage of theta_pop is a statement about the estimator
                and should be nominal whatever the model does.

  theta_latent  the generative parameter that produced the data. Coverage of
                theta_latent is only expected when the fitted W actually spans
                the truth; where it does not, the gap is misspecification, not
                a software defect.

Reporting only one of them is how a suite convinces itself of something false.
Both are returned.

DESIGN
------
The phenotype is defined by its CONDITIONAL QUANTILE FUNCTION given genotype,

    Q(tau | g) = mu0 + a.g + s0 exp(b.g/s0) * S(z_tau ; c0 + c.g, exp(d0 + d.g))
    S(z; e, d)  = sinh((asinh(z) + e) / d)          [sinh-arcsinh]

S is strictly increasing in z for every finite e and every d > 0, so every
effect size produces a genuine distribution -- which a Cornish-Fisher cubic
does not: measured, its monotone range caps a kurtosis effect at |Z| ~ 1.5 at
N = 50,000, and that cap would be indistinguishable in the results from a
property of the method. The scale enters through an exponential for the same
reason: an additive scale coefficient large enough to be interesting drives the
conditional standard deviation NEGATIVE for some genotypes. The exponential is
what sim_meanvar/src used (archived at tag archive/sim-campaigns-20260917).

The four derivatives at the baseline are closed form (`cf_weight_matrix`), so

because the derivative of exactly that expression is the weight matrix the
suite fits (tests/lib/models.R, gmoq_cf_columns). So for small effects

    theta_latent = (a, b, c, d)

with no algebra in between, and the second-order departure from it is the
genotype-skewness residue that docs/FINDINGS_pooled_dispersion.md quantifies --
which the suite wants PRESENT, not removed.

`c0` and `d0` set the BASELINE skewness and kurtosis of the phenotype, so a
skewed trait is the same machinery and not a second code path.

The mixture families use the analytic two-Normal mixture instead, because
their weight columns are the shipped make_weights_normal_mixture_vqtl().

Stage 1 is the SHIPPED code: fungwas_stage1.core.run_stage1 and
fungwas_stage1.pooled_cov.compute_sigma0, exactly as sim_meanvar/src and
sim/density/src drove it (both archived at tag
archive/sim-campaigns-20260917). Nothing about the estimator is reimplemented
here.
"""

from __future__ import annotations

import gzip
import json
import os
import struct
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.special import ndtr, ndtri, expit, logit


# --------------------------------------------------------------------------
# locating the shipped Stage 1
# --------------------------------------------------------------------------

def add_stage1_to_path(repo_root: Optional[str] = None) -> str:
    """Put the Stage 1 package on sys.path and return the directory used.

    GMOQ_STAGE1_SRC wins, then an explicit repo_root, then the repo this file
    lives in. Never a hard-coded absolute path: the existing sim drivers pin
    one, and it silently analyses a different tree than the one under test.
    """
    src = os.environ.get("GMOQ_STAGE1_SRC")
    if not src:
        root = repo_root or os.environ.get("GMOQ_REPO")
        if not root:
            root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
        src = os.path.join(root, "stage1")
    src = os.path.abspath(src)
    if src not in sys.path:
        sys.path.insert(0, src)
    return src


PROD_TAUS = np.round(np.linspace(0.02, 0.98, 45), 3)


# --------------------------------------------------------------------------
# the Cornish-Fisher conditional law
# --------------------------------------------------------------------------

def sas(z, e, d):
    """sinh-arcsinh shape transform: S(z; e, d) = sinh((asinh(z) + e)/d).

    e skews, d controls tail weight, and S is STRICTLY INCREASING in z for
    every finite e and every d > 0. That is why it is used here rather than a
    Cornish-Fisher polynomial: a cubic quantile expansion stops being monotone
    once the shape coefficients grow, which caps the largest shape effect that
    can be simulated at all (measured: |Z| ~ 1.5 on a kurtosis parameter at
    N = 50,000) and makes the cap look like a property of the method rather
    than of the simulation. The sinh-arcsinh family has no such cap, and its
    inverse is closed form, so no root-finding enters the truth at all.
    """
    return np.sinh((np.arcsinh(z) + e) / d)


def sas_inverse(w, e, d):
    """Exact inverse of `sas`: z = sinh(d asinh(w) - e)."""
    return np.sinh(d * np.arcsinh(w) - e)


def sas_dz(z, e, d):
    """dS/dz, needed for the baseline density."""
    return np.cosh((np.arcsinh(z) + e) / d) / (d * np.sqrt(1.0 + z * z))


# --------------------------------------------------------------------------
# the eight families
# --------------------------------------------------------------------------

FAMILIES = ("null", "location", "scale", "locscale", "skew", "kurtosis",
            "mix_weight", "mix_mean")

#: Which parameter each family perturbs, and which model in tests/lib/models.R
#: is the parsimonious truth.
FAMILY_PARAM = {
    "null":       (),
    "location":   ("a",),
    "scale":      ("b",),
    "locscale":   ("a", "b"),
    "skew":       ("a", "b", "c"),
    "kurtosis":   ("a", "b", "d"),
    "mix_weight": ("e",),
    "mix_mean":   ("m",),
}


@dataclass
class Design:
    """One simulated replicate's variant configuration."""
    family: str
    maf: np.ndarray                 # (M,)
    causal: np.ndarray              # indices of effect-carrying variants
    a: np.ndarray = field(default=None)   # location coefficient per variant
    b: np.ndarray = field(default=None)   # scale coefficient: sd = s0 exp(b.g/s0)
    c: np.ndarray = field(default=None)   # skew coefficient: e = c0 + c.g
    d: np.ndarray = field(default=None)   # tail-weight coefficient: log d = d0 + d.g
    e: np.ndarray = field(default=None)   # mixture logit coefficient
    m: np.ndarray = field(default=None)   # mixture component-1 mean coefficient
    mu0: float = 0.0
    s0: float = 1.0
    c0: float = 0.0                 # baseline sinh-arcsinh skewness of the phenotype
    d0: float = 0.0                 # baseline log tail-weight
    mix: Optional[dict] = None      # p1, mu1, sd1, mu2, sd2 for the mixture families
    # Post-hoc discretisation of the SIMULATED phenotype, used only by the
    # tier 2 failure-mode pins. One of
    #   {"kind": "round",   "h": <step>}       ties at multiples of h
    #   {"kind": "ceiling", "tau": <level>}    everything above q_tau piled on it
    #   {"kind": "floor",   "tau": <level>}    everything below q_tau piled on it
    # These change the ESTIMAND as well as the estimator's conditioning (see
    # docs/BP2_lumpy_phenotypes.md section 4), so Truth does NOT describe the
    # observed scale for them and the pins use them only under the null.
    discretise: Optional[dict] = None
    #: The |Z| the design was calibrated for, recorded so a replicate can be
    #: checked against what was asked for.
    z_requested: float = float("nan")
    #: Factor the SHAPE coefficients were shrunk by to keep the simulated
    #: distribution within a plausible range. 1.0 means the design was used as
    #: asked for; below 1.0 the realised |Z| is shrink x z_requested.
    shrink: float = 1.0

    def __post_init__(self):
        M = len(self.maf)
        for nm in ("a", "b", "c", "d", "e", "m"):
            if getattr(self, nm) is None:
                setattr(self, nm, np.zeros(M))

    @property
    def is_mixture(self) -> bool:
        return self.family.startswith("mix")


def _lin(vec, G):
    """G (n, M) times a per-variant coefficient vector."""
    return G @ vec


class Truth:
    """Exact population quantities for one Design.

    The marginal law of Y is the average of the conditional law over the
    genotypes of the effect-carrying variants, taken on a large common
    background sample (common random numbers, so the CONTRAST across a single
    variant's genotype -- which is all b_l(tau) depends on -- is estimated far
    more precisely than the level).
    """

    def __init__(self, design: Design, nbg: int = 20_000, seed: int = 7):
        self.D = design
        rng = np.random.default_rng(seed)
        idx = design.causal
        self.idx = idx
        self.gb = rng.binomial(2, design.maf[idx], size=(nbg, len(idx))).astype(float)
        self.nbg = nbg

    # ---- conditional CDF -------------------------------------------------
    def _cdf_cf(self, y, A, B, C, D):
        """P(Y <= y | linear predictors), y scalar or (ny,), A..D (n,).

        Closed form: Y = mu + sigma sinh((asinh(Z) + e)/d) with Z standard
        normal, so F(y) = Phi(sinh(d asinh((y - mu)/sigma) - e)) exactly.
        """
        y = np.atleast_1d(np.asarray(y, dtype=float))[:, None]
        s = self.D.s0 * np.exp(B[None, :] / self.D.s0)
        w = (y - (self.D.mu0 + A[None, :])) / s
        e = self.D.c0 + C[None, :]
        d = np.exp(self.D.d0 + D[None, :])
        return ndtr(sas_inverse(w, e, d))

    def _cdf_mix(self, y, E, Mm):
        y = np.atleast_1d(np.asarray(y, dtype=float))[:, None]
        mx = self.D.mix
        p1 = expit(logit(mx["p1"]) + E[None, :])
        mu1 = mx["mu1"] + Mm[None, :]
        return (p1 * ndtr((y - mu1) / mx["sd1"])
                + (1 - p1) * ndtr((y - mx["mu2"]) / mx["sd2"]))

    def _predictors(self, gb):
        D = self.D
        i = self.idx
        return (gb @ D.a[i], gb @ D.b[i], gb @ D.c[i], gb @ D.d[i],
                gb @ D.e[i], gb @ D.m[i])

    def cdf(self, y, gb=None):
        """Marginal CDF of Y at y (averaged over the background)."""
        gb = self.gb if gb is None else gb
        A, B, C, Dd, E, Mm = self._predictors(gb)
        if self.D.is_mixture:
            return self._cdf_mix(y, E, Mm).mean(axis=1)
        return self._cdf_cf(y, A, B, C, Dd).mean(axis=1)

    def pdf(self, y, h=None):
        y = np.atleast_1d(np.asarray(y, dtype=float))
        if h is None:
            h = 1e-4 * max(self.D.s0, 1.0)
        return (self.cdf(y + h) - self.cdf(y - h)) / (2 * h)

    def quantiles(self, taus, coarse: int = 4000):
        """Marginal quantiles: a coarse grid pass to bracket, then Newton.

        The grid pass is the expensive one (a grid x background outer product),
        so it uses a SUBSAMPLE of the background. The Newton refinement that
        follows uses the full background and is what sets the accuracy, so the
        subsample costs nothing but is worth an order of magnitude in time.
        """
        taus = np.asarray(taus, dtype=float)
        if self.D.is_mixture:
            mx = self.D.mix
            lo = min(mx["mu1"], mx["mu2"]) - 12 * max(mx["sd1"], mx["sd2"])
            hi = max(mx["mu1"], mx["mu2"]) + 12 * max(mx["sd1"], mx["sd2"])
        else:
            # Bracket by the CONDITIONAL quantile function's own range at
            # |z| = 6 over the background rather than by a fixed multiple of
            # s0, which a strongly shaped cell would overflow.
            A, B, C, Dd, _, _ = self._predictors(self.gb)
            sc = self.D.s0 * np.exp(B / self.D.s0)
            e = self.D.c0 + C; d = np.exp(self.D.d0 + Dd)
            lo = float(np.min(self.D.mu0 + A + sc * sas(-6.0, e, d)))
            hi = float(np.max(self.D.mu0 + A + sc * sas(6.0, e, d)))
        grid = np.linspace(lo, hi, 1201)
        sub = self.gb if self.gb.shape[0] <= coarse else self.gb[:coarse]
        F = self.cdf(grid, sub)
        q = np.interp(taus, F, grid)
        for _ in range(8):
            f = np.maximum(self.pdf(q), 1e-12)
            q = q - (self.cdf(q) - taus) / f
        return q

    # ---- the exact Stage 1 estimand --------------------------------------
    def b_true(self, taus):
        """b_l(tau) for every effect-carrying variant, EXACTLY.

        b_l(tau) = -Cov_g(F(q_tau | g_l = g), g_l) / (f(q_tau) Var(g_l))
        with the conditional CDF averaged over the other variants' genotypes.
        """
        D = self.D
        q = self.quantiles(taus)
        f = np.maximum(self.pdf(q), 1e-300)
        out = np.zeros((len(self.idx), len(taus)))
        for r, l in enumerate(self.idx):
            p = D.maf[l]
            pg = np.array([(1 - p) ** 2, 2 * p * (1 - p), p ** 2])
            vg = 2 * p * (1 - p)
            gb_others = self.gb.copy()
            cov = np.zeros(len(taus))
            for k, gval in enumerate((0.0, 1.0, 2.0)):
                gb_others[:, r] = gval
                Fk = self.cdf(q, gb_others)
                cov += pg[k] * (gval - 2 * p) * Fk
            out[r] = -cov / (f * vg)
        return q, f, out

    def theta_latent(self, taus, sigma):
        """The generative parameter, expressed on the suite's W columns."""
        D = self.D
        K = {"beta_mu": D.a, "beta_sigma": D.b, "beta_skew": D.c, "beta_kurt": D.d,
             "gamma": D.e, "beta_mu1": D.m}
        return {k: v[self.idx] for k, v in K.items()}


# --------------------------------------------------------------------------
# simulating a phenotype
# --------------------------------------------------------------------------

def apply_discretisation(y: np.ndarray, spec: Optional[dict]) -> np.ndarray:
    """Round or censor a simulated phenotype, for the tier 2 pins."""
    if not spec:
        return y
    kind = spec["kind"]
    if kind == "round":
        h = float(spec["h"]) * float(np.std(y))
        return np.round(y / h) * h
    if kind == "ceiling":
        lim = float(np.quantile(y, float(spec["tau"]), method="median_unbiased"))
        return np.minimum(y, lim)
    if kind == "floor":
        lim = float(np.quantile(y, float(spec["tau"]), method="median_unbiased"))
        return np.maximum(y, lim)
    raise ValueError(f"unknown discretisation {kind!r}")


def simulate_phenotype(design: Design, G: np.ndarray, rng) -> np.ndarray:
    """Draw Y from the conditional law defined by the design."""
    n = G.shape[0]
    u = rng.random(n)
    if design.is_mixture:
        mx = design.mix
        p1 = expit(logit(mx["p1"]) + _lin(design.e, G))
        mu1 = mx["mu1"] + _lin(design.m, G)
        pick = rng.random(n) < p1
        y = np.where(pick,
                     rng.normal(mu1, mx["sd1"]),
                     rng.normal(mx["mu2"], mx["sd2"], size=n))
        return apply_discretisation(y, design.discretise)
    z = ndtri(u)
    A = _lin(design.a, G)
    s = design.s0 * np.exp(_lin(design.b, G) / design.s0)
    e = design.c0 + _lin(design.c, G)
    d = np.exp(design.d0 + _lin(design.d, G))
    y = design.mu0 + A + s * sas(z, e, d)
    return apply_discretisation(y, design.discretise)


# --------------------------------------------------------------------------
# building a design for one cell
# --------------------------------------------------------------------------

def bridge_sigma0(taus, f):
    """Population Sigma_0: (min(s,t) - s t) / (f_s f_t)."""
    taus = np.asarray(taus, dtype=float); f = np.asarray(f, dtype=float)
    return (np.minimum.outer(taus, taus) - np.outer(taus, taus)) / np.outer(f, f)


def cf_weight_matrix(taus, sigma, e0=0.0, d0=0.0):
    """The four weight columns, matching gmoq_cf_columns() in tests/lib/models.R.

    Each is d q_tau / d theta_j for the sinh-arcsinh conditional quantile
    function at the baseline (e = e0, log d = d0):

        mu    : 1
        sigma : S(z; e0, exp(d0))
        skew  : sigma dS/de
        tail  : sigma dS/dlog(d)
    """
    z = ndtri(np.asarray(taus, dtype=float))
    d = np.exp(d0)
    u = (np.arcsinh(z) + e0) / d
    return np.column_stack([
        np.ones_like(z),
        np.sinh(u),
        sigma * np.cosh(u) / d,
        -sigma * np.cosh(u) * u])


def cf_baseline_density(taus, sigma, e0=0.0, d0=0.0):
    """f(q_tau) for Y = mu + sigma S(Z; e0, exp(d0)), analytically."""
    z = ndtri(np.asarray(taus, dtype=float))
    dq = sigma * sas_dz(z, e0, np.exp(d0))
    return np.exp(-0.5 * z * z) / np.sqrt(2 * np.pi) / dq


def mix_quantiles_density(taus, mx, grid=40001):
    """Quantiles and density of a two-Normal mixture, by inversion on a grid."""
    lo = min(mx["mu1"], mx["mu2"]) - 12 * max(mx["sd1"], mx["sd2"])
    hi = max(mx["mu1"], mx["mu2"]) + 12 * max(mx["sd1"], mx["sd2"])
    yy = np.linspace(lo, hi, grid)
    F = mx["p1"] * ndtr((yy - mx["mu1"]) / mx["sd1"]) + \
        (1 - mx["p1"]) * ndtr((yy - mx["mu2"]) / mx["sd2"])
    q = np.interp(np.asarray(taus, dtype=float), F, yy)
    f = (mx["p1"] * np.exp(-0.5 * ((q - mx["mu1"]) / mx["sd1"]) ** 2)
         / (mx["sd1"] * np.sqrt(2 * np.pi))
         + (1 - mx["p1"]) * np.exp(-0.5 * ((q - mx["mu2"]) / mx["sd2"]) ** 2)
         / (mx["sd2"] * np.sqrt(2 * np.pi)))
    return q, f


def mix_weight_matrix(taus, q, mx):
    """gamma, mu1, mu2, sigma1, sigma2 columns -- the shipped mixture weights."""
    f1 = np.exp(-0.5 * ((q - mx["mu1"]) / mx["sd1"]) ** 2) / (mx["sd1"] * np.sqrt(2 * np.pi))
    F1 = ndtr((q - mx["mu1"]) / mx["sd1"])
    f2 = np.exp(-0.5 * ((q - mx["mu2"]) / mx["sd2"]) ** 2) / (mx["sd2"] * np.sqrt(2 * np.pi))
    F2 = ndtr((q - mx["mu2"]) / mx["sd2"])
    p1 = mx["p1"]
    f = np.maximum(p1 * f1 + (1 - p1) * f2, 1e-12)
    return np.column_stack([
        -(p1 * (1 - p1) / f) * (F1 - F2),
        (p1 * f1) / f,
        ((1 - p1) * f2) / f,
        (p1 * f1 * (q - mx["mu1"]) / mx["sd1"]) / f,
        ((1 - p1) * f2 * (q - mx["mu2"]) / mx["sd2"]) / f])


def unit_parameter_se(W, Sigma0, N, vg):
    """sqrt(diag(A Sigma0 A')) / sqrt(N vg): the OLS parameter SE for one variant.

    This is the analytic standard error Stage 2 will report, so setting an
    effect to z_target times it makes the realised |Z| land near z_target
    instead of wherever a hand-picked multiplier happens to put it.
    """
    A = np.linalg.solve(W.T @ W, W.T)
    return np.sqrt(np.diag(A @ Sigma0 @ A.T)) / np.sqrt(N * vg)


def build_design(family: str, maf: float, n_causal: int, z_target: float,
                 N: int, taus, rng,
                 baseline_skew: float = 0.0, baseline_kurt: float = 0.0,
                 s0: float = 1.0, n_null: int = 200,
                 discretise: Optional[dict] = None) -> Design:
    """One replicate's design: n_null null variants and n_causal causal ones.

    The causal effect on the family's leading parameter is set to z_target
    times that parameter's ANALYTIC standard error, so the realised |Z| is
    close to the target by construction rather than by a tuned multiplier.
    Signs alternate so the causal variants do not shift the marginal mean.

    Nothing downstream depends on the target being hit: the truth is computed
    exactly from the design that results.
    """
    if family not in FAMILIES:
        raise ValueError(f"unknown family {family!r}; expected one of {FAMILIES}")
    M = n_null + n_causal
    mafv = np.full(M, float(maf))
    causal = np.arange(n_null, M)
    vg = 2 * maf * (1 - maf)
    sgn = np.where(np.arange(n_causal) % 2 == 0, 1.0, -1.0)
    taus = np.asarray(taus, dtype=float)

    mix = None
    if family.startswith("mix"):
        mix = dict(p1=0.45, mu1=-1.0, sd1=0.8, mu2=1.6, sd2=1.0)

    D = Design(family=family, maf=mafv, causal=causal, mu0=0.0, s0=s0,
               c0=baseline_skew, d0=baseline_kurt, mix=mix, discretise=discretise)
    if not n_causal or family == "null":
        return D

    if D.is_mixture:
        q, f = mix_quantiles_density(taus, mix)
        W = mix_weight_matrix(taus, q, mix)
        se = unit_parameter_se(W, bridge_sigma0(taus, f), N, vg)
        col = {"mix_weight": 0, "mix_mean": 1}[family]
        val = z_target * se[col]
        getattr(D, {"mix_weight": "e", "mix_mean": "m"}[family])[causal] = sgn * val
        return D

    f = cf_baseline_density(taus, s0, baseline_skew, baseline_kurt)
    W4 = cf_weight_matrix(taus, s0)
    S0 = bridge_sigma0(taus, f)
    # The standard error must come from the model that will ACTUALLY be fitted
    # for this family, not from the full four-column design: dropping a column
    # changes the remaining rows of (W'W)^-1 W' and therefore their standard
    # errors, and calibrating against the wrong one puts the realised |Z| a
    # factor of three away from the target.
    cols = {"location": [0, 1], "scale": [0, 1], "locscale": [0, 1],
            "skew": [0, 1, 2], "kurtosis": [0, 1, 3]}[family]
    se_fit = unit_parameter_se(W4[:, cols], S0, N, vg)
    se = np.zeros(4)
    for k, cidx in enumerate(cols):
        se[cidx] = se_fit[k]

    if family == "location":
        D.a[causal] = sgn * z_target * se[0]
    elif family == "scale":
        D.b[causal] = sgn * z_target * se[1]
    elif family == "locscale":
        D.a[causal] = sgn * z_target * se[0]
        D.b[causal] = sgn * z_target * se[1]
    elif family == "skew":
        D.a[causal] = sgn * 0.5 * z_target * se[0]
        D.b[causal] = sgn * 0.5 * z_target * se[1]
        D.c[causal] = sgn * z_target * se[2]
    elif family == "kurtosis":
        D.a[causal] = sgn * 0.5 * z_target * se[0]
        D.b[causal] = sgn * 0.5 * z_target * se[1]
        D.d[causal] = sgn * z_target * se[3]

    # No MONOTONICITY clamp is needed -- the sinh-arcsinh quantile function is
    # increasing for every finite e and every d > 0. A RANGE cap still is: the
    # shape columns are weakly identified (they retain about 1-5 per cent of
    # their norm after projecting out location and scale, measured in tier 0),
    # so their standard errors are large and a |Z| = 12 design asks for a
    # skewness or a tail weight far outside anything a real trait shows. The
    # cap is on the parameter itself, at 5 sd of its linear predictor:
    #   |e| <= 1.5   (sinh-arcsinh skewness; 1.5 is already extreme)
    #   |log d| <= 0.8  (tail weight between 0.45x and 2.2x)
    # When it binds, the realised |Z| is shrink x z_target and `shrink` says so,
    # which is reported rather than hidden.
    D.z_requested = float(z_target)
    shrink = 1.0
    for coef, cap in ((D.c, 1.5), (D.d, 0.8)):
        sd_lin = float(np.sqrt(np.sum(coef[causal] ** 2 * vg)))
        if 5.0 * sd_lin > cap:
            shrink = min(shrink, cap / (5.0 * sd_lin))
    if shrink < 1.0:
        D.c[causal] *= shrink
        D.d[causal] *= shrink
    D.shrink = shrink
    return D


# --------------------------------------------------------------------------
# running the SHIPPED Stage 1
# --------------------------------------------------------------------------

def run_stage1_on(design: Design, N: int, taus, seed: int,
                  n_blocks: int = 1000, n_covar: int = 0, n_threads: int = 1):
    """Simulate genotypes and a phenotype, then run the shipped Stage 1."""
    add_stage1_to_path()
    from fungwas_stage1 import core, pooled_cov  # noqa: E402

    rng = np.random.default_rng(seed)
    M = len(design.maf)
    G = rng.binomial(2, np.broadcast_to(design.maf, (N, M))).astype(np.float64)
    y = simulate_phenotype(design, G, rng)

    C = rng.normal(size=(N, n_covar)) if n_covar else None
    res = core.run_stage1(G, y, np.asarray(taus, dtype=float),
                          covariates=C, n_blocks=n_blocks, seed=seed,
                          n_threads=n_threads)

    # Sigma_0, exactly as Stage 1's own --sigma0 path forms it.
    rif, q_hat, fhat = core.compute_rif(y, np.asarray(taus, dtype=float))
    X = np.ones((N, 1)) if C is None else np.column_stack([np.ones(N), C])
    rif_r, Qb = core.residualize(rif, X)
    sigma0 = pooled_cov.compute_sigma0(rif_r)

    res["sigma0"] = sigma0
    res["fhat"] = fhat
    res["y_mean"] = float(y.mean())
    res["y_sd"] = float(y.std(ddof=1))
    res["n_samples"] = int(N)
    res["n_covariates"] = int(n_covar)
    res["covariate_names"] = [f"C{i+1}" for i in range(n_covar)]
    return res


# --------------------------------------------------------------------------
# writing a VALID Stage 1 artefact set
# --------------------------------------------------------------------------

def write_stage1_files(prefix: str, taus, res: dict, phenotype: str = "simpheno",
                       cov_dtype: str = "float32", chrom: str = "1",
                       snp_ids: Optional[Sequence[str]] = None,
                       write_sigma0_files: bool = True,
                       corrupt: Optional[str] = None) -> Dict[str, str]:
    """Write the exact file set the shipped Stage 1 CLI writes.

    Column names, number formats, the row-major packed upper triangle, the int8
    scale sidecar and both JSON sidecars all follow
    stage1/fungwas_stage1/cli.py and cov_meta.py / pooled_cov.py, and the JSON
    is built by those modules rather than hand-written, so a schema change in
    Stage 1 breaks this loudly instead of silently.

    `corrupt` deliberately damages one part of the artefact set, for the tier 1
    refusal tests: one of
        "cov_ids_shuffled", "no_cov_ids", "no_d_resid", "no_cov_meta",
        "sigma0_wrong_taus", "sigma0_wrong_covariates", "sigma0_wrong_n",
        "sigma0_wrong_loco_scale", "future_schema", "duplicate_snp_ids".
    """
    add_stage1_to_path()
    from fungwas_stage1 import cov_meta as cm
    from fungwas_stage1 import pooled_cov as pc

    taus = np.asarray(taus, dtype=float)
    T = len(taus)
    betas = np.asarray(res["betas"])
    ses = np.asarray(res["se"])
    covs = res["covariances"]
    d_resid = np.asarray(res["d_resid"])
    Mv = betas.shape[0]
    if snp_ids is None:
        snp_ids = [f"rs{j+1:07d}" for j in range(Mv)]

    n_elem = T * (T + 1) // 2
    triu = np.triu_indices(T)

    paths = {}
    p_tsv = f"{prefix}.stage1.tsv.gz"
    p_cov = f"{prefix}.cov.gz"
    p_ids = f"{prefix}.cov.ids.tsv.gz"
    p_scale = f"{prefix}.cov.scale.gz"
    p_meta = f"{prefix}.cov.meta.json"

    cols = ["snp_id", "chr", "bp", "effect_allele", "other_allele"]
    for t in taus:
        cols += [f"beta_tau_{t:.3f}", f"se_tau_{t:.3f}"]
    if corrupt != "no_d_resid":
        cols.append("d_resid")

    ids = list(snp_ids)
    if corrupt == "duplicate_snp_ids" and Mv >= 2:
        ids[1] = ids[0]

    with gzip.open(p_tsv, "wt") as f_tsv, gzip.open(p_cov, "wb") as f_cov, \
         gzip.open(p_ids, "wt") as f_ids:
        f_tsv.write("\t".join(cols) + "\n")
        f_ids.write("snp_id\tchr\tbp\teffect_allele\tother_allele\n")
        f_scale = gzip.open(p_scale, "wb") if cov_dtype == "int8" else None
        try:
            id_rows = []
            for j in range(Mv):
                row = [ids[j], chrom, j + 1, "A", "G"]
                for bb, ss in zip(betas[j], ses[j]):
                    row += [f"{bb:.6g}", f"{ss:.6g}"]
                if corrupt != "no_d_resid":
                    row.append(f"{d_resid[j]:.9g}")
                f_tsv.write("\t".join(map(str, row)) + "\n")
                id_rows.append([ids[j], chrom, j + 1, "A", "G"])

                cu = np.asarray(covs[j])[triu]
                if cov_dtype == "int8":
                    sc = float(np.max(np.abs(cu)) / 127.0)
                    if not np.isfinite(sc) or sc <= 0:
                        sc = 1.0
                    f_cov.write(np.clip(np.rint(cu / sc), -127, 127)
                                .astype(np.int8).tobytes(order="C"))
                    f_scale.write(struct.pack("f", sc))
                else:
                    f_cov.write(struct.pack(f"{n_elem}f", *cu.astype(np.float64)))

            if corrupt == "cov_ids_shuffled" and Mv >= 2:
                id_rows[0], id_rows[1] = id_rows[1], id_rows[0]
            for r in id_rows:
                f_ids.write("\t".join(map(str, r)) + "\n")
        finally:
            if f_scale is not None:
                f_scale.close()

    paths.update(stage1=p_tsv, cov=p_cov, cov_ids=p_ids, cov_meta=p_meta)
    if cov_dtype == "int8":
        paths["cov_scale"] = p_scale

    if corrupt == "no_cov_ids":
        os.remove(p_ids)
        paths.pop("cov_ids")

    # --- the Sigma_0 sidecar, and the analysis block both files share -----
    covars = list(res.get("covariate_names", []))
    n_samples = int(res.get("n_samples", 0))
    s0_taus = list(np.round(taus, 3))
    loco_mode = None
    if corrupt == "sigma0_wrong_taus":
        s0_taus = list(np.round(taus + 0.001, 3))
    if corrupt == "sigma0_wrong_covariates":
        covars = covars + ["an_extra_covariate"]
    if corrupt == "sigma0_wrong_n":
        n_samples = n_samples + 1
    if corrupt == "sigma0_wrong_loco_scale":
        loco_mode = "none"

    s0_meta = pc.build_sigma0_meta(
        phenotype=phenotype, taus=s0_taus, n_samples=n_samples,
        covariates=covars, covariate_basis_rank=len(covars) + 1,
        loco_predictor=None, loco_scale_mode=loco_mode, loco_scale_factor=None,
        chromosomes=[chrom])

    # The cov.meta analysis block must describe THE SCAN, so it is built from
    # the true values even when the sigma0 sidecar has been corrupted. That is
    # the whole point: Stage 2 compares the two and refuses a disagreement.
    true_meta = pc.build_sigma0_meta(
        phenotype=phenotype, taus=list(np.round(taus, 3)),
        n_samples=int(res.get("n_samples", 0)),
        covariates=list(res.get("covariate_names", [])),
        covariate_basis_rank=len(res.get("covariate_names", [])) + 1,
        loco_predictor=None, loco_scale_mode=None, loco_scale_factor=None,
        chromosomes=[chrom])

    meta = cm.build_cov_meta(
        cov_dtype=cov_dtype, n_taus=T, n_blocks=int(res["n_blocks"]),
        seed=42, n_blocks_nonempty=int(res["n_blocks_nonempty"]),
        analysis=pc.sigma0_analysis_block(true_meta))
    if corrupt == "future_schema":
        meta["cov_meta_version"] = 99
    if corrupt != "no_cov_meta":
        cm.write_cov_meta(p_meta, meta)
    else:
        paths.pop("cov_meta")

    if write_sigma0_files:
        pc.write_sigma0(prefix, {chrom: np.asarray(res["sigma0"])}, s0_meta)
        sp, mp = pc.sigma0_paths(prefix)
        paths.update(sigma0=sp, sigma0_meta=mp)

    return paths


# --------------------------------------------------------------------------
# the standard cell driver used by tiers 1-3
# --------------------------------------------------------------------------

def run_replicate(family: str, maf: float, N: int, taus, seed: int,
                  z_target: float = 12.0, n_null: int = 200, n_causal: int = 8,
                  baseline_skew: float = 0.0, baseline_kurt: float = 0.0,
                  n_blocks: int = 1000, n_covar: int = 0, nbg: int = 20_000,
                  n_threads: int = 1, want_truth: bool = True) -> dict:
    """One replicate: a design, the shipped Stage 1, and the exact truth."""
    rng = np.random.default_rng(seed)
    design = build_design(family, maf, n_causal, z_target, N, taus, rng,
                          baseline_skew=baseline_skew, baseline_kurt=baseline_kurt,
                          n_null=n_null)
    res = run_stage1_on(design, N, taus, seed, n_blocks=n_blocks,
                        n_covar=n_covar, n_threads=n_threads)
    out = dict(betas=res["betas"], se=res["se"], d_resid=res["d_resid"],
               sigma0=res["sigma0"], q_tau=res["q_tau"], fhat=res["fhat"],
               y_mean=res["y_mean"], y_sd=res["y_sd"],
               covariances=res["covariances"],
               n_blocks=res["n_blocks"], n_blocks_nonempty=res["n_blocks_nonempty"],
               maf=design.maf, causal=design.causal, family=family)
    if want_truth:
        tr = Truth(design, nbg=nbg, seed=seed + 991)
        q_pop, f_pop, b_c = tr.b_true(taus)
        b_true = np.zeros_like(res["betas"])
        b_true[design.causal] = b_c
        out.update(b_true=b_true, q_pop=q_pop, f_pop=f_pop,
                   theta_latent=tr.theta_latent(taus, res["y_sd"]))
    return out
