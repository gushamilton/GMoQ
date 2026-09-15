#!/usr/bin/env python3
"""Bounded reproduction of the testosterone two-density Stage 1/Stage 2 gap.

This intentionally uses the production Stage 1 RIF calculation with an oracle
two-normal mixture nuisance fit.  It isolates the density interface: Stage 1
uses its KDE at empirical quantiles, while the historical mixture weights use
the fitted mixture density at those same quantiles.  No UK Biobank inputs are
read and no paper outputs are touched.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import ndtr


TAUS = np.round(np.linspace(0.02, 0.98, 45), 3)
PARAMETERS = ("beta_lower", "beta_upper", "gamma_upper")


def type8_quantile(y: np.ndarray, probs: np.ndarray) -> np.ndarray:
    """R type-8 / NumPy median-unbiased quantiles without version dependence."""
    z = np.sort(np.asarray(y, dtype=float))
    n = z.size
    h = (n + 1.0 / 3.0) * probs + 1.0 / 3.0
    j = np.floor(h).astype(int)
    g = h - j
    j0 = np.clip(j - 1, 0, n - 1)
    j1 = np.clip(j, 0, n - 1)
    return (1.0 - g) * z[j0] + g * z[j1]


def production_rif(y: np.ndarray, taus: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The calculation in stage1/fungwas_stage1/core.py::compute_rif."""
    n = y.size
    q = type8_quantile(y, taus)
    h = 1.06 * np.std(y, ddof=0) * n ** (-0.2)
    fhat = np.mean(
        np.exp(-0.5 * ((y[:, None] - q[None, :]) / h) ** 2), axis=0
    ) / (h * np.sqrt(2.0 * np.pi))
    fhat = np.maximum(fhat, 1e-8)
    indicators = (y[:, None] <= q[None, :]).astype(float)
    rif = q[None, :] + (taus[None, :] - indicators) / fhat[None, :]
    return rif, q, fhat


def mixture_design(q: np.ndarray, separation: float = 8.0, p_upper: float = 0.5):
    """CDF derivatives and fitted-mixture density for the oracle nuisance fit."""
    mu_lower, mu_upper = -separation / 2.0, separation / 2.0
    f_lower = np.exp(-0.5 * (q - mu_lower) ** 2) / np.sqrt(2.0 * np.pi)
    f_upper = np.exp(-0.5 * (q - mu_upper) ** 2) / np.sqrt(2.0 * np.pi)
    F_lower, F_upper = ndtr(q - mu_lower), ndtr(q - mu_upper)
    dF = np.column_stack(
        (
            -(1.0 - p_upper) * f_lower,
            -p_upper * f_upper,
            p_upper * (1.0 - p_upper) * (F_upper - F_lower),
        )
    )
    f_model = (1.0 - p_upper) * f_lower + p_upper * f_upper
    return dF, f_model


def run_one(seed: int, n: int, separation: float, gamma: float):
    rng = np.random.default_rng(seed)
    maf, p_upper = 0.3, 0.5
    g = rng.binomial(2, maf, size=n).astype(float)
    # Stage 1 includes an intercept in its covariate projection, so use the
    # realized sample-centered dosage rather than only the population mean.
    x = g - np.mean(g)
    upper = rng.random(n) < (1.0 / (1.0 + np.exp(-gamma * x)))
    y = np.where(upper, separation / 2.0, -separation / 2.0) + rng.normal(size=n)

    rif, q, f_stage1 = production_rif(y, TAUS)
    s2 = float(x @ x)
    beta_rif = (x @ rif) / s2
    a_cdf = (x @ (y[:, None] <= q[None, :])) / s2
    dF, f_model = mixture_design(q, separation, p_upper)

    # Historical construction: beta_RIF = (-dF/f_model) theta.
    w_mismatch = -dF / f_model[:, None]
    # Proposed repair: use exactly Stage 1's fhat in the same mapping.
    w_matched = -dF / f_stage1[:, None]
    theta_mismatch = np.linalg.lstsq(w_mismatch, beta_rif, rcond=1e-12)[0]
    theta_matched = np.linalg.lstsq(w_matched, beta_rif, rcond=1e-12)[0]
    theta_cdf = np.linalg.lstsq(dF, a_cdf, rcond=1e-12)[0]

    # The identity is the diagnostic's key invariant; it should hold up to
    # floating-point error for every replicate.
    identity_error = np.max(np.abs(beta_rif + a_cdf / f_stage1))
    return {
        "mismatch": theta_mismatch,
        "matched": theta_matched,
        "cdf": theta_cdf,
        "identity_error": identity_error,
        "density_ratio_median": float(np.median(f_model / f_stage1)),
        "density_ratio_min": float(np.min(f_model / f_stage1)),
        "density_ratio_max": float(np.max(f_model / f_stage1)),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=60)
    ap.add_argument("--n", type=int, default=10_000)
    ap.add_argument("--seed", type=int, default=20260915)
    ap.add_argument("--separation", type=float, default=8.0)
    ap.add_argument("--gamma", type=float, default=0.08)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    if a.reps < 1 or a.n < 100 or a.separation <= 0:
        ap.error("require reps >= 1, n >= 100 and separation > 0")

    rows = [run_one(a.seed + i, a.n, a.separation, a.gamma) for i in range(a.reps)]
    result = {
        "design": {
            "reps": a.reps,
            "n": a.n,
            "seed": a.seed,
            "separation": a.separation,
            "gamma_truth": a.gamma,
            "maf": 0.3,
            "taus": TAUS.tolist(),
            "parameters": list(PARAMETERS),
            "stage1_density": "1.06*population_sd*n^(-0.2) Gaussian KDE at type-8 empirical quantiles",
            "stage2_density": "oracle two-normal mixture PDF at the same quantiles",
        },
        "mean_estimate": {
            "mismatched": np.mean([r["mismatch"] for r in rows], axis=0).tolist(),
            "matched": np.mean([r["matched"] for r in rows], axis=0).tolist(),
            "cdf": np.mean([r["cdf"] for r in rows], axis=0).tolist(),
        },
        "sd_estimate": {
            "mismatched": np.std([r["mismatch"] for r in rows], axis=0, ddof=1).tolist(),
            "matched": np.std([r["matched"] for r in rows], axis=0, ddof=1).tolist(),
            "cdf": np.std([r["cdf"] for r in rows], axis=0, ddof=1).tolist(),
        },
        "diagnostic": {
            "max_rif_cdf_identity_error": float(max(r["identity_error"] for r in rows)),
            "density_ratio_median_mean": float(np.mean([r["density_ratio_median"] for r in rows])),
            "density_ratio_min_mean": float(np.mean([r["density_ratio_min"] for r in rows])),
            "density_ratio_max_mean": float(np.mean([r["density_ratio_max"] for r in rows])),
            "max_matched_cdf_estimate_error": float(
                max(np.max(np.abs(r["matched"] - r["cdf"])) for r in rows)
            ),
        },
    }
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
