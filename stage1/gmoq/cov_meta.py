"""Covariance metadata (``<prefix>.cov.meta.json``) written alongside Stage 1 output.

Stage 1 estimates the T x T tau-covariance ``Sigma`` of the quantile-regression
coefficients by a delete-one-block jackknife with ``B`` blocks.  Because
``Sigma`` is *estimated* rather than known, the Stage 2 Wald statistic
``z = theta / SE`` is **t** distributed, not standard normal: for a fixed
contrast ``a``,

    a Sigma_jk a / a Sigma a  ~  chi2_nu / nu

for a Wishart-like ``Sigma_jk`` on ``nu`` degrees of freedom, which makes
``z`` exactly ``t_nu``.  Reading a normal reference instead inflates the
Type-I error rate in the far tail when the number of blocks is small (see
``docs/METHOD.md``).

Stage 2 cannot apply the correct reference unless Stage 1 tells it ``nu``.
That is what this metadata carries.

Schema (version 2)::

    {
      "cov_meta_version": 2,
      "cov_dtype": "int8",
      "n_taus": 45,
      "n_cov_elements": 1035,
      "record_bytes": 1035,
      "scale_file": true,
      "jackknife": {
        "method": "delete_one_block",
        "n_blocks": 1000,
        "n_blocks_nonempty": 1000,
        "seed": 42,
        "df_model": {
          "form": "1/nu = 1/(n_blocks_nonempty - 1) + c / n_eff",
          "n_eff_column": "d_resid",
          "c_location": 0.88,
          "c_scale": 2.62,
          "fallback_n_eff": 9134,
          "source": "simulation calibration (docs/METHOD.md)"
        }
      }
    }

**There is deliberately no single ``nu`` field.**  In simulation the effective
degrees of freedom fall well below ``B - 1``, vary by SNP with allele
frequency, and differ by fitted parameter (lower for a scale-type contrast than
for a location-type one), so no single number is right.

What is stable is the law, fitted to simulations::

    1/nu = 1/(B - 1) + c / N_eff        N_eff = D^2 / sum_i g_ri^4

``c = 0.88`` for a location-type contrast and ``2.62`` for a scale-type one.
``N_eff`` collapses to ``D = sum_i g_ri^2`` because a genotype's fourth central
moment equals its second, and ``D`` is the denominator Stage 1 already forms
for every SNP.  So the whole model costs one float per SNP, emitted as the
``d_resid`` column of the TSV, plus the constants recorded here.

``fallback_n_eff`` is ``D`` for a variant of MAF 0.01 in a biobank-scale sample
(about 460,000 people), used only when ``d_resid`` is unavailable; it is
conservative for commoner variants.

Files written before the ``jackknife`` block existed have none of this.  They
remain readable: ``resolve_nu`` returns ``None``, and Stage 2 refuses to emit a
parameter p-value rather than silently falling back to a normal reference.
"""

from __future__ import annotations

import json
from typing import Any, Optional, Tuple

import numpy as np

#: Bumped whenever the on-disk schema gains or changes a field.
#:
#: Version 4 adds the ``analysis`` block.  It records what the scan was: the
#: phenotype, the covariate columns, the analysed sample count, the LOCO
#: predictor, and a fingerprint over all of those plus the tau grid.  None of
#: it is needed to read the jackknife covariance; it exists so that Stage 2 can
#: check a *pooled* covariance sidecar (``.sigma0.meta.json``) against the scan
#: it is being applied to and refuse on any disagreement.  ``Sigma_0`` built on
#: a different covariate set is wrong in a way that leaves the standard errors
#: plausible and the point estimates untouched, so nothing downstream can
#: detect it.
COV_META_VERSION = 4

#: ``nu_source`` value for the analytic delete-one-block jackknife df.
NU_SOURCE_BLOCKS_MINUS_ONE = "blocks_minus_one"

JACKKNIFE_METHOD = "delete_one_block"

#: Stage 1 TSV column holding ``D = sum_i g_resid^2``, which stands in for N_eff.
N_EFF_COLUMN = "d_resid"

#: Satterthwaite coefficient for a flat (location-type) contrast over taus.
C_LOCATION = 0.88

#: Satterthwaite coefficient for a tail-loaded (scale-type) contrast over taus.
C_SCALE = 2.62

#: N_eff to assume when ``d_resid`` is unavailable: ``D`` at MAF 0.01 in a
#: biobank-scale sample, so conservative at every commoner variant.
FALLBACK_N_EFF = 9134

DF_MODEL_FORM = "1/nu = 1/(n_blocks_nonempty - 1) + c / n_eff"

DF_MODEL_SOURCE = "simulation calibration (docs/METHOD.md)"


def df_model() -> dict:
    """The degrees-of-freedom law Stage 2 should apply, as recorded on disk."""
    return {
        "form": DF_MODEL_FORM,
        "n_eff_column": N_EFF_COLUMN,
        "c_location": C_LOCATION,
        "c_scale": C_SCALE,
        "fallback_n_eff": FALLBACK_N_EFF,
        "source": DF_MODEL_SOURCE,
    }


def effective_nu(n_blocks: int, n_eff: float, c: float) -> float:
    """Effective degrees of freedom of the jackknife covariance for one contrast.

    ``1/nu = 1/(B - 1) + c/N_eff``.  The first term is the nominal jackknife df;
    the second is the Satterthwaite penalty for unequal block contributions,
    which is what makes ``nu`` fall short of ``B - 1`` and makes it depend on
    the SNP through ``N_eff``.
    """
    inv = 1.0 / (int(n_blocks) - 1) + float(c) / float(n_eff)
    if not np.isfinite(inv) or inv <= 0:
        raise ValueError(f"non-positive 1/nu from n_blocks={n_blocks}, n_eff={n_eff}, c={c}")
    return 1.0 / inv


def jackknife_nu(n_blocks: int) -> int:
    """Degrees of freedom of a delete-one-block jackknife covariance.

    The estimator is ``Sigma = ((B-1)/B) * sum_k (theta_k - theta_bar)^2``
    over the ``B`` leave-one-block-out estimates, i.e. a sample covariance of
    ``B`` centred vectors, so ``nu = B - 1``.
    """
    n_blocks = int(n_blocks)
    if n_blocks < 2:
        raise ValueError(f"n_blocks must be at least 2 to have any df, got {n_blocks}")
    return n_blocks - 1


def build_cov_meta(
    cov_dtype: str,
    n_taus: int,
    n_blocks: int,
    seed: Optional[int] = None,
    nu: Optional[float] = None,
    nu_source: Optional[str] = None,
    n_blocks_nonempty: Optional[int] = None,
    analysis: Optional[dict] = None,
    se_source: Optional[str] = None,
) -> dict:
    """Build the ``.cov.meta.json`` payload.

    Parameters
    ----------
    cov_dtype : {"int8", "float32"}
        Storage dtype of the covariance payload.
    n_taus : int
        Number of tau levels (T).
    n_blocks : int
        Number of jackknife blocks actually used (B).
    seed : int, optional
        Seed used to assign samples to blocks; recorded for reproducibility.
    nu : float, optional
        A constant ``nu`` to advertise to Stage 2.  Normally omitted: a single
        ``nu`` is not defensible, because the effective degrees of freedom vary
        by SNP and by fitted parameter.  Stage 2 computes them from
        ``df_model`` instead.  Kept only so a caller can pin a value
        deliberately, and it must then say where the value came from.
    nu_source : str, optional
        Required whenever ``nu`` is supplied.
    n_blocks_nonempty : int, optional
        Number of blocks that actually received a sample.  Deleting an empty
        block reproduces the full-sample estimate, so empty blocks contribute
        no replicate and the degrees of freedom follow the populated count.
        Equal to ``n_blocks`` whenever N is comfortably larger than B, which
        is every production run.
    analysis : dict, optional
        Provenance of the scan, from ``pooled_cov.sigma0_analysis_block``.
        Recorded verbatim so Stage 2 can compare it with the matching field
        of a ``.sigma0.meta.json`` sidecar and refuse a mismatch.
    se_source : str, optional
        What the ``se_tau_*`` columns of the Stage 1 TSV actually contain, when
        it is not the jackknife.  ``"pooled_sigma0"`` means
        ``sqrt(Sigma_0[t,t] / d_resid)``, written by ``--cov-mode pooled-only``,
        which does not accumulate the per-block statistics a jackknife standard
        error is a function of.  Omitted entirely when the jackknife ran, so
        every file written before this field existed is unchanged and a reader
        may treat "absent" as "jackknife".
    """
    if cov_dtype not in ("int8", "float32"):
        raise ValueError(f"cov_dtype must be int8 or float32, got {cov_dtype!r}")
    n_taus = int(n_taus)
    n_blocks = int(n_blocks)
    n_cov_elements = n_taus * (n_taus + 1) // 2

    if n_blocks_nonempty is None:
        n_blocks_nonempty = n_blocks
    n_blocks_nonempty = int(n_blocks_nonempty)
    if n_blocks_nonempty > n_blocks:
        raise ValueError(
            f"n_blocks_nonempty ({n_blocks_nonempty}) cannot exceed n_blocks ({n_blocks})"
        )

    jackknife: dict[str, Any] = {
        "method": JACKKNIFE_METHOD,
        "n_blocks": n_blocks,
        "n_blocks_nonempty": n_blocks_nonempty,
        "df_model": df_model(),
    }
    if nu is not None:
        if nu_source is None:
            raise ValueError("nu_source must be given whenever nu is supplied explicitly")
        jackknife["nu"] = nu
        jackknife["nu_source"] = nu_source
    if seed is not None:
        jackknife["seed"] = int(seed)

    payload = {
        "cov_meta_version": COV_META_VERSION,
        "cov_dtype": cov_dtype,
        "n_taus": n_taus,
        "n_cov_elements": int(n_cov_elements),
        "record_bytes": int(n_cov_elements if cov_dtype == "int8" else n_cov_elements * 4),
        "scale_file": bool(cov_dtype == "int8"),
        "jackknife": jackknife,
    }
    if analysis is not None:
        payload["analysis"] = dict(analysis)
    if se_source is not None:
        payload["se_source"] = str(se_source)
    return payload


def write_cov_meta(path: str, meta: dict) -> None:
    """Write a metadata payload as JSON."""
    with open(path, "w") as handle:
        json.dump(meta, handle, indent=2)


def read_cov_meta(path: str) -> dict:
    """Read a ``.cov.meta.json`` file."""
    with open(path) as handle:
        return json.load(handle)


def _jackknife(meta: dict) -> Optional[dict]:
    jk = meta.get("jackknife") if isinstance(meta, dict) else None
    return jk if isinstance(jk, dict) else None


def resolve_nu(meta: dict) -> Tuple[Optional[float], Optional[str]]:
    """Read a constant ``nu`` if one was deliberately recorded.

    Returns ``(None, None)`` otherwise -- including when the block count is
    known.  There is no ``B - 1`` fallback here on purpose: the effective
    degrees of freedom are measurably below ``B - 1`` and vary by SNP and by
    fitted parameter, so deriving a constant from the block count would be
    manufacturing a number that measurement contradicts.  Use
    ``resolve_block_count`` plus ``resolve_df_model`` instead.
    """
    jk = _jackknife(meta)
    if jk is None:
        return None, None
    nu = jk.get("nu")
    if isinstance(nu, (int, float)) and not isinstance(nu, bool) and nu > 0:
        source = jk.get("nu_source") or "unrecorded"
        return float(nu), str(source)
    return None, None


def resolve_block_count(meta: dict) -> Optional[int]:
    """Number of jackknife blocks that carried a replicate.

    Prefers ``n_blocks_nonempty`` over the requested ``n_blocks``: deleting an
    empty block reproduces the full-sample estimate and contributes nothing.
    """
    jk = _jackknife(meta)
    if jk is None:
        return None
    for key in ("n_blocks_nonempty", "n_blocks"):
        value = jk.get(key)
        if isinstance(value, int) and not isinstance(value, bool) and value >= 2:
            return int(value)
    return None


def resolve_analysis(meta: dict) -> Optional[dict]:
    """The ``analysis`` provenance block, if the writer supplied one."""
    block = meta.get("analysis") if isinstance(meta, dict) else None
    return block if isinstance(block, dict) else None


def resolve_df_model(meta: dict) -> Optional[dict]:
    """The recorded degrees-of-freedom law, if the writer supplied one."""
    jk = _jackknife(meta)
    if jk is None:
        return None
    model = jk.get("df_model")
    return model if isinstance(model, dict) else None
