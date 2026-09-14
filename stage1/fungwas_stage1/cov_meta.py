"""Covariance metadata (``<prefix>.cov.meta.json``) written alongside Stage 1 output.

Stage 1 estimates the T x T tau-covariance ``Sigma`` of the quantile-regression
coefficients by a delete-one-block jackknife with ``B`` blocks.  Because
``Sigma`` is *estimated* rather than known, the Stage 2 Wald statistic
``z = theta / SE`` is **t** distributed, not standard normal: for a fixed
contrast ``a``,

    a Sigma_jk a / a Sigma a  ~  chi2_nu / nu

for a Wishart-like ``Sigma_jk`` on ``nu`` degrees of freedom, which makes
``z`` exactly ``t_nu``.  Reading a normal reference instead inflates Type-I
error by 17x at p = 5e-8 when B = 65 (see ``docs/KNOWN_ISSUES.md``).

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
        "nu": 999,
        "nu_source": "blocks_minus_one"
      }
    }

``nu`` is stored explicitly and separately from ``n_blocks`` on purpose.  The
analytic value for a delete-one-block jackknife is ``B - 1`` (the covariance
is a sample covariance of ``B`` leave-one-out estimates, so its rank and its
Wishart degrees of freedom are ``B - 1``), and that is what Stage 1 writes,
tagged ``nu_source = "blocks_minus_one"``.  A separate workstream is measuring
whether the *effective* degrees of freedom are materially lower than ``B - 1``;
when it reports, a better number can be written into the same ``nu`` field with
a different ``nu_source`` and Stage 2 will pick it up with no format change.

Files written before this schema existed have no ``jackknife`` block at all.
They remain readable: ``resolve_nu`` returns ``None`` for them, and Stage 2
refuses to emit a parameter p-value rather than silently falling back to a
normal reference.
"""

from __future__ import annotations

import json
from typing import Any, Optional, Tuple

#: Bumped whenever the on-disk schema gains or changes a field.
COV_META_VERSION = 2

#: ``nu_source`` value for the analytic delete-one-block jackknife df.
NU_SOURCE_BLOCKS_MINUS_ONE = "blocks_minus_one"

JACKKNIFE_METHOD = "delete_one_block"


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
        Degrees of freedom to advertise to Stage 2.  Defaults to ``B - 1``.
    nu_source : str, optional
        How ``nu`` was arrived at.  Defaults to ``"blocks_minus_one"``.
    n_blocks_nonempty : int, optional
        Number of blocks that actually received a sample.  Deleting an empty
        block reproduces the full-sample estimate, so empty blocks contribute
        no replicate and the degrees of freedom follow the populated count.
        Equal to ``n_blocks`` whenever N is comfortably larger than B, which
        is every production run.
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

    if nu is None:
        nu = jackknife_nu(n_blocks_nonempty)
        nu_source = NU_SOURCE_BLOCKS_MINUS_ONE
    elif nu_source is None:
        raise ValueError("nu_source must be given whenever nu is supplied explicitly")

    jackknife: dict[str, Any] = {
        "method": JACKKNIFE_METHOD,
        "n_blocks": n_blocks,
        "n_blocks_nonempty": n_blocks_nonempty,
        "nu": nu,
        "nu_source": nu_source,
    }
    if seed is not None:
        jackknife["seed"] = int(seed)

    return {
        "cov_meta_version": COV_META_VERSION,
        "cov_dtype": cov_dtype,
        "n_taus": n_taus,
        "n_cov_elements": int(n_cov_elements),
        "record_bytes": int(n_cov_elements if cov_dtype == "int8" else n_cov_elements * 4),
        "scale_file": bool(cov_dtype == "int8"),
        "jackknife": jackknife,
    }


def write_cov_meta(path: str, meta: dict) -> None:
    """Write a metadata payload as JSON."""
    with open(path, "w") as handle:
        json.dump(meta, handle, indent=2)


def read_cov_meta(path: str) -> dict:
    """Read a ``.cov.meta.json`` file."""
    with open(path) as handle:
        return json.load(handle)


def resolve_nu(meta: dict) -> Tuple[Optional[float], Optional[str]]:
    """Resolve the degrees of freedom to use for the Stage 2 t reference.

    Returns ``(nu, nu_source)``.  Resolution order:

    1. an explicit ``jackknife.nu`` (whatever wrote it);
    2. otherwise ``jackknife.n_blocks - 1``;
    3. otherwise ``(None, None)`` -- a pre-v2 file that does not record the
       block count, for which the right reference is simply not knowable.
    """
    jk = meta.get("jackknife") if isinstance(meta, dict) else None
    if not isinstance(jk, dict):
        return None, None

    nu = jk.get("nu")
    if isinstance(nu, (int, float)) and not isinstance(nu, bool) and nu > 0:
        source = jk.get("nu_source") or NU_SOURCE_BLOCKS_MINUS_ONE
        return float(nu), str(source)

    n_blocks = jk.get("n_blocks_nonempty", jk.get("n_blocks"))
    if isinstance(n_blocks, int) and not isinstance(n_blocks, bool) and n_blocks >= 2:
        return float(jackknife_nu(n_blocks)), NU_SOURCE_BLOCKS_MINUS_ONE

    return None, None
