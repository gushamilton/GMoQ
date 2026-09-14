"""Cell definitions: phenotype shape x discretisation."""
import numpy as np
from shapes import Base, Dmap

TAUS = np.array([0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50,
                 0.60, 0.70, 0.80, 0.85, 0.90, 0.95, 0.98])

# floor thresholds are the latent standard-normal quantiles giving the stated
# censored fraction
from scipy import stats as sps

CELLS = {
    # ---- smooth controls
    "normal_smooth":   dict(base="normal",   D=("identity", {})),
    "lognorm_smooth":  dict(base="lognormal", D=("identity", {})),
    "bimodal_smooth":  dict(base="bimodal",  D=("identity", {})),
    "bimodal_uneven":  dict(base="bimodal_uneven", D=("identity", {})),
    # ---- rounding / ties  (step in units of sd(Y))
    "round_0p10":      dict(base="normal", D=("round", dict(step=0.10))),
    "round_0p25":      dict(base="normal", D=("round", dict(step=0.25))),
    "round_0p50":      dict(base="normal", D=("round", dict(step=0.50))),
    "round_1p00":      dict(base="normal", D=("round", dict(step=1.00))),
    # ---- floor censoring
    "floor_05":        dict(base="normal", D=("floor", dict(at=float(sps.norm.ppf(0.05))))),
    "floor_20":        dict(base="normal", D=("floor", dict(at=float(sps.norm.ppf(0.20))))),
    "floor_40":        dict(base="normal", D=("floor", dict(at=float(sps.norm.ppf(0.40))))),
    "ceil_20":         dict(base="normal", D=("ceiling", dict(at=float(sps.norm.ppf(0.80))))),
    # ---- count-like: few distinct values
    "count_8":         dict(base="normal", D=("clip_round", dict(step=0.60, lo=-1.8, hi=1.8))),
    "count_4":         dict(base="normal", D=("clip_round", dict(step=1.00, lo=-1.0, hi=2.0))),
    # ---- lumpy + bimodal together
    "bimodal_round":   dict(base="bimodal", D=("round", dict(step=0.50))),
}


def build(name):
    c = CELLS[name]
    b = Base(c["base"])
    D = Dmap(c["D"][0], **c["D"][1])
    return b, D
