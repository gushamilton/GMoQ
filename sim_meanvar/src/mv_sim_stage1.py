"""Stage 1 driver for the mean-variance (dispersion) Monte Carlo.

One replicate = ONE phenotype vector + M = 200 independent variants, run through
the SHIPPED Stage 1 (core.run_stage1 / core.compute_rif /
pooled_cov.compute_sigma0).

It also computes the EXACT population truth b_l(tau) for every effect-carrying
variant, so Stage 2's projected theta can be compared with a known population
value.

Population truth
----------------
Latent scale:  X = 0.5 + sum_l delta_l g_l + 0.5 * exp(sum_l gamma_l g_l) * e.
Stage 1's estimand for variant l is the population OLS slope of
RIF(Y;tau) = q_tau + (tau - 1{Y<=q_tau})/f_Y(q_tau) on g_l, i.e.

    b_l(tau) = -Cov_g( F_Y(q_tau | g_l = g), g ) / ( f_Y(q_tau) * Var(g_l) ).

Every transform used (identity, exp, rank-inverse-normal) is MONOTONE, so
F_Y(q_tau | g) = F_X(x_tau | g) with x_tau the latent tau-quantile; only the
marginal density f_Y(q_tau) changes.  F_X(x | g_l = g) is the average of
Phi((x - 0.5 - delta_l g - S_-l) / (0.5 exp(gamma_l g + T_-l))) over a large
COMMON-RANDOM-NUMBER background sample (S_-l, T_-l), so the contrast across g --
all that b_l(tau) depends on -- is estimated very precisely.

Population Sigma_0 is exact: RIF depends on Y only through 1{Y<=q_tau}, so
Sigma_0[i,j] = (min(t_i,t_j) - t_i t_j) / (f_i f_j), matching Stage 1's
(1/n) Psi' Psi.
"""
import argparse, json, os, sys, time
import numpy as np
from scipy.special import ndtr, ndtri

STAGE1_SRC = "/user/work/fh6520/wt_pooled-guards/stage1"
sys.path.insert(0, STAGE1_SRC)
from fungwas_stage1 import core, pooled_cov            # noqa: E402

TAUS = np.linspace(0.02, 0.98, 45)
T = len(TAUS)
M = 200
N_NULL, N_MEAN, N_DISP = 170, 20, 10
MAFSET = np.array([0.05, 0.10, 0.20, 0.30, 0.50])

# 10 levels x 2 variants (one +, one -) = 20 pure-MEAN variants.
# (The brief said "5,10,15,...,40, two variants per level", which is 16 variants,
#  not the 20 the class counts require; seq(4,40,4) is the minimal change that
#  satisfies both the class count and "two variants per level".)
MEAN_Z = np.repeat(np.arange(4.0, 40.1, 4.0), 2) * np.tile([1.0, -1.0], 10)
DISP_Z = np.repeat(np.array([3.0, 5.0, 8.0, 12.0, 16.0]), 2) * np.tile([1.0, -1.0], 5)

CLS = np.array([0] * N_NULL + [1] * N_MEAN + [2] * N_DISP)   # 0 null 1 mean 2 disp
SIG_E = 0.5            # latent residual SD at g = 0
GAMMA_REF = 0.005      # reference gamma for the linear calibration of beta_sigma2


class Truth:
    """Population machinery for one replicate's variant configuration."""

    def __init__(self, maf, delta, gamma, nbg, rng):
        self.maf, self.delta, self.gamma = maf, delta, gamma
        self.eff = np.where((delta != 0) | (gamma != 0))[0]
        gb = rng.binomial(2, maf[self.eff], size=(nbg, len(self.eff))).astype(np.float64)
        self.gb = gb
        self.S = gb @ delta[self.eff]
        self.Tt = gb @ gamma[self.eff]
        self.s = SIG_E * np.exp(self.Tt)
        self.mu = 0.5 + self.S
        self.nbg = nbg

    def F(self, x):
        return ndtr((np.asarray(x, dtype=float)[..., None] - self.mu) / self.s).mean(-1)

    def f(self, x):
        z = (np.asarray(x, dtype=float)[..., None] - self.mu) / self.s
        return (np.exp(-0.5 * z * z) / (np.sqrt(2 * np.pi) * self.s)).mean(-1)

    def quantiles(self, taus):
        lo = (self.mu - 10 * self.s).min()
        hi = (self.mu + 10 * self.s).max()
        grid = np.linspace(lo, hi, 601)
        x = np.interp(taus, self.F(grid), grid)
        for _ in range(6):
            x = x - (self.F(x) - taus) / np.maximum(self.f(x), 1e-12)
        return x

    def moments(self, scale):
        if scale == "log":
            return float(self.mu.mean()), float(np.sqrt(self.mu.var() + (self.s ** 2).mean()))
        if scale == "raw":
            ey = np.exp(self.mu + 0.5 * self.s ** 2)
            ey2 = np.exp(2 * self.mu + 2 * self.s ** 2)
            return float(ey.mean()), float(np.sqrt(ey2.mean() - ey.mean() ** 2))
        return 0.0, 1.0

    def b_curves(self, x_tau, f_y, which=None):
        """b_l(tau); rows follow `which` (indices into self.eff) or all of eff."""
        idx = np.arange(len(self.eff)) if which is None else np.asarray(which)
        out = np.zeros((len(idx), len(x_tau)))
        xt = np.asarray(x_tau, dtype=float)[:, None]
        for r, j in enumerate(idx):
            l = self.eff[j]
            p = self.maf[l]
            pg = np.array([(1 - p) ** 2, 2 * p * (1 - p), p ** 2])
            vg = 2 * p * (1 - p)
            gl = self.gb[:, j]
            S_ml = self.S - self.delta[l] * gl
            T_ml = self.Tt - self.gamma[l] * gl
            cov = np.zeros(len(x_tau))
            for k, g in enumerate((0.0, 1.0, 2.0)):
                mu = 0.5 + self.delta[l] * g + S_ml
                sd = SIG_E * np.exp(self.gamma[l] * g + T_ml)
                cov += pg[k] * (g - 2 * p) * ndtr((xt - mu) / sd).mean(1)
            out[r] = -cov / (f_y * vg)
        return out


def scale_density(scale, x_tau, f_x, taus):
    if scale == "log":
        return f_x, x_tau
    if scale == "raw":
        return f_x / np.exp(x_tau), np.exp(x_tau)
    if scale == "rint":
        z = ndtri(taus)
        return np.exp(-0.5 * z * z) / np.sqrt(2 * np.pi), z
    raise ValueError(scale)


def transform(scale, X):
    if scale == "log":
        return X
    if scale == "raw":
        return np.exp(X)
    if scale == "rint":
        n = len(X)
        r = np.empty(n)
        r[np.argsort(X, kind="stable")] = np.arange(1, n + 1)
        return ndtri((r - 0.5) / n)
    raise ValueError(scale)


def sigma0_pop(taus, f_y):
    return (np.minimum.outer(taus, taus) - np.outer(taus, taus)) / np.outer(f_y, f_y)


def make_W(taus, sd):
    return np.column_stack([np.ones(len(taus)), ndtri(taus) / (2.0 * sd)])


def one_rep(scale, N, seed, nbg):
    rng = np.random.default_rng(seed)
    maf = rng.choice(MAFSET, size=M)
    vg = 2 * maf * (1 - maf)

    # --- mean effects: fixed point so the realised |Z_mu| tracks the target ---
    sd_tot = SIG_E
    delta = np.zeros(M)
    sl = slice(N_NULL, N_NULL + N_MEAN)
    for _ in range(50):
        delta[sl] = MEAN_Z * sd_tot / np.sqrt(N * vg[sl])
        new = float(np.sqrt(SIG_E ** 2 + np.sum(delta ** 2 * vg)))
        if abs(new - sd_tot) < 1e-13:
            sd_tot = new
            break
        sd_tot = new

    # --- gamma calibration at a reference gamma ------------------------------
    gamma = np.zeros(M)
    gprobe = np.zeros(M)
    gprobe[N_NULL + N_MEAN:] = GAMMA_REF
    rp = np.random.default_rng(seed + 11)
    trp = Truth(maf, delta, gprobe, nbg, rp)
    x_p = trp.quantiles(TAUS)
    f_yp, _ = scale_density(scale, x_p, trp.f(x_p), TAUS)
    mY0, sY0 = trp.moments(scale)
    W0 = make_W(TAUS, sY0)
    A0 = np.linalg.solve(W0.T @ W0, W0.T)
    S0p0 = sigma0_pop(TAUS, f_yp)
    se_unit = np.sqrt(np.diag(A0 @ S0p0 @ A0.T))           # SE * sqrt(d_resid)
    jdisp = np.array([int(np.where(trp.eff == l)[0][0])
                      for l in range(N_NULL + N_MEAN, M)])
    bp = trp.b_curves(x_p, f_yp, which=jdisp)
    thp = np.linalg.solve(W0.T @ W0, W0.T @ bp.T)
    slope = thp[1] / GAMMA_REF
    se_sig2 = se_unit[1] / np.sqrt(N * vg[N_NULL + N_MEAN:])
    gamma[N_NULL + N_MEAN:] = np.clip(DISP_Z * se_sig2 / slope, -0.5, 0.5)

    # --- exact truth with the final effects ----------------------------------
    tr = Truth(maf, delta, gamma, nbg, np.random.default_rng(seed + 11))
    x_tau = tr.quantiles(TAUS)
    f_y, q_y = scale_density(scale, x_tau, tr.f(x_tau), TAUS)
    b_true = np.zeros((M, T))
    b_true[tr.eff] = tr.b_curves(x_tau, f_y)
    S0_pop = sigma0_pop(TAUS, f_y)
    mY_p, sY_p = tr.moments(scale)

    # --- simulate and run the SHIPPED Stage 1 --------------------------------
    G = rng.binomial(2, np.broadcast_to(maf, (N, M))).astype(np.float64)
    e = rng.normal(size=N)
    X = 0.5 + G @ delta + SIG_E * np.exp(G @ gamma) * e
    y = transform(scale, X)

    res = core.run_stage1(G, y, TAUS, n_blocks=1000)
    rif, q_tau_hat, fhat = core.compute_rif(y, TAUS)
    Qb, _ = np.linalg.qr(np.ones((N, 1)))
    S0 = pooled_cov.compute_sigma0(rif - Qb @ (Qb.T @ rif))

    return dict(betas=res["betas"], d_resid=res["d_resid"], sigma0=S0,
                q_tau=q_tau_hat, fhat=fhat, mean_y=float(y.mean()),
                sd_y=float(y.std(ddof=1)), maf=maf, delta=delta, gamma=gamma,
                b_true=b_true, sigma0_pop=S0_pop, f_pop=f_y, q_pop=q_y,
                mean_y_pop=mY_p, sd_y_pop=sY_p, sd_latent=float(sd_tot))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scale", required=True, choices=["raw", "log", "rint"])
    ap.add_argument("--N", type=int, default=50000)
    ap.add_argument("--chunk", type=int, required=True)
    ap.add_argument("--reps", type=int, required=True)
    ap.add_argument("--nbg", type=int, default=200000)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    tag = f"{args.scale}_chunk{args.chunk:04d}"
    if os.path.exists(os.path.join(args.out, tag + ".json")):
        print(tag, "already done", flush=True)
        return

    R = args.reps
    acc = {k: [] for k in ("beta", "dresid", "sigma0", "qtau", "fhat", "fpop",
                           "btrue", "sigma0pop", "maf", "delta", "gamma", "mom")}
    t0 = time.time()
    for r in range(R):
        seed = int(1_000_000 * (["raw", "log", "rint"].index(args.scale) + 1)
                   + args.chunk * 10000 + r)
        o = one_rep(args.scale, args.N, seed, args.nbg)
        acc["beta"].append(o["betas"]); acc["dresid"].append(o["d_resid"])
        acc["sigma0"].append(o["sigma0"]); acc["qtau"].append(o["q_tau"])
        acc["fhat"].append(o["fhat"]); acc["fpop"].append(o["f_pop"])
        acc["btrue"].append(o["b_true"]); acc["sigma0pop"].append(o["sigma0_pop"])
        acc["maf"].append(o["maf"]); acc["delta"].append(o["delta"])
        acc["gamma"].append(o["gamma"])
        acc["mom"].append([o["mean_y"], o["sd_y"], o["mean_y_pop"],
                           o["sd_y_pop"], o["sd_latent"]])
        print(f"  rep {r+1}/{R} {time.time()-t0:.1f}s", flush=True)

    for k, v in acc.items():
        np.asarray(v, dtype="<f8").tofile(os.path.join(args.out, f"{tag}.{k}.bin"))
    with open(os.path.join(args.out, tag + ".json"), "w") as fh:
        json.dump(dict(scale=args.scale, N=args.N, chunk=args.chunk, reps=R,
                       M=M, T=T, nbg=args.nbg, cls=CLS.tolist(),
                       secs=time.time() - t0), fh)
    print(f"{tag}: {R} reps in {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
