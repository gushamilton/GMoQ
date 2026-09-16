## in-place patch: teach draw_eps the standardised 2-normal mixture law
import io, re, sys
p = "/user/work/fh6520/gmoq_covmixW/src/sim_stage1_mix.py"
s = io.open(p).read()
old = '''def draw_eps(rng, n, shape):
    if shape == "normal":
        return rng.normal(size=n)
    k = GAMMA_SHAPE[shape]'''
new = '''def draw_eps(rng, n, shape):
    if shape == "normal":
        return rng.normal(size=n)
    if shape.startswith("mix"):
        # standardised 50/50 two-normal mixture, component separation S in
        # units of the component sd:  X ~ .5 N(-S/2,1) + .5 N(S/2,1),
        # Var(X) = 1 + S^2/4, eps = X / sqrt(1 + S^2/4).
        # This is exactly eps_cdf/eps_pdf in src/truth_mix.R.
        S = float(shape[3:])
        kap = np.sqrt(1.0 + S * S / 4.0)
        z = rng.normal(size=n)
        sgn = rng.integers(0, 2, size=n) * 2.0 - 1.0
        return (z + sgn * S / 2.0) / kap
    k = GAMMA_SHAPE[shape]'''
assert old in s, "anchor not found"
s = s.replace(old, new)
io.open(p, "w").write(s)
print("patched", p)
