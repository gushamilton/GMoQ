## Derive summarise_mix.R from the PUBLISHED sim/covsource summariser.
## Only the paths and the power block change; every rate, its Monte Carlo
## error, the paired differences and the 5e-8 tail extrapolation are the
## published code verbatim.
import io
src = "/user/work/fh6520/covsource_run/src/summarise_cs.R"
dst = "/user/work/fh6520/gmoq_covmixW/src/summarise_mix.R"
s = io.open(src).read()

s = s.replace('WS <- "/user/work/fh6520/covsource_run"',
              'WS <- "/user/work/fh6520/gmoq_covmixW"')
assert 'gmoq_covmixW' in s

# record the column count of each W on every row
old = 'fit = fn, est = est, cov_source = cs, param = p, nrep = n, nfail = nfail,'
new = 'fit = fn, K = ncol(x$W[[fn]]), est = est, cov_source = cs, param = p, nrep = n, nfail = nfail,'
assert s.count(old) == 1
s = s.replace(old, new)

# power at matched ACTUAL Type-I: generalise from the hardcoded beta_sigma2 to
# whichever parameter this W carries.
old = '''  if (rr$dgp != "disp" || rr$param != "beta_sigma2") next
  nc <- nullmap[nullmap$maf == rr$maf & nullmap$shape == rr$shape &
                nullmap$fit == rr$fit & nullmap$est == rr$est &
                nullmap$cov_source == rr$cov_source & nullmap$param == "beta_sigma2", ]
  if (!nrow(nc)) next
  Tn <- stats_store[[paste(nc$cell[1], paste(rr$fit, rr$est, rr$cov_source, sep = "|"), "beta_sigma2", sep = "@@")]]$TT
  Td <- stats_store[[paste(rr$cell, paste(rr$fit, rr$est, rr$cov_source, sep = "|"), "beta_sigma2", sep = "@@")]]$TT'''
new = '''  if (rr$dgp == "null") next
  nc <- nullmap[nullmap$maf == rr$maf & nullmap$shape == rr$shape &
                nullmap$fit == rr$fit & nullmap$est == rr$est &
                nullmap$cov_source == rr$cov_source & nullmap$param == rr$param, ]
  if (!nrow(nc)) next
  Tn <- stats_store[[paste(nc$cell[1], paste(rr$fit, rr$est, rr$cov_source, sep = "|"), rr$param, sep = "@@")]]$TT
  Td <- stats_store[[paste(rr$cell, paste(rr$fit, rr$est, rr$cov_source, sep = "|"), rr$param, sep = "@@")]]$TT
  if (is.null(Tn) || is.null(Td)) next'''
assert s.count(old) == 1
s = s.replace(old, new)

old = '''      cell = rr$cell, maf = rr$maf, shape = rr$shape, Ztarget = rr$Ztarget,
      fit = rr$fit, est = rr$est, cov_source = rr$cov_source, alpha_actual = al,'''
new = '''      cell = rr$cell, dgp = rr$dgp, maf = rr$maf, shape = rr$shape, Ztarget = rr$Ztarget,
      fit = rr$fit, K = rr$K, param = rr$param, est = rr$est,
      cov_source = rr$cov_source, alpha_actual = al,'''
assert s.count(old) == 1
s = s.replace(old, new)

s = s.replace('"summary_cs.tsv"', '"summary_mix.tsv"')
s = s.replace('"paired_cs.tsv"', '"paired_mix.tsv"')
s = s.replace('"power_cs.tsv"', '"power_mix.tsv"')
io.open(dst, "w").write(s)
print("wrote", dst)
