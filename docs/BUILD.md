# Building and testing

## Requirements

- **Stage 1:** Python >= 3.10, a C++ compiler with OpenMP, Armadillo and zlib (headers
  included: `zlib1g-dev` on Debian/Ubuntu). libdeflate is used if present, for faster BGEN
  decoding. The conda environment (`stage1/environment.yml`) provides all of them, on Linux and
  macOS.
- **Stage 2:** R >= 4.1 with `Rcpp`, `RcppArmadillo` and `data.table`.

## Stage 1 (Python package `gmoq`, command `gmoq-stage1`)

```bash
conda env create -f stage1/environment.yml
conda activate gmoq
pip install ./stage1
python -c "from gmoq import core; print(core.HAVE_CPP)"   # must print True
```

**Without conda:**

```bash
brew install armadillo libomp          # macOS
# sudo apt install libarmadillo-dev    # Debian/Ubuntu
python3 -m venv venv && . venv/bin/activate
pip install ./stage1
```

If Armadillo is somewhere else, set `ARMADILLO_DIR` to its prefix (the folder holding
`include/armadillo` and `lib/`). If it is missing, the build stops and says how to install it.

**OpenMP.** The build uses the first OpenMP setup that works: `-fopenmp` (GCC, conda's clang),
or Apple clang with a separately installed libomp. If none works it still builds, single-threaded,
and Stage 1 warns at run time that `--threads` has no effect. Set `GMOQ_STAGE1_OPENMP=0` to
build without it on purpose.

**Instruction set.** The default build targets the portable baseline (x86-64 SSE2 or arm64), so
one build runs on any machine of that architecture; the fast kernel selects AVX2 code paths at
run time where available. `GMOQ_STAGE1_ARCH` adds compiler flags for the exact kernel only, for
example `GMOQ_STAGE1_ARCH="-march=native" pip install ./stage1`.

**If `HAVE_CPP` is False,** the import error is printed when Stage 1 starts. Reinstall with
`pip install ./stage1`. After changing `stage1/src/`, reinstall too: a stale extension fails with
missing-symbol errors.

The command `gmoq-stage1` and `python -m gmoq.cli` are equivalent.

## Stage 2 (R package `gmoq`)

```bash
R CMD INSTALL stage2                    # or R CMD INSTALL -l /path/to/lib stage2
```

On macOS, Homebrew R can ship a `data.table` without zlib, which breaks reading and writing
`.gz` files. If `data.table::fwrite(data.frame(a = 1), tempfile(fileext = ".gz"))` fails,
reinstall it: `install.packages("data.table", type = "source")`.

## Tests

Each package has one compact test suite. Neither needs data files.

```bash
# Stage 1: run from outside stage1/ so the installed package is imported
# (with a plain venv, pip install "./stage1[dev]" first)
python -m pytest stage1/tests -q

# Stage 2
Rscript -e 'testthat::test_dir("stage2/tests/testthat", package = "gmoq", load_package = "installed")'
R CMD build stage2 && R CMD check --no-manual gmoq_*.tar.gz
```

The end-to-end check is the 1000 Genomes example ([examples/1kg_chr22/](examples/1kg_chr22/README.md)).
It downloads about 200 MB and needs `plink2`, so continuous integration runs it only on request.
