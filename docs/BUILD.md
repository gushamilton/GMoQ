# Building and testing

## Requirements

- **Stage 1:** Python >= 3.10, a C++ compiler with OpenMP, and Armadillo. The conda
  environment (`stage1/environment.yml`) provides all of them, on Linux and macOS.
- **Stage 2:** R >= 4.1 with `Rcpp`, `RcppArmadillo` and `data.table`.

## Stage 1

```bash
conda env create -f stage1/environment.yml
conda activate fungwas-stage1
pip install ./stage1
python -c "from fungwas_stage1 import core; print(core.HAVE_CPP)"   # must print True
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
and Stage 1 warns at run time that `--threads` has no effect. Set `FUNGWAS_STAGE1_OPENMP=0` to
build without it on purpose.

**If `HAVE_CPP` is False,** the import error is printed when Stage 1 starts. Reinstall with
`pip install ./stage1`. After pulling changes to `stage1/src/`, reinstall too: a stale extension
fails with missing-symbol errors.

**Threads.** Set `--threads` to the number of *physical* cores.

## Stage 2

```bash
R CMD INSTALL stage2                    # or R CMD INSTALL -l /path/to/lib stage2
```

On macOS, Homebrew R can ship a `data.table` without zlib, which breaks reading and writing
`.gz` files. If `data.table::fwrite(data.frame(a = 1), tempfile(fileext = ".gz"))` fails,
reinstall it: `install.packages("data.table", type = "source")`.

## Tests

```bash
# Stage 1 (with a plain venv: pip install "./stage1[dev]" first)
python -m pytest stage1/tests -q

# Stage 2
for t in stage2/tests/test-*.R; do Rscript --vanilla "$t" || exit 1; done
Rscript -e 'testthat::test_dir("stage2/tests/testthat", package = "fungwasStage2", load_package = "installed")'
```
