"""
Build script for Stage 1 C++ extension using pybind11.

Usage:
    python setup.py build_ext --inplace
    
Or for development:
    pip install -e .
"""

import os
import shutil
import subprocess
import sys
import tempfile
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

class get_pybind_include:
    """Helper class to determine pybind11 include path."""
    def __str__(self):
        import pybind11
        return pybind11.get_include()

def get_armadillo_include():
    """Find Armadillo include path from conda."""
    conda_prefix = os.environ.get('CONDA_PREFIX', '')
    if conda_prefix:
        arma_path = os.path.join(conda_prefix, 'include')
        if os.path.exists(os.path.join(arma_path, 'armadillo')):
            return arma_path
    # Fallback to common locations
    for path in ['/usr/include', '/usr/local/include']:
        if os.path.exists(os.path.join(path, 'armadillo')):
            return path
    return None

def get_armadillo_lib():
    """Find Armadillo library path from conda."""
    conda_prefix = os.environ.get('CONDA_PREFIX', '')
    if conda_prefix:
        lib_path = os.path.join(conda_prefix, 'lib')
        if os.path.exists(lib_path):
            return lib_path
    return None

# Vectorisation flags for the C++ kernel.
#
# The default is EMPTY, i.e. the x86-64 SSE2 baseline, and that is now a
# deliberate correctness choice rather than a performance one.  Two kernels live
# in this extension and they want different things:
#
#   * The exact kernel (stage1_pybind.cpp) must stay bit-identical to the
#     historical output that `--exact-kernel` exists to reproduce.  Raising the
#     baseline lets GCC contract a*b+c into an FMA, which changes float64 results
#     in the last bits (measured max absolute difference 2.2e-16, max
#     |delta beta|/SE 4.6e-15).  So its ISA is pinned here.
#
#   * The fast kernel (stage1_fast.hpp) needs AVX2 + FMA to be worth anything -
#     measured 4x on its streaming sweep - so it carries its own ISA: the hot
#     leaves are compiled twice, once at this baseline and once under
#     `#pragma GCC target("avx,avx2,fma")`, and dispatched at run time on
#     __builtin_cpu_supports.  That keeps the binary portable across the mix of
#     Skylake-SP, Cascade Lake and AMD nodes a BluePebble `short` allocation can
#     land on, without a -march=native build that could SIGILL on the next one.
#
# What raising the baseline is worth FOR THE EXACT KERNEL, measured rather than
# assumed:
#
# An older revision of this comment said raising the ISA "bought nothing", on a
# microbenchmark of a 500-SNP batch against the then-current kernel.  That measurement no
# longer holds.  Re-measured end to end on BluePebble (2026-09-16, Xeon Gold 5120, one
# chr22 testosterone scan = 17,550 variants, N=399,710, T=45, --blocks 1000), building the
# SAME sources three ways and interleaving the builds in rotating order so that a drift in
# node load cannot land preferentially on one of them:
#
#   kernel seconds, 4 threads (one per physical core), all three builds on one node:
#     --cov-mode full          SSE2 215.4   AVX2 193.0 (1.12x)   AVX2 -ffp-contract=off 198.5 (1.09x)
#     --cov-mode pooled-only   SSE2 160.6   AVX2 141.2 (1.14x)   AVX2 -ffp-contract=off 147.9 (1.09x)
#   wall clock, 4 threads:
#     --cov-mode full          SSE2 258.9   AVX2 235.2 (1.10x)
#     --cov-mode pooled-only   SSE2 187.5   AVX2 167.5 (1.12x)
#   kernel seconds, 8 threads (4 physical cores + SMT siblings), two independent passes:
#     --cov-mode pooled-only   SSE2 128.2/128.1   AVX2 102.5/105.1 (1.23x)
#
# So the gain is real but modest, and it is larger for the pooled-only kernel - which has
# no block accumulator and is therefore compute-bound rather than L3-bound - than for the
# block kernel.  GCC 8 still emits mostly SCALAR FMAs (158 of 167) and no zmm registers at
# all, so the win is the wider ymm loads and the contraction, not vectorised reductions.
#
# The old comment's *other* claim was that FMA contraction changes float64 results in the
# last bits.  That is true of the arithmetic and false of the artefacts: measured on the
# full chr22 scan above, the AVX2 build's .stage1.tsv.gz, .sigma0.tsv.gz, .cov.ids.tsv.gz
# and the whole 72,657,000-byte .cov.gz payload (18,164,250 float32 values) are byte-for-byte
# IDENTICAL to the SSE2 build's once decompressed.  The last-bit float64 differences are
# some nine orders of magnitude below the %.6g text and float32 binary resolution the
# outputs are written at, so they cannot survive serialisation.  -ffp-contract=off removes
# every FMA (objdump: 167 -> 0) while keeping the AVX2 width, and is likewise byte-identical;
# it is not needed to preserve the outputs, and it gives up about a third of the speedup.
#
# The default is nevertheless still EMPTY, i.e. the portable SSE2 baseline, and the reason
# is now PORTABILITY ALONE rather than results or verification.  Raising the baseline is a
# build-time choice but the binary outlives the build host: the extension is built once into
# a package directory on a shared filesystem and then imported by whatever node an
# allocation lands on.  A -march=haswell binary SIGILLs, with no graceful diagnostic, on any
# host without AVX2.  Every current BluePebble node has it, but "every current node" is not
# a property setup.py can check at import time on the machine that will run the code, and a
# 10-14% saving does not pay for a hard crash on the one node that does not.
#
# Turn it on deliberately, per site, when the run hosts are known:
#   FUNGWAS_STAGE1_ARCH="-march=haswell -mtune=skylake-avx512"  # AVX2+FMA; measured above
#   FUNGWAS_STAGE1_ARCH="-march=skylake-avx512"                 # Intel-only; GCC 8 emits no zmm
#   FUNGWAS_STAGE1_ARCH="-march=native"                         # NOT portable, see above
# Setting the variable to the empty string is an explicit request for the baseline.
#
# The fast kernel does not need any of that, because it dispatches its own ISA at
# run time; the paragraphs above are what the override buys the exact kernel, and
# the reason the default stays at the portable baseline for both.  Setting the
# override changes the exact kernel's arithmetic, so a run that must reproduce
# historical output byte for byte should leave it alone.
#
# -ffast-math is deliberately never used.  It is the one flag that WOULD change the
# artefacts, because it licenses reassociating the reductions rather than contracting them.
_DEFAULT_ARCH_FLAGS = ''


def resolve_arch_flags(env=None):
    """Resolve the ISA flags for the kernel build from the environment.

    Unset means the portable baseline.  Set means exactly what the caller asked for,
    split on whitespace, including the empty string as an explicit request for the
    baseline.  Kept as a function so the policy is testable without a compiler.
    """
    if env is None:
        env = os.environ
    arch_env = env.get('FUNGWAS_STAGE1_ARCH')
    return (_DEFAULT_ARCH_FLAGS if arch_env is None else arch_env).split()


_arch_flags = resolve_arch_flags()

ext_modules = [
    Extension(
        'fungwas_stage1._stage1_cpp',
        sources=['src/stage1_pybind.cpp'],
        depends=['src/stage1_fast.hpp', 'src/stage1_fast_hot.inc'],
        include_dirs=[
            get_pybind_include(),
            get_armadillo_include(),
        ],
        library_dirs=[get_armadillo_lib()] if get_armadillo_lib() else [],
        runtime_library_dirs=[get_armadillo_lib()] if get_armadillo_lib() else [],
        libraries=['armadillo'],
        # OpenMP flags are added by BuildExt once a probe compile shows which work.
        extra_compile_args=['-std=c++14', '-O3', '-fPIC'] + _arch_flags,
        language='c++'
    ),
]


def _libomp_prefixes():
    """Where a separately installed libomp may live (macOS: conda, Homebrew)."""
    out = []
    if os.environ.get('CONDA_PREFIX'):
        out.append(os.environ['CONDA_PREFIX'])
    brew = shutil.which('brew')
    if brew:
        try:
            out.append(subprocess.run([brew, '--prefix', 'libomp'], capture_output=True,
                                      text=True, check=True).stdout.strip())
        except (OSError, subprocess.CalledProcessError):
            pass
    out += ['/opt/homebrew/opt/libomp', '/usr/local/opt/libomp']
    return [p for p in dict.fromkeys(out) if p and os.path.exists(os.path.join(p, 'include', 'omp.h'))]


def openmp_candidates(platform=sys.platform):
    """(compile_args, link_args) pairs to try, in order.

    GCC and LLVM clang (including conda-forge's) take -fopenmp. Apple clang does
    not: it needs -Xpreprocessor -fopenmp plus a separately installed libomp.
    """
    cands = [(['-fopenmp'], ['-fopenmp'])]
    if platform == 'darwin':
        for p in _libomp_prefixes():
            inc, lib = os.path.join(p, 'include'), os.path.join(p, 'lib')
            cands.append((['-Xpreprocessor', '-fopenmp', '-I' + inc],
                          ['-L' + lib, '-lomp', '-Wl,-rpath,' + lib]))
    return cands


class BuildExt(build_ext):
    """Custom build_ext command with compiler customization."""

    def _openmp_works(self, cflags, lflags):
        from distutils.errors import CompileError, LinkError
        tmp = tempfile.mkdtemp()
        try:
            src = os.path.join(tmp, 'omp_probe.cpp')
            with open(src, 'w') as fh:
                fh.write('#include <omp.h>\nint main() { return omp_get_max_threads() > 0 ? 0 : 1; }\n')
            objs = self.compiler.compile([src], output_dir=tmp, extra_postargs=cflags)
            self.compiler.link_executable(objs, 'omp_probe', output_dir=tmp, extra_postargs=lflags,
                                          target_lang='c++')
            return True
        except (CompileError, LinkError):
            return False
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    def build_extensions(self):
        # Remove -Wstrict-prototypes (not valid for C++)
        if '-Wstrict-prototypes' in self.compiler.compiler_so:
            self.compiler.compiler_so.remove('-Wstrict-prototypes')
        flags = None
        if os.environ.get('FUNGWAS_STAGE1_OPENMP', '1') != '0':
            flags = next((c for c in openmp_candidates() if self._openmp_works(*c)), None)
        if flags is None:
            sys.stderr.write(
                '\n*** fungwas-stage1: building WITHOUT OpenMP; the C++ kernel will run on one '
                'thread. On macOS: conda install -c conda-forge cxx-compiler llvm-openmp '
                '(or brew install libomp), then reinstall. ***\n\n')
        for ext in self.extensions:
            if flags is not None:
                ext.extra_compile_args = ext.extra_compile_args + flags[0]
                ext.extra_link_args = ext.extra_link_args + flags[1]
        super().build_extensions()

# Guarded so that a test can import this file and check the ISA policy without a
# compiler and without running a build.  Both `python setup.py build_ext --inplace`
# and pip's PEP 517 backend execute setup.py with __name__ == "__main__", so the
# guard changes nothing about how the extension is actually built.
# Package metadata (name, version, dependencies, entry point) lives in
# pyproject.toml only; this file just builds the extension.
if __name__ == '__main__':
    setup(
        ext_modules=ext_modules,
        cmdclass={'build_ext': BuildExt},
    )
