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

def _brew_prefix(formula):
    """Homebrew prefix of a formula, or None."""
    brew = shutil.which('brew')
    if not brew:
        return None
    try:
        return subprocess.run([brew, '--prefix', formula], capture_output=True,
                              text=True, check=True).stdout.strip() or None
    except (OSError, subprocess.CalledProcessError):
        return None


def find_armadillo():
    """(include_dir, lib_dir) of an Armadillo install, or None.

    Looks in $ARMADILLO_DIR, the active conda environment, Homebrew and the system
    prefixes, in that order.
    """
    prefixes = [os.environ.get('ARMADILLO_DIR'), os.environ.get('CONDA_PREFIX'),
                _brew_prefix('armadillo'), '/opt/homebrew', '/usr/local', '/usr']
    for p in dict.fromkeys(p for p in prefixes if p):
        inc = os.path.join(p, 'include')
        if not os.path.exists(os.path.join(inc, 'armadillo')):
            continue
        for sub in ('lib', 'lib64', 'lib/x86_64-linux-gnu', 'lib/aarch64-linux-gnu'):
            d = os.path.join(p, sub)
            if os.path.isdir(d) and any(f.startswith('libarmadillo') for f in os.listdir(d)):
                return inc, d
    return None

ARMADILLO_MISSING = (
    "Armadillo (C++ linear algebra library) was not found. Install it with one of:\n"
    "  conda env create -f environment.yml   (recommended; see README)\n"
    "  brew install armadillo                (macOS)\n"
    "  apt install libarmadillo-dev          (Debian/Ubuntu)\n"
    "or point ARMADILLO_DIR at a prefix containing include/armadillo and lib/libarmadillo*."
)

# ISA flags for the C++ kernel. The default is EMPTY: the portable x86-64 SSE2
# (or arm64) baseline, chosen for portability. The extension is often built once
# on a shared filesystem and imported on whatever node a job lands on, and a
# -march=haswell binary SIGILLs on a host without AVX2.
#
# The fast kernel does not need a flag: its hot loops are compiled twice (baseline
# and `#pragma GCC target("avx,avx2,fma")`, x86-64 only) and dispatched at run time.
# GMOQ_STAGE1_ARCH raises the baseline for the --exact-kernel path only, for
# 10-14% (measured; outputs byte-identical at the precision written), e.g.
#   GMOQ_STAGE1_ARCH="-march=haswell -mtune=skylake-avx512" pip install .
# An empty value is an explicit request for the baseline. -ffast-math is never
# used: it would reassociate reductions and change the outputs.
# See docs/METHOD.md.
_DEFAULT_ARCH_FLAGS = ''


def resolve_arch_flags(env=None):
    """Resolve the ISA flags for the kernel build from the environment.

    Unset means the portable baseline.  Set means exactly what the caller asked for,
    split on whitespace, including the empty string as an explicit request for the
    baseline.  Kept as a function so the policy is testable without a compiler.
    """
    if env is None:
        env = os.environ
    arch_env = env.get('GMOQ_STAGE1_ARCH')
    return (_DEFAULT_ARCH_FLAGS if arch_env is None else arch_env).split()


_arch_flags = resolve_arch_flags()
_arma = find_armadillo()

ext_modules = [
    Extension(
        'gmoq._stage1_cpp',
        sources=['src/stage1_pybind.cpp'],
        depends=['src/stage1_fast.hpp', 'src/stage1_fast_hot.inc', 'src/stage1_pooled_fast.hpp'],
        include_dirs=[get_pybind_include()] + ([_arma[0]] if _arma else []),
        library_dirs=[_arma[1]] if _arma else [],
        # An explicit rpath: on macOS setuptools' runtime_library_dirs emits -L only,
        # which leaves @rpath/libarmadillo unresolvable outside conda.
        extra_link_args=['-Wl,-rpath,' + _arma[1]] if _arma else [],
        libraries=['armadillo'],
        # OpenMP flags are added by BuildExt once a probe compile shows which work.
        extra_compile_args=['-std=c++14', '-O3', '-fPIC'] + _arch_flags,
        language='c++'
    ),
]


def _libomp_prefixes():
    """Where a separately installed libomp may live (macOS: conda, Homebrew)."""
    out = [os.environ.get('CONDA_PREFIX'), _brew_prefix('libomp'),
           '/opt/homebrew/opt/libomp', '/usr/local/opt/libomp']
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
        if _arma is None:
            raise SystemExit('\n' + ARMADILLO_MISSING + '\n')
        flags = None
        if os.environ.get('GMOQ_STAGE1_OPENMP', '1') != '0':
            flags = next((c for c in openmp_candidates() if self._openmp_works(*c)), None)
        if flags is None:
            sys.stderr.write(
                '\n*** gmoq-stage1: building WITHOUT OpenMP; the C++ kernel will run on one '
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
