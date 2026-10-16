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

# ISA flags for the C++ kernel. The default is EMPTY: the portable x86-64 SSE2
# (or arm64) baseline, chosen for portability. The extension is often built once
# on a shared filesystem and imported on whatever node a job lands on, and a
# -march=haswell binary SIGILLs on a host without AVX2.
#
# The fast kernel does not need a flag: its hot loops are compiled twice (baseline
# and `#pragma GCC target("avx,avx2,fma")`, x86-64 only) and dispatched at run time.
# FUNGWAS_STAGE1_ARCH raises the baseline for the --exact-kernel path only, for
# 10-14% (measured; outputs byte-identical at the precision written), e.g.
#   FUNGWAS_STAGE1_ARCH="-march=haswell -mtune=skylake-avx512" pip install .
# An empty value is an explicit request for the baseline. -ffast-math is never
# used: it would reassociate reductions and change the outputs.
# See docs/STAGE1_FAST_KERNEL.md.
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
