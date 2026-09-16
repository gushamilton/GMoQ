"""
Build script for Stage 1 C++ extension using pybind11.

Usage:
    python setup.py build_ext --inplace
    
Or for development:
    pip install -e .
"""

import os
import sys
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
# The kernel was historically built at the x86-64 SSE2 baseline - two-wide
# double-precision, no FMA - even on AVX-512 hardware.  Raising the baseline is a
# build-time choice, not a run-time one, and it has to be made conservatively: a
# BluePebble `short` allocation can land on Skylake-SP, Cascade Lake *or* an AMD node,
# and a binary built with -march=native on one of those can SIGILL on another.
#
# The default is left EMPTY, i.e. the SSE2 baseline is kept, for PORTABILITY -- not
# because the ISA is worthless. That second claim was true of the pre-2026-09-16 kernel
# and is now stale; it is kept here only as history.
#
#   HISTORICAL (old kernel, superseded): interleaved on one node, 15 timings per build of
#   a production-shaped 500-SNP batch, the median kernel time with -march=haswell was
#   0.9997 of the SSE2 baseline against a 9-22% run-to-run IQR. The kernel was
#   memory-bound (0.38 flop/byte, IPC 0.7) with inner loops striding down the wrong axis,
#   so GCC contracted FMAs -- 131 of the 136 it emitted were *scalar* -- but could not
#   widen anything.
#
#   CURRENT (2026-09-16, after the panel kernels were reworked): -march=haswell is worth
#   about 13% of wall clock. Full chr22 testosterone with LOCO on 4 physical cores,
#   --cov-mode pooled-only: 172.4s -> 142.2s, i.e. 101.8 -> 123.4 SNPs/s. Decompressed
#   output is BYTE-IDENTICAL to the SSE2 build in both full and pooled-only modes, so
#   unlike the historical note there is now measured gain at no measured cost.
#
# It is still not the default because this is a build-time choice for a package other
# people install. -march=haswell needs AVX2, i.e. Intel 2013+ or AMD Zen+; that covers
# every current BluePebble node, but a stale binary shipped to older hardware SIGILLs.
# Sites that control their own hardware should turn it on:
#
#   FUNGWAS_STAGE1_ARCH="-march=haswell -mtune=skylake-avx512" pip install -e stage1
#
# -march=skylake-avx512 is no better: GCC 8 emits zero zmm registers for it and produces
# bitwise identical results to -march=haswell. -march=native is NOT portable -- a
# BluePebble `short` allocation can land on Skylake-SP, Cascade Lake or an AMD node.
#
# -ffast-math is deliberately never used.
_DEFAULT_ARCH_FLAGS = ''
_arch_env = os.environ.get('FUNGWAS_STAGE1_ARCH')
_arch_flags = (_DEFAULT_ARCH_FLAGS if _arch_env is None else _arch_env).split()

ext_modules = [
    Extension(
        'fungwas_stage1._stage1_cpp',
        sources=['src/stage1_pybind.cpp'],
        include_dirs=[
            get_pybind_include(),
            get_armadillo_include(),
        ],
        library_dirs=[get_armadillo_lib()] if get_armadillo_lib() else [],
        runtime_library_dirs=[get_armadillo_lib()] if get_armadillo_lib() else [],
        libraries=['armadillo'],
        extra_compile_args=['-std=c++14', '-O3', '-fopenmp', '-fPIC'] + _arch_flags,
        extra_link_args=['-fopenmp'],
        language='c++'
    ),
]

class BuildExt(build_ext):
    """Custom build_ext command with compiler customization."""
    def build_extensions(self):
        # Remove -Wstrict-prototypes (not valid for C++)
        if '-Wstrict-prototypes' in self.compiler.compiler_so:
            self.compiler.compiler_so.remove('-Wstrict-prototypes')
        super().build_extensions()

setup(
    name='fungwas-stage1',
    version='0.1.0',
    ext_modules=ext_modules,
    cmdclass={'build_ext': BuildExt},
    packages=['fungwas_stage1'],
    python_requires='>=3.8',
    install_requires=[
        'numpy>=1.20',
        'scipy',
    ],
    extras_require={
        'bgen': ['bgen>=1.9.7'],
    },
    entry_points={
        'console_scripts': [
            'fungwas-stage1=fungwas_stage1.cli:main',
        ],
    },
)
