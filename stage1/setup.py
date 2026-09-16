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
# The old note here said arch flags "bought nothing".  That was measured on the
# exact kernel, whose inner loops stride down the wrong axis so GCC can contract
# FMAs but cannot widen anything - 131 of the 136 it emits are scalar.  It does
# not carry over to the fast kernel, which is why that one dispatches its own.
#
# The override still works and still applies to the whole translation unit:
#   FUNGWAS_STAGE1_ARCH="-march=haswell -mtune=skylake-avx512"  # AVX2+FMA, portable
#   FUNGWAS_STAGE1_ARCH="-march=native"                         # NOT portable
# Setting it changes the exact kernel's arithmetic, so a run that must reproduce
# historical output byte for byte should leave it alone.
#
# -ffast-math is deliberately never used.
_DEFAULT_ARCH_FLAGS = ''
_arch_env = os.environ.get('FUNGWAS_STAGE1_ARCH')
_arch_flags = (_DEFAULT_ARCH_FLAGS if _arch_env is None else _arch_env).split()

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
