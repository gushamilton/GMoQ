"""The kernel's ISA baseline is an opt-in, and these pin why.

Measured end to end on BluePebble (2026-09-16, Xeon Gold 5120, one chr22 testosterone
scan: 17,550 variants, N=399,710, T=45, --blocks 1000), building the same sources three
ways and rotating the build order between passes so node-load drift cannot land
preferentially on one build:

* ``-march=haswell -mtune=skylake-avx512`` is worth about 1.12x on the block kernel and
  1.14x on the pooled-only kernel at 4 threads, and 1.23x on the pooled-only kernel at 8.
* It is **byte-identical**: the full chr22 ``.stage1.tsv.gz``, ``.sigma0.tsv.gz``,
  ``.cov.ids.tsv.gz`` and the entire 72,657,000-byte ``.cov.gz`` (18,164,250 float32
  values) all decompress to the same bytes as the SSE2 build's.  FMA contraction moves
  float64 sums in the last bits; the artefacts are written at ``%.6g`` text and float32
  binary, which is some nine orders of magnitude coarser, so nothing survives to the file.
* ``-ffp-contract=off`` keeps the AVX2 width and removes every FMA, and is byte-identical
  too -- so it is not needed to protect the outputs, and it costs about a third of the gain.

So the reason the default stays at the portable baseline is NOT results and NOT the
project's byte-comparison verification, both of which survive AVX2 untouched.  It is that
the built extension outlives its build host: it is compiled once into a package directory
on a shared filesystem and then imported by whatever node an allocation lands on, and a
``-march=haswell`` binary SIGILLs without a diagnostic on a host without AVX2.

These tests pin the policy, not the speed.  If someone flips ``_DEFAULT_ARCH_FLAGS``, the
first test fails and this docstring is what they have to argue with.
"""

import importlib.util
import os
from pathlib import Path

import pytest

SETUP_PY = Path(__file__).resolve().parents[1] / "setup.py"


@pytest.fixture(scope="module")
def setup_module_obj():
    """Import stage1/setup.py without running setup().

    The ``if __name__ == "__main__"`` guard in setup.py is what makes this possible;
    importing under any other name must not build, install or write anything.
    """
    spec = importlib.util.spec_from_file_location("fungwas_stage1_setup", SETUP_PY)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_setup_py_is_importable_without_running_setup(setup_module_obj):
    """Guard is in place: importing setup.py defines the policy and builds nothing."""
    assert setup_module_obj.__name__ != "__main__"
    assert hasattr(setup_module_obj, "resolve_arch_flags")


def test_the_default_isa_is_the_portable_baseline(setup_module_obj):
    """AVX2 is opt-in.  See the module docstring for why, and why it is not about FP."""
    assert setup_module_obj._DEFAULT_ARCH_FLAGS == ''
    assert setup_module_obj.resolve_arch_flags({}) == []


def test_unset_and_explicitly_empty_both_mean_the_baseline(setup_module_obj):
    resolve = setup_module_obj.resolve_arch_flags
    assert resolve({}) == []
    assert resolve({'FUNGWAS_STAGE1_ARCH': ''}) == []
    assert resolve({'FUNGWAS_STAGE1_ARCH': '   '}) == []


def test_the_documented_opt_in_parses_to_the_flags_it_names(setup_module_obj):
    """The exact string in setup.py's comment is the string that was measured."""
    resolve = setup_module_obj.resolve_arch_flags
    assert resolve({'FUNGWAS_STAGE1_ARCH': '-march=haswell -mtune=skylake-avx512'}) == [
        '-march=haswell', '-mtune=skylake-avx512',
    ]
    # The bit-identity-belt-and-braces variant, also measured byte-identical.
    assert resolve({
        'FUNGWAS_STAGE1_ARCH': '-march=haswell -mtune=skylake-avx512 -ffp-contract=off'
    }) == ['-march=haswell', '-mtune=skylake-avx512', '-ffp-contract=off']
    assert resolve({'FUNGWAS_STAGE1_ARCH': '-march=native'}) == ['-march=native']


def test_resolve_reads_the_real_environment_when_not_given_one(setup_module_obj, monkeypatch):
    resolve = setup_module_obj.resolve_arch_flags
    monkeypatch.setenv('FUNGWAS_STAGE1_ARCH', '-march=haswell')
    assert resolve() == ['-march=haswell']
    monkeypatch.delenv('FUNGWAS_STAGE1_ARCH')
    assert resolve() == []


def test_fast_math_is_never_compiled_in(setup_module_obj):
    """-ffast-math licenses reassociating the reductions, which WOULD move the outputs.

    Unlike -mfma, whose effect is nine orders of magnitude below the serialisation
    resolution, reassociation changes the summation itself.  It must never reach the
    baseline compile line, whatever the environment asks for.
    """
    base = [a for a in setup_module_obj.ext_modules[0].extra_compile_args
            if a not in setup_module_obj.resolve_arch_flags()]
    assert '-ffast-math' not in base
    assert '-Ofast' not in base          # -Ofast implies -ffast-math
    assert '-O3' in base
    assert '-fopenmp' in base


def test_the_arch_flags_land_on_the_compile_line(setup_module_obj):
    """Whatever the policy resolves to is what the extension is actually built with."""
    args = setup_module_obj.ext_modules[0].extra_compile_args
    for flag in setup_module_obj._arch_flags:
        assert flag in args
    assert os.environ.get('FUNGWAS_STAGE1_ARCH') is not None or \
        setup_module_obj._arch_flags == []
