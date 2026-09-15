import pathlib
import subprocess
import sys
import tempfile

root = pathlib.Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix='gmoq-output-guard-') as directory:
    out = pathlib.Path(directory)
    sentinel = out / 'completed-result.txt'
    sentinel.write_text('preserve this result\n')
    cases = [
        ('tests/validity_pilot/pilot_stage1.py', []),
        ('tests/validity_pilot/covariate_basis_probe.py', ['--repo', str(root)]),
        ('scripts/validation/real_chr22_smoke.py', ['--bgen', 'unused', '--sample', 'unused', '--snps', 'unused']),
        ('scripts/validation/build_logbmi_matched_baseline.py',
         [item for flag in ('--bgen', '--sample', '--snps', '--pheno', '--covar', '--keep', '--loco', '--stage1')
          for item in (flag, 'unused')]),
    ]
    for script, extra in cases:
        result = subprocess.run([sys.executable, str(root / script), '--out', directory, *extra],
                                text=True, capture_output=True)
        assert result.returncode == 2, (script, result.returncode, result.stderr)
        assert 'Refusing to overwrite' in result.stderr, result.stderr
        assert sentinel.read_text() == 'preserve this result\n'
        assert list(out.iterdir()) == [sentinel]
        print('PASS existing-output refusal:', script)
    script = 'scripts/validation/run_real_chr22_empirical.sh'
    result = subprocess.run(['bash', str(root / script), directory], text=True, capture_output=True)
    assert result.returncode == 2, result.stderr
    assert 'Refusing to overwrite experiment directory' in result.stderr, result.stderr
    assert sentinel.read_text() == 'preserve this result\n'
    assert list(out.iterdir()) == [sentinel]
    print('PASS existing-output refusal:', script)
