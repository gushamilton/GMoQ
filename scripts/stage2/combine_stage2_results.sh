#!/usr/bin/env bash
set -euo pipefail

case "$(hostname)" in
  bp1-compute*) ;;
  *) echo "Run Stage 2 concatenation on a BP compute node (bp1-compute*)." >&2; exit 2 ;;
esac

## Concatenate already-computed Stage 2 chromosome streams.  This intentionally
## does not combine Stage 1 covariance: each chromosome was fitted against its
## own LOCO-adjusted baseline and Sigma_0 sidecar.

input_dir=""
out_prefix=""
trait=""
model=""
source=""

usage() {
  cat >&2 <<'EOF'
Usage: bash combine_stage2_results.sh --input-dir DIR --out PREFIX \
  --trait TRAIT --model MODEL --cov-source pooled|jackknife
EOF
  exit 2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-dir) input_dir="${2:?missing value for --input-dir}"; shift 2 ;;
    --out) out_prefix="${2:?missing value for --out}"; shift 2 ;;
    --trait) trait="${2:?missing value for --trait}"; shift 2 ;;
    --model) model="${2:?missing value for --model}"; shift 2 ;;
    --cov-source) source="${2:?missing value for --cov-source}"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "unknown option: $1" >&2; usage ;;
  esac
done

[[ -n "$input_dir" && -n "$out_prefix" && -n "$trait" && -n "$model" && -n "$source" ]] || usage
[[ "$source" == "pooled" || "$source" == "jackknife" ]] || {
  echo "--cov-source must be pooled or jackknife" >&2
  exit 2
}

mkdir -p "$(dirname "$out_prefix")"
out_tsv="${out_prefix}.tsv.gz"
out_cov="${out_prefix}.param_cov.gz"
out_manifest="${out_prefix}.manifest.tsv"
for f in "$out_tsv" "$out_cov" "$out_manifest"; do
  [[ ! -e "$f" ]] || { echo "refusing to overwrite $f; use a new exclusive run root" >&2; exit 3; }
done

files=()
cov_files=()
manifest_files=()
for chr in $(seq 1 22); do
  pad=$(printf '%02d' "$chr")
  f="${input_dir}/chr${pad}.tsv.gz"
  m="${input_dir}/chr${pad}.manifest.tsv"
  [[ -s "$f" ]] || { echo "missing Stage 2 chromosome output: $f" >&2; exit 4; }
  gzip -t "$f" >/dev/null
  [[ -s "$m" ]] || { echo "missing Stage 2 chromosome manifest: $m" >&2; exit 4; }
  files+=("$f")
  manifest_files+=("$m")
  c="${input_dir}/chr${pad}.param_cov.gz"
  if [[ -e "$c" ]]; then
    [[ -s "$c" ]] || { echo "empty parameter covariance: $c" >&2; exit 4; }
    gzip -t "$c" >/dev/null
    cov_files+=("$c")
  fi
done

if [[ "${#cov_files[@]}" -ne 0 && "${#cov_files[@]}" -ne 22 ]]; then
  echo "parameter covariance is present for ${#cov_files[@]} of 22 chromosomes" >&2
  exit 4
fi

first_header=$(gzip -cd "${files[0]}" | sed -n '1p')
[[ -n "$first_header" ]] || { echo "empty Stage 2 header" >&2; exit 4; }
for f in "${files[@]}"; do
  h=$(gzip -cd "$f" | sed -n '1p')
  [[ "$h" == "$first_header" ]] || {
    echo "Stage 2 header mismatch in $f" >&2
    exit 5
  }
done

tmp_tsv="${out_tsv}.tmp.$$"
{
  gzip -cd "${files[0]}" | sed -n '1p'
  for f in "${files[@]}"; do
    gzip -cd "$f" | tail -n +2
  done
} | gzip > "$tmp_tsv"
mv "$tmp_tsv" "$out_tsv"

if [[ "${#cov_files[@]}" -eq 22 ]]; then
  tmp_cov="${out_cov}.tmp.$$"
  cat "${cov_files[@]}" > "$tmp_cov"
  mv "$tmp_cov" "$out_cov"
else
  out_cov=""
fi

tmp_manifest="${out_manifest}.tmp.$$"
has_cov=0
[[ "${#cov_files[@]}" -eq 22 ]] && has_cov=1
{
  echo -e "trait\tmodel\tcovariance_source\tchromosome\tstage2_file\tparameter_covariance_file"
  for m in "${manifest_files[@]}"; do
    tail -n +2 "$m" | awk -F '\t' -v trait="$trait" -v model="$model" -v source="$source" \
      -v input="$input_dir" -v has_cov="$has_cov" \
      'BEGIN { OFS="\t" } { cov = has_cov ? input "/chr" sprintf("%02d", $3) ".param_cov.gz" : ""; print trait, model, source, $3, input "/chr" sprintf("%02d", $3) ".tsv.gz", cov }'
  done
} > "$tmp_manifest"
mv "$tmp_manifest" "$out_manifest"

echo "combined_trait=${trait} model=${model} cov_source=${source}"
echo "stage2=${out_tsv}"
if [[ -n "$out_cov" ]]; then echo "param_cov=${out_cov}"; fi
echo "manifest=${out_manifest}"
