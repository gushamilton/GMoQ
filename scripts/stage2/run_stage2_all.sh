#!/usr/bin/env bash
set -euo pipefail

case "$(hostname)" in
  bp1-compute*) ;;
  *) echo "Run Stage 2 on a BP compute node (bp1-compute*)." >&2; exit 2 ;;
esac

## Prepare matched W (optional), run one model for one or both covariance
## sources chromosome-by-chromosome, then concatenate result streams.  This
## wrapper submits nothing and allocates nothing; invoke it inside the one
## coordinated compute allocation.

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
trait=""
model=""
sources_raw="pooled"
stage1_dir=""
stage1_root=""
loco_file=""
pheno_file=""
covar_file=""
keep_file=""
bgen_sample_file="/bp1/mrcieu1/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample"
weights_root="results/paper/fungwas/weights_current"
out_root="results/paper/fungwas/stage2_current"
r_lib=""
threads="${SLURM_CPUS_PER_TASK:-1}"
prepare=0

usage() {
  cat >&2 <<'EOF'
Usage: bash run_stage2_all.sh --trait TRAIT --model MODEL \
  (--stage1-dir TRAIT_DIR | --stage1-root ARRAY_ROOT) \
  --loco FILE --pheno FILE --covar FILE --keep FILE \
  [--cov-source pooled[,jackknife]] [--weights-root DIR] [--out-root DIR] \
  [--bgen-sample FILE] [--r-lib DIR] [--threads N] [--prepare-weights]

Models: bmi {additive,shifted,both,logscale_vqtl}; log_bmi/weight
{additive,shifted,both}; testosterone and shrink traits {skew};
jcv_log_all/mcv_log_all {mixture}.
EOF
  exit 2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --trait) trait="${2:?missing value for --trait}"; shift 2 ;;
    --model) model="${2:?missing value for --model}"; shift 2 ;;
    --cov-source) sources_raw="${2:?missing value for --cov-source}"; shift 2 ;;
    --stage1-dir) stage1_dir="${2:?missing value for --stage1-dir}"; shift 2 ;;
    --stage1-root) stage1_root="${2:?missing value for --stage1-root}"; shift 2 ;;
    --loco) loco_file="${2:?missing value for --loco}"; shift 2 ;;
    --pheno) pheno_file="${2:?missing value for --pheno}"; shift 2 ;;
    --covar) covar_file="${2:?missing value for --covar}"; shift 2 ;;
    --keep) keep_file="${2:?missing value for --keep}"; shift 2 ;;
    --bgen-sample) bgen_sample_file="${2:?missing value for --bgen-sample}"; shift 2 ;;
    --weights-root) weights_root="${2:?missing value for --weights-root}"; shift 2 ;;
    --out-root) out_root="${2:?missing value for --out-root}"; shift 2 ;;
    --r-lib) r_lib="${2:?missing value for --r-lib}"; shift 2 ;;
    --threads) threads="${2:?missing value for --threads}"; shift 2 ;;
    --prepare-weights) prepare=1; shift ;;
    -h|--help) usage ;;
    *) echo "unknown option: $1" >&2; usage ;;
  esac
done

[[ -n "$trait" && -n "$model" && -n "$loco_file" && \
  -n "$pheno_file" && -n "$covar_file" && -n "$keep_file" ]] || usage
if [[ -n "$stage1_dir" && -n "$stage1_root" ]]; then
  echo "pass only one of --stage1-dir or --stage1-root" >&2
  exit 2
fi
if [[ -z "$stage1_dir" && -z "$stage1_root" ]]; then
  echo "one of --stage1-dir or --stage1-root is required" >&2
  exit 2
fi

case "$trait:$model" in
  bmi:additive|bmi:shifted|bmi:both|bmi:logscale_vqtl) family=shifted_lognormal ;;
  log_bmi:additive|log_bmi:shifted|log_bmi:both) family=shifted_lognormal ;;
  weight:additive|weight:shifted|weight:both) family=shifted_lognormal ;;
  testosterone:skew|testosterone_shrink_0.2:skew|testosterone_shrink_0.5:skew|testosterone_shrink_0.7:skew) family=skew ;;
  jcv_log_all:mixture|mcv_log_all:mixture) family=mixture ;;
  *) echo "unsupported trait/model: ${trait}/${model}" >&2; exit 2 ;;
esac

stage1_input_root="$stage1_dir"
[[ -n "$stage1_input_root" ]] || stage1_input_root="$stage1_root"
for f in "$stage1_input_root" "$loco_file" "$pheno_file" "$covar_file" "$keep_file" "$bgen_sample_file"; do
  [[ -e "$f" ]] || { echo "missing required path: $f" >&2; exit 3; }
done
[[ "$threads" =~ ^[0-9]+$ && "$threads" -ge 1 ]] || { echo "--threads must be positive" >&2; exit 2; }

trait_tag=${trait//./_}
weights_trait_dir="${weights_root}/${trait_tag}"
IFS=',' read -r -a sources <<< "$sources_raw"
[[ "${#sources[@]}" -gt 0 ]] || { echo "no covariance source requested" >&2; exit 2; }
for source in "${sources[@]}"; do
  [[ "$source" == "pooled" || "$source" == "jackknife" ]] || {
    echo "unsupported covariance source: $source" >&2; exit 2;
  }
done

resolve_stage1() {
  local pad="$1"
  local base="${trait_tag}_chr${pad}.stage1.tsv.gz"
  local candidates=()
  if [[ -n "$stage1_root" ]]; then
    # Current array layout: ROOT/<trait-tag>/chrNN/<trait-tag>_chrNN.*.
    candidates+=("${stage1_root}/${trait_tag}/chr${pad}/${base}")
    # Also accept a literal trait directory where punctuation was retained.
    [[ "$trait" == "$trait_tag" ]] || candidates+=("${stage1_root}/${trait}/chr${pad}/${base}")
    # A trait-specific root with chrNN subdirectories is unambiguous too.
    candidates+=("${stage1_root}/chr${pad}/${base}")
  else
    # Legacy flat layout under a trait-specific --stage1-dir.
    candidates+=("${stage1_dir}/chr${pad}/${base}" "${stage1_dir}/${base}")
  fi
  local found=()
  local candidate
  for candidate in "${candidates[@]}"; do
    [[ -s "$candidate" ]] && found+=("$candidate")
  done
  if [[ "${#found[@]}" -ne 1 ]]; then
    echo "expected exactly one Stage 1 file for chr${pad}; found ${#found[@]}" >&2
    printf '  candidate: %s\n' "${candidates[@]}" >&2
    printf '  existing: %s\n' "${found[@]}" >&2
    exit 4
  fi
  printf '%s\n' "${found[0]}"
}

for chr_num in $(seq 1 22); do
  chr_pad=$(printf '%02d' "$chr_num")
  stage1=$(resolve_stage1 "$chr_pad")
  weight_dir="${weights_trait_dir}/chr${chr_pad}"

  if [[ "$prepare" -eq 1 ]]; then
    prep_args=(
      --trait "$trait" --chr "$chr_num" --stage1 "$stage1"
      --pheno "$pheno_file" --covar "$covar_file" --keep "$keep_file"
      --loco "$loco_file" --bgen-sample "$bgen_sample_file"
      --family "$family" --out-dir "$weight_dir" --loco-scale regression
    )
    [[ -n "$r_lib" ]] && prep_args+=(--r-lib "$r_lib")
    Rscript "${script_dir}/prepare_matched_weights.R" "${prep_args[@]}"
  fi

  weight_file="${weight_dir}/weights_${model}.rds"
  [[ -s "$weight_file" ]] || {
    echo "missing matched weight file: $weight_file (use --prepare-weights)" >&2
    exit 4
  }

  for source in "${sources[@]}"; do
    source_dir="${out_root}/${trait_tag}/${model}/${source}"
    mkdir -p "$source_dir"
    out_prefix="${source_dir}/chr${chr_pad}"
    run_args=(
      --trait "$trait" --model "$model" --chr "$chr_num"
      --stage1 "$stage1" --weights "$weight_file" --out "$out_prefix"
      --cov-source "$source" --threads "$threads"
    )
    if [[ "$source" == "pooled" ]]; then
      run_args+=(--sigma0-file "${stage1%.stage1.tsv.gz}.sigma0.tsv.gz")
    else
      run_args+=(--cov-file "${stage1%.stage1.tsv.gz}.cov.gz")
    fi
    [[ -n "$r_lib" ]] && run_args+=(--r-lib "$r_lib")
    Rscript "${script_dir}/run_stage2_chr.R" "${run_args[@]}"
  done
done

for source in "${sources[@]}"; do
  input_dir="${out_root}/${trait_tag}/${model}/${source}"
  combine_prefix="${out_root}/${trait_tag}/${model}/${trait_tag}_${model}_${source}"
  combine_args=(
    --input-dir "$input_dir" --out "$combine_prefix"
    --trait "$trait" --model "$model" --cov-source "$source"
  )
  bash "${script_dir}/combine_stage2_results.sh" "${combine_args[@]}"
done
