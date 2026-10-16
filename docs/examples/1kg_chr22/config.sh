# Settings for the 1000 Genomes example. Override any of them from the environment.
# Data and outputs go OUTSIDE the repository by default.
GMOQ_1KG_DIR="${GMOQ_1KG_DIR:-$HOME/gmoq_1kg_example}"
VCF="${VCF:-$GMOQ_1KG_DIR/data/chr22.vcf.gz}"   # 1000 Genomes phase 3 chr22 (GRCh37)
WORK="${WORK:-$GMOQ_1KG_DIR/work}"               # generated files (safe to delete)
PLINK2="${PLINK2:-plink2}"                       # plink2 on PATH, or a full path
THREADS="${THREADS:-4}"
