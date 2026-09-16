#!/bin/bash
source /etc/profile.d/modules.sh
module load languages/R/4.5.1
cd /user/work/fh6520/gmoq_covmixW
for s in verify_W mechanism robust onset mixdgp; do
  echo "=== running $s ==="
  Rscript --vanilla "src/$s.R" > "out/$s.txt" 2>&1
  echo "  exit $? -> out/$s.txt ($(wc -l < out/$s.txt) lines)"
done
echo DIAGDONE
