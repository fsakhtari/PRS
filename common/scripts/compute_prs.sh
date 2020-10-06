#!/bin/bash
#SBATCH --mem=81920
#SBATCH --partition=discovery
#SBATCH --mail-user=akhtarifs
#SBATCH --mail-type=END, FAIL

Rscript --verbose $PRS_COMMON_DIR/scripts/compute_prs.R $*
