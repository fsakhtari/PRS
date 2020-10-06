#!/bin/bash

#SBATCH --partition=discovery
#SBATCH --mail-user=akhtarifs
#SBATCH --mail-type=END, FAIL

Rscript --verbose $PRS_COMMON_DIR/scripts/co_prs_e.R $*
