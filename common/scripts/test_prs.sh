#!/bin/bash

#SBATCH --partition=discovery
#SBATCH --mail-user=akhtarifs
#SBATCH --mail-type=END, FAIL

Rscript --verbose $PRS_COMMON_DIR/scripts/test_prs.R