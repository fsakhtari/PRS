#!/bin/bash

#SBATCH --partition=discovery
#SBATCH --mail-user=akhtarifs
#SBATCH --mail-type=END, FAIL

Rscript --verbose $PRS_COMMON_DIR/scripts/setup_LDrefpanel.R $*