#!/bin/bash

#SBATCH --partition=discovery
#SBATCH --mem=500g
#SBATCH --exclusive
#SBATCH --mail-user=akhtarifs@nih.gov
#SBATCH --mail-type=END, FAIL

Rscript --verbose $PRS_COMMON_DIR/scripts/setup_validn_data.R $*