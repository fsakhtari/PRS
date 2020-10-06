#!/bin/bash

#SBATCH --mem=400g
#SBATCH --partition=discovery
#SBATCH --mail-user=akhtarifs@nih.gov
#SBATCH --mail-type=END, FAIL

Rscript --verbose $PRS_COMMON_DIR/scripts/setup_validn_data.R $*