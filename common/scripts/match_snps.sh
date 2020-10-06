#!/bin/bash

#SBATCH --partition=discovery
#SBATCH --mem=400g
#SBATCH --exclusive
#SBATCH --mail-user=akhtarifs@nih.gov
#SBATCH --mail-type=END, FAIL

Rscript --verbose $PRS_COMMON_DIR/scripts/match_snps.R $*
