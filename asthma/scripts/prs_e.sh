#!/bin/bash

### Cases-only PRSxE test ###


## init() ##

# The common PRS directory for all diseases
export PRS_COMMON_DIR='../../common'
# The PRS directory of this specific phenotype for which PRS is being computed
export PRS_PHENO_DIR='..'

sbatch --output=$PRS_PHENO_DIR/results/tmp/co_prs_e.log $PRS_COMMON_DIR/scripts/co_prs_e.sh --pheno he_d030_asthma_PARQ
