#!/bin/bash

### Risk score analysis ###


## init() ##

# The common PRS directory for all diseases
export PRS_COMMON_DIR='../../common'
# The PRS directory of this specific phenotype for which PRS is being computed
export PRS_PHENO_DIR='..'


sbatch --output=$PRS_PHENO_DIR/risk_scores_tt50/risk_scores.log $PRS_COMMON_DIR/scripts/risk_scores.sh

sbatch --output=$PRS_PHENO_DIR/risk_scores_tt50_nopd/risk_scores.log $PRS_COMMON_DIR/scripts/risk_scores.sh --prediabetes FALSE