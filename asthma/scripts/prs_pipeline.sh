#!/bin/bash

### PRS analysis pipeline ###


## init() ##

# The common PRS directory for all diseases
export PRS_COMMON_DIR='../../common'
# The PRS directory of this specific phenotype for which PRS is being computed
export PRS_PHENO_DIR='..'


## execute() ##
# This is the ordered list of commands/scripts to be executed for this analysis.

echo "PRS pipeline for phenotype 'Asthma' initiated"

# QC and setup the LD reference panel data.
sbatch --output=$PRS_PHENO_DIR/logs/setup_LDrefpanel.log $PRS_COMMON_DIR/scripts/setup_LDrefpanel.sh

# QC and setup the validation dataset.
# asthma = f.22127.0.0 in the UK Biobank data
sbatch --output=$PRS_PHENO_DIR/logs/setup_validn_data.log $PRS_COMMON_DIR/scripts/setup_validn_data.sh --pheno "f.22127.0.0"

# QC and setup the test dataset.
sbatch --output=$PRS_PHENO_DIR/logs/setup_test_data.log $PRS_COMMON_DIR/scripts/setup_test_data.sh --pheno he_d030_asthma_PARQ

# QC and setup the summary statistics data. There is a different script for each summary statistic dataset
sbatch --wait --output=$PRS_PHENO_DIR/logs/setup_summstats.log $PRS_PHENO_DIR/scripts/setup_summstats.sh

# Match SNPs + QC across datasets.
sbatch --wait --output=$PRS_PHENO_DIR/logs/match_snps.log $PRS_COMMON_DIR/scripts/match_snps.sh

# For chr = 1 to 22
  # Compute PRS per chromosome. Runs in parallel across chromosomes.
  sbatch --wait --output=$PRS_PHENO_DIR/logs/compute_prs_batch.log $PRS_COMMON_DIR/scripts/compute_prs_batch.sh

# Compute PRS across all chromosomes for test data. Test combined PRS on test (EPR) data.
sbatch --wait --output=$PRS_PHENO_DIR/logs/test_prs.log $PRS_COMMON_DIR/scripts/test_prs.sh


## uninit() ##

# Wait for all jobs/commands to finish before uninitializing
wait

unset PRS_COMMON_DIR
unset PRS_PHENO_DIR

echo $PRS_COMMON_DIR
echo $PRS_PHENO_DIR

echo "PRS pipeline for phenotype 'Asthma' complete"

### End ###
