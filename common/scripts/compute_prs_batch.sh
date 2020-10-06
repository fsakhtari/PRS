#!/bin/bash

#SBATCH --partition=discovery
#SBATCH --mail-user=akhtarifs
#SBATCH --mail-type=END, FAIL
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=6

for (( i=1; i<=22; i++ ));
do
    echo "Computing PRS for chromosome $i"
    sbatch --output $PRS_PHENO_DIR/logs/compute_prs.chr$i.log $PRS_COMMON_DIR/scripts/compute_prs.sh --chr=$i
done