#!/bin/bash

#SBATCH --partition=discovery
#SBATCH --mail-user=akhtarifs
#SBATCH --mail-type=END, FAIL

Rscript --verbose setup_summstats.R $*