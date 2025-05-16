#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=120
#SBATCH --time=48:00:00
#SBATCH --partition=standardqueue
#SBATCH --mem=512G
#SBATCH --job-name=gwas-e2g-enrichment          # Job name
#SBATCH --output=e2g-enrichment_colata.out  # Standard output file
#SBATCH --output=e2g-enrichment_colata.err  # Error log file

source /opt/software/anaconda3/2023.03-py3.10/etc/profile.d/conda.sh
conda activate /projects/cbmr_shared/people/wkq953/non-GDPR/.conda/envs/run_snakemake

# config_file="/projects/hh_liver_cohorts-AUDIT/data/FLINC/multiome/processed/pipeline_output_freezes/FLINC/2.10.0/objects/real/GWAS_E2G_benchmarking/config/config_FLINC.yml"
config_file="/projects/colata-AUDIT/data/multiome/processed/pipeline_output_freezes/colata/2.11.0/objects/real/GWAS_E2G_benchmarking/config/config_colata.yml"

snakemake -j1 --use-conda --configfile $config_file 