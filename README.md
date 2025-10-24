# GWAS E2G benchmarking
Benchmark enhancer-gene predictions against fine-mapped GWAS variants

Code adapted from the following pipelines:
- [ABC-Max](https://github.com/EngreitzLab/ABC-GWAS-Paper/blob/main/ABC-Max)
- [eQTLEnrichment](https://github.com/EngreitzLab/eQTLEnrichment/tree/integrated)

This pipeline performs two analyses to benchmark enhancer-gene predictions against noncoding fine-mapped GWAS variants:
1. **Variant enrichment and recall:** How many fine-mapped GWAS variants for a given trait overlap predicted enhancers in a given cell type? How enriched are the GWAS variants compared to 1000G SNPs? Predictions are evaluated in both threshold-dependent and -independent analyses. 
2. **Linking variants to causal genes:** What are the precision and recall of enhancer-gene predictions in a given cell type at identifying causal genes for credible sets for a given trait? This analysis uses a silver-standard set of noncoding credible set–causal gene links (inferred using independent coding variants as described in [Weeks et al., 2023](https://doi.org/10.1038/s41588-023-01443-6)). Predictions are evaluted independently and in conjunction with the orthogonal Polygenic Priority Score (PoPS).

<hr>

## System requirements

### Hardware requirements

For optimal performance, we suggest using a computer equipped with 32+ GB RAM to run this pipeline.

### Software requirements

This pipeline has been successfully tested with Linux systems.

The software dependencies and versions to run the pipeline are: `python<=3.11`, `mamba=1.5.11`, `snakemake=7`, `polars>1.0`
All software dependencies and versions used within the pipeline are listed in `workflow/envs/GWAS_env.yml`. 

<hr>

## Running the pipeline
1. Clone this repository
```
git@github.com:EngreitzLab/GWAS_E2G_benchmarking.git
```
2. Edit configuration files, including downloading necessary reference files (see below)
3. Activate a conda environment with [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html), [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html), and [polars](https://docs.pola.rs/) installed.
4. Run the pipeline
```
snakemake -j1 --configfile config/config_example.yml --use-conda
```

<hr>

## Configuration files
This pipeline requires five configuration files to provide flexibility in evaluating a enhancer predictions from many models in many cell types against GWAS variants in many traits in any desired combination while accounting for redundancy in related cell types and traits. These files are explained below as clearly as possible.

### Main configuration file
An example of the main config file is included at `config/config.yml`. The following fields are allowed and required unless otherwise noted:
- **results:** output directory name
- **methodsTable:** configuration file with information about enhancer-gene predictions methods (see below for specifications). An example is included at `config/config_methods.tsv`. 
- **predictionsTable:** configuration file for predictions for each cell type and method. An example is included at `config/config_predictions.tsv`. This is a .tsv file with the column `biosample` (referring to cell types) and additional columns for each predictive method, titled with the identifier in the `methodsTable`. Entries are file paths for non-thresholdded predictions, which must include the following columns: `chr`, `start`,`end`, `TargetGene` (gene symbol), and the designated score column. If predictions by a method do not exist for an included biosample, leave the entry blank.
- **methods:** list of methods from `methodsTable` to be benchmarked
- **comparisonsTable:** configuration file delineating which groups of biosample–trait pairs to analyze together. An example is included at `config/config_comparisons.tsv`. This is a .tsv file with the columns `name`, the identifier for a set of biosample–trait pairs, `biosample`, corresponding to a biosample in the `predictionsTable` *or* a biosample group (defined below), and `trait`, corresponding to a GWAS trait included in the `variantKey` *or* a trait group (both defined below). Each row corresponds to a biosample–trait pair in the given comparison. The aggregated variant overlap and gene-linking results for each comparison will be computed.
- **biosampleGroups:** (optional) lists indicating groups of related biosamples to for which predictions will be aggregated, then treated as an additional biosample in all analyses. Names of biosample groups must be distinct from those of biosamples in the `predictionsTable`. 
- **tratiGroups:** (optional) lists indicating groups of related traits for which variants will be merged and deduplicated, then treated as an additional trait in all analyses. Names of trait groups must be distinct from those in the `variantKey`. 
- **baseDir:** full file path of pipeline directory 
- **scratchDir:** file path to write temporary files to during processing; can be the same as **baseDir**
- **nThresholdSteps:**  number of thresholds to use to compute enrichment-recall curves; we recommend using 25 
- **thresholdPIP:** minimum PIP value for variants included in analysis; we recommend using 0.1
- **thresholdPval:** threshold p-value to use in statistical tests
- **numPoPSGenes:** when evaluating gene-linking, maximum rank of PoPS score to consider a positive prediction for a credible set; we recommend using 2
- **numPredGenes:** when evaluating gene-linking, maximum rank of prediction score to consider a positive prediction for a credible set; we recommend using 2
- **plotFixedScale:** boolean indicating whether heatmaps should be plotted with the same color scale across all predictors, or with color scales fit to the range of data displayed
- **variantKey:** configuration file for variant files for each trait. An example of this file is included at `resources/UKBB_variant_key.tsv`, which also points to variant files from fine-mapping of GWAS data for 94 traits from the UK Biobank, obtained from https://www.finucanelab.org/data. The variant key contains columns `trait` and `variant_file`. Each `variant_file` must contain columns  `chr`,`start`,`end`,`rsid` (can be any unique identifier), `pip`,`CredibleSet` (in the format `chr:start-end`),`trait`.
- **genePrioritizationTable:** file with deduplicated silver-standard variant–gene links and annotations of all proximal genes to each credible set. This file compatible with the UK Biobank traits is provided at `resources/UKBiobank.ABCGene.anyabc.tsv`
- **bgVariants:** a bed file of 1000G SNPs  with columns (no header) `chr`, `start`, `end`, `rsid`. The list of background variants we use can be downloaded on Synapse [here](https://www.synapse.org/#!Synapse:syn52264319) and was collated from https://alkesgroup.broadinstitute.org/LDSCORE/baseline_v1.1_hg38_annots/
- Genome annotations **chrSizes**, **partition**, and **TSS**, all provided here in the `resources/genome_annotation` directory.

### Methods configuration file
The following columns describing each enhancer-gene prediction method are required for the `methodsTable` configuraiton file:
- **method:** predictive method name, to match listed methods in the main config file outlined above
- **pred_name_long:** a string with the full method name to be used in plots
- **threshold:** score threshold value to be used to generate boxplots stratified by distance and enrichment heat maps (outputs (1) and (4)). We used the threshold values corresponding to 70% recall from our CRISPR benchmarking pipeline for each model.
- **score_col:** name of the column in prediction files with the prediction score
- **color:** hex code specifying what color to use in plots for this method. May be left blank and will be automatically filled.
- **inverse_predictor:** `TRUE` if high scores correspond to lower prediction confidence (*e.g.* distance to TSS), otherwise `FALSE`
- **boolean:** `TRUE` if this is a binary 0 or 1 predictor, otherwise `FALSE`
