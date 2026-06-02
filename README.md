# CRASH

## Cancer-associated gut micRobiome Analysis of Strains in Haematological patients

CRASH is a reproducible bioinformatics workflow developed for my MSc Bioinformatics thesis at Lund University. The project analyses publicly available whole-genome shotgun metagenomic datasets from haematological cancer patients and healthy controls.

The workflow was designed to study gut microbiome dysbiosis across multiple analysis layers:

- species-level taxonomic composition
- strain-level phylogenetic structure
- functional potential
- antimicrobial resistance gene profiles
- cohort-wise and cross-cohort statistical patterns

The main biological question is whether haematological cancer is associated with reproducible gut microbiome changes across independent cohorts, and whether these changes are reflected at the taxonomic, functional, AMR, and strain levels.

## Project overview

The workflow was applied to four analytical cohorts:

- Belgian AML cohort
- Polish AML cohort
- Polish lymphoma cohort
- Chinese NKTCL cohort

The pipeline combines several complementary tools:

- **KneadData** for read trimming and host decontamination
- **MetaPhlAn4** for marker-based species-level profiling
- **StrainPhlAn4** for strain-level phylogenetic analysis
- **HUMAnN4** for gene-family, KO, and pathway-level functional profiling
- **deepARG** for antimicrobial resistance gene profiling
- **Kraken2 + Bracken** as an independent taxonomic profiling check
- **R-based analysis** for diversity analysis, differential abundance testing, D-statistics, and figure generation

MetaPhlAn4 was used as the primary taxonomic profiling method. Kraken2/Bracken outputs were used as complementary checks and were not the main basis for the biological interpretation.

## Repository structure

```text
CRASH/
│
├── scripts/
│   ├── 00_setup/
│   │   ├── parallel_downloader.sh
│   │   ├── download_fastq.sh
│   │   ├── metaphlan_installer.sh
│   │   └── kraken_installer.sh
│   │
│   ├── 01_kneaddata/
│   │   └── QC, adapter trimming, and host decontamination scripts
│   │
│   ├── 02_metaphlan4/
│   │   ├── metaphlan.sh
│   │   └── humann_metaphlan.sh
│   │
│   ├── 03_strainphlan4/
│   │   ├── sample2markers.sh
│   │   ├── strainphlan_master.sh
│   │   ├── strainphlan_worker.sh
│   │   └── strainphlan.py
│   │
│   ├── 04_humann4/
│   │   └── humann.sh
│   │
│   ├── 05_kraken2_bracken/
│   │   ├── kraken.sh
│   │   └── kraken2mpa.sh
│   │
│   └── 06_deeparg/
│       └── deepARG profiling and post-processing scripts
│
├── R/
│   ├── alpha_diversity/
│   ├── beta_diversity/
│   ├── differential_abundance/
│   ├── functional_analysis/
│   ├── amr_analysis/
│   ├── strainphlan/
│   └── plotting/
│
├── environments/
│   └── Conda environment YAML files
│
├── metadata/
│   └── sample metadata and accession lists
│
├── docs/
│   └── pipeline notes, troubleshooting, and known issues
│
├── results/
│   └── output tables and figures
│
├── logs/
│   └── SLURM log files
│
└── README.md
```

Note: large intermediate files are not tracked by Git. This includes cleaned FASTQ files, Bowtie2 outputs, HUMAnN intermediate files, Kraken reports, and large temporary mapping files.

## Installation

Clone the repository:

```bash
git clone https://github.com/vikkxtar123/CRASH.git
cd CRASH
```

## Conda environments

The workflow uses separate Conda environments for different pipeline modules. Environment YAML files are stored in the `environments/` directory.

```bash
conda env create -f environments/kneaddata.yml
conda env create -f environments/metaphlan.yml
conda env create -f environments/krakentools.yml
conda env create -f environments/humann4.yml
conda env create -f environments/deeparg.yml
conda env create -f environments/r_analysis.yml
```

Activate the required environment before running each pipeline step.

Example:

```bash
conda activate metaphlan
```

## Required databases

The pipeline requires databases for MetaPhlAn4, HUMAnN4, Kraken2, host decontamination, and deepARG.

### MetaPhlAn4 / StrainPhlAn4

This analysis used the `mpa_vJun23_CHOCOPhlAnSGB_202403` database. Install that specific version so results are reproducible:

```bash
metaphlan --install \
    --index mpa_vJun23_CHOCOPhlAnSGB_202403 \
    --bowtie2db <METAPHLAN_DB_DIR>
```

StrainPhlAn4 uses the same database; no separate download is required.

### HUMAnN4

The full ChocoPhlAn nucleotide database and the UniRef90 DIAMOND protein database were used:

```bash
humann_databases --download chocophlan full <HUMANN_DB_DIR>
humann_databases --download uniref uniref90_diamond <HUMANN_DB_DIR>
```

KO regrouping uses the utility mapping bundled with HUMAnN4 (`humann_regroup_table`); no extra download is required.

### Human host database for KneadData

```bash
kneaddata_database --download human_genome bowtie2 <KNEADDATA_DB_DIR>
```

### Kraken2

This analysis used the prebuilt **PlusPF** index (archaea, bacteria, viral, plasmid, human, UniVec_Core, protozoa, and fungi). Download the prebuilt index rather than building a custom database, so the Kraken2/Bracken validation results are reproducible:

```bash
mkdir -p <KRAKEN_DB_DIR> && cd <KRAKEN_DB_DIR>
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_<YYYYMMDD>.tar.gz
tar -xzf k2_pluspf_<YYYYMMDD>.tar.gz
```

Replace `<YYYYMMDD>` with the index build date used in your analysis. Available versions are listed at https://benlangmead.github.io/aws-indexes/k2.

### deepARG

Install or download the deepARG database according to the deepARG documentation. The path to the database should be updated in the deepARG scripts before running.

## Input data

Raw FASTQ files were downloaded from public repositories using accession lists stored under `metadata/` or `scripts/00_setup/`.

Example:

```bash
sbatch scripts/00_setup/parallel_downloader.sh
```

The download script uses GNU parallel and is designed to download multiple SRA/ENA accessions efficiently.

Before running, edit the script to point to the correct accession list and output directory.

## Pipeline usage

Most scripts were written for the LUNARC COSMOS cluster and use SLURM. The general design is a worker-farm pattern, usually using one 48-core node split into multiple parallel workers.

Many scripts include skip logic, so interrupted jobs can usually be resubmitted without rerunning completed samples.

## 1. Quality control

Quality control was performed using KneadData.

Main steps:

- adapter trimming
- low-quality read removal
- human host read removal
- output of cleaned paired-end FASTQ files

Run:

```bash
sbatch scripts/01_kneaddata/kneaddata.sh
```

Expected outputs:

```text
cleaned_fastq/
kneaddata_logs/
fastqc_reports/
```

## 2. Taxonomic profiling with MetaPhlAn4

MetaPhlAn4 was used for the primary species-level taxonomic profiling.

```bash
sbatch scripts/02_metaphlan4/metaphlan.sh
```

This step produces:

- per-sample MetaPhlAn profiles
- Bowtie2 mapping outputs
- marker files for StrainPhlAn
- merged species-level abundance tables

The script also runs `sample2markers.py` where applicable, so the MetaPhlAn mapping output can be reused for strain-level analysis.

## 3. Strain-level profiling with StrainPhlAn4

If marker files were already generated during MetaPhlAn profiling, run the StrainPhlAn workflow directly:

```bash
sbatch scripts/03_strainphlan4/strainphlan_master.sh
```

If starting from existing Bowtie2 mapouts, first generate marker files:

```bash
sbatch scripts/03_strainphlan4/sample2markers.sh
```

The StrainPhlAn master script:

- extracts eligible clades
- maps SGB names to GTDB/species names where possible
- dispatches one worker per clade
- builds species-level strain phylogenies
- logs failed clades

Summarise tree outputs:

```bash
python scripts/03_strainphlan4/strainphlan.py \
    --root /path/to/strainphlan/output \
    --names /path/to/output_sgb_names.tsv \
    --out strainphlan_summary.tsv
```

Expected outputs:

```text
strainphlan_trees/
strainphlan_summary.tsv
failure_logs/
```

## 4. Functional profiling with HUMAnN4

HUMAnN4 was used to profile gene families, KO features, and MetaCyc pathways.

```bash
sbatch scripts/04_humann4/humann.sh
```

The workflow uses precomputed MetaPhlAn profiles as taxonomic input where possible.

Expected outputs:

```text
genefamilies.tsv
pathabundance.tsv
pathcoverage.tsv
KO_tables/
MetaCyc_pathway_tables/
```

Post-processing scripts were used to:

- join per-sample tables
- normalise gene-family and KO tables to CPM
- convert pathway outputs to relative abundance
- apply prevalence filtering
- prepare input tables for R analysis

## 5. Antimicrobial resistance profiling with deepARG

deepARG was used to profile antimicrobial resistance gene classes from metagenomic data.

Run:

```bash
sbatch scripts/06_deeparg/deeparg.sh
```

Expected outputs:

```text
deeparg_raw_outputs/
ARG_type_tables/
ARG_subtype_tables/
deeparg_summary_tables/
```

The final ARG tables were used for:

- ARG alpha diversity
- ARG beta diversity
- differential abundance testing
- cross-cohort consensus analysis

## 6. Independent taxonomic profiling with Kraken2 + Bracken

Kraken2 and Bracken were used as an independent taxonomic profiling check.

```bash
sbatch scripts/05_kraken2_bracken/kraken.sh
```

Parse outputs into MetaPhlAn-style tables:

```bash
sbatch scripts/05_kraken2_bracken/kraken2mpa.sh
```

These outputs were used for comparison and quality checking, but MetaPhlAn4 was used as the primary taxonomic profiling method in the thesis.

## 7. Statistical analysis in R

R scripts are stored in the `R/` directory.

The main statistical analyses include:

- alpha diversity testing
- Bray-Curtis and Jaccard beta diversity
- PERMANOVA using `adonis2`
- dispersion testing using `betadisper`
- differential abundance analysis using MaAsLin2
- supporting Wilcoxon rank-sum testing
- supporting ALDEx2 analysis
- Fritz and Purvis D-statistics for strain-level clustering
- cross-cohort consensus analysis
- figure generation

The main differential abundance framework used three methods:

1. MaAsLin2 as the anchor model
2. Wilcoxon rank-sum test as a non-parametric support method
3. ALDEx2 as a compositional support method

High-confidence features were defined based on MaAsLin2 significance and directional support from at least one supporting method.

## Main output files

The expected final outputs include:

```text
results/
├── taxonomic_profiles/
├── functional_profiles/
├── amr_profiles/
├── strainphlan/
├── differential_abundance/
├── alpha_diversity/
├── beta_diversity/
├── consensus_tables/
└── figures/
```

Main figure types generated by the workflow:

- cohort overview plots
- alpha diversity boxplots
- beta diversity PCoA plots
- taxonomic differential abundance barplots
- consensus heatmaps
- functional KO/pathway plots
- ARG class enrichment plots
- strain-level D-statistic summaries

## Reproducibility notes

This repository is intended to make the analysis workflow transparent and reproducible. The scripts are organised by analysis step and can be rerun from public sequencing accessions.

However, some paths are specific to the LUNARC COSMOS cluster and may need to be edited before running elsewhere. These include:

- project storage paths
- database paths
- SLURM account names
- temporary directory paths
- Conda environment paths

Large intermediate files are excluded from Git because of storage size. These files can be regenerated using the provided scripts, accession lists, and database setup instructions.

## Data availability

The sequencing datasets used in this project are publicly available under the following accessions:

- **PRJNA813705** — Belgian AML cohort (NCBI SRA / ENA)
- **PRJNA1116523 and PRJNA885289** — Polish AML and lymphoma cohorts (NCBI SRA / ENA)
- **CRA007433** — Chinese NKTCL cohort (NGDC Genome Sequence Archive)

The repository contains:

- scripts
- environment files
- post-processing scripts
- R analysis scripts
- figure-generation scripts

The repository does not contain:

- raw FASTQ files
- cleaned FASTQ files
- large Bowtie2 mapout files
- HUMAnN intermediate outputs
- Kraken2 databases
- MetaPhlAn/HUMAnN/deepARG databases

## Thesis context

This repository was developed as part of the MSc thesis:

**Multi-cohort metagenomic profiling of gut microbiome dysbiosis in haematological cancer patients**

MSc Bioinformatics  
Lund University

The workflow was developed to support reproducible analysis of taxonomic, functional, antimicrobial resistance, and strain-level patterns in gut metagenomes from haematological cancer cohorts.

## Notes

This pipeline was developed for research and thesis purposes. It is not intended as a clinical diagnostic tool.

Some scripts may require manual path editing before use on a different system. The workflow is most directly reproducible on an HPC system with SLURM and Conda available.

## Citation

If using or adapting this workflow, please cite the GitHub repository and the associated MSc thesis.

```text
Natarajan J. (2026). Multi-cohort metagenomic profiling of gut microbiome dysbiosis in haematological cancer patients. MSc Thesis, Lund University.
```

## Contact

For questions about the workflow, open an issue on the GitHub repository or contact the repository owner.
