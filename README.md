# CRASH

## Cancer-associated gut micRobiome Analysis of Strains in Haematological patients

CRASH is a reproducible bioinformatics workflow developed for the MSc Bioinformatics thesis at Lund University. The project analyses publicly available whole-genome shotgun metagenomic datasets from haematological cancer patients and healthy controls.

The workflow studies gut microbiome dysbiosis across multiple analysis layers:

- species-level taxonomic composition
- strain-level phylogenetic structure
- functional potential
- antimicrobial resistance gene profiles
- cohort-wise and cross-cohort statistical patterns

The main biological question is if haematological cancer is associated with reproducible gut microbiome changes across independent cohorts, and whether these changes are reflected at the taxonomic, functional, AMR, and strain levels.

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
- **Kraken2 + Bracken** as an independent, complementary taxonomic profiling check
- **R-based analysis** for diversity analysis, differential abundance testing, D-statistics, and figure generation

MetaPhlAn4 was used as the primary taxonomic profiling method. Kraken2/Bracken outputs were used as complementary checks and were not the main basis for the biological interpretation.

## Repository structure

All HPC scripts are in `scripts/` and source the shared `config.sh`. The R analysis lives flat in `R/`. Large intermediate files are not tracked by Git.

```text
CRASH/
│
├── scripts/                    # HPC pipeline (SLURM); each script sources config.sh
│   ├── config.sh               # dataset registry, DB paths, module/scratch helpers
│   ├── kneaddata.sh            # QC, trimming, human host decontamination
│   ├── metaphlan_farm.sh       # MetaPhlAn4 profiling + sample2markers (strain input)
│   ├── humann_meta.sh          # MetaPhlAn (Oct22 index) profiles as HUMAnN input
│   ├── humann4.sh              # HUMAnN4 functional profiling
│   ├── humann_utility.sh       # HUMAnN4 post-processing (join, regroup KO, normalise)
│   ├── strainphlan_all.sh      # multi-cohort StrainPhlAn run (4 cohorts pooled)
│   ├── strainphlan_master.sh   # generic StrainPhlAn master/worker dispatcher
│   ├── strainphlan_worker.sh   # per-clade StrainPhlAn worker
│   ├── strainphlan.py          # summarise strain trees into a summary table
│   ├── kraken.sh               # Kraken2 + Bracken farm (validation)
│   ├── kraken2mpa.sh           # convert Bracken output to MetaPhlAn-style tables
│   ├── deeparg.sh              # deepARG farm (AMR profiling)
│   └── deeparg_worker.sh       # per-sample deepARG worker
│
├── R/                          # statistical analysis and figures (R 4.4.0)
│   ├── MetaAnalysis_1.R        # MetaPhlAn per-cohort: MaAsLin2 + Wilcoxon + ALDEx2 + diversity
│   ├── MetaAnalysis_2.R        # MetaPhlAn cross-cohort consensus overlap
│   ├── MetaAnalysis_3.R        # MetaPhlAn cross-cohort summary figures
│   ├── MetaAnalysis_Tables.R   # thesis tables (cohort summary, consensus, diversity, DA counts)
│   ├── KrakenAnalysis_1.R      # Kraken2/Bracken per-cohort DA + diversity
│   ├── KrakenAnalysis_2.R      # Kraken cross-cohort overlap
│   ├── KrakenAnalysis_3.R      # Kraken figures
│   ├── DeepARG_1.R             # deepARG per-cohort DA (type/subtype) + Fisher + diversity
│   ├── DeepARG_2.R             # deepARG cross-cohort overlap
│   ├── DeepARG_3.R             # deepARG figures
│   ├── humann_figures.R        # KO heatmap with KEGG enzyme annotations
│   └── Strainphlan.R           # StrainPhlAn D-test (MODE: AML / Lymphoma / all)
│
├── environments/               # Conda environment YAML files
├── metadata/                   # Metadata for all cohorts
└── README.md
```

Excluded from Git (regenerable from the scripts and databases): cleaned FASTQ files, Bowtie2 mapouts, HUMAnN intermediate files, Kraken reports, and large temporary mapping files.

## Installation

```bash
git clone https://github.com/vikkxtar123/CRASH.git
cd CRASH
```

## Conda environments

Each pipeline module uses its own Conda environment. The environment names referenced by the scripts are `kneaddata`, `metaphlan`, `humann4`, `kraken2`, `krakentools`, and `deeparg_env`.

```bash
conda env create -f environments/kneaddata.yml       # kneaddata
conda env create -f environments/metaphlan.yml       # metaphlan  (MetaPhlAn4 + StrainPhlAn4)
conda env create -f environments/humann4.yml         # humann4
conda env create -f environments/kraken2.yml         # kraken2
conda env create -f environments/krakentools.yml     # krakentools (kraken2mpa)
conda env create -f environments/deeparg.yml         # deeparg_env (Python 2.7)
```

The R analysis is **not** run through Conda. R scripts were run locally on Windows (R 4.4.0) and on the cluster via `module load R/4.4.0`.

## Dataset keys

The pipeline scripts take a dataset key, resolved to a storage path by `resolve_dataset()` in `config.sh`. The four keys are:

| Key | Cohort |
| --- | --- |
| `prjna813705` | Belgian AML |
| `leukemia` | Polish AML (Kulecka) |
| `lymphoma` | Polish lymphoma (Kulecka) |
| `cra007433` | Chinese NKTCL |

`config.sh` also holds the SLURM account, shared database paths, and the `load_miniforge` / `make_scratch` / `register_outputs` helpers. Edit these paths before running on a different system.

## Required databases

The pipeline requires databases for MetaPhlAn4 (two versions), HUMAnN4, Kraken2, host decontamination, and deepARG.

### MetaPhlAn4 / StrainPhlAn4 (primary taxonomy + strain)

Primary profiling and StrainPhlAn4 use the `mpa_vJun23_CHOCOPhlAnSGB_202403` database:

```bash
metaphlan --install \
    --index mpa_vJun23_CHOCOPhlAnSGB_202403 \
    --bowtie2db <METAPHLAN_DB_DIR>
```

### MetaPhlAn4 (HUMAnN4 taxonomic prescreen)

`humann_meta.sh` regenerates MetaPhlAn profiles with the HUMAnN4-compatible `mpa_vOct22_CHOCOPhlAnSGB_202403`, which HUMAnN4 then uses:

```bash
metaphlan --install \
    --index mpa_vOct22_CHOCOPhlAnSGB_202403 \
    --bowtie2db <METAPHLAN_OCT22_DB_DIR>
```

### HUMAnN4

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

Install or download the deepARG database according to the deepARG documentation, and set the `DEEPARG_DB` path in `deeparg.sh` before running.

## Input data

Raw FASTQ files were downloaded from public repositories using the accession lists from the projects. The download step is independent of the analysis pipeline; any standard fetcher (SRA Toolkit `prefetch`/`fasterq-dump`, ENA, or a GNU parallel batch wrapper) can be used. Place paired-end FASTQ files in each cohort's `FASTQ/` directory before running QC.

Note the naming convention: NCBI cohorts use `*_1.fastq.gz` / `*_2.fastq.gz`; the CRA cohort uses `*_f1.fq.gz` / `*_r2.fq.gz`. `kneaddata.sh` handles both automatically.

## Pipeline usage

Scripts were written for the LUNARC COSMOS cluster and use SLURM, typically a worker-farm pattern on one 48-core node split into parallel workers. Most scripts include skip logic, so interrupted jobs can usually be resubmitted without rerunning completed samples.

### 1. Quality control (KneadData)

```bash
sbatch scripts/kneaddata.sh <dataset> [SRR1 SRR2 ...]
```

Adapter trimming, low-quality read removal, and human host removal. Optional sample IDs restrict the run to a subset. Outputs go to `<dataset>/KNEADDATA/`.

### 2. Taxonomic profiling (MetaPhlAn4)

```bash
sbatch scripts/metaphlan_farm.sh <dataset> [SRR1 SRR2 ...]
```

Runs MetaPhlAn4 (`mpa_vJun23` index) and `sample2markers` in one farm, producing per-sample profiles, Bowtie2 outputs, and consensus marker files. The consensus markers are reused directly by StrainPhlAn.

### 3. Strain-level profiling (StrainPhlAn4)

Run the multi-cohort StrainPhlAn analysis, which pools the consensus markers from all four cohorts:

```bash
sbatch scripts/strainphlan_all.sh
```

`strainphlan_master.sh` is the generic master/worker variant if you want to dispatch a custom marker set. Either way, workers (`strainphlan_worker.sh`) build per-clade phylogenies and log failed clades.

Summarise the trees:

```bash
python scripts/strainphlan.py \
    --root /path/to/strainphlan/output \
    --names /path/to/output_sgb_names.tsv \
    --out strainphlan_summary.tsv
```

The downstream Fritz & Purvis D-test is run in `R/Strainphlan.R`, which separates cohorts by disease via `MODE` (`AML` = Belgian + Polish leukemia; `Lymphoma` = Polish lymphoma + NKTCL) and includes a dataset-of-origin D-test to flag geographic confounding.

### 4. Functional profiling (HUMAnN4)

HUMAnN4 needs HUMAnN-compatible MetaPhlAn profiles first, so run the prescreen, then HUMAnN4, then post-processing:

```bash
sbatch scripts/humann_meta.sh <dataset>     # MetaPhlAn Oct22 profiles for HUMAnN
sbatch scripts/humann4.sh <dataset>         # gene families, KO, pathways
sbatch scripts/humann_utility.sh <dataset>  # join, regroup to KO, normalise, prevalence filter
```

Post-processing flattens per-sample subfolders, joins tables, regroups gene families to KO, and normalises gene-family/KO tables to CPM and pathway tables to relative abundance, producing the input tables for R analysis.

### 5. Antimicrobial resistance profiling (deepARG)

deepARG takes a sample list (TSV) rather than a dataset key:

```bash
sbatch scripts/deeparg.sh <sample_list.tsv>
```

Workers (`deeparg_worker.sh`) profile each sample; outputs feed the ARG type/subtype tables used for ARG diversity, differential abundance, and cross-cohort consensus.

### 6. Independent taxonomic profiling (Kraken2 + Bracken)

```bash
sbatch scripts/kraken.sh <dataset>       # Kraken2 + Bracken farm
sbatch scripts/kraken2mpa.sh <dataset>   # parse into MetaPhlAn-style tables
```

Used for comparison and quality checking; MetaPhlAn4 remains the primary method.

### 7. Statistical analysis (R)

The R scripts follow a consistent three-script pattern per data layer:

- **Script 1** — per-cohort differential abundance and diversity. Set `DATASET` (and `LEVEL` for deepARG) at the top and source the file.
- **Script 2** — cross-cohort consensus overlap, fed by the high-confidence hits from Script 1.
- **Script 3** — cross-cohort summary figures.

`MetaAnalysis_*.R` covers the primary MetaPhlAn layer, `KrakenAnalysis_*.R` the Kraken validation layer, and `DeepARG_*.R` the AMR layer. `MetaAnalysis_Tables.R` builds the thesis tables (cohort summary, consensus species, diversity, DA triangulation counts, curated functional groups) and writes a combined `.xlsx`. `humann_figures.R` produces the functional figures. `Strainphlan.R` runs the strain-level D-test.

The main differential abundance framework uses three methods:

1. MaAsLin2 as the anchor model
2. Wilcoxon rank-sum test as a non-parametric support method
3. ALDEx2 as a compositional support method

High-confidence features require MaAsLin2 significance plus directional support from at least one other method.

## Reproducibility notes

The scripts are organised by analysis step and can be rerun from public sequencing accessions. Some settings are environment-specific and must be edited before running elsewhere:

- HPC paths, database paths, SLURM account, and scratch directories live in `scripts/config.sh`. Replace the paths on SLURM with your usage case, they have been commented or edited out from all the scripts.
- R scripts use local Windows paths and a `DATASET` / `LEVEL` / `MODE` variable set at the top of each file — edit these before sourcing.
- After uploading any Windows-edited shell script to the cluster, run `dos2unix` on it.

## Data availability

The sequencing datasets used in this project are publicly available under the following accessions:

- **PRJNA813705** — Belgian AML cohort (NCBI SRA / ENA)
- **PRJNA1116523** — Polish AML and lymphoma cohorts (NCBI SRA / ENA)
- **CRA007433** — Chinese NKTCL cohort (NGDC Genome Sequence Archive)

The repository contains scripts, environment files, post-processing scripts, R analysis scripts, and figure-generation scripts. It does not contain raw or cleaned FASTQ files, Bowtie2 mapouts, HUMAnN intermediate outputs, or any of the reference databases.

## Thesis context

This repository was developed as part of the MSc thesis:

**Taxonomic, functional, resistome and strain-level profiling of gut microbiomes in haematological cancer patients**

MSc Bioinformatics  
Lund University

## Notes

This pipeline was developed for research and thesis purposes. It is not intended as a clinical diagnostic tool. The workflow is most directly reproducible on an HPC system with SLURM and Conda available.

## Citation

If using or adapting this workflow, please cite the GitHub repository and the associated MSc thesis.

```text
Natarajan J. (2026). Multi-cohort metagenomic profiling of gut microbiome dysbiosis in haematological cancer patients. MSc Thesis, Lund University.
```

## Contact

For questions about the workflow, open an issue on the GitHub repository or contact the repository owner.
