# CRASH
## Strain-Level Profiling of Gut Microbiomes in Hematological Cancer - CRASH (Cancer-associated gut micRobiome Analysis of Strains in Haematological patients)

This is a multilayered metagenomics pipeline for taxonomic, strain-level, and functional profiling of the gut microbiome across haemotological cancer datasets. Developed as part of a MSc Bioinformatics thesis at Lund University.

This pipeline makes use of three complementary profiling approaches:
- MetaPhlAn4: Marker-based species-level taxonomic profiling
- StrainPhlAn4: Strain-level phylogenetic profiling
- HUMAnN4: Functional profiling (gene families and metabolic pathways)

An independent taxonomic profiling approach for cross-validation and depth: 
- Kraken2 + Bracken: k-mer based taxonomic classification.

GitHub repo can be found here.

## Installation

### Clone the repository:

```
git clone https://github.com/vikkxtar123/CRASH.git
cd CRASH
```
### Installing Anaconda
If you don't have Anaconda installed, you can download and install it by following these steps:

- Visit the official Anaconda website.
- Download the installer for your operating system (Windows, macOS, or Linux).
- Follow the installation instructions provided for your platform.
- After installation, you can verify that Anaconda is installed correctly by running:

```
conda --version
```
### Setting up the Conda Environment

To run the pipeline, multiple conda environments are required. Profiles of these environments are found under the ` environments/ ` directory. They can be installed using: 
```
conda env create -f environments/kneaddata.yml # For KneadData (QC)
conda env create -f environments/metaphlan.yml # For MetaPhlAn4 and StrainPhlAn4
conda env create -f environments/krakentools.yml # For Kraken2 and Bracken
conda env create -f environments/humann4.yml # For HUMAnN4
```

## Usage

At this point, the directory structure should be as follows:
 
```
CRASH/
│
├── scripts/
│   ├── 00_setup/
│   │   ├── parallel_downloader.sh       # Parallel FASTQ download via ENA (GNU parallel)
│   │   ├── download_fastq.sh            # Per-accession download worker (called by parallel_downloader.sh)
│   │   ├── metaphlan_installer.sh       # Install MetaPhlAn mpa_vJun23 database
│   │   └── kraken_installer.sh          # Download and install k2_pluspf Kraken2 database
│   ├── 01_kneaddata/                    # QC — adapter trimming and host decontamination
│   ├── 02_metaphlan4/
│   │   ├── metaphlan.sh                 # MetaPhlAn4 farm (Jun23 DB) + inline sample2markers
│   │   └── humann_metaphlan.sh          # MetaPhlAn4 farm (Oct22 DB, profiling only for HUMAnN4)
│   ├── 03_strainphlan4/
│   │   ├── sample2markers.sh            # Standalone sample2markers farm (from bowtie2 mapouts)
│   │   ├── strainphlan_master.sh        # Master farm — clade extraction + parallel tree building
│   │   ├── strainphlan_worker.sh        # Per-clade worker (called by master via srun)
│   │   └── strainphlan.py     # Summarise all clade outputs into a single TSV
│   ├── 04_humann4/
│   │   └── humann.sh                    # HUMAnN4 parallel farm (2 workers × 24 threads)
│   └── 05_kraken2_bracken/
│       ├── kraken.sh                    # Kraken2 + Bracken farm (6 workers × 8 threads)
│       └── kraken2mpa.sh               # Parse Bracken outputs into MetaPhlAn-style tables
├── R/
│   ├── beta_diversity/                  # PCoA, PERMANOVA, betadisper
│   ├── maaslin2/                        # Differential abundance analysis
│   └── strainphlan/                     # Phylogenetic D-test (caper::phylo.d)
├── envs/                                # Conda environment YAML files
├── docs/                                # Pipeline notes and known gotchas
├── results/                             # Output tables and figures (not tracked by git)
├── logs/                                # SLURM log files (not tracked by git)
└── README.md
```
Download the required databases:
 
```
# MetaPhlAn / StrainPhlAn
metaphlan --install --index mpa_vJun23_CHOCOPhlAnSGB_202403 --bowtie2db <DB_DIR>
 
# HUMAnN
humann_databases --download chocophlan full <DB_DIR>
humann_databases --download uniref uniref90_diamond <DB_DIR>
 
# Kraken2
kraken2-build --download-library ... --db <DB_DIR>/k2_pluspf
```

Download raw FASTQ data:
 
```
# Edit parallel_downloader.sh to point at your SRA ID list, then:
sbatch scripts/00_setup/parallel_downloader.sh
```
 
SRA ID lists for each dataset are provided in `scripts/00_setup/` as `<ACCESSION>_sra_ids.txt`.

### Usage
 
All SLURM scripts use a 6-worker farm pattern on a single 48-core COSMOS node unless otherwise noted. Per-sample skip logic is included — safe to resubmit if a job is interrupted.
 
### 1. Quality control
 
Quality control is performed with KneadData (Trimmomatic + human host removal). Scripts are in `scripts/01_kneaddata/`.
 
### 2. Taxonomic profiling
 
```
sbatch scripts/02_metaphlan4/metaphlan.sh
```
 
Runs MetaPhlAn4 against the `mpa_vJun23_CHOCOPhlAnSGB_202403` database using a 6-worker farm (6 × 8 threads). Saves `--bowtie2out` and runs `sample2markers.py` inline on the local SAM before cleanup — this avoids writing large SAM files to shared storage. Profiles, bowtie2 mapouts, and consensus marker `.json.bz2` files are all written to storage.
 
`humann_metaphlan.sh` is an alternative run against the Oct22 database (`mpa_vOct22_CHOCOPhlAnSGB_202403`) for cross-version comparison. It uses `--no_map` and does not produce alignment files.
 
### 3. Strain-level profiling
 
If `metaphlan.sh` was used, consensus markers are already produced. If starting from existing bowtie2 mapouts:
 
```
sbatch scripts/03_strainphlan4/sample2markers.sh
```
 
Then run the StrainPhlAn farm:
 
```
sbatch scripts/03_strainphlan4/strainphlan_master.sh
```
 
The master script extracts the clade list using `--print_clades_only`, assigns GTDB names via the SGB2GTDB mapping file, then dispatches one `strainphlan_worker.sh` per clade via `srun --exclusive`. Workers run `extract_markers.py` and `strainphlan` sequentially per clade. Failed clades are logged to `failures.log`.
 
After all trees are built, summarise results:
 
```
python scripts/03_strainphlan4/strainphlan.py \
    --root  /path/to/Output \
    --names /path/to/output_sgb_names.tsv \
    --out   strainphlan_summary.tsv
```
 
This produces a single TSV with per-clade stats: samples in tree, markers used, polymorphic rates, and species/genus names.
 
### 4. Functional profiling
 
```
sbatch scripts/04_humann4/humann.sh
```
 
Runs HUMAnN4 with 2 concurrent workers (24 threads each) on a single 48-core node. Samples are processed using pre-computed MetaPhlAn profiles (`--taxonomic-profile`) to skip re-profiling. Uses `--prescreen-threshold 1.0` to prevent parallel bowtie2-build workers from exhausting `/local/slurmtmp`. Intermediate files are written to `SNIC_TMP` and cleaned up after each sample.
 
### 5. Independent taxonomic profiling
 
```
sbatch scripts/05_kraken2_bracken/kraken.sh
```
 
Runs Kraken2 + Bracken (species and genus level) using a 6-worker farm. The Kraken2 database is staged once to local disk at job start to avoid repeated Lustre I/O across workers.
 
Parse outputs into MetaPhlAn-style tables:
 
```
sbatch scripts/05_kraken2_bracken/kraken2mpa.sh
```
 
This produces both Bracken-derived relative abundance tables and `kreport2mpa.py`-derived full-taxonomy MPA tables for direct comparison with MetaPhlAn4 output.
