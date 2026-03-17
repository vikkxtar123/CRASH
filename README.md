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


