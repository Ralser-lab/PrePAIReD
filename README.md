# Preprocessing Pipeline --- Nextflow

This repository contains a **Nextflow-based preprocessing workflow** for
proteomics data.

The pipeline is designed to run reproducibly using **Docker containers**
and requires only a small number of dependencies.

------------------------------------------------------------------------

# Table of Contents

-   [Requirements](#requirements)
-   [Installation](#installation)
    -   [Docker Installation](#docker-installation)
    -   [Java Installation](#java-installation)
    -   [Nextflow Installation](#nextflow-installation)
    -   [Git Installation](#git-installation)
-   [Download the Pipeline](#download-the-pipeline)
-   [Starting Docker](#starting-docker)
-   [Working Directory Setup](#working-directory-setup)
-   [Configuration Files](#configuration-files)
-   [Running the Pipeline](#running-the-pipeline)
-   [Resuming a Pipeline](#resuming-a-pipeline)
-   [Folder Structure Example](#folder-structure-example)
-   [Troubleshooting](#troubleshooting)
-   [Minimal Requirements Summary](#minimal-requirements-summary)
-   [Support](#support)

------------------------------------------------------------------------

# Requirements

You need:

-   Java (JDK 11 or newer recommended)
-   Docker
-   Nextflow
-   Git (optional but recommended)

[⬆ Back to top](#top)

------------------------------------------------------------------------

# Installation

## Docker Installation

Docker is required to run the workflow in a reproducible environment.

Official installation page:

https://docs.docker.com/get-docker/

You can install Docker via:

-   Docker Desktop (recommended for macOS and Windows)
-   Terminal installation (Linux)

------------------------------------------------------------------------

## macOS --- Docker Desktop

Installation guide:

https://docs.docker.com/desktop/setup/install/mac-install/

### Option 1 — Install via GUI

Steps:

1.  Download Docker Desktop for Mac.
2.  Install the `.dmg` file.
3.  Start Docker Desktop.

Verify:

``` bash
docker --version
```
### Option 2 — Install via Terminal (Homebrew)

```bash
brew install --cask docker
open /Applications/Docker.app
docker --version
```

[⬆ Back to top](#top)

------------------------------------------------------------------------
## Windows --- Docker Desktop (Recommended with WSL2)

Installation guide:

https://docs.docker.com/desktop/setup/install/windows-install/

### Step 1 --- Install WSL2

Open PowerShell as Administrator:

``` powershell
wsl --install
```

Reboot if requested.

### Step 2 --- Install Docker Desktop

1.  Install Docker Desktop.
2.  Open Settings → enable:

-   Use WSL 2 based engine
-   Integration with your Linux distribution

Verify inside WSL terminal:

``` bash
docker --version
```

[⬆ Back to top](#top)

------------------------------------------------------------------------

## Linux --- Docker via Terminal

Example for Ubuntu / Debian:

``` bash
sudo apt-get update
sudo apt-get install -y docker.io
sudo systemctl enable --now docker
```

Add user permissions:

``` bash
sudo usermod -aG docker $USER
newgrp docker
```

Verify:

``` bash
docker --version
```

[⬆ Back to top](#top)

------------------------------------------------------------------------

# Java Installation

Nextflow requires Java.

## macOS (Homebrew)

``` bash
brew install openjdk@17
```

## Linux

``` bash
sudo apt-get install -y openjdk-17-jre
```

## Windows (WSL)

Install inside WSL terminal:

``` bash
sudo apt-get install -y openjdk-17-jre
```

Verify:

``` bash
java -version
```

[⬆ Back to top](#top)

------------------------------------------------------------------------
# Nextflow Installation

Official documentation:  
https://www.nextflow.io/docs/latest/getstarted.html

Nextflow requires **Java (JDK 11 or newer)**.

---

# macOS Installation

Download Nextflow using curl:

```bash
curl -s https://get.nextflow.io | bash
```

Move the executable to your system path:

```bash
sudo mv nextflow /usr/local/bin/
```

Verify installation:

```bash
nextflow -version
```

[⬆ Back to top](#top)

---

# Linux Installation

Download Nextflow:

```bash
curl -s https://get.nextflow.io | bash
```

Move Nextflow to your system path:

```bash
sudo mv nextflow /usr/local/bin/
```

Verify installation:

```bash
nextflow -version
```

[⬆ Back to top](#top)

---

# Windows Installation (WSL2 Recommended)

Nextflow runs best in **Windows Subsystem for Linux (WSL)**.

---

## 1. Install WSL2

Open **PowerShell as Administrator** and run:

```powershell
wsl --install
```

Restart your computer if required.

---

## 2. Open the WSL Terminal (Ubuntu)

Install Java:

```bash
sudo apt-get update
sudo apt-get install -y openjdk-17-jre
```

Verify installation:

```bash
java -version
```

---

## 3. Install Nextflow

Download Nextflow:

```bash
curl -s https://get.nextflow.io | bash
```

Move Nextflow to the system path:

```bash
sudo mv nextflow /usr/local/bin/
```

Verify installation:

```bash
nextflow -version
```

---

# Verify Installation

Run the following command:

```bash
nextflow -version
```

If the installation was successful, the Nextflow version will be printed.

[⬆ Back to top](#top)

------------------------------------------------------------------------

# Git Installation

Optional but recommended.

``` bash
# macOS
brew install git

# Linux / WSL
sudo apt-get install -y git
```

[⬆ Back to top](#top)

------------------------------------------------------------------------

# Download the Pipeline

Clone the repository:

``` bash
git clone https://github.com/Ralser-lab/Preprocessing.git
cd Preprocessing
```

Or download manually via GitHub **Code → Download ZIP**.

[⬆ Back to top](#top)

------------------------------------------------------------------------

# Starting Docker

Docker must be running before starting the pipeline.

### macOS / Windows

Start Docker Desktop.

### Linux

``` bash
sudo systemctl start docker
```

Verify:

``` bash
docker info
```

### Docker Images

You **do not need to manually download images**.\
Nextflow automatically pulls the required Docker image when the pipeline
starts.

[⬆ Back to top](#top)

------------------------------------------------------------------------

# Working Directory Setup

Create a working directory where configuration files and results will be
stored.

Example:

``` bash
mkdir -p "/Users/**username**/Desktop/test pipeline"
```

⚠️ If your path contains spaces, always use quotes.

[⬆ Back to top](#top)

------------------------------------------------------------------------

# Configuration Files

Your working directory must contain the following files:

    test pipeline/
    ├── parameters.txt
    ├── config_evaluation.txt
    ├── Config_data.txt
    └── Config.txt

These files define:

  File                    Purpose
  ----------------------- --------------------------
  parameters.txt          Pipeline parameters
  config_evaluation.txt   Evaluation settings
  Config_data.txt         Input/output paths
  Config.txt              Preprocessing step order

[⬆ Back to top](#top)

------------------------------------------------------------------------

## Configure Data Paths

Edit `Config_data.txt`:

``` txt
WORKING_PATH = "/Users/**username**/Desktop/test pipeline"
ORIGINALDATA = "/Users/**username**/Desktop/Data/diann.tsv"
METADATA     = "/Users/**username**/Desktop/Data/Metadata.tsv"
```

Make sure:

-   Paths exist on your system
-   Use quotes if paths contain spaces

[⬆ Back to top](#top)

------------------------------------------------------------------------

## Configure Preprocessing Steps

Example `config.txt`:

``` txt
Run = new
order = 01_FILTERING_LONG_py, 02_FILTERING_WIDE_py, 03_NORMALISATION_py, 04_IMPUTATION_py, 05_CORRECTION_py, 06_SUMMARISATION_py
```
or in case of rerun
``` txt
Run = 02_FILTERING_WIDE_py
order = 07_NORMALISATION_py, 08_IMPUTATION_py, 09_CORRECTION_py, 10_SUMMARISATION_py
```

[⬆ Back to top](#top)

------------------------------------------------------------------------

# Running the Pipeline

Navigate to the working directory:

``` bash
cd "/Users/username/Desktop/test_pipeline"
```

Run pipeline:

``` bash
nextflow run /path/to/Preprocessing/nextflow.nf -profile docker
```

Example:

``` bash
nextflow run "/Users/username/.../workflow/nextflow.nf"
```

[⬆ Back to top](#top)

------------------------------------------------------------------------

# Resuming a Pipeline

Resume previous run:

``` bash
nextflow run nextflow.nf -resume
```

Resume with aggressive cleanup:

``` bash
nextflow run nextflow.nf -resume --cleanup_aggressive true
```

[⬆ Back to top](#top)

------------------------------------------------------------------------

# Folder Structure Example

    project/
    │
    ├── Preprocessing/        		# cloned repository
    │   └── nextflow.nf
    │
    └── test pipeline/        		# working directory
        ├── 01_FILTERING_LONG		# result folder for each module
		│   ├── data_01.tsv			# intermediate data
		│   ├── report.pdf			# report file include all report till this step
		│   ├── report_folder		# folder includes evaluation fiigure and the corresponding figure captions (all reported in report.pdf)
		│   └── step_log.txt		# log file for the module to trace errors
		├── 02_FILTERING_WIDE
		├── 03_IMPUTATION
		├── 04_SUMMARISATION
		├── Parameter_Files   		# parameter direcory
		├── parameters.txt
        ├── config_evaluation.txt
        ├── Config_data.txt
        └── Config.txt

[⬆ Back to top](#top)

------------------------------------------------------------------------

# Troubleshooting

## Paths with Spaces

Always wrap paths in quotes.

    "/Users/name/My Folder/data.tsv"

## Docker Permission Issues

    sudo usermod -aG docker $USER

Log out and log back in.

[⬆ Back to top](#top)

------------------------------------------------------------------------

## Windows Recommendation

Use WSL paths instead of Windows paths.

Good:

    /home/user/project/

Avoid:

    C:\Users\user\project\

------------------------------------------------------------------------

# Minimal Requirements Summary

You only need:

-   Docker
-   Java
-   Nextflow

Nothing else is required.

[⬆ Back to top](#top)

------------------------------------------------------------------------

# Support

If you encounter issues, please open a GitHub issue in this repository.

[⬆ Back to top](#top)

------------------------------------------------------------------------

IMPUTATION documentation

- loading data, initialize module
- log parameters
- replace zeros with NaNs in data, if switched on
- if no values are missing, skip imputation
- if hierarchical_imputation is activated (hierarchical_imputation==1), load (if imputation_table_fp existent and overwrite_imputation_procedures==0) or create and save the imputation procedures table (imputation_table_fp not existent or overwrite_imputation_procedures==1) specifying imputation strategies over hierarchies:
	- one can specify different hierarchies (e.g. superbatch, batch, plate), sorted from highest to lowest hierarchy and/or missingness thresholds
    - hierarchy instances need to be unique,  e.g. there cannot be plate "1" in several batches, or name them as "batch1_plate1" and "batch2_plate1" etc. by concatenating hierarchies
    - an imputation method can be specified for every missingness type (MAR, MNAR, ProtAbsence, ProtAbsence/MAR, ProtAbsence/MNAR), missingness threshold and hierarchy combination
    - e.g. select as hierarchies superbatch;superbatch_batch;superbatch_batch_plate and missingness_thresh 1_0.33_0.25_0.3 to impute for every plate in every batch in every superbatch with missingness up to 25% with one method (e.g. knn with k=15) and with missingness up to 0.3 with a different (or the same) method (e.g knn with k=5), impute over the whole batch with missingness up to 0.33 with a different method and in every superbatch all remaining missing values, again with a different method. Within every imputation step, one can additionally distinguish different imputation methods for all the different missingness types. Also, the parameters for the imputation methods can be adjusted in the imputation procedures table 
    - the table can be created before calling the module_imputation and manually adjusted with different parameters for advanced application, e.g. setting different missingness thresholds within a hierarchy (like batch1 with a differnt threshold than batch2) or e.g. different parameter settings like k-values for knn. The parameters missingness_imputation_thresholds, imputation_sample_hierarchies and imputation_method_* are to generate a first template and in most cases sufficient for the table generation. 
	 
- if eval_misingness_on_protein_presence create a mask specifying whether any precursor of a protein is present in a sample or not
- if eval_missingness_by_detection_proportion identifies hierarchy to evaluate missingness (e.g. over whole dataset or per batch) evaluates missingness (MAR or MNAR) by modelling precursor completeness over abundance and stores missingness type in another mask
- combines mask specifying missingness by protein presence with mask specifying missingness by detection proportion (two masks from previous steps); a precursor in a sample can fall into several categories and thus be assigned with several missingness types (e.g. missing by protein presence and MAR); if save_intermediate_output_for_QA this mask is saved
- if missingness types (by Protein presence, MAR, MNAR) should be saved into filter mask, assigns their values
- if hierarchical imputation, do hierarchical imputation:
	- hierarchies could be something like batch - superbatch or different conditions (healthy/diseased, growing conditions) or covariates (age, sex)
    - loop over all hierarchies (select samples of interest)
        - log number of samples in hierarchy level
        - loop over all missingness thresholds (select precursors of interest)
            - log number of precursors in hierarchy and missingness level
            - extract all missingness types present in precursor and sample subset
            - loop over missingness types
                - extract imputation method and imputation parameters for this hierarchy, missingness threshold and missingness type
                - remember missing values
                - apply imputation on precursor and sample subset, precursors and samples in given subset without any measured values not imputed
                - update data with imputed values for given missingness type, missingness threshold and hierachy samples
                - count number of still missing values
                - update mask, log number of imputed values

- if no hierarchical imputation or imputation over whole dataset (e.g. after hierarchical imputation) is switched on, apply imputation over complete dataset:
    - extract missingness types
    - loop over missingness types
        - define imputation method from parameters dict
        - remember which and how many values were missing before imputation
        - apply imputation
        - remember which and how many values are missing after imputation
        - log them and update mask
- count final number of missing and imputed values and log them
- close module
