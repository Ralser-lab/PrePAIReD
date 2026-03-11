# Preprocessing Pipeline --- Nextflow

![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A522.10-brightgreen)
![Docker](https://img.shields.io/badge/Container-Docker-blue)
![License](https://img.shields.io/badge/license-MIT-lightgrey)

This repository contains a **Nextflow-based preprocessing workflow for
proteomics data**.

The pipeline performs standardized preprocessing of proteomics results
including:

-   filtering
-   normalization
-   missing value imputation
-   batch correction
-   protein summarisation
-   evaluation reporting

The workflow runs reproducibly using **Docker containers**.

------------------------------------------------------------------------

# Table of Contents

-   Requirements
-   Installation
-   Download the Pipeline
-   Workflow Overview
-   Working Directory Setup
-   Configuration Files
-   Container Environment Setup
-   Running the Pipeline
-   Resuming a Pipeline
-   Example Commands
-   Parameter Overview
-   Folder Structure Example
-   Troubleshooting
-   Minimal Requirements Summary
-   Support
-   Imputation Module Documentation

------------------------------------------------------------------------

# Requirements

Required software:

-   Java (JDK ≥11)
-   Docker
-   Nextflow

Optional:

-   Git

------------------------------------------------------------------------
# Installation

## Install Docker

Official installation guide:\
https://docs.docker.com/get-docker/

### macOS

``` bash
brew install --cask docker
open /Applications/Docker.app
docker --version
```

### Windows (WSL2 recommended)

``` powershell
wsl --install
```

Then install **Docker Desktop**.

Verify:

``` bash
docker --version
```

### Linux

``` bash
sudo apt-get update
sudo apt-get install -y docker.io
sudo systemctl enable --now docker
docker --version
```

------------------------------------------------------------------------
## Install Java

Nextflow requires Java.

### macOS (Homebrew)

``` bash
brew install openjdk@17
```

### Linux

``` bash
sudo apt-get install -y openjdk-17-jre
```

### Windows (WSL)

Install inside WSL terminal:

``` bash
sudo apt-get install -y openjdk-17-jre
```

Verify:

``` bash
java -version
```

------------------------------------------------------------------------

## Install Nextflow

Official documentation:  
https://www.nextflow.io/docs/latest/getstarted.html

Nextflow requires **Java (JDK 11 or newer)**.

---

### macOS Installation

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
---

### Linux Installation

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
---

#### Windows Installation (WSL2 Recommended)

Nextflow runs best in **Windows Subsystem for Linux (WSL)**.

---

#### 1. Install WSL2

Open **PowerShell as Administrator** and run:

```powershell
wsl --install
```

Restart your computer if required.

---

#### 2. Open the WSL Terminal (Ubuntu)

Install Java:

```bash
sudo apt-get update
sudo apt-get install -y openjdk-17-jre
```

Verify installation:

```bash
java -version
```
------------------------------------------------------------------------

## Git Installation

Optional but recommended.

```bash
# macOS
brew install git

# Linux / WSL
sudo apt-get install -y git

# Windows (winget)
winget install --id Git.Git -e --source winget

Verify the installation by opening PowerShell, Command Prompt, or Git Bash and running:
git --version

If `winget` is not available, install Git manually:

1. Download Git for Windows: https://git-scm.com/download/win  
2. Run the downloaded `.exe` installer.  
3. Follow the setup wizard and keep the default settings unless you have a specific preference.  
4. Finish the installation.

------------------------------------------------------------------------

# Download the Pipeline

Clone the repository:

``` bash
git clone https://github.com/Ralser-lab/Preprocessing.git
cd Preprocessing
```

Or download manually from GitHub.

------------------------------------------------------------------------
# Workflow Overview

The pipeline executes various preprocessing modules, the following gives you one possible example of preprocessing modules:

    Raw Proteomics Table
            │
            ▼
    01_FILTERING_LONG
            │
            ▼
    02_FILTERING_WIDE
            │
            ▼
    03_NORMALISATION
            │
            ▼
    04_IMPUTATION
            │
            ▼
    05_CORRECTION
            │
            ▼
    06_SUMMARISATION
            │
            ▼
    Evaluation Reports

Each module produces:

-   processed data tables
-   evaluation figures
-   logs
-   summary reports

------------------------------------------------------------------------

# Working Directory Setup

Create a directory where the pipeline will run.

Example:

``` bash
mkdir -p "/Users/username/Desktop/test_pipeline"
```

------------------------------------------------------------------------

# Configuration Files

The working directory must contain:

    test_pipeline/
    ├── parameters.txt
    ├── config_evaluation.txt
    ├── Config_data.txt
    └── Config.txt

  File                    Description
  ----------------------- ---------------------
  parameters.txt          module parameters
  config_evaluation.txt   evaluation settings
  Config_data.txt         input/output paths
  Config.txt              preprocessing order

------------------------------------------------------------------------

## Configure Data Paths

Example `Config_data.txt`:

``` txt
WORKING_PATH = "/Users/username/Desktop/test_pipeline"
ORIGINALDATA = "/Users/username/Desktop/Data/diann.tsv"
METADATA     = "/Users/username/Desktop/Data/Metadata.tsv"
```

------------------------------------------------------------------------
