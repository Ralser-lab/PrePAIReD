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
-   [Docker Environment Setup](#docker-environment-setup)
-   [Running the Pipeline](#running-the-pipeline)
-   [Resuming a Pipeline](#resuming-a-pipeline)
-   [Example Commands](#example-commands)
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

## macOS Installation

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

## Linux Installation

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

## Windows Installation (WSL2 Recommended)

Nextflow runs best in **Windows Subsystem for Linux (WSL)**.

---

### 1. Install WSL2

Open **PowerShell as Administrator** and run:

```powershell
wsl --install
```

Restart your computer if required.

---

### 2. Open the WSL Terminal (Ubuntu)

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

### 3. Install Nextflow

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

### Verify Installation

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

For an automatic setup of the remaining files run Initialization/parameter_generation.py (Only requirement: Python 3.9+).

An example `config.txt` output will look like:

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

# Docker Environment Setup

Before running the pipeline, the Docker image must be built locally from the `Dockerfile` in the repository.

First ensure that **Docker is installed and running**.

---
## Start Docker

### macOS

If you are using **Docker Desktop**, start it from Applications:

```bash
open /Applications/Docker.app
```

Verify Docker is running:

```bash
docker info
```

If you are using **Colima**, start it with:

```bash
colima start
```
Verify:

```bash
docker info
```

---

### Windows

Start **Docker Desktop** from the Start Menu.

Alternatively using PowerShell:

```powershell
Start-Process "Docker Desktop"
```

Verify Docker is running:

```bash
docker info
```

---

### Linux

Start the Docker service:

```bash
sudo systemctl start docker
```

Enable Docker to start automatically:

```bash
sudo systemctl enable docker
```

Verify Docker is running:

```bash
docker info
```

---

## Build the Image
Navigate to the repository root where the `Dockerfile` is located.

```bash
docker build --no-cache -t unified_env_all:1.0 .
```

Check that the image was created:

```bash
docker images
```

You should see something similar to:

```
REPOSITORY          TAG       IMAGE ID
unified_env_all     1.0       <image_id>
```

---

## Notes

- Docker must be running before starting the pipeline.
- The image only needs to be built once unless the `Dockerfile` changes.
- If the image is published to a container registry, users will not need to build it locally.

---

## Container Configuration in Nextflow

Local Docker image in **Dockerfile**:

``` groovy
process.container = 'unified_env_all:1.0'
docker.enabled = true
```

Remote container example:

``` groovy
process.container = 'ghcr.io/ralser-lab/preprocessing:1.0'
docker.enabled = true
```

Nextflow will automatically download the container if needed.

[⬆ Back to top](#top)

------------------------------------------------------------------------

# Running the Pipeline

Navigate to the working directory:

``` bash
cd "/Users/username/Desktop/test_pipeline"
```

Run pipeline:

``` bash
nextflow run /path/to/Preprocessing/nextflow.nf 
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


------------------------------------------------------------------------

# Example Commands

## Option 1: Use the prebuilt GitHub container image

Clone the repository and enter the directory:

```bash
git clone https://github.com/Ralser-lab/Preprocessing.git
cd Preprocessing
```

Check that `nextflow.config` contains:

```groovy
process.container = 'ghcr.io/ralser-lab/preprocessing:1.0'
```

Move to your working data directory:

```bash
cd /Users/username/Desktop/test_data
```

Prepare the required input files in this directory:

- `config_data.txt`  
  Include:
  - path to your input data
  - path to your metadata
  - working directory path

- `config.txt`
- `config_evaluation.txt`
- `Parameters.txt`

Run the pipeline using the absolute path to the workflow:

```bash
nextflow run "/Users/username/Desktop/Preprocessing/nextflow.nf" 
```

---

## Option 2: Build the Docker image locally

Clone the repository and enter the directory:

```bash
git clone https://github.com/Ralser-lab/Preprocessing.git
cd Preprocessing
```

Make sure [Docker is running](#starting-docker) before building the image. 

Build the Docker image locally:

```bash
docker build -t unified_env_all:1.0 .
```

Check that `nextflow.config` contains:

```groovy
process.container = 'unified_env_all:1.0'
```

Move to your working data directory:

```bash
cd /Users/username/Desktop/test_data
```

Prepare the required input files in this directory:

- `config_data.txt`  
  Include:
  - path to your input data
  - path to your metadata
  - working directory path

- `config.txt`
- `config_evaluation.txt`
- `Parameters.txt`

Run the pipeline using the absolute path to the workflow:

```bash
nextflow run "/Users/username/Desktop/Preprocessing/nextflow.nf" 
```

---

## Before running

Make sure the following are installed and available:

- Nextflow
- Docker
- Git

For the local-image option, also make sure the Docker application is started before running:

- **macOS / Windows:** open Docker Desktop and wait until it shows Docker is running
- **Linux:** start the Docker service if needed

Example check commands:

```bash
git --version
docker --version
nextflow -version
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
