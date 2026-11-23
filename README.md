![License](https://img.shields.io/badge/License-MIT-blue.svg)
![Python](https://img.shields.io/badge/Python-3.12-green.svg)
![R](https://img.shields.io/badge/R-4.4-276DC3.svg?logo=r&logoColor=white&v=2)
![Nextflow](https://img.shields.io/badge/Nextflow-26.0-3ac486.svg?logo=nextflow&logoColor=white)
![Version](https://img.shields.io/badge/Version-0.1.3-orange.svg)
[![Docs](https://img.shields.io/badge/Docs-PDF-blue.svg)](docs/Manual_PrePAIReD.pdf)
![Last Commit](https://img.shields.io/github/last-commit/Ralser-lab/Preprocessing)
![Status](https://img.shields.io/badge/Status-Active-brightgreen.svg)

# Preprocessing Pipeline

This repository contains a **Nextflow-based preprocessing workflow** for MS-based bottom-up proteomics data. It is designed to run reproducibly using **Docker containers** and requires only a small number of dependencies (see [Requirements](#requirements)). 

The provided modules also allow for manual execution of every module or the use of single functions within a module (see [Environment setup](environment.yml)).

For a more detailed description and other use cases please check out the [Manual](docs/Manual_alpha_version.pdf).

------------------------------------------------------------------------

# Table of Contents

-   [Requirements](#requirements)
-   [Download the Pipeline](#download-the-pipeline)
-   [Installation](#installation)
    -   [Docker Installation](#docker-installation)
    -   [Java Installation](#java-installation)
    -   [Nextflow Installation](#nextflow-installation)
    -   [Git Installation](#git-installation)
-   [Docker Environment Setup](#docker-environment-setup)
-   [Running the Pipeline](#running-the-pipeline)
-   [Rerunning from a Step](#rerunning-from-a-step)
-   [Example Commands](#example-commands)
-   [Folder Structure Example](#folder-structure-example)
-   [Run History and Reproducibility](#run-history-and-reproducibility)
-   [Troubleshooting](#troubleshooting)
-   [Minimal Requirements Summary](#minimal-requirements-summary)
-   [Support](#support)

------------------------------------------------------------------------

# Requirements

For a fully automated run you need:

-   Java (JDK 17 or later)
-   Docker
-   Nextflow
-   Git (optional but recommended)

For manual use no further software is required (see [Environment setup](environment.yml) for the required packages)

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

Nextflow requires **Java (JDK 17 or later)**.

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

# Configuration and Parameter Files

Your working directory must contain the following files:

    test pipeline/
    ├── Parameters.txt
    ├── Config_evaluation.txt
    ├── Config_data.txt
    └── Config.txt

All files are automatically generated by Initialization/gui.py. A detailed description can be found in the [Manual](docs/Manual_alpha_version.pdf).

------------------------------------------------------------------------

# Docker Environment Setup

First ensure that **Docker is installed and running**

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

## Notes

- Docker must be running before starting the pipeline.
- The image only needs to be built once unless the `Dockerfile` changes.
- If the image is published to a container registry, users will not need to build it locally.

---

## Container Configuration in Nextflow

Before running the pipeline, choose the prebuilt registry image or build a local image from the repository’s `Dockerfile`. 

Local image setting in `nextflow.config`:

``` groovy
params.container_image = 'unified_env_all:1.0'
docker.enabled = true
```

Remote container example:

``` groovy
params.container_image = 'ghcr.io/ralser-lab/preprocessing:1.0'
docker.enabled = true
```

Nextflow can pull the registry image if it is missing locally. On that first run, the startup log may show an unknown image ID because the workflow checks the local image before the pull.

### Updating the environment

Update the build inputs that changed, such as `environment.yml` or `Dockerfile`, then build from the repository root under a new tag.

```bash
docker build -t unified_env_all:1.1 .
```

Update `params.container_image`, or pass `--container_image unified_env_all:1.1` for one run. If distributing the image through a registry, tag and push it there and use the registry reference. Keep the exact image or its registry digest for reproducibility.

[⬆ Back to top](#top)

------------------------------------------------------------------------

# Running the Pipeline

Launch from the working directory and point Nextflow to the workflow in the repository.

```bash
cd "PATH/TO/WORKDIR/"
nextflow run "PATH/TO/PrePAIReD/nextflow.nf"
```

To launch elsewhere, pass `--config_data /absolute/path/to/Config_data.txt` and set an absolute `WORKING_PATH` in that file. A relative `WORKING_PATH` resolves from the launch directory. Relative data, metadata, and optional imputation procedures paths resolve from `WORKING_PATH`.

Docker is enabled in the supplied configuration. The optional `-profile linux_docker_user` runs the container under your host user and group IDs.

### Optional resource settings

```bash
nextflow run "PATH/TO/PrePAIReD/nextflow.nf" \
  --memory '8 GB' --cpus 4 --max_memory '24 GB'
```

- `--memory` is the initial memory request and defaults to 4 GB. `--cpus` is the common CPU request and defaults to 4.
- `--module_memory` and `--module_cpus` override the requests for the main processing task. `--init_memory` and `--init_cpus` apply to initialisation, which defaults to one CPU. The environment-log task uses 1 GB and one CPU.
- `--max_memory` sets the retry ceiling. By default, the ceiling is detected Docker memory minus `--memory_margin` (2 GB by default), subject to availability and the minimum in the config. If Docker memory cannot be detected and no ceiling is given, the initial request is used as the ceiling, subject to that minimum.

### Freeing disk space

`--cleanup_aggressive true` frees disk space by removing Nextflow’s `work/` directory, the launch directory’s `.nextflow/` and `.nextflow.log` when the run ends, including after a failure. The numbered step folders remain, so reruns from a step still work. Leave it off if you want to inspect the Nextflow files of a failed task.

### Retries

With the supplied `nextflow.config`, a task that exits with 137, 139 or 143 gets up to three retries with more memory, within the memory ceiling. These exit codes can have causes other than memory exhaustion. A retry of the main processing task starts again at the first entry in `Order`, not at the step that failed.

> **Developer note:** the `RUN_ALL_MODULES` selector in `nextflow.config` overrides the inline 137/143 rule in `nextflow.nf`. The initialisation task uses the same retry codes, but a nonzero exit from `Initialize.py` is converted to task exit 100, so it is not retried.

[⬆ Back to top](#top)

------------------------------------------------------------------------

# Rerunning from a Step

To try other settings without repeating earlier steps, start a new branch from a saved step folder with `Run = <seed step>`.

List the new steps in `Order` and use distinct step numbers to keep both sets of step outputs.

The workflow writes a `step_status.txt` into each step folder when it can. It is set to `RUNNING` when a step starts and changes to `COMPLETED` only after the workflow has checked the seed files it published: the numbered table, `id_mapping.tsv`, `step_log.txt`, the cumulative `report.pdf`, the masks the step produced, and that no leftover masks from an older run remain. It also records the Nextflow run name and the time. When a seed folder has a status file, it must say `COMPLETED`; the rerun prints the seed's run name and time when it starts. It also requires the seed folder, its numbered data table and `id_mapping.tsv`.

- If the start marker was written, a step that failed or was interrupted stays `RUNNING` and is refused, even if older files from a previous run remain in the folder.
- If a step finished but its published files could not be verified, its status is set to `INCOMPLETE` (if the status update succeeds) and `pipeline_log.txt` shows a `[STEP STATUS WARNING]`. The run itself continues.
- If `RUNNING` cannot be written at the start of a step, that step publishes nothing into its folder and a warning appears in `pipeline_log.txt`. The overall run continues, but that step folder still holds a previous run’s outputs (with their own status and run name), or nothing.
- A folder whose log was written by this version but whose `step_status.txt` is missing is refused.

**Folders from older pipeline versions** have no `step_status.txt`. For those, the workflow checks `step_log.txt` if it exists. It refuses the seed if the log shows a failed module or report (QA problems do not count, because they are non-fatal). Otherwise, or if there is no log at all, it continues with a warning, so such a folder may be accepted without proof that it is complete.

Before branching, confirm that the seed is from the run you intend to reuse. `COMPLETED` shows the step finished, but in a reused working folder it can come from an earlier run if a later run stopped before reaching that step. Check the run name and time against the run you want, then review `qa_status.tsv` for warnings. If you can’t establish which run produced the files, use a seed from a separate run folder.

```bash
cat "02_FILTERING_WIDE/step_status.txt"
```

[⬆ Back to top](#top)

------------------------------------------------------------------------

# Example Commands

## Option 1: Use the prebuilt GitHub container image

Clone the repository and enter the directory:

```bash
git clone https://github.com/Ralser-lab/PrePAIReD.git
cd PrePAIReD
```

Check that `nextflow.config` contains:

```groovy
params.container_image = 'ghcr.io/ralser-lab/preprocessing:1.0'
```
## Option 2: Build the Docker image locally

Clone the repository and enter the directory:

```bash
git clone https://github.com/Ralser-lab/PrePAIReD.git
cd PrePAIReD
```

Make sure [Docker is running](#docker-environment-setup) before building the image. 

Build the Docker image locally:

```bash
docker build -t unified_env_all:1.0 .
```

Check that `nextflow.config` contains:

```groovy
params.container_image = 'unified_env_all:1.0'
```

## In both cases

Move to your working data directory:

```bash
cd "PATH/TO/WORKDIR/"
```

Prepare the required input files in this directory running:

```bash
python gui.py
```

Note: As of now, the GUI requires a minimal environment setup as described in [Manual](docs/Manual_alpha_version.pdf).
Run the pipeline using the absolute path to the workflow:

```bash
nextflow run "PATH/TO/PrePAIReD/nextflow.nf"
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

This example launches from `WORKING_PATH`, uses the default `--outdir`, and runs the six steps shown above.

```text
test pipeline/
├── Config_data.txt, Config.txt, Config_evaluation.txt, Parameters.txt
├── 01_FILTERING_LONG/
│   ├── data_01.tsv
│   ├── data_01_mask.tsv, data_01_mask_description.tsv
│   ├── id_mapping.tsv
│   ├── report.pdf, report_folder/
│   ├── step_log.txt, log.txt, qa_status.tsv
│   └── step_status.txt
├── 02_FILTERING_WIDE/
├── 03_NORMALISATION/
├── 04_IMPUTATION/
├── 05_CORRECTION/
├── 06_SUMMARISATION/
├── Parameter_Files/, logfile_setup.txt
├── report.pdf, pipeline_log.txt
├── Metadata_generated.tsv            (only if no METADATA was given)
├── env_log.txt, env_log.pdf
├── pipeline.done
├── run_history/, previous_reports/
├── nextflow_report.html, nextflow_timeline.html, nextflow_trace.txt
└── work/, .nextflow/
```

Some step folders hold one extra file: `Normalization_factors.tsv` (normalisation), `imputation_procedures.tsv` (imputation) and `covariate_associations.tsv` (correction, when confounders are adjusted).

Step folders contain the files produced and copied for that step. `qa_status.tsv` is cumulative through that step. The figures folder is per step, and the PDF is carried forward. The environment PDF is best effort. `work/` and the three Nextflow reports go under `--outdir`, which defaults to the launch directory. `.nextflow/` is in the launch directory.

**Before starting another run in the same folder:** a new run can replace files in matching step folders and leave older files behind. The workflow moves the existing nonempty top-level `report.pdf` into `previous_reports/`. Use a fresh analysis folder or save the prior outputs separately to preserve the complete result.

[⬆ Back to top](#top)

------------------------------------------------------------------------

# Run History and Reproducibility

After a successful run, look for its record in `run_history/<timestamp>_newrun/` or `run_history/<timestamp>_rerun_from_<seed-folder>/`. A run can complete with QA warnings.

The record can contain startup copies of the four settings files and any supplied imputation procedures file, plus `run_info.txt`, environment logs, generated parameters, setup and pipeline logs, and the available final PDF and trace. If you see `[RUN HISTORY WARNING]`, the record may be missing or incomplete.

The history records SHA-256 checksums of the original data and metadata, and of the seed table for a rerun. Keep the input data, output tables, figures, and complete rerun seed bundle separately. To repeat a run, restore the applicable settings and inputs, confirm the available checksums, and use the recorded command, code, and image. Avoid editing inputs and settings while the run is active.

[⬆ Back to top](#top)

------------------------------------------------------------------------

# Troubleshooting

| Message | Fix |
|---------|-----|
| `Missing config_data file` | Launch from the working folder, or pass `--config_data` with the file's path |
| `Config file missing 'Run =' entry` / `'Order =' entry` | Regenerate `Config.txt` or add the line. Keep `Order = ...` on one line |
| `Missing step script for 03_...` | Check the spelling. Each entry must match `module_<name>.py` or `.R` |
| `Rerun folder '...' not found` | `Run =` must name an existing step folder |
| `Cannot rerun from '...': step_status.txt says 'RUNNING'` or `'INCOMPLETE'` | Run the seed step again, or choose another seed. For `INCOMPLETE`, look for `[STEP STATUS WARNING]` in `pipeline_log.txt`; see [Rerunning from a Step](#rerunning-from-a-step) |
| `Cannot rerun from '...': step_status.txt is missing` | The step was written by this version but lost its status file. Run the seed step again, or choose another seed |
| `Cannot read input file` / `Incomplete read` | Read the reason printed and check permissions and host access. For an online-only file, make it available offline or move it to a local folder |
| Initialisation failure | Start with `pipeline_log.txt`, then open `logfile_setup.txt` if present |
| `[MEMORY ERROR] ... exit 137` | Check Docker's memory, then raise `--memory` and `--max_memory` |
| Docker permission denied (Linux) | `sudo usermod -aG docker $USER`, then log out and in |

- **Module, QA or report errors.** Search `step_log.txt` for `###MOD_ERROR###`, `###QA_ERROR###` or `###REP_ERROR###`. A module or report failure stops processing. A QA failure is recorded and processing continues. Check the current run status and QA summary before using the results.

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

If you generate settings with the repository’s Python generator, you also need Python on the host.

[⬆ Back to top](#top)

------------------------------------------------------------------------

# Support

If you encounter issues, please open a GitHub issue in this repository.

[⬆ Back to top](#top)

------------------------------------------------------------------------
