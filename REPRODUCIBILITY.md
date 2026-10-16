{\rtf1\ansi\ansicpg1252\cocoartf2907
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 # Reproducibility\
\
The pipeline automatically creates two PDF reports:\
\
| PDF | When | Purpose |\
|---|---|---|\
| `image_versions.pdf` | Docker image build | Software and package versions in the image |\
| `env_log.pdf` | Every pipeline run | Image, pipeline version and settings actually used |\
\
For a release, the image PDF should show **3/3 checks passed** and the run PDF **5/5 checks passed**.\
\
---\
\
## New release\
\
Run everything from the `Preprocessing` folder.\
\
### 1. Choose the image version\
\
If only pipeline code changed and the software environment is unchanged, keep the existing Docker image.\
\
If `Dockerfile`, `environment.yml`, R packages or the base image changes, use a new version, for example:\
\
`1.0 \uc0\u8594  1.1`\
\
A released image version should never be rebuilt with different contents.\
\
---\
\
### 2. Pin the base image\
\
For every new image version:\
\
```bash\
docker pull mambaorg/micromamba:latest\
docker image inspect --format '\{\{index .RepoDigests 0\}\}' mambaorg/micromamba:latest\
```\
\
Copy the returned `mambaorg/micromamba@sha256:...` value into the Dockerfile:\
\
```dockerfile\
ARG BASE_IMAGE=mambaorg/micromamba@sha256:...\
```\
\
This keeps the base image reproducible.\
\
---\
\
### 3. Set and build the image version\
\
Set the version in `nextflow.config`:\
\
```groovy\
params.container_image = params.container_image ?: 'unified_env_all:1.1'\
```\
\
Build:\
\
```bash\
docker build -t unified_env_all:1.1 .\
```\
\
The build automatically creates `image_versions.pdf` inside the image.\
\
Save the PDF:\
\
```bash\
mkdir -p image_versions\
\
docker run --rm unified_env_all:1.1 \\\
  cat /workflow/image_versions.pdf \\\
  > image_versions/image_versions_unified_env_all_1.1.pdf\
```\
\
Open it and check:\
\
**3/3 reproducibility checks passed**\
\
Also verify that the important Python and R package versions look correct.\
\
---\
\
### 4. Commit\
\
```bash\
git add -A\
git commit -m "Release: unified_env_all:1.1"\
git status\
```\
\
Before the final pipeline run, Git should report:\
\
```text\
nothing to commit, working tree clean\
```\
\
---\
\
### 5. Run the pipeline\
\
Run the pipeline normally, for example:\
\
```bash\
nextflow run nextflow.nf \\\
  -profile linux_docker_user \\\
  --memory '5 GB' \\\
  --cpus 1\
```\
\
The pipeline automatically creates:\
\
```text\
env_log.pdf\
```\
\
in the working folder before processing starts.\
\
Open it and check:\
\
**5/5 reproducibility checks passed**\
\
Also confirm that the expected Docker image and Git commit are shown.\
\
---\
\
### 6. Release\
\
If both PDFs are correct:\
\
```bash\
git push\
```\
\
Optionally create a Git release tag:\
\
```bash\
git tag v1.1\
git push --tags\
```\
\
---\
\
## What should I keep?\
\
For a publication, keep:\
\
- `image_versions_unified_env_all_<version>.pdf`\
- `env_log.pdf` from the analysis\
- the Git commit or release tag\
\
That's enough to identify the software environment and pipeline code used for the analysis.\
\
---\
\
## If something says ATTENTION\
\
| Problem | Solution |\
|---|---|\
| Base image is a floating `:latest` tag | Pin the digest and rebuild |\
| Package/environment check fails | Check the Docker build and rebuild |\
| Container images differ | Check `params.container_image` |\
| Pipeline has uncommitted changes | Commit the changes and rerun |\
| `env_log.pdf` is missing | Rebuild the Docker image |\
\
---\
\
## In short\
\
For a new software environment:\
\
```text\
Pin base image\
     \uc0\u8595 \
Set new image version\
     \uc0\u8595 \
docker build\
     \uc0\u8595 \
Check image_versions.pdf (3/3)\
     \uc0\u8595 \
Commit\
     \uc0\u8595 \
Run pipeline\
     \uc0\u8595 \
Check env_log.pdf (5/5)\
     \uc0\u8595 \
git push\
```\
\
For normal analysis runs with an existing image:\
\
```text\
Run pipeline \uc0\u8594  check env_log.pdf \u8594  done\
```}