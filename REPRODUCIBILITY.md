# Reproducibility

The pipeline automatically creates two PDF reports:

| PDF | When | Purpose |
|---|---|---|
| `image_versions.pdf` | Docker image build | Software and package versions in the image |
| `env_log.pdf` | Every pipeline run | Image, pipeline version and settings actually used |

For a release, the image PDF should show **3/3 checks passed** and the run PDF **5/5 checks passed**.

---

## New release

Run everything from the `Preprocessing` folder.

### 1. Choose the image version

If only pipeline code changed and the software environment is unchanged, keep the existing Docker image.

If `Dockerfile`, `environment.yml`, R packages or the base image changes, use a new version, for example:

`1.0 → 1.1`

A released image version should never be rebuilt with different contents.

---

### 2. Pin the base image

For every new image version:

```bash
docker pull mambaorg/micromamba:latest
docker image inspect --format '{{index .RepoDigests 0}}' mambaorg/micromamba:latest
```

Copy the returned `mambaorg/micromamba@sha256:...` value into the Dockerfile:

```dockerfile
ARG BASE_IMAGE=mambaorg/micromamba@sha256:...
```

This keeps the base image reproducible.

---

### 3. Set and build the image version

Set the version in `nextflow.config`:

```groovy
params.container_image = params.container_image ?: 'unified_env_all:1.1'
```

Build:

```bash
docker build -t unified_env_all:1.1 .
```

The build automatically creates `image_versions.pdf` inside the image.

Save the PDF:

```bash
mkdir -p image_versions

docker run --rm unified_env_all:1.1 \
  cat /workflow/image_versions.pdf \
  > image_versions/image_versions_unified_env_all_1.1.pdf
```

Open it and check:

**3/3 reproducibility checks passed**

Also verify that the important Python and R package versions look correct.

---

### 4. Commit

```bash
git add -A
git commit -m "Release: unified_env_all:1.1"
git status
```

Before the final pipeline run, Git should report:

```text
nothing to commit, working tree clean
```

---

### 5. Run the pipeline

Run the pipeline normally, for example:

```bash
nextflow run nextflow.nf \
  -profile linux_docker_user \
  --memory '5 GB' \
  --cpus 1
```

The pipeline automatically creates:

```text
env_log.pdf
```

in the working folder before processing starts.

Open it and check:

**5/5 reproducibility checks passed**

