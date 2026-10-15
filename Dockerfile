# syntax=docker/dockerfile:1
# (The line above enables the "heredoc" block at the end of this file, used to
#  write a script directly from the Dockerfile. Docker Desktop supports it.)

# -----------------------------------------------------------------------------
# Base image.
#
# A tag such as ":latest" can point to a different image over time, so it does
# NOT make the build reproducible. Only a digest does. Pin it once (before a
# release / paper submission):
#   docker pull mambaorg/micromamba:latest
#   docker image inspect --format '{{index .RepoDigests 0}}' mambaorg/micromamba:latest
# and put the printed value here, e.g.
#   ARG BASE_IMAGE=mambaorg/micromamba@sha256:<digest>
# The build command stays the same: docker build -t unified_env_all:1.0 .
#
# The reference given here (tag or digest) is stored in the image as
# UNIFIED_BASE_IMAGE and shown in image_versions.txt and env_log.txt. A tag
# like ":latest" is recorded as written; the digest it resolved to is NOT
# discovered. The "Reproducibility status" therefore shows a WARNING for a
# floating tag and "OK - pinned to an immutable digest" once a digest is set.
# -----------------------------------------------------------------------------
ARG BASE_IMAGE=mambaorg/micromamba@sha256:e0a99b0f17a759e14c2f967dc0ca2d3a3c1ca3c62955f4d20bba770eaaf0184d
FROM ${BASE_IMAGE}
ARG BASE_IMAGE
ENV UNIFIED_BASE_IMAGE=${BASE_IMAGE}

WORKDIR /workflow

COPY environment.yml /tmp/environment.yml
RUN micromamba env create -f /tmp/environment.yml -n unifiedenv && \
    micromamba clean --all --yes

# CRAN packages not taken from conda. 'iq' provides the MaxLFQ implementation
# that replaces the former C++ code. install.packages() only warns on failure,
# so the build is stopped explicitly if any package is still missing.
RUN micromamba run -n unifiedenv R -e "\
  pkgs <- c('iq','mice','missForest'); \
  to_install <- pkgs[!pkgs %in% rownames(installed.packages())]; \
  if(length(to_install)) install.packages(to_install, repos='https://cloud.r-project.org'); \
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]; \
  if(length(missing)) stop('R packages failed to install: ', paste(missing, collapse=', ')); \
  cat('R packages OK:', paste(pkgs, collapse=', '), '\n')"

ENV PATH="/opt/conda/envs/unifiedenv/bin:$PATH"
ENV TZ=Europe/Berlin

COPY Initialization/ ./Initialization/
COPY modules_python/ ./modules_python/
COPY modules_R/ ./modules_R/
COPY docs/ ./docs/

# Give the (non-root) micromamba user ownership of /workflow
USER root
RUN chown -R $MAMBA_USER:$MAMBA_USER /workflow
USER $MAMBA_USER

# =============================================================================
# REPRODUCIBILITY RECORD  (last step of the build)
#
# Installs the command "environment_versions" and runs it once, writing
#   /workflow/image_versions.pdf   the report ("What did we build?"): summary,
#                                  software snapshot, provenance, all packages
#   /workflow/image_versions.txt   the same information as plain text (source
#                                  of the PDF; required, the PDF is best effort)
#
# The pipeline runs the SAME command at the start of every run and writes
# env_log.txt / env_log.pdf into the working directory ("What did this run use?").
#
# Copy the PDF of a built image to your computer:
#   docker run --rm unified_env_all:1.0 cat /workflow/image_versions.pdf > image_versions_unified_env_all_1.0.pdf
# =============================================================================
USER root

COPY --chmod=755 <<'EOF_SCRIPT' /usr/local/bin/environment_versions
#!/usr/bin/env bash
# =============================================================================
# environment_versions
#
# Writes the reproducibility record of the software environment:
#   <name>.txt  authoritative plain-text record (grep / diff / parse)
#   <name>.pdf  polished report: summary, software snapshot, provenance and
#               complete package tables (if ReportLab is available)
#
# Used for both records, so they share one layout and can be compared:
#   IMAGE mode  docker build    -> /workflow/image_versions.txt (+ .pdf)
#   RUN mode    pipeline start  -> <WORKING_PATH>/env_log.txt   (+ .pdf)
#
# Usage:
#   environment_versions <output.txt> IMAGE|RUN [header_file] [changes_file]
#
#   header_file   optional text inserted after the title, already formatted
#                 as sections (the pipeline passes RUN / PIPELINE / CONTAINER /
#                 COMPUTE details collected on the host).
#   changes_file  optional "git diff" of uncommitted pipeline changes; it is
#                 stored in the appendix of the record (tracked files only;
#                 untracked files are listed by name).
#
# The command never fails: anything missing is written as "not available".
# The .txt is always written first; the .pdf is best effort and a PDF problem
# never affects the .txt or the exit status.
# =============================================================================

set +e

OUT="${1:-environment_versions.txt}"
MODE="$(printf '%s' "${2:-IMAGE}" | tr '[:lower:]' '[:upper:]')"
HEADER_FILE="${3:-}"
CHANGES_FILE="${4:-}"
# Ignore a changes file that is empty or only whitespace.
if [ -n "$CHANGES_FILE" ] && ! grep -q '[^[:space:]]' "$CHANGES_FILE" 2>/dev/null; then
    CHANGES_FILE=""
fi

# The pipeline's conda environment. CONDA_PREFIX is NOT trusted first: when
# the pipeline runs, the base image's entrypoint activates the (empty) base
# environment and sets CONDA_PREFIX=/opt/conda, which would record the wrong
# environment. Order: the known environment path, then the environment of
# the python on PATH, then CONDA_PREFIX.
ENV_PREFIX="/opt/conda/envs/unifiedenv"
if [ ! -d "$ENV_PREFIX/conda-meta" ]; then
    py_prefix="$(command -v python >/dev/null 2>&1 && python -c 'import sys; print(sys.prefix)' 2>/dev/null)"
    if [ -n "$py_prefix" ] && [ -d "$py_prefix/conda-meta" ]; then
        ENV_PREFIX="$py_prefix"
    elif [ -n "${CONDA_PREFIX:-}" ]; then
        ENV_PREFIX="$CONDA_PREFIX"
    fi
fi

EQ="==============================================================================="
DASH="-------------------------------------------------------------------------------"

center() { printf '%*s\n' $(( (${#1} + 79) / 2 )) "$1"; }
banner() { printf '\n%s\n' "$EQ"; for t in "$@"; do center "$t"; done; printf '%s\n' "$EQ"; }
section() { printf '\n%s\n%s\n' "$1" "$DASH"; }
kv() { printf '%-23s%s\n' "$1" "${2:-not available}"; }

have() { command -v "$1" >/dev/null 2>&1; }

os_name="not available"
[ -f /etc/os-release ] && os_name=$(. /etc/os-release && printf '%s' "${PRETTY_NAME:-not available}")
arch=$(uname -m 2>/dev/null)
case "$arch" in
    x86_64)  platform="linux/amd64" ;;
    aarch64) platform="linux/arm64" ;;
    *)       platform="linux/${arch:-unknown}" ;;
esac
python_v=$(have python && python -c 'import platform; print(platform.python_version())' 2>/dev/null)
r_v=$(have Rscript && Rscript --vanilla -e 'cat(R.version$major, ".", R.version$minor, sep = "")' 2>/dev/null)
mamba_v=$(have micromamba && micromamba --version 2>/dev/null | head -n 1)

# Explicit conda lock (collected first: needed for the status block).
LOCK_TMP="$(mktemp 2>/dev/null || echo /tmp/environment_versions_lock.$$)"
if have micromamba; then
    micromamba env export -p "$ENV_PREFIX" --explicit > "$LOCK_TMP" 2>&1
fi

# ---------------------------------------------------------------------------
# Reproducibility status (shown at the top of the record)
# ---------------------------------------------------------------------------
missing=""
[ -n "$python_v" ] || missing="$missing Python"
[ -n "$r_v" ] || missing="$missing R"
[ -n "$mamba_v" ] || missing="$missing micromamba"
if [ -z "$missing" ]; then
    st_env="OK - Python, R and all conda/pip/R packages recorded"
else
    st_env="WARNING - not found:$missing"
fi

# "OK" only if the lock really lists packages, not just if the command ran.
lock_count=$(grep -c '^https\?://' "$LOCK_TMP" 2>/dev/null)
lock_count=${lock_count:-0}
if grep -q '^@EXPLICIT' "$LOCK_TMP" 2>/dev/null && [ "$lock_count" -gt 0 ]; then
    st_lock="OK - exact conda lock list recorded, $lock_count packages (see Rebuild information)"
else
    st_lock="WARNING - conda lock list could not be created or contains no packages"
fi

base_ref="${UNIFIED_BASE_IMAGE:-unknown}"
case "$base_ref" in
    *@sha256:*) st_base="OK - pinned to an immutable digest" ;;
    unknown)    st_base="WARNING - base image reference unknown" ;;
    *)          st_base="WARNING - floating tag ($base_ref); pin it to a digest in the Dockerfile" ;;
esac

header_value() {
    [ -n "$HEADER_FILE" ] && [ -f "$HEADER_FILE" ] || return 0
    grep -m1 "^$1 " "$HEADER_FILE" | cut -c24- | sed 's/^ *//; s/ *$//'
}
img_id="$(header_value 'Image ID')"
img_name="$(header_value 'Image')"
img_task="$(header_value 'Image run by ENV_LOG')"
case "$img_id" in
    sha256:*) st_img="OK - immutable image ID recorded" ;;
    *)        st_img="WARNING - image ID could not be determined" ;;
esac
if [ -n "$img_task" ] && [ -n "$img_name" ] && [ "$img_task" != "$img_name" ]; then
    st_img="WARNING - recorded image ($img_name) differs from the image the ENV_LOG task ran in ($img_task)"
fi
chg="$(header_value 'Uncommitted changes')"
case "$chg" in
    NONE*) st_code="OK - clean git checkout, the commit fully identifies the code" ;;
    YES*)
        if [ -n "$CHANGES_FILE" ]; then
            st_code="WARNING - uncommitted changes; diff of tracked files recorded"
            if grep -qE '^\[\.\.\. diff truncated after [0-9]+ characters \.\.\.\]$' "$CHANGES_FILE" 2>/dev/null; then
                st_code="$st_code (TRUNCATED, over size limit)"
            fi
            if ! grep -A1 -F '[Untracked files' "$CHANGES_FILE" 2>/dev/null | tail -n 1 | grep -q '^(none)'; then
                st_code="$st_code; untracked files listed by name only"
            fi
        else
            st_code="WARNING - uncommitted changes; diff not available"
        fi ;;
    *)     st_code="WARNING - pipeline code not under git; version cannot be verified" ;;
esac

{
    if [ "$MODE" = "RUN" ]; then
        printf '%s\n' "$EQ"
        center "PIPELINE RUN"
        center "REPRODUCIBILITY RECORD"
        printf '%s\n' "$EQ"
    else
        printf '%s\n' "$EQ"
        center "UNIFIED ENVIRONMENT"
        center "REPRODUCIBILITY RECORD"
        printf '%s\n' "$EQ"
        section "IMAGE"
        kv "Generated" "$(date '+%Y-%m-%d %H:%M:%S %Z')"
        kv "Platform" "$platform"
        kv "Base image" "${UNIFIED_BASE_IMAGE:-mambaorg/micromamba (not pinned)}"
        kv "Image ID" "not known during the build (the ID exists only after the last layer); recorded in every env_log.txt and via: docker image inspect"
        kv "Conda environment" "$ENV_PREFIX"
    fi

    section "REPRODUCIBILITY STATUS"
    kv "Software environment" "$st_env"
    kv "Package lock" "$st_lock"
    kv "Base image" "$st_base"
    if [ "$MODE" = "RUN" ]; then
        kv "Container image" "$st_img"
        kv "Pipeline code" "$st_code"
    fi

    if [ -n "$HEADER_FILE" ] && [ -f "$HEADER_FILE" ]; then
        cat "$HEADER_FILE"
    fi

    section "CORE ENVIRONMENT"
    kv "Operating system" "$os_name"
    kv "Platform" "$platform"
    kv "Python" "$python_v"
    kv "R" "$r_v"
    kv "micromamba" "$mamba_v"
    if [ "$MODE" = "RUN" ] && [ -f /workflow/image_versions.txt ]; then
        built=$(grep -m1 '^Generated' /workflow/image_versions.txt | sed 's/^Generated *//')
        kv "Image built" "$built"
    fi

    section "KEY PYTHON PACKAGES"
    if have python; then
        ENV_PREFIX="$ENV_PREFIX" python - <<'EOF_PY' 2>/dev/null
# Version source: the conda installation record (conda-meta) first, because
# it is what the environment was built from; pip metadata only as fallback
# (pip-only packages, or packages whose pip metadata is wrong, e.g. 0.0.0).
import glob
import json
import os
from importlib import metadata as md

conda = {}
for f in glob.glob(os.path.join(os.environ.get("ENV_PREFIX", ""), "conda-meta", "*.json")):
    try:
        with open(f) as fh:
            meta = json.load(fh)
        conda[meta["name"].lower()] = meta["version"]
    except Exception:
        pass

packages = [
    # label             conda names                      pip distribution names
    ("NumPy",           ["numpy"],                       ["numpy"]),
    ("pandas",          ["pandas"],                      ["pandas"]),
    ("SciPy",           ["scipy"],                       ["scipy"]),
    ("scikit-learn",    ["scikit-learn"],                ["scikit-learn"]),
    ("statsmodels",     ["statsmodels"],                 ["statsmodels"]),
    ("matplotlib",      ["matplotlib", "matplotlib-base"], ["matplotlib"]),
    ("seaborn",         ["seaborn", "seaborn-base"],     ["seaborn"]),
    ("PyTorch",         ["pytorch"],                     ["torch"]),
    ("torchvision",     ["torchvision"],                 ["torchvision"]),
    ("torchaudio",      ["torchaudio"],                  ["torchaudio"]),
    ("PyArrow",         ["pyarrow", "pyarrow-core"],     ["pyarrow"]),
    ("fastparquet",     ["fastparquet"],                 ["fastparquet"]),
    ("openpyxl",        ["openpyxl"],                    ["openpyxl"]),
    ("ReportLab",       ["reportlab"],                   ["reportlab"]),
    ("pypdf",           ["pypdf"],                       ["pypdf"]),
    ("python-igraph",   ["python-igraph"],               ["igraph", "python-igraph"]),
    ("scikit-misc",     ["scikit-misc"],                 ["scikit-misc"]),
    ("tqdm",            ["tqdm"],                        ["tqdm"]),
    ("harmony-pytorch", ["harmony-pytorch"],             ["harmony-pytorch", "harmony_pytorch"]),
    ("pyComBat",        ["pycombat"],                    ["pycombat"]),
    ("pingouin",        ["pingouin"],                    ["pingouin"]),
]
for label, conda_names, dists in packages:
    version = next((conda[c] for c in conda_names if c in conda), None)
    if version is None:
        for dist in dists:
            try:
                version = md.version(dist)
                break
            except Exception:
                pass
    print(f"{label:<23}{version or 'not installed'}")
EOF_PY
    else
        echo "Python not available"
    fi

    section "KEY R PACKAGES"
    if have Rscript; then
        Rscript --vanilla -e '
            pkgs <- c("iq", "mice", "missForest", "imputeLCMD", "limma", "sva", "vsn",
                      "preprocessCore", "pcaMethods", "DESeq2", "edgeR", "genefilter",
                      "SummarizedExperiment", "Biobase", "BiocParallel", "GenomicRanges",
                      "data.table", "tidyverse", "pheatmap", "optparse")
            for (p in pkgs) {
                v <- tryCatch(as.character(packageVersion(p)), error = function(e) "not installed")
                cat(sprintf("%-23s%s\n", p, v))
            }
        ' 2>/dev/null
    else
        echo "R not available"
    fi

    banner "COMPLETE ENVIRONMENT"

    printf '\n[Conda packages: name, version, build, channel]\n\n'
    if have micromamba; then
        micromamba list -p "$ENV_PREFIX" 2>&1
    else
        echo "micromamba not available"
    fi

    printf '\n[Python packages (pip)]\n\n'
    if have python; then
        python -m pip list --format=freeze --disable-pip-version-check 2>&1
    else
        echo "Python not available"
    fi

    printf '\n[R packages: name, version, library]\n\n'
    if have Rscript; then
        Rscript --vanilla -e '
            ip <- installed.packages()[, c("Package", "Version", "LibPath")]
            ip <- ip[order(tolower(ip[, "Package"])), , drop = FALSE]
            for (i in seq_len(nrow(ip))) {
                cat(sprintf("%-32s%-16s%s\n", ip[i, "Package"], ip[i, "Version"], ip[i, "LibPath"]))
            }
        ' 2>&1
    else
        echo "R not available"
    fi

    banner "REBUILD INFORMATION"
    printf '\nThe Docker image (identified by its image ID / digest) is the unit of\n'
    printf 'reproducibility. The information below documents how it was assembled:\n\n'
    printf '  Conda explicit lock   reconstructs the exact conda environment\n'
    printf '  pip freeze            records the Python packages as seen by pip\n'
    printf '  R package inventory   records the R package versions (incl. CRAN installs)\n'
    printf '  Base image            %s\n' "$base_ref"
    printf '\nTo rebuild only the conda part: save the lock list below as unifiedenv.lock and run\n'
    printf '    micromamba create -n unifiedenv --file unifiedenv.lock\n'
    printf 'R packages installed from CRAN in the Dockerfile (iq, mice, missForest) are not\n'
    printf 'part of the conda lock; their versions are listed in the R package inventory.\n'
    if [ -f /tmp/environment.yml ]; then
        printf '\n[environment.yml used for the image build]\n\n'
        cat /tmp/environment.yml
    fi
    printf '\n[Explicit conda lock list]\n\n'
    if [ -s "$LOCK_TMP" ]; then
        cat "$LOCK_TMP"
    else
        echo "not available"
    fi

    if [ -n "$CHANGES_FILE" ] && [ -s "$CHANGES_FILE" ]; then
        banner "UNCOMMITTED PIPELINE CHANGES"
        printf '\nThe pipeline code differed from the recorded git commit.\n'
        printf 'Tracked modifications are preserved below as a git diff against that commit.\n'
        printf 'Untracked files are listed by name, but their contents are not captured.\n'
        if grep -qE '^\[\.\.\. diff truncated after [0-9]+ characters \.\.\.\]$' "$CHANGES_FILE" 2>/dev/null; then
            printf '\nWARNING: the git diff exceeded the size limit and was TRUNCATED.\n'
        fi
        printf '\n'
        cat "$CHANGES_FILE"
    fi

    printf '\n%s\n' "$EQ"
    center "END OF RECORD"
    printf '%s\n' "$EQ"
} > "$OUT" 2>&1

rm -f "$LOCK_TMP"
echo "Reproducibility record written: $OUT"

# -----------------------------------------------------------------------------
# PDF version: summary pages + full appendix. Optional; the .txt stays the
# authoritative record. Skipped quietly if Python/ReportLab are missing.
# -----------------------------------------------------------------------------
PDF="${OUT%.txt}.pdf"
rm -f "$PDF"
if ! { have python && python -c 'import reportlab' >/dev/null 2>&1; }; then
    echo "WARNING: ReportLab not available, $PDF not created. $OUT was written successfully."
else
    python - "$OUT" "$PDF" <<'EOF_PDF'
"""Render the reproducibility record (.txt) as a polished PDF report.

Page 1   summary: versions, reproducibility checks, identifiers
Page 2   software snapshot (key Python / R packages, core environment)
Page 3   provenance and rebuild information
Appendix complete package inventory as compact tables, uncommitted
         pipeline changes (run records), exact conda lock (file names)
"""
import re
import sys
from datetime import datetime

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import (CondPageBreak, KeepTogether, PageBreak, Paragraph,
                                Preformatted, SimpleDocTemplate, Spacer, Table,
                                TableStyle)

txt_path, pdf_path = sys.argv[1], sys.argv[2]
lines = open(txt_path, encoding="utf-8", errors="replace").read().splitlines()

EQ_RE = re.compile(r"^=+$")
DASH_RE = re.compile(r"^-+$")

# --------------------------------------------------------------------------
# Parse the text record
# --------------------------------------------------------------------------
n = len(lines)
i = 0
title_lines = []
if n and EQ_RE.match(lines[0]):
    i = 1
    while i < n and not EQ_RE.match(lines[i]):
        title_lines.append(lines[i].strip())
        i += 1
    i += 1

# summary sections ("NAME" + dashed line + "label   value" rows)
sections, current = [], None
while i < n and not EQ_RE.match(lines[i]):
    line = lines[i]
    if i + 1 < n and DASH_RE.match(lines[i + 1]) and line.strip():
        current = {"name": line.strip(), "rows": []}
        sections.append(current)
        i += 2
        continue
    if current is not None and line.strip():
        label, value = line[:23].strip(), line[23:].strip()
        if not value:
            label, value = "", line.strip()
        current["rows"].append((label, value))
    i += 1

# banner blocks ("=====" / TITLE / "=====" / content)
blocks = {}
while i < n:
    if EQ_RE.match(lines[i]):
        j = i + 1
        names = []
        while j < n and not EQ_RE.match(lines[j]):
            names.append(lines[j].strip())
            j += 1
        k = j + 1
        while k < n and not EQ_RE.match(lines[k]):
            k += 1
        blocks[" ".join(names)] = lines[j + 1:k]
        i = k
    else:
        i += 1


def subsections(block):
    """Split a block into {"[name]": lines}; text before the first is 'intro'."""
    out, name = {"intro": []}, "intro"
    for line in block:
        s = line.strip()
        if re.match(r"^\[[A-Za-z]", s) and s.endswith("]"):
            name = s[1:-1]
            out[name] = []
        else:
            out[name].append(line)
    return out


def find(subs, prefix):
    for key, value in subs.items():
        if key.lower().startswith(prefix.lower()):
            return value
    return []


sec = {s["name"]: s["rows"] for s in sections}
lookup = {}
for s in sections:
    for label, value in s["rows"]:
        lookup.setdefault(label, value)

is_run = any("PIPELINE RUN" in t for t in title_lines)

env_subs = subsections(blocks.get("COMPLETE ENVIRONMENT", []))
reb_subs = subsections(blocks.get("REBUILD INFORMATION", []))
chg_subs = subsections(blocks.get("UNCOMMITTED PIPELINE CHANGES", []))

# conda packages
conda_rows = []
for line in find(env_subs, "Conda packages"):
    s = line.strip()
    if not s or s.startswith("List of packages") or set(s) <= set("─-━ "):
        continue
    parts = s.split()
    if parts[:2] == ["Name", "Version"] or len(parts) < 3:
        continue
    conda_rows.append((parts[0], parts[1], parts[2], parts[3] if len(parts) > 3 else ""))

# pip packages
pip_rows = []
for line in find(env_subs, "Python packages"):
    s = line.strip()
    if "==" in s:
        name, version = s.split("==", 1)
        pip_rows.append((name, version))
    elif " @ " in s:
        name, src = s.split(" @ ", 1)
        pip_rows.append((name, "from " + src.rsplit("/", 1)[-1][:40]))

# R packages
r_rows = []
for line in find(env_subs, "R packages"):
    parts = line.split()
    if len(parts) >= 2 and not line.startswith("R not"):
        r_rows.append((parts[0], parts[1]))

# exact conda lock (file names + channel)
lock_rows = []
for line in find(reb_subs, "Explicit conda lock"):
    s = line.strip()
    if s.startswith("http"):
        url = s.split("#", 1)[0]
        segs = url.split("/")
        fname = segs[-1]
        channel = "/".join(segs[3:-1]) if len(segs) > 4 else ""
        lock_rows.append((fname, channel))

env_yml = [l for l in find(reb_subs, "environment.yml") if l.strip() or l == ""]
while env_yml and not env_yml[-1].strip():
    env_yml.pop()
while env_yml and not env_yml[0].strip():
    env_yml.pop(0)

# --------------------------------------------------------------------------
# Styles
# --------------------------------------------------------------------------
INK = colors.HexColor("#17202a")
MUTED = colors.HexColor("#5f6b7a")
ACCENT = colors.HexColor("#1d4e89")
ACCENT_LIGHT = colors.HexColor("#e8eef7")
LINE = colors.HexColor("#d9dee5")
PANEL = colors.HexColor("#f4f6f9")
OK = colors.HexColor("#1e7b3a")
OK_LIGHT = colors.HexColor("#e7f4ec")
WARN = colors.HexColor("#b45309")
WARN_LIGHT = colors.HexColor("#fdf1e3")
WHITE = colors.white

W = A4[0] - 36 * mm  # usable width


def st(name, **kw):
    base = dict(fontName="Helvetica", fontSize=8.5, leading=11, textColor=INK)
    base.update(kw)
    return ParagraphStyle(name, **base)


s_kicker = st("kicker", fontName="Helvetica-Bold", fontSize=8, leading=10, textColor=colors.HexColor("#c9d6ea"))
s_title = st("title", fontName="Helvetica-Bold", fontSize=22, leading=26, textColor=WHITE)
s_meta = st("meta", fontSize=9, leading=12, textColor=colors.HexColor("#dbe4f1"))
s_h1 = st("h1", fontName="Helvetica-Bold", fontSize=15, leading=19, spaceAfter=2)
s_h2 = st("h2", fontName="Helvetica-Bold", fontSize=10, leading=13, textColor=ACCENT, spaceBefore=9, spaceAfter=4)
s_note = st("note", fontSize=8, leading=10.5, textColor=MUTED)
s_lab = st("lab", fontSize=8, leading=10.5, textColor=MUTED)
s_val = st("val")
s_mono = st("mono", fontName="Courier", fontSize=7.8, leading=10)
s_card_l = st("cardl", fontName="Helvetica-Bold", fontSize=6.8, leading=8.5, textColor=MUTED)
s_card_v = st("cardv", fontName="Helvetica-Bold", fontSize=13, leading=16)
s_tab = st("tab", fontSize=6.8, leading=8.2)
s_tabh = st("tabh", fontName="Helvetica-Bold", fontSize=6.8, leading=8.2, textColor=MUTED)
s_code = st("code", fontName="Courier", fontSize=6.2, leading=7.4)


def esc(text):
    return (text or "").replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def P(text, style=s_val):
    return Paragraph(esc(text), style)


def kv_table(rows, label_w=44 * mm):
    data = [[P(l, s_lab), P(v, s_val)] for l, v in rows]
    t = Table(data, colWidths=[label_w, W - label_w])
    t.setStyle(TableStyle([
        ("LINEBELOW", (0, 0), (-1, -2), 0.3, LINE),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("TOPPADDING", (0, 0), (-1, -1), 3), ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
        ("LEFTPADDING", (0, 0), (-1, -1), 0),
    ]))
    return t


def data_table(header, rows, widths, font=s_tab, gap_col=None):
    data = [[P(h, s_tabh) for h in header]] + [[P(c, font) for c in r] for r in rows]
    t = Table(data, colWidths=widths, repeatRows=1)
    style = [
        ("LINEBELOW", (0, 0), (-1, 0), 0.6, ACCENT),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("TOPPADDING", (0, 0), (-1, -1), 1.0), ("BOTTOMPADDING", (0, 0), (-1, -1), 1.0),
        ("LEFTPADDING", (0, 0), (-1, -1), 3), ("RIGHTPADDING", (0, 0), (-1, -1), 3),
    ]
    if gap_col is not None:
        style.append(("LINEAFTER", (gap_col, 0), (gap_col, -1), 0.6, LINE))
    for r in range(2, len(data), 2):
        style.append(("BACKGROUND", (0, r), (-1, r), PANEL))
    t.setStyle(TableStyle(style))
    return t


def paired(rows):
    """Two entries per table row (left, right), to halve the page count."""
    out = []
    width = len(rows[0]) if rows else 0
    for k in range(0, len(rows), 2):
        a = tuple(rows[k])
        b = tuple(rows[k + 1]) if k + 1 < len(rows) else ("",) * width
        out.append(a + b)
    return out


def single_channel(values):
    vals = {v for v in values if v}
    return vals.pop() if len(vals) == 1 else None


# --------------------------------------------------------------------------
# Content
# --------------------------------------------------------------------------
story = []
kind = "PIPELINE RUN" if is_run else "IMAGE BUILD"
title = "Pipeline Run" if is_run else "Unified Environment"
when = lookup.get("Run started") or lookup.get("Generated") or datetime.now().strftime("%Y-%m-%d %H:%M")
when = when.replace("T", " ")[:19]
platform = lookup.get("Platform", "")
meta = " · ".join(x for x in [when, platform, lookup.get("Image", "") if is_run else lookup.get("Operating system", "")] if x)

header = Table([[P(kind + " · REPRODUCIBILITY REPORT", s_kicker)],
                [P(title, s_title)],
                [P(meta, s_meta)]], colWidths=[W])
header.setStyle(TableStyle([
    ("BACKGROUND", (0, 0), (-1, -1), ACCENT),
    ("LEFTPADDING", (0, 0), (-1, -1), 12), ("RIGHTPADDING", (0, 0), (-1, -1), 12),
    ("TOPPADDING", (0, 0), (-1, 0), 12), ("BOTTOMPADDING", (0, -1), (-1, -1), 12),
    ("TOPPADDING", (0, 1), (-1, -1), 2), ("BOTTOMPADDING", (0, 0), (-1, -2), 2),
]))
story += [header, Spacer(1, 6 * mm)]

# version cards
if is_run:
    cards = [("NEXTFLOW", (lookup.get("Nextflow") or "n/a").split(" (")[0]),
             ("PYTHON", lookup.get("Python")), ("R", lookup.get("R")),
             ("MICROMAMBA", lookup.get("micromamba"))]
else:
    cards = [("PYTHON", lookup.get("Python")), ("R", lookup.get("R")),
             ("MICROMAMBA", lookup.get("micromamba")), ("PLATFORM", platform)]
cw = W / 4
card_t = Table([[[P(l, s_card_l), Spacer(1, 2), P(v or "n/a", s_card_v)] for l, v in cards]],
               colWidths=[cw] * 4)
card_t.setStyle(TableStyle([
    ("BACKGROUND", (0, 0), (-1, -1), PANEL),
    ("LINEAFTER", (0, 0), (-2, -1), 1.5, WHITE),
    ("LEFTPADDING", (0, 0), (-1, -1), 9), ("TOPPADDING", (0, 0), (-1, -1), 8),
    ("BOTTOMPADDING", (0, 0), (-1, -1), 9), ("VALIGN", (0, 0), (-1, -1), "TOP"),
]))
story += [card_t, Spacer(1, 6 * mm)]

# reproducibility checks
checks = sec.get("REPRODUCIBILITY STATUS", [])
n_ok = sum(1 for _, v in checks if v.startswith("OK"))
n_all = len(checks)
if checks:
    good = n_ok == n_all
    verdict = (f"All {n_all} reproducibility checks passed" if good
               else f"{n_all - n_ok} of {n_all} reproducibility checks need attention")
    v_t = Table([[P(verdict, st("verdict", fontName="Helvetica-Bold", fontSize=10.5, leading=13,
                                 textColor=OK if good else WARN))]], colWidths=[W])
    v_t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), OK_LIGHT if good else WARN_LIGHT),
        ("LINEBEFORE", (0, 0), (0, -1), 3, OK if good else WARN),
        ("LEFTPADDING", (0, 0), (-1, -1), 10), ("TOPPADDING", (0, 0), (-1, -1), 7),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 7),
    ]))
    rows = []
    for label, value in checks:
        head, _, detail = value.partition(" - ")
        ok = head == "OK"
        color = OK if ok else WARN
        rows.append([
            P("●", st("dot", fontSize=11, leading=12, textColor=color)),
            P(label, st("cl", fontName="Helvetica-Bold", fontSize=8.8, leading=11)),
            P("READY" if ok else "ATTENTION", st("cs", fontName="Helvetica-Bold", fontSize=7.5, leading=11, textColor=color)),
            P(detail or value, s_val),
        ])
    c_t = Table(rows, colWidths=[6 * mm, 40 * mm, 22 * mm, W - 68 * mm])
    c_t.setStyle(TableStyle([
        ("LINEBELOW", (0, 0), (-1, -2), 0.3, LINE),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("TOPPADDING", (0, 0), (-1, -1), 4), ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ("LEFTPADDING", (0, 0), (-1, -1), 2),
    ]))
    story += [P("Reproducibility", s_h2), v_t, Spacer(1, 2 * mm), c_t]

# immutable identifiers
ids = []
if lookup.get("Image ID", "").startswith("sha256:"):
    ids.append(("Image ID", lookup["Image ID"]))
dig = lookup.get("Image digest", "")
if dig and not dig.lower().startswith("none"):
    ids.append(("Image digest", dig))
if lookup.get("Git commit"):
    ids.append(("Pipeline git commit", lookup["Git commit"]))
base = lookup.get("Base image", "") if not is_run else ""
if ids:
    id_t = Table([[P(l, s_lab), P(v, s_mono)] for l, v in ids], colWidths=[36 * mm, W - 36 * mm])
    id_t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), ACCENT_LIGHT),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LEFTPADDING", (0, 0), (-1, -1), 8), ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]))
    story += [P("Identifiers", s_h2), id_t]

# key facts of page 1
if is_run:
    facts = [r for r in sec.get("RUN", []) if r[0] in
             ("Run name", "Run started", "Command line", "Run mode", "Steps", "Working path")]
    facts += [r for r in sec.get("CONTAINER", []) if r[0] in ("Image", "Image run by ENV_LOG")]
    if facts:
        story += [P("Run", s_h2), kv_table(facts)]
else:
    facts = [r for r in sec.get("IMAGE", []) if r[0] in ("Generated", "Base image", "Image ID", "Conda environment")]
    facts += [r for r in sec.get("CORE ENVIRONMENT", []) if r[0] in ("Operating system",)]
    if facts:
        story += [P("Image", s_h2), kv_table(facts)]

# ---- page 2: software snapshot ----
story.append(PageBreak())
story.append(P("Software snapshot", s_h1))
story.append(P("Key packages of the analysis environment. The complete inventory is in the appendix.", s_note))
story.append(Spacer(1, 4 * mm))
py_rows = sec.get("KEY PYTHON PACKAGES", [])
r_key = sec.get("KEY R PACKAGES", [])
half = (W - 8 * mm) / 2


def mini(header, rows):
    if not rows:
        rows = [("not available", "")]
    return data_table([header, "Version"], rows, [half * 0.62, half * 0.38], font=st("m8", fontSize=8, leading=10.5))


snap = Table([[mini("Python package", py_rows), mini("R / Bioconductor package", r_key)]],
             colWidths=[half + 4 * mm, half + 4 * mm])
snap.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP"), ("LEFTPADDING", (0, 0), (-1, -1), 0)]))
story.append(snap)
core = sec.get("CORE ENVIRONMENT", [])
if core:
    story += [P("Core environment", s_h2), kv_table(core)]

# ---- page 3: provenance & rebuild ----
story.append(PageBreak())
story.append(P("Provenance and rebuild", s_h1))
story.append(P("Where this environment comes from and how it can be reconstructed.", s_note))
for name in (("RUN", "PIPELINE", "CONTAINER", "COMPUTE") if is_run else ("IMAGE",)):
    if sec.get(name):
        story += [P(name.title(), s_h2), kv_table(sec[name])]

intro = [l.strip() for l in reb_subs.get("intro", []) if l.strip()]
if intro:
    story.append(P("Rebuild information", s_h2))
    explain = [l for l in intro if not re.match(r"^(Conda explicit lock|pip freeze|R package inventory|Base image)\b", l)]
    mapping = [l for l in intro if re.match(r"^(Conda explicit lock|pip freeze|R package inventory|Base image)\b", l)]
    if mapping:
        rows = []
        for l in mapping:
            parts = re.split(r"\s{2,}", l, maxsplit=1)
            rows.append((parts[0], parts[1] if len(parts) > 1 else ""))
        story.append(kv_table(rows))
        story.append(Spacer(1, 2 * mm))
    para = []
    for l in explain + [""]:
        if l.startswith("micromamba ") or not l:
            if para:
                story.append(P(" ".join(para), s_note))
                story.append(Spacer(1, 1.5 * mm))
                para = []
            if l:
                story.append(Preformatted("    " + l, s_mono))
                story.append(Spacer(1, 1.5 * mm))
        else:
            para.append(l)
    story.append(P(f"Exact conda lock: {len(lock_rows)} package files (Appendix D). "
                   "The environment.yml used for the build is in Appendix E.", s_note))

# ---- appendix ----
story.append(PageBreak())
story.append(P("Appendix A — Conda packages", s_h1))
conda_channel = single_channel(r[3] for r in conda_rows)
story.append(P(f"{len(conda_rows)} packages in the conda environment"
               + (f", all from channel {conda_channel}." if conda_channel else " (name, version, build, channel)."), s_note))
story.append(Spacer(1, 3 * mm))
if conda_rows and conda_channel:
    h = W / 2
    story.append(data_table(["Package", "Version", "Build"] * 2,
                            paired([r[:3] for r in conda_rows]),
                            [h * 0.42, h * 0.24, h * 0.34] * 2, gap_col=2))
elif conda_rows:
    story.append(data_table(["Package", "Version", "Build", "Channel"], conda_rows,
                            [W * 0.34, W * 0.18, W * 0.30, W * 0.18]))
else:
    story.append(P("not available", s_val))

pw = [W * 0.32, W * 0.18] * 2
story.append(CondPageBreak(60 * mm))
story.append(P("Appendix B — Python packages (pip)", s_h1))
story.append(P(f"{len(pip_rows)} packages as reported by pip.", s_note))
story.append(Spacer(1, 3 * mm))
story.append(data_table(["Package", "Version", "Package", "Version"], paired(pip_rows), pw, gap_col=1)
             if pip_rows else P("not available", s_val))

story.append(CondPageBreak(60 * mm))
story.append(P("Appendix C — R packages", s_h1))
story.append(P(f"{len(r_rows)} installed R packages, including CRAN installs.", s_note))
story.append(Spacer(1, 3 * mm))
story.append(data_table(["Package", "Version", "Package", "Version"], paired(r_rows), pw, gap_col=1)
             if r_rows else P("not available", s_val))

story.append(CondPageBreak(60 * mm))
story.append(P("Appendix D — Exact conda lock", s_h1))
lock_channel = single_channel(r[1] for r in lock_rows)
story.append(P("Every conda package pinned to an exact file (name-version-build)"
               + (f", channel {lock_channel}" if lock_channel else "")
               + ". The full download URLs with checksums are in the .txt record.", s_note))
story.append(Spacer(1, 3 * mm))
if lock_rows and lock_channel:
    story.append(data_table(["Package file"] * 2, paired([(r[0],) for r in lock_rows]), [W / 2] * 2, gap_col=0))
elif lock_rows:
    story.append(data_table(["Package file", "Channel / platform"], lock_rows, [W * 0.7, W * 0.3]))
else:
    story.append(P("not available", s_val))

if env_yml:
    story.append(CondPageBreak(80 * mm))
    story.append(P("Appendix E — environment.yml used for the image build", s_h1))
    story.append(Spacer(1, 3 * mm))
    story.append(Preformatted("\n".join(env_yml), s_code))

if chg_subs.get("intro") or len(chg_subs) > 1:
    story.append(PageBreak())
    story.append(P("Appendix F — Uncommitted pipeline changes", s_h1))
    for l in [x.strip() for x in chg_subs.get("intro", []) if x.strip()]:
        story.append(P(l, st("w", fontSize=8, leading=10.5, textColor=WARN if "WARNING" in l else MUTED)))
    files = [l.rstrip() for l in find(chg_subs, "git status") if l.strip()]
    if files:
        story += [P("Changed files", s_h2), Preformatted("\n".join(files), s_mono)]
    untracked = [l.strip() for l in find(chg_subs, "Untracked") if l.strip()]
    if untracked:
        story += [P("Untracked files (names only)", s_h2), Preformatted("\n".join(untracked), s_mono)]
    diff = find(chg_subs, "git diff")
    if diff:
        limit = 2500
        shown = diff[:limit]
        story.append(P("Diff of tracked files", s_h2))
        if len(diff) > limit:
            story.append(P(f"Showing the first {limit} of {len(diff)} lines; the complete diff is in the .txt record.",
                           st("w2", fontSize=8, leading=10.5, textColor=WARN)))
        story.append(Preformatted("\n".join(shown), s_code, maxLineLength=150, newLineChars="  "))


def decorate(canvas, doc):
    canvas.saveState()
    canvas.setFillColor(ACCENT)
    canvas.rect(0, A4[1] - 4, A4[0], 4, stroke=0, fill=1)
    canvas.setFont("Helvetica", 7)
    canvas.setFillColor(MUTED)
    canvas.drawString(18 * mm, 10 * mm, f"{title} · reproducibility report · {when}")
    canvas.drawRightString(A4[0] - 18 * mm, 10 * mm, f"Page {doc.page}")
    canvas.restoreState()


doc = SimpleDocTemplate(pdf_path, pagesize=A4, leftMargin=18 * mm, rightMargin=18 * mm,
                        topMargin=15 * mm, bottomMargin=17 * mm,
                        title=f"{title} reproducibility report", author="environment_versions")
doc.build(story, onFirstPage=decorate, onLaterPages=decorate)
print(f"PDF written: {pdf_path}")
EOF_PDF
    if [ "$?" -ne 0 ]; then
        rm -f "$PDF"
        echo "WARNING: Could not generate $PDF. $OUT was written successfully; continuing."
    fi
fi

exit 0
EOF_SCRIPT

RUN environment_versions /workflow/image_versions.txt IMAGE && \
    chown $MAMBA_USER:$MAMBA_USER /workflow/image_versions.txt /workflow/image_versions.pdf 2>/dev/null; \
    test -s /workflow/image_versions.txt

USER $MAMBA_USER

CMD ["bash"]
