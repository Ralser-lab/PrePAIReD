#!/usr/bin/env bash
set -uo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE="${SCRIPT_DIR}/nextflow.nf"

pipeline_args=()
color_mode=0

for argument in "$@"; do
    if [ "$argument" = "--color" ]; then
        color_mode=1
    elif [ "$argument" = "--clean" ]; then
        :
    else
        pipeline_args+=("$argument")
    fi
done

if [ "$color_mode" -eq 1 ]; then
    exec nextflow run "$PIPELINE" \
        -ansi-log true \
        "${pipeline_args[@]}"
fi

FILTER='
BEGIN {
    reset = "\033[0m"
    cyan = "\033[36m"
    green = "\033[32m"
    red = "\033[31m"
    yellow = "\033[33m"
    blue = "\033[34m"
    dim = "\033[2m"
}

{
    clean = $0
    gsub(/\033\[[0-9;?]*[ -\/]*[@-~]/, "", clean)
}

clean ~ /ERROR ~/ {
    skip = 1
    next
}

skip && clean ~ /Check .*\.nextflow\.log.*for details/ {
    skip = 0
    next
}

skip { next }

{
    output = $0

    if (clean ~ /^executor >/) {
        output = cyan output reset
    } else if (clean ~ /^\[[^]]+\].*(✔|COMPLETED|PASSED)/) {
        output = green output reset
    } else if (clean ~ /^\[[^]]+\].*(✘|FAILED|FATAL)/ ||
               clean ~ /PIPELINE FAILED/) {
        output = red output reset
    } else if (clean ~ /WARN/) {
        output = yellow output reset
    } else if (clean ~ /CLEANUP/) {
        output = blue output reset
    } else if (clean ~ /^={10,}$/) {
        output = dim output reset
    }

    print output
    fflush()
}
'

set -o pipefail

nextflow run "$PIPELINE" \
    -ansi-log false \
    "${pipeline_args[@]}" 2>&1 |
    awk "$FILTER"

nextflow_exit=${PIPESTATUS[0]}
exit "$nextflow_exit"
