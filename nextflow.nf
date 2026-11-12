nextflow.enable.dsl=2

def pget(String key, Object defaultValue) {
    return params.containsKey(key) ? params[key] : defaultValue
}

def die(String msg) {
    throw new RuntimeException(msg)
}

Map parseKeyValueFile(String path) {
    def m = [:]
    def f = new File(path)
    if( !f.exists() ) die("Missing config_data file: ${path}")

    f.eachLine { line ->
        def s = line.trim()
        if( !s || s.startsWith('#') ) return
        def parts = s.split('=', 2)
        if( parts.size() == 2 ) {
            m[parts[0].trim()] = parts[1].trim().replaceAll(/^["']|["']$/, '')
        }
    }

    return m
}

String resolveLaunch(String p) {
    def f = new File(p)
    return f.isAbsolute() ? f.absolutePath : new File(launchDir.toString(), p).absolutePath
}

String resolveFrom(String base, String p) {
    def f = new File(p)
    return f.isAbsolute() ? f.absolutePath : new File(base, p).absolutePath
}

String resolveCode(String p) {
    def f = new File(p)
    return f.isAbsolute() ? f.absolutePath : new File(projectDir.toString(), p).absolutePath
}

File mustExist(String path, String label) {
    def f = new File(path)
    if( !f.exists() ) die("Missing ${label}: ${path}")
    return f
}

File ensureExistsOrCreateEmpty(String path, String label) {
    def f = new File(path)
    if( !f.exists() ) {
        f.parentFile?.mkdirs()
        f.text = ''
        println "Created empty ${label}: ${path}"
    }
    return f
}

/*
 * CLOUD-SYNC SAFEGUARD (iCloud, OneDrive, Dropbox, Google Drive, Box, ...)
 *
 * Cloud-sync clients can keep files "online-only": the file is listed and
 * has a size, but its contents are not on the disk. The host OS downloads
 * such a file automatically when a program on the host reads it, but a
 * Docker container cannot trigger that download and fails with
 * "OSError: [Errno 5] Input/output error" (or similar).
 *
 * Nextflow runs on the host, so every input is read once here, before any
 * container starts. This forces the sync client to download online-only
 * files. Files that are already local are just read quickly.
 * If a file cannot be downloaded (offline, sync paused, signed out), the
 * pipeline stops here with a clear message instead of failing in Docker.
 */
def materialize(File f, String label, Map stats) {
    if( f == null || !f.exists() ) return

    if( f.isDirectory() ) {
        f.listFiles()?.each { child -> materialize(child, label, stats) }
        return
    }

    def t0 = System.currentTimeMillis()
    // Count bytes in a map: closures should not reassign outer local variables.
    def counter = [bytes: 0L]
    try {
        f.eachByte(1048576) { buf, len -> counter.bytes += len }
    }
    catch( Exception e ) {
        die("""
Cannot read input file: ${f.absolutePath} (${label})
Reason: ${e.message}

This file is probably stored online-only by a cloud-sync service
(iCloud Drive, OneDrive, Dropbox, Google Drive, ...) and could not be downloaded.
Fix: make sure you are online and the sync client is running and signed in,
or mark the folder as "Keep Downloaded" / "Always keep on this device",
or move the project to a folder that is not synced.
""")
    }

    def ms = System.currentTimeMillis() - t0
    def n = counter.bytes
    if( n != f.length() ) {
        die("Incomplete read of ${f.absolutePath}: got ${n} of ${f.length()} bytes. " +
            "The file may be online-only in a cloud-sync folder and could not be fully downloaded.")
    }

    stats.files += 1
    stats.bytes += n
    if( ms > 500 ) {
        stats.slow << "${f.absolutePath} (${ms} ms)"
    }
}

/*
 * Run an external command WITHOUT a shell.
 *
 * Each list element is passed to the program as exactly one argument,
 * so paths containing spaces (or quotes) are never split or misread.
 * (A plain "command string".execute() splits on every space and does not
 * understand quotes, which breaks for any path with spaces.)
 * Output is printed; failures only produce a warning.
 */
def runCommand(List cmd) {
    try {
        def proc = cmd.collect { it.toString() }.execute()
        def out = new StringBuilder()
        def err = new StringBuilder()
        proc.waitForProcessOutput(out, err)
        if( out.toString().trim() ) println out.toString().trim()
        if( proc.exitValue() != 0 && err.toString().trim() ) {
            println "[CLEANUP WARNING] ${err.toString().trim()}"
        }
    }
    catch( Exception e ) {
        println "[CLEANUP WARNING] Could not run ${cmd[0]}: ${e.message}"
    }
}

/*
 * Run an external command WITHOUT a shell and return its output (trimmed).
 * Returns an empty string if the command fails or does not exist.
 * Used for the version information in env_log.txt.
 */
String captureCommand(List cmd) {
    try {
        def proc = cmd.collect { it.toString() }.execute()
        def out = new StringBuilder()
        def err = new StringBuilder()
        proc.waitForProcessOutput(out, err)
        return proc.exitValue() == 0 ? out.toString().trim() : ''
    }
    catch( Exception e ) {
        return ''
    }
}

/*
 * env_log.txt layout helpers (same layout as the in-image command
 * "environment_versions"): a section title with a dashed line, then
 * "label   value" lines with the value starting at column 24.
 */
String envSection(String title) {
    return "\n${title}\n" + ('-' * 79)
}

String envLine(String label, Object value) {
    def shown = (value == null || value.toString().trim() == '') ? 'unknown' : value.toString().trim()
    return "${label.padRight(23)}${shown}"
}

/*
 * Quote a value for the shell: wraps it in single quotes, so spaces and
 * special characters stay inside one argument. An apostrophe inside the
 * value is written as '\'' (close quote, escaped apostrophe, reopen quote).
 */
String shellQuote(String s) {
    return "'" + s.replace("'", "'\\''") + "'"
}

/*
 * Extra Docker volume options for folders the processes write to directly.
 *
 * nextflow.config only mounts params.docker_mount_root (the parent of the
 * launch folder). If WORKING_PATH (or another written/read folder) lies
 * outside it, Docker would write into the container's own temporary disk,
 * and the files would silently disappear when the container stops.
 * Every folder not already covered by that mount is mounted here.
 */
String dockerMountOpts(List folders, String mountRoot) {
    def roots = []
    folders.findAll { it }.each { folder ->
        def abs = new File(folder.toString()).absolutePath
        def covered = (mountRoot && (abs == mountRoot || abs.startsWith(mountRoot + '/'))) ||
                      roots.any { r -> abs == r || abs.startsWith(r + '/') }
        if( !covered ) {
            roots << abs
        }
    }
    return roots.collect { r -> "-v " + shellQuote("${r}:${r}") }.join(' ')
}

/*
 * RUN HISTORY (traceability of every successful run)
 *
 * At the START of every run the settings files are read into memory
 * (captureSettings). They are NOT moved: the latest Config.txt,
 * Config_data.txt, Config_evaluation.txt and Parameters.txt always stay in
 * their normal place (working folder) and the GUI keeps working with them.
 *
 * Only if the run COMPLETES SUCCESSFULLY, writeRunHistory() saves them to
 *
 *   WORKING_PATH/run_history/<yyyyMMdd_HHmmss>_<newrun | rerun_from_<step>>/
 *       Config.txt, Config_data.txt, Config_evaluation.txt, Parameters.txt
 *       (+ imputation_procedures.tsv if one is used)
 *       run_info.txt   when, which run, code version, image, SHA-256 of the data
 *       env_log.txt / env_log.pdf
 *       Parameter_Files/, logfile_setup.txt   parameters the modules used
 *       pipeline_log.txt     only the part written by this run
 *       nextflow_trace.txt   resources, retries and exit code per task
 *       report.pdf           final report of this run
 *
 * All of these are COPIES: the latest version always stays in the working
 * folder as well.
 *
 * The timestamp is the START time of the run, and the files are the
 * versions read at the start, so later changes in the GUI during the run do
 * not end up in the history. A failed or interrupted run writes nothing.
 *
 * Nothing in run_history/ is ever moved, overwritten or deleted by the
 * pipeline. A problem here never stops the pipeline; it prints a warning.
 */
/*
 * SHA-256 checksum of a file (hex), or a short reason if it cannot be read.
 * Proves later which exact data a run used, without copying the data.
 */
String sha256(File f) {
    if( f == null || !f.isFile() ) return 'not available (file not found)'
    try {
        def md = java.security.MessageDigest.getInstance('SHA-256')
        f.eachByte(1048576) { buf, len -> md.update(buf, 0, len) }
        return md.digest().encodeHex().toString()
    }
    catch( Exception e ) {
        return "not available (${e.message})"
    }
}

/* Copy a file or a whole folder (recursively) to dst. */
def copyTree(File src, File dst) {
    if( src.isDirectory() ) {
        dst.mkdirs()
        src.listFiles()?.each { child -> copyTree(child, new File(dst, child.name)) }
    }
    else if( src.isFile() ) {
        dst.parentFile?.mkdirs()
        dst.bytes = src.bytes
        dst.setLastModified(src.lastModified())
    }
}

List captureSettings(Map files) {
    def captured = []
    files.each { label, src ->
        def f = src == null ? null : (src instanceof File ? src : new File(src.toString()))
        if( f == null ) return
        def entry = [label: label.toString(), name: f.name, path: f.absolutePath, bytes: null, mtime: 0L]
        try {
            if( f.isFile() ) {
                entry.bytes = f.bytes
                entry.mtime = f.lastModified()
            }
        }
        catch( Exception e ) {
            println "[RUN HISTORY WARNING] Could not read ${f.absolutePath}: ${e.message}"
        }
        captured << entry
    }
    return captured
}

String writeRunHistory(String workdir, String folderName, List captured, List infoLines, List extraFiles) {
    try {
        def root = new File(workdir, 'run_history')
        def base = folderName.replaceAll(/[^A-Za-z0-9._-]+/, '_')
        // Two runs started within the same second get _1, _2, ... (never overwrite).
        def k = (0..999).find { n -> !new File(root, n == 0 ? base : "${base}_${n}").exists() }
        if( k == null ) {
            println "[RUN HISTORY WARNING] No free folder name for ${base}; run history not saved."
            return ''
        }
        def dir = new File(root, k == 0 ? base : "${base}_${k}")
        dir.mkdirs()
        if( !dir.isDirectory() ) {
            println "[RUN HISTORY WARNING] Could not create ${dir}; run history not saved."
            return ''
        }

        def fileLines = []
        captured.each { c ->
            if( c.bytes != null ) {
                def dst = new File(dir, c.name.toString())
                dst.bytes = c.bytes
                if( c.mtime ) dst.setLastModified(c.mtime as long)
                fileLines << envLine(c.label.toString(), c.path)
            }
            else {
                fileLines << envLine(c.label.toString(), "not found at the start of the run (${c.path})")
            }
        }
        // Extra files/folders of this run (copies; the originals stay).
        extraFiles.findAll { it != null && it.exists() }.each { f ->
            try {
                copyTree(f, new File(dir, f.name))
            }
            catch( Exception e ) {
                println "[RUN HISTORY WARNING] Could not copy ${f}: ${e.message}"
            }
        }

        def info = ['=' * 79, 'SETTINGS OF ONE SUCCESSFUL PIPELINE RUN', '=' * 79, ''] +
                   infoLines +
                   ['', 'Saved files (as read at the start of the run; original location):'] + fileLines +
                   ['',
                    'To repeat this run, or restart from the same step with the same settings:',
                    'copy Config.txt and Parameters.txt (and, if they differ, Config_data.txt,',
                    'Config_evaluation.txt) from this folder back to their original location',
                    '(or load them in the GUI), then start the pipeline again.',
                    '']
        new File(dir, 'run_info.txt').text = info.join('\n') + '\n'
        return dir.absolutePath
    }
    catch( Exception e ) {
        println "[RUN HISTORY WARNING] Run history could not be saved: ${e.message}"
        return ''
    }
}

/*
 * ENV_LOG
 *
 * Runs first, inside the same Docker image as all other steps, and writes
 * env_log.txt (+ env_log.pdf) into WORKING_PATH: "What did this particular run use?"
 *
 * - The run details (Nextflow version, pipeline git commit, Docker image
 *   name AND unique image ID, Docker version, resources, ...) are collected
 *   on the host by the workflow and passed in as run_header.
 * - The package versions are collected inside the container by the
 *   command "environment_versions", which the Dockerfile installs in the
 *   image and also uses to write image_versions.txt at build time, so both
 *   files have the same layout and can be compared directly.
 *
 * env_log.txt / env_log.pdf in WORKING_PATH always belong to the latest run.
 * A successful run also saves both in its run_history/ folder.
 */
process ENV_LOG {

    tag "ENV_LOG"

    maxForks 1

    publishDir { publish_root }, mode: 'copy', overwrite: true

    input:
    tuple val(publish_root),
          val(run_header),
          val(run_changes)

    output:
    path("env_log.*")

    script:
    """
    cat > run_header.txt <<'EOF_RUN_HEADER'
${run_header}
EOF_RUN_HEADER

    cat > run_changes.txt <<'EOF_RUN_CHANGES'
${run_changes}
EOF_RUN_CHANGES

    # The container this logging task actually runs in. All processes take
    # their image from the single setting process.container = params.container_image
    # in nextflow.config (no process overrides it), and the record compares
    # this value with the recorded image name, so a mismatch is visible.
    sed -i 's|@@TASK_CONTAINER@@|${task.container}|' run_header.txt || true

    # Writes env_log.txt (authoritative) and env_log.pdf (best effort).
    # Nothing here may stop the pipeline: every command is allowed to fail,
    # and env_log.txt is guaranteed to exist afterwards.
    if command -v environment_versions >/dev/null 2>&1; then
        environment_versions env_log.txt RUN run_header.txt run_changes.txt || true
    fi

    if [ ! -s env_log.txt ]; then
        # Image built before the reproducibility record was added to the Dockerfile.
        {
            printf '%s\n' "==============================================================================="
            printf '%s\n' "                                 PIPELINE RUN"
            printf '%s\n' "                            REPRODUCIBILITY RECORD"
            printf '%s\n' "==============================================================================="
            cat run_header.txt
            printf '\n%s\n' "WARNING: the package versions could not be recorded (this Docker image"
            printf '%s\n'   "has no working 'environment_versions' command). Rebuild the image with the"
            printf '%s\n'   "current Dockerfile. Run details above were recorded on the host."
            if grep -q '[^[:space:]]' run_changes.txt 2>/dev/null; then
                printf '\n%s\n\n' "UNCOMMITTED PIPELINE CHANGES"
                cat run_changes.txt
            fi
        } > env_log.txt || true
    fi
    """
}

process INIT_PARAMS {

    tag "INIT_PARAMS | attempt=${task.attempt} | mem=${task.memory} | cpus=${task.cpus}"

    maxForks 1

    // Mount WORKING_PATH if it lies outside the default Docker mount (see dockerMountOpts).
    containerOptions { container_mounts }

    publishDir { publish_root }, mode: 'copy', overwrite: true

    input:
    tuple val(publish_root),
          val(pipeline_log_path),
          path(init_script),
          path(flags_pkl),
          path(params_pkl),
          path(params_txt),
          path(original_data),
          path(meta_file),
          val(container_mounts)

    output:
    tuple path("Parameter_Files"), path("logfile_setup.txt"), path("init.done")

    script:
    """
    set -uo pipefail

    mkdir -p "${publish_root}"
    touch "${pipeline_log_path}"

    persistent_setup_log="${publish_root}/logfile_setup.txt"

    (
        while [ ! -f logfile_setup.txt ]; do
            sleep 0.1
        done

        tail -n 0 -F logfile_setup.txt
    ) >> "\$persistent_setup_log" 2>/dev/null &

    setup_tail_pid=\$!

    stop_setup_tail() {
        kill "\$setup_tail_pid" 2>/dev/null || true
        wait "\$setup_tail_pid" 2>/dev/null || true

        if [ -f logfile_setup.txt ]; then
            cp -f logfile_setup.txt "\$persistent_setup_log"
        fi
    }

    {
        printf '%s\n' ""
        printf '%s\n' "============================================================"
        printf '%s\n' "[INITIALIZATION START]"
        printf '%s\n' "[RESOURCE INFO] Attempt=${task.attempt} | CPUs=${task.cpus} | Memory=${task.memory}"
        printf '%s\n' "============================================================"
    } >> "${pipeline_log_path}"

    set +e

    # Staged input file names are inserted by Nextflow already shell-escaped
    # (e.g. my\\ metadata.xlsx), so they must NOT be wrapped in double quotes,
    # otherwise the backslash becomes part of the file name.
    python ${init_script} \
        -w "." \
        -m ${meta_file} \
        -d ${original_data} 2>&1 |
        tee -a "${pipeline_log_path}"

    py_exit=\${PIPESTATUS[0]}

    set -e

    if [ "\$py_exit" -ne 0 ]; then
        stop_setup_tail

        printf '%s\n' \
            "[FATAL] Initialize.py failed with exit code \$py_exit" \
            >> "${pipeline_log_path}"

        printf '%s\n' \
            "[FATAL] Initialize.py failed with exit code \$py_exit" >&2

        exit 100
    fi

    if [ ! -f logfile_setup.txt ]; then
        stop_setup_tail

        printf '%s\n' \
            "[FATAL] logfile_setup.txt was not created by Initialize.py" \
            >> "${pipeline_log_path}"

        printf '%s\n' \
            "[FATAL] logfile_setup.txt was not created by Initialize.py" >&2

        exit 100
    fi

    stop_setup_tail

    cat logfile_setup.txt >> "${pipeline_log_path}"

    if grep -Fq '###ERROR###' logfile_setup.txt; then
        printf '%s\n' \
            "[FATAL] Initialize.py produced ###ERROR### in logfile_setup.txt" \
            >> "${pipeline_log_path}"

        printf '%s\n' \
            "[FATAL] Initialize.py produced ###ERROR### in logfile_setup.txt" >&2

        exit 100
    fi

    if [ ! -d Parameter_Files ]; then
        printf '%s\n' \
            "[FATAL] Parameter_Files was not created by Initialize.py" \
            >> "${pipeline_log_path}"

        printf '%s\n' \
            "[FATAL] Parameter_Files was not created by Initialize.py" >&2

        exit 100
    fi

    printf '%s\n' "[INITIALIZATION STATUS] PASSED" |
        tee -a "${pipeline_log_path}"

    printf '%s\n' "[INITIALIZATION END]" |
        tee -a "${pipeline_log_path}"

    touch init.done
    """
}

process RUN_ALL_MODULES {

    tag "RUN_ALL_MODULES | attempt=${task.attempt} | mem=${task.memory} | cpus=${task.cpus}"

    maxForks 1

    // Mount WORKING_PATH if it lies outside the default Docker mount (see dockerMountOpts).
    containerOptions { container_mounts }

    /*
     * Retry only memory/preemption-related termination.
     *
     * All other nonzero exits terminate the workflow.
     * QA failures are handled internally and do not cause
     * a nonzero process exit.
     */
    errorStrategy {
        task.exitStatus in [137, 143] ? 'retry' : 'terminate'
    }

    beforeScript """
    echo "[RESOURCE INFO] Attempt ${task.attempt} | CPUs=${task.cpus} | Memory=${task.memory}"
    """

    input:
    tuple val(publish_root),
          val(steps_tsv),
          path(code_py_dir),
          path(code_r_dir),
          path(report_script),
          path(qa_script),
          path(eval_cfg),
          path(figscptn),
          path(seed_data),
          path(meta_file),
          path(id_file),
          path(mask_file),
          path(mask_dp_file),
          val(pipeline_log_path),
          path(prev_report),
          path(doc_tsv),
          path(paramdir),
          val(imputation_procedures_input),
          val(container_mounts)

    output:
    path("pipeline.done")
    path("qa_status.tsv")

    script:
    """
    set -euo pipefail

    BASE="\$PWD"

    cat > steps.tsv <<'EOF_STEPS'
${steps_tsv}
EOF_STEPS

    # Staged input file names are inserted by Nextflow already shell-escaped
    # (e.g. my\\ data.tsv), so they are placed OUTSIDE the double quotes:
    # "\$BASE/"${seed_data} gives one argument, even with spaces in the name.
    current_data="\$BASE/"${seed_data}
    current_id="\$BASE/"${id_file}
    current_mask="\$BASE/"${mask_file}
    current_mask_dp="\$BASE/"${mask_dp_file}
    current_report="\$BASE/"${prev_report}


    # The one cumulative log for the whole workflow.
    # pipeline_log_path is the absolute WORKING_PATH/pipeline_log.txt path.
    persistent_log="${pipeline_log_path}"

    mkdir -p "\$(dirname "\$persistent_log")"
    touch "\$persistent_log"

    # Append messages to the current module log and the one pipeline log.
    log_both() {
        printf '%s\\n' "\$*" |
            tee -a \
                step_log.txt \
                "\$persistent_log" >/dev/null
    }

    # Publish currently available artifacts for the active module.

    #

    # This is called after successful execution and before fatal exits,

    # so diagnostics remain available when the process fails.
    publish_current_step() {
        local publish_step="\$1"
        local counter_value="\$2"

        mkdir -p "\$publish_step"

        cp -f step_log.txt "\$publish_step/" || true
        cp -f log.txt "\$publish_step/" || true
        cp -f report.pdf "\$publish_step/" || true
        cp -rf report_folder "\$publish_step/" || true
        cp -f "data_\${counter_value}.tsv" "\$publish_step/" || true
        cp -f id_mapping.tsv "\$publish_step/" || true
        cp -f "data_\${counter_value}_mask.tsv" "\$publish_step/" || true
        cp -f "data_\${counter_value}_mask_description.tsv" \
            "\$publish_step/" || true
        cp -f qa_status.tsv "\$publish_step/" || true

        # Normalization factors: only published by the step that created
        # them (normalization). The file is removed from the task folder at
        # the end of that step, so later steps do not carry a copy.
        # An empty file (no factors provided/created) is not published.
        if [ -s Normalization_factors.tsv ]; then
            cp -f Normalization_factors.tsv "\$publish_step/" || true
        fi

        # Covariate associations: written occasionally by the correction
        # module (when adjust_confounders is used). Same rules as the
        # normalization factors: published only in the step that created
        # it, not if empty, and removed afterwards (see end of the step).
        if [ -s covariate_associations.tsv ]; then
            cp -f covariate_associations.tsv "\$publish_step/" || true
        fi

        case "\$publish_step" in
            *IMPUTATION*)
                if [ -f imputation_procedures.tsv ]; then
                    cp -f imputation_procedures.tsv "\$publish_step/" || true
                fi
                ;;
        esac
    }

    # Pipeline-level QA status table.
    printf 'module\\tqa_status\\texit_code\\n' > qa_status.tsv

    data_two_back=""
    mask_two_back=""
    mask_dp_two_back=""
    data_one_back="\$current_data"
    mask_one_back="\$current_mask"
    mask_dp_one_back="\$current_mask_dp"

    cleanup_data_bundle() {
        for input_path in "\$1" "\$2" "\$3"; do
            if [ -n "\$input_path" ] &&
               [ "\$input_path" != "\$current_data" ] &&
               [ "\$input_path" != "\$current_mask" ] &&
               [ "\$input_path" != "\$current_mask_dp" ] &&
               [ -f "\$input_path" ]; then
                rm -f -- "\$input_path"
                log_both "[CLEANUP] Removed \$(basename "\$input_path")"
            fi
        done
    }

    log_both ""
    log_both "============================================================"
    log_both "[NEXTFLOW PROCESS] RUN_ALL_MODULES"
    log_both "[NEXTFLOW ATTEMPT] ${task.attempt}"
    log_both "[RESOURCE INFO] CPUs=${task.cpus} | Memory=${task.memory}"
    log_both "============================================================"

    rm -f "\$BASE/imputation_procedures.tsv"
    if [ -n "${imputation_procedures_input}" ]; then
        cp -f "${imputation_procedures_input}" "\$BASE/imputation_procedures.tsv"
        log_both "[IMPUTATION] Using provided imputation_procedures.tsv"
    else
        log_both "[IMPUTATION] No imputation_procedures.tsv provided; generation is enabled."
    fi

    while IFS=\$'\\t' read -r counter base engine script_name param_name; do

        # Ignore completely empty rows.
        if [ -z "\${counter}" ]; then
            continue
        fi

        step_name="\${base}"
        module_name="\${counter}_\${step_name}"

        cd "\$BASE"

        # Remove marker files from the previous module.
        rm -f MOD.ok MOD.fail
        rm -f QA.ok QA.fail QA.skip
        rm -f REP.ok REP.fail

        rm -rf report_folder
        mkdir -p report_folder

        # Clear only the per-module log.

        #

        # The cumulative pipeline logs are never cleared.
        : > step_log.txt

        publish_step="${publish_root}/\${module_name}"

        log_both ""
        log_both "============================================================"
        log_both "===== \${module_name} (\${engine}) ====="
        log_both "============================================================"
        log_both "[MODULE INFO] Running \${module_name}"
        log_both "----- MODULE START (\${module_name}) -----"

        # Resolve module runner and script.
        if [ "\$engine" = "py" ]; then
            script_path="\$BASE/"${code_py_dir}"/\${script_name}"
            runner="python"

        elif [ "\$engine" = "r" ] ||
             [ "\$engine" = "R" ]; then

            script_path="\$BASE/"${code_r_dir}"/\${script_name}"
            runner="Rscript"

        else
            log_both "###MOD_ERROR### Unsupported engine '\${engine}' for \${module_name}"
            log_both "[MODULE STATUS] FAILED"
            log_both "[QA STATUS] SKIPPED because the engine is unsupported."
            log_both "----- MODULE END -----"

            touch MOD.fail
            touch QA.skip

            printf '%s\\t%s\\t%s\\n' \
                "\$module_name" "SKIPPED" "NA" >> qa_status.tsv

            cp -f step_log.txt log.txt || true
            publish_current_step "\$publish_step" "\$counter"
            printf '%s\\n' "[FATAL] \${module_name}: unsupported engine '\${engine}'" >> "\$persistent_log"
            exit 100
        fi

        param_file="\$BASE/"${paramdir}"/\${param_name}"

        # Validate module script.
        if [ ! -f "\$script_path" ]; then
            log_both "###MOD_ERROR### Missing module script: \$script_path"
            log_both "[MODULE STATUS] FAILED"
            log_both "[QA STATUS] SKIPPED because the module script is missing."
            log_both "----- MODULE END -----"

            touch MOD.fail
            touch QA.skip

            printf '%s\\t%s\\t%s\\n' \
                "\$module_name" "SKIPPED" "NA" >> qa_status.tsv

            cp -f step_log.txt log.txt || true
            publish_current_step "\$publish_step" "\$counter"
            printf '%s\\n' "[FATAL] \${module_name}: module script is missing" | tee -a "\$persistent_log" >&2
            exit 100
        fi

        # Validate module parameter file.
        if [ ! -f "\$param_file" ]; then
            log_both "###MOD_ERROR### Missing parameter file: \$param_file"
            log_both "[MODULE STATUS] FAILED"
            log_both "[QA STATUS] SKIPPED because the parameter file is missing."
            log_both "----- MODULE END -----"

            touch MOD.fail
            touch QA.skip

            printf '%s\\t%s\\t%s\\n' \
                "\$module_name" "SKIPPED" "NA" >> qa_status.tsv

            cp -f step_log.txt log.txt || true
            publish_current_step "\$publish_step" "\$counter"
            printf '%s\\n' "[FATAL] \${module_name}: parameter file is missing" | tee -a "\$persistent_log" >&2
            exit 100
        fi

        rm -f ./*_Parameters.json
        cp -f "\$param_file" .

        # Record where module-specific output begins.

        #

        # This allows generic ###ERROR### markers to be attributed

        # specifically to the module phase.
        module_log_start=\$(wc -l < step_log.txt)

        # Execute module.

        #

        # Output is streamed immediately to:

        #

        # - step_log.txt

        # - work/.../pipeline_log.txt

        # - WORKDIR_ABS/pipeline_log.txt
        set +e

        "\$runner" "\$script_path" \
            -d "\$current_data" \
            -p "\$param_file" \
            -o . \
            -m "\$BASE/"${meta_file} \
            -l step_log.txt \
            -i "\$current_id" \
            -c "\$counter" 2>&1 |
            tee -a \
                step_log.txt \
                "\$persistent_log"

        mod_rc=\${PIPESTATUS[0]}

        set -e

        # Inspect only lines produced during module execution.
        module_marker_found=0

        if tail -n "+\$((module_log_start + 1))" step_log.txt |
            grep -Eq "###MOD_ERROR###|###ERROR###"; then

            module_marker_found=1
        fi

        if [ "\$mod_rc" -ne 0 ] ||
           [ "\$module_marker_found" -eq 1 ]; then

            if [ "\$mod_rc" -ne 0 ]; then
                log_both "###MOD_ERROR### module failed exit=\$mod_rc for \${module_name}"
            else
                log_both "###MOD_ERROR### module produced an error marker for \${module_name}"
            fi

            log_both "[MODULE STATUS] FAILED"
            log_both "[QA STATUS] SKIPPED because module output may be incomplete."
            log_both "----- MODULE END -----"

            touch MOD.fail
            touch QA.skip

            printf '%s\\t%s\\t%s\\n' \
                "\$module_name" "SKIPPED" "NA" >> qa_status.tsv

            if [ "\$mod_rc" -eq 137 ] ||
               [ "\$mod_rc" -eq 143 ]; then

                log_both "[MEMORY ERROR] \${module_name} was killed with exit \$mod_rc."
                log_both "[MEMORY ERROR] The assigned memory may be insufficient."
                log_both "[MEMORY ERROR] Nextflow will retry if retries remain."
            fi

            cp -f step_log.txt log.txt || true
            publish_current_step "\$publish_step" "\$counter"
            if [ "\$mod_rc" -ne 0 ]; then
                printf '%s\\n' "[FATAL] \${module_name}: module failed with exit \$mod_rc" | tee -a "\$persistent_log" >&2
                exit "\$mod_rc"
            else
                printf '%s\\n' "[FATAL] \${module_name}: module error marker found" | tee -a "\$persistent_log" >&2
                exit 100
            fi
        fi

        touch MOD.ok

        log_both "[MODULE STATUS] PASSED"
        log_both "----- MODULE END -----"

        # QA phase.

        #

        # QA failures are always non-fatal.

        #

        # A QA failure may be detected by:

        #

        # - a nonzero QA exit status

        # - ###QA_ERROR### in QA output

        # - generic ###ERROR### in QA output

        #

        # None of these conditions causes an exit.
        log_both ""
        log_both "----- QA START (\${module_name}) -----"

        qa_log_start=\$(wc -l < step_log.txt)

        set +e

        python "\$BASE/"${qa_script} \
            -w . \
            -m "\$BASE/"${meta_file} \
            -s "\${module_name}" \
            -e "\$BASE/"${eval_cfg} \
            --docker 2>&1 |
            tee -a \
                step_log.txt \
                "\$persistent_log"

        qa_rc=\${PIPESTATUS[0]}

        set -e

        # Inspect only output generated during QA.
        qa_marker_found=0

        if tail -n "+\$((qa_log_start + 1))" step_log.txt |
            grep -Eq "###QA_ERROR###|###ERROR###"; then

            qa_marker_found=1
        fi

        if [ "\$qa_rc" -ne 0 ] ||
           [ "\$qa_marker_found" -eq 1 ]; then

            log_both "###QA_ERROR### QA failed for \${module_name}; exit=\$qa_rc"
            log_both "[QA STATUS] FAILED for module \${module_name}"
            log_both "[QA WARNING] QA failure is non-fatal."
            log_both "[QA WARNING] Report generation will continue."
            log_both "[QA WARNING] Pipeline execution will continue with the next module."

            touch QA.fail
            rm -f QA.ok QA.skip

            printf '%s\\t%s\\t%s\\n' \
                "\$module_name" "FAILED" "\$qa_rc" >> qa_status.tsv

        else
            touch QA.ok
            rm -f QA.fail QA.skip

            log_both "[QA STATUS] PASSED for module \${module_name}"

            printf '%s\\t%s\\t%s\\n' \
                "\$module_name" "PASSED" "0" >> qa_status.tsv
        fi

        log_both "----- QA END -----"

        # Report generation.

        #

        # Report failures remain fatal.

        #

        # Record the report start line so generic ###ERROR### markers

        # produced during QA cannot be mistaken for report errors.
        log_both ""
        log_both "----- REPORT START (\${module_name}) -----"

        report_log_start=\$(wc -l < step_log.txt)

        set +e

        if [ -f "\$current_report" ] &&
           [ -s "\$current_report" ]; then

            python "\$BASE/"${report_script} \
                --log step_log.txt \
                --in_pdf "\$current_report" \
                --out report_new.pdf \
                --doc_tsv "\$BASE/"${doc_tsv} 2>&1 |
                tee -a \
                    step_log.txt \
                    "\$persistent_log"

            rep_rc=\${PIPESTATUS[0]}

        else
            python "\$BASE/"${report_script} \
                --log step_log.txt \
                --out report_new.pdf \
                --doc_tsv "\$BASE/"${doc_tsv} 2>&1 |
                tee -a \
                    step_log.txt \
                    "\$persistent_log"

            rep_rc=\${PIPESTATUS[0]}
        fi

        set -e

        # Inspect only output generated during report creation.
        report_marker_found=0

        if tail -n "+\$((report_log_start + 1))" step_log.txt |
            grep -Eq "###REP_ERROR###|###ERROR###"; then

            report_marker_found=1
        fi

        if [ "\$rep_rc" -ne 0 ] ||
           [ "\$report_marker_found" -eq 1 ]; then

            if [ "\$rep_rc" -ne 0 ]; then
                log_both "###REP_ERROR### report failed exit=\$rep_rc for \${module_name}"
            else
                log_both "###REP_ERROR### report produced an error marker for \${module_name}"
            fi

            log_both "[REPORT STATUS] FAILED"
            log_both "----- REPORT END -----"

            touch REP.fail
            rm -f REP.ok

            cp -f step_log.txt log.txt || true
            publish_current_step "\$publish_step" "\$counter"
            if [ "\$rep_rc" -ne 0 ]; then
                printf '%s\\n' "[FATAL] \${module_name}: report failed with exit \$rep_rc" | tee -a "\$persistent_log" >&2
                exit "\$rep_rc"
            else
                printf '%s\\n' "[FATAL] \${module_name}: report error marker found" | tee -a "\$persistent_log" >&2
                exit 100
            fi
        fi

        if [ ! -f report_new.pdf ] ||
           [ ! -s report_new.pdf ]; then

            log_both "###REP_ERROR### report_new.pdf was not created or is empty for \${module_name}"
            log_both "[REPORT STATUS] FAILED"
            log_both "----- REPORT END -----"

            touch REP.fail
            rm -f REP.ok

            cp -f step_log.txt log.txt || true
            publish_current_step "\$publish_step" "\$counter"
            printf '%s\\n' "[FATAL] \${module_name}: generated report is missing or empty" | tee -a "\$persistent_log" >&2
            exit 100
        fi

        mv -f report_new.pdf report.pdf

        touch REP.ok
        rm -f REP.fail

        log_both "[REPORT STATUS] PASSED"
        log_both "----- REPORT END -----"

        # Keep the cumulative report in the working directory up to date.
        cp -f report.pdf "${publish_root}/report.pdf.tmp"
        mv -f "${publish_root}/report.pdf.tmp" "${publish_root}/report.pdf"

        # Keep a module-local copy of the complete step log.
        cp -f step_log.txt log.txt

        # If QA failed, report the affected module as a warning only.

        #

        # There is intentionally no exit here.
        if [ -f QA.fail ]; then
            log_both "[PIPELINE WARNING] \${module_name} completed with a QA problem."
            log_both "[PIPELINE WARNING] The QA problem did not stop the pipeline."
            log_both "[PIPELINE WARNING] See QA.fail, step_log.txt and qa_status.tsv."
        fi

        # Publish the module artifacts after all warning messages

        # have been added to step_log.txt.
        cp -f step_log.txt log.txt
        publish_current_step "\$publish_step" "\$counter"

        # Promote this module's outputs to become inputs for

        # the next module.
        current_data="\$BASE/data_\${counter}.tsv"
        current_id="\$BASE/id_mapping.tsv"
        current_report="\$BASE/report.pdf"

        if [ -f "\$BASE/data_\${counter}_mask.tsv" ]; then
            current_mask="\$BASE/data_\${counter}_mask.tsv"
        fi

        if [ -f "\$BASE/data_\${counter}_mask_description.tsv" ]; then
            current_mask_dp="\$BASE/data_\${counter}_mask_description.tsv"
        fi

        # Normalization_factors.tsv is needed only by the step that created it
        # and its QA (normalization evaluation). It is already published in
        # that step's folder, so remove it here: no copy in later steps.
        rm -f "\$BASE/Normalization_factors.tsv"
        # Same for covariate_associations.tsv (correction step only).
        rm -f "\$BASE/covariate_associations.tsv"

        # Safety synchronization at the end of every module iteration.
        cd "\$BASE"

        cleanup_data_bundle \
            "\$data_two_back" \
            "\$mask_two_back" \
            "\$mask_dp_two_back"

        data_two_back="\$data_one_back"
        mask_two_back="\$mask_one_back"
        mask_dp_two_back="\$mask_dp_one_back"
        data_one_back="\$current_data"
        mask_one_back="\$current_mask"
        mask_dp_one_back="\$current_mask_dp"

    done < steps.tsv

    # Remove the last retained previous-generation data files.
    cleanup_data_bundle \
        "\$data_two_back" \
        "\$mask_two_back" \
        "\$mask_dp_two_back"

    # Pipeline-level QA summary.

    #

    # QA failures remain non-fatal.
    qa_failure_count=\$(awk -F '\\t' '
        NR > 1 && \$2 == "FAILED" {
            count++
        }
        END {
            print count + 0
        }
    ' qa_status.tsv)

    qa_skip_count=\$(awk -F '\\t' '
        NR > 1 && \$2 == "SKIPPED" {
            count++
        }
        END {
            print count + 0
        }
    ' qa_status.tsv)

    qa_pass_count=\$(awk -F '\\t' '
        NR > 1 && \$2 == "PASSED" {
            count++
        }
        END {
            print count + 0
        }
    ' qa_status.tsv)

    printf '%s\\n' "" |
        tee -a "\$persistent_log"

    printf '%s\\n' "============================================================" |
        tee -a "\$persistent_log"

    printf '%s\\n' "PIPELINE QA SUMMARY" |
        tee -a "\$persistent_log"

    printf '%s\\n' "QA passed:  \$qa_pass_count" |
        tee -a "\$persistent_log"

    printf '%s\\n' "QA failed:  \$qa_failure_count" |
        tee -a "\$persistent_log"

    printf '%s\\n' "QA skipped: \$qa_skip_count" |
        tee -a "\$persistent_log"

    if [ "\$qa_failure_count" -gt 0 ]; then

        printf '%s\\n' \
            "###QA_WARNING_SUMMARY### One or more modules failed QA." |
            tee -a "\$persistent_log"

        printf '%s\\n' \
            "[QA WARNING] The pipeline completed, but QA problems were detected in:" |
            tee -a "\$persistent_log"

        awk -F '\\t' '
            NR > 1 && \$2 == "FAILED" {
                printf "[QA WARNING]   - %s, exit code %s\\n", \$1, \$3
            }
        ' qa_status.tsv |
            tee -a "\$persistent_log"

        printf '%s\\n' \
            "[QA WARNING] QA failures were non-fatal and did not stop the pipeline." |
            tee -a "\$persistent_log"

    else
        printf '%s\\n' \
            "[QA SUMMARY] No QA failures were detected." |
            tee -a "\$persistent_log"
    fi

    printf '%s\\n' "============================================================" |
        tee -a "\$persistent_log"

    # Final synchronization before successful completion.
    # Created only when:
    #
    # - all modules completed
    # - all reports completed
    # - QA failures, if present, were recorded as non-fatal
    printf '%s\\n' "[PIPELINE STATUS] COMPLETED" |
        tee -a "\$persistent_log"

    # Required by the declared Nextflow output.
    printf '%s\\n' "completed" > pipeline.done

    # Persistent success marker used by onComplete and the launcher.
    printf '%s\\n' "completed" > "${publish_root}/pipeline.done"
    """
}

workflow {

    main:

    def config_data = pget('config_data', 'Config_data.txt')
    def code_py = pget('code_py', 'modules_python')
    def code_r = pget('code_r', 'modules_R')
    def method_doc = pget('method_doc', 'docs/Method_description.tsv')
    def qa_script_param = pget('qa_script', 'modules_python/module_quality_assessment.py')
    def eval_config_param = pget('eval_config', 'Config_evaluation.txt')
    def figs_cptn_param = pget('figs_cptn', 'modules_python/Figures_Captions.tsv')
    def cleanup_aggressive = pget('cleanup_aggressive', false)

    cfgDataPath = resolveLaunch(config_data as String)
    cfgData = parseKeyValueFile(cfgDataPath)

    WORKING_PATH = cfgData.get('WORKING_PATH') ?: '.'
    WORKDIR_ABS = resolveFrom(launchDir.toString(), WORKING_PATH)
    new File(WORKDIR_ABS, 'pipeline.done').delete()

    CONFIG_TXT_ABS = new File(WORKDIR_ABS, 'Config.txt').absolutePath

    ORIGINALDATA_ABS = cfgData.get('ORIGINALDATA') ? resolveFrom(WORKDIR_ABS, cfgData['ORIGINALDATA']) : null
    METADATA_ABS = cfgData.get('METADATA') ? resolveFrom(WORKDIR_ABS, cfgData['METADATA']) : null

    def imputationProceduresValue = cfgData.get('IMPUTATION_PROCEDURES')?.toString()?.trim() ?: ''
    IMPUTATION_PROCEDURES_ABS = imputationProceduresValue ?
        resolveFrom(WORKDIR_ABS, imputationProceduresValue) : null

    imputationProceduresFile = IMPUTATION_PROCEDURES_ABS ?
        mustExist(IMPUTATION_PROCEDURES_ABS, 'imputation_procedures.tsv') : null

    cfgFile = mustExist(CONFIG_TXT_ABS, 'Config.txt')

    lines = cfgFile.readLines()

    // Match the keys exactly ("Run =" / "Order ="), not just lines that
    // start with these letters (e.g. "run_error_check = 1" or "Runtime = ...").
    run_value = lines.find {
        it ==~ /(?i)^\s*run\s*=.*/
    }?.split('=', 2)?.getAt(1)?.trim()

    if( !run_value ) {
        die("Config file missing 'Run =' entry")
    }

    order_line = lines.find {
        it ==~ /(?i)^\s*order\s*=.*/
    }

    if( !order_line ) {
        die("Config file missing 'Order =' entry")
    }

    ordered = order_line.split('=', 2)[1]
        .replaceAll(/\s+/, '')
        .split(',')
        .findAll { it }
        .toList()

    println "Run mode          : ${run_value}"
    println "Ordered processes : ${ordered}"

    if( run_value.toLowerCase() == 'new' ) {
        start_dir = WORKDIR_ABS
        precounter = '00'
    }
    else {
        rerun_name = run_value.replaceFirst(/_(?i)(py|r)$/, '')
        rerun_folder = new File(WORKDIR_ABS, rerun_name)

        if( !rerun_folder.isDirectory() ) {
            die("Rerun folder '${rerun_name}' not found: ${rerun_folder}")
        }

        stepLog = new File(rerun_folder, 'step_log.txt')
        if( stepLog.exists() ) {
            bad = stepLog.readLines().any { it.contains('###ERROR###') }
            if( bad ) {
                die("Cannot rerun from '${rerun_name}' because step_log.txt contains ###ERROR###")
            }
        }
        else {
            println "[WARN] step_log.txt not found in ${rerun_folder}"
        }

        start_dir = rerun_folder.absolutePath
        precounter = rerun_name.split('_', 2)[0]
    }

    println "Start dir         : ${start_dir}"
    println "Precounter        : ${precounter}"

    seed_data_path = run_value.equalsIgnoreCase('new')
        ? ORIGINALDATA_ABS
        : new File(start_dir, "data_${precounter}.tsv").absolutePath

    seed_data = mustExist(seed_data_path, 'seed data')
    original_data = mustExist(ORIGINALDATA_ABS, 'ORIGINALDATA')

    seed_idmap_path = new File(start_dir, 'id_mapping.tsv').absolutePath
    seed_mask_path = new File(start_dir, "data_${precounter}_mask.tsv").absolutePath
    seed_maskdp_path = new File(start_dir, "data_${precounter}_mask_description.tsv").absolutePath

    seed_mask = ensureExistsOrCreateEmpty(seed_mask_path, 'seed mask')
    seed_maskdp = ensureExistsOrCreateEmpty(seed_maskdp_path, 'seed mask description')

    seed_idmap = run_value.equalsIgnoreCase('new')
        ? ensureExistsOrCreateEmpty(seed_idmap_path, 'seed idmap')
        : mustExist(seed_idmap_path, 'seed idmap')

    PIPELINE_LOG_ABS = new File(WORKDIR_ABS, 'pipeline_log.txt').absolutePath
    ensureExistsOrCreateEmpty(PIPELINE_LOG_ABS, 'pipeline_log.txt')

    // Extra Docker mounts for folders outside params.docker_mount_root.
    CONTAINER_MOUNTS = dockerMountOpts(
        [
            WORKDIR_ABS,
            IMPUTATION_PROCEDURES_ABS ? new File(IMPUTATION_PROCEDURES_ABS).parent : null
        ],
        params.docker_mount_root ? params.docker_mount_root.toString() : ''
    )

    if( CONTAINER_MOUNTS ) {
        println "Extra Docker mounts: ${CONTAINER_MOUNTS}"
    }

    def reportSeedPath = new File(start_dir, 'report.pdf').absolutePath

    if( run_value.equalsIgnoreCase('new') ) {
        def oldReport = new File(reportSeedPath)

        if( oldReport.exists() && oldReport.length() > 0 ) {

            def archiveDir = new File(start_dir, 'previous_reports')
            archiveDir.mkdirs()

            def timestamp = new Date().format('yyyyMMdd_HHmmss')
            def archivedReport = new File(
                archiveDir,
                "report_${timestamp}.pdf"
            )

            println "Archiving previous complete report to: ${archivedReport}"

            if( !oldReport.renameTo(archivedReport) ) {
                die("Could not archive previous report.pdf")
            }
        }
    }

    reportSeed = ensureExistsOrCreateEmpty(
        reportSeedPath,
        'report seed'
    )

    metaFile = mustExist(METADATA_ABS, 'Metadata.tsv')

    methodDoc = mustExist(resolveCode(method_doc as String), 'method doc TSV')

    reportPy = mustExist(resolveCode("${code_py}/report.py"), 'modules_python/report.py')
    qaScript = mustExist(resolveCode(qa_script_param as String), 'QA script')
    evalConfig = mustExist(resolveFrom(WORKDIR_ABS, eval_config_param as String), 'Config_evaluation.txt')

    figscptn = mustExist(resolveCode(figs_cptn_param as String), 'Figures Captions')
    codePyDir = mustExist(resolveCode(code_py as String), 'modules_python')
    codeRDir = mustExist(resolveCode(code_r as String), 'modules_R')

    // Cloud-sync safeguard: download any online-only inputs before Docker starts.
    println "Checking inputs for online-only cloud files ..."
    def cloudStats = [files: 0, bytes: 0L, slow: []]
    [
        'ORIGINALDATA'          : original_data,
        'seed data'             : seed_data,
        'METADATA'              : metaFile,
        'seed idmap'            : seed_idmap,
        'seed mask'             : seed_mask,
        'seed mask description' : seed_maskdp,
        'report seed'           : reportSeed,
        'method doc'            : methodDoc,
        'report.py'             : reportPy,
        'QA script'             : qaScript,
        'Config_evaluation.txt' : evalConfig,
        'Figures Captions'      : figscptn,
        'IMPUTATION_PROCEDURES' : imputationProceduresFile,
        'Parameters.txt'        : new File("${WORKDIR_ABS}/Parameters.txt"),
        'modules_python'        : codePyDir,
        'modules_R'             : codeRDir,
        'Initialization'        : new File("${projectDir}/Initialization")
    ].each { label, f -> materialize(f, label, cloudStats) }

    println "Inputs checked    : ${cloudStats.files} files, ${cloudStats.bytes.intdiv(1048576)} MB"
    if( cloudStats.slow ) {
        println "[CLOUD] Downloaded / slow-read files (were likely online-only):"
        cloudStats.slow.each { println "[CLOUD]   ${it}" }
    }


    stepsResolved = []

    ordered.each { item ->

        def parts = item.split('_', 2)
        def counter = parts[0]
        def rest = parts[1]

        def suffix = rest.tokenize('_')[-1].toLowerCase()
        def engine = suffix == 'r' ? 'r' : 'py'
        def base = rest.replaceAll(/_(?i)(py|r)$/, '')

        def codeDir = engine == 'py' ? codePyDir : codeRDir
        def ext = engine == 'py' ? 'py' : 'R'

        def targetName = "module_${base}.${ext}".toLowerCase()

        def scriptFile = codeDir.listFiles()?.find {
            it.name.toLowerCase() == targetName
        }

        if( !scriptFile ) {
            die("Missing step script for ${item} in ${codeDir}")
        }

        def paramName = "${counter}_${base}_Parameters.json"

        stepsResolved << [
            counter: counter,
            base: base,
            engine: engine,
            scriptName: scriptFile.name,
            paramName: paramName
        ]
    }

    steps_tsv = stepsResolved.collect {
        "${it.counter}\t${it.base}\t${it.engine}\t${it.scriptName}\t${it.paramName}"
    }.join('\n')

    println "Working path      : ${WORKDIR_ABS}"
    println "Original data     : ${ORIGINALDATA_ABS}"
    println "Metadata          : ${METADATA_ABS}"
    println "Steps             : ${stepsResolved.size()}"

    // ------------------------------------------------------------------
    // Environment log (env_log.txt): run details collected on the host.
    // ------------------------------------------------------------------
    def imageName = params.container_image ? params.container_image.toString() : 'not set (params.container_image missing in nextflow.config)'
    def gitDir = projectDir.toString()
    def gitCommit = captureCommand(['git', '-C', gitDir, 'rev-parse', 'HEAD'])
    def gitBranch = captureCommand(['git', '-C', gitDir, 'rev-parse', '--abbrev-ref', 'HEAD'])
    def gitRemote = captureCommand(['git', '-C', gitDir, 'config', '--get', 'remote.origin.url'])
    def gitChanges = captureCommand(['git', '-C', gitDir, 'status', '--porcelain'])

    // If the code has uncommitted changes, preserve the tracked git diff
    // (subject to the size limit) and list untracked files by name.
    def runChanges = ''
    if( gitChanges ) {
        def diffText = captureCommand(['git', '-C', gitDir, 'diff', 'HEAD', '--no-color', '--no-ext-diff'])
        def untracked = captureCommand(['git', '-C', gitDir, 'ls-files', '--others', '--exclude-standard'])
        def maxChars = 5000000
        if( diffText.length() > maxChars ) {
            diffText = diffText.substring(0, maxChars) +
                "\n[... diff truncated after ${maxChars} characters ...]"
        }
        runChanges = "[git status --porcelain]\n${gitChanges}\n\n" +
                     "[git diff HEAD]\n${diffText ?: '(no changes to tracked files)'}\n\n" +
                     "[Untracked files (content not included)]\n${untracked ?: '(none)'}"
    }
    def imageId = captureCommand(['docker', 'image', 'inspect', '--format', '{{.Id}}', imageName])
    def imageDigests = captureCommand(['docker', 'image', 'inspect', '--format', '{{join .RepoDigests ", "}}', imageName])
    def imageCreated = captureCommand(['docker', 'image', 'inspect', '--format', '{{.Created}}', imageName])
    def imagePlatform = captureCommand(['docker', 'image', 'inspect', '--format', '{{.Os}}/{{.Architecture}}', imageName])
    def dockerVersion = captureCommand(['docker', 'version', '--format', '{{.Server.Version}}'])
    def dockerSystem = captureCommand(['docker', 'info', '--format', '{{.OperatingSystem}} | {{.NCPU}} CPUs | {{.MemTotal}} bytes memory'])

    // ------------------------------------------------------------------
    // Run history: read the settings of THIS run now (start); they are saved
    // to run_history/<start time>_<newrun|rerun_from_step>/ only if the run
    // completes successfully (see workflow.onComplete). Nothing is moved.
    // ------------------------------------------------------------------
    def historyLabel = run_value.equalsIgnoreCase('new') ? 'newrun' : "rerun_from_${new File(start_dir.toString()).name}"
    def historyFolderName = new Date().format('yyyyMMdd_HHmmss') + '_' + historyLabel
    def capturedSettings = captureSettings([
        'Config.txt'            : cfgFile,
        'Config_data.txt'       : new File(cfgDataPath),
        'Config_evaluation.txt' : evalConfig,
        'Parameters.txt'        : new File(WORKDIR_ABS, 'Parameters.txt'),
        'Imputation procedures' : imputationProceduresFile
    ])
    // Checksums of the data this run starts from (as read at the start).
    def isNewRun = run_value.equalsIgnoreCase('new')
    def dataChecksums = [
        envLine('SHA-256 ORIGINALDATA', sha256(original_data)),
        envLine('SHA-256 METADATA', sha256(metaFile))
    ]
    if( !isNewRun ) {
        dataChecksums << envLine('SHA-256 seed data', "${sha256(seed_data)} (${seed_data.name})")
    }
    // pipeline_log.txt grows over all runs: remember where this run begins.
    def pipelineLogStart = new File(PIPELINE_LOG_ABS).length()

    def historyInfo = [
        envLine('Run started', workflow.start),
        envLine('Run name', workflow.runName),
        envLine('Run mode', run_value),
        envLine('Start folder', start_dir),
        envLine('Steps', ordered.join(', ')),
        envLine('Working path', WORKDIR_ABS),
        envLine('Command line', workflow.commandLine),
        envLine('Git commit', gitCommit ?: 'unknown (not a git checkout)'),
        envLine('Uncommitted changes', gitChanges ? 'YES (the diff is in env_log.txt)' : 'NONE'),
        envLine('Image', imageName),
        envLine('Image ID', imageId ?: 'unknown (image not found locally)'),
        envLine('Original data', ORIGINALDATA_ABS),
        envLine('Metadata', METADATA_ABS)
    ] + dataChecksums

    def headerLines = [
        envSection('RUN'),
        envLine('Run name', workflow.runName),
        envLine('Run started', workflow.start),
        envLine('Nextflow', "${workflow.nextflow.version} (build ${workflow.nextflow.build})"),
        envLine('Session', workflow.sessionId),
        envLine('Command line', workflow.commandLine),
        envLine('Profile', workflow.profile),
        envLine('Run mode', run_value),
        envLine('Steps', ordered.join(', ')),
        envLine('Working path', WORKDIR_ABS),
        envLine('Original data', ORIGINALDATA_ABS),
        envLine('Metadata', METADATA_ABS),
        envLine('Run history', "run_history/${historyFolderName} (saved only if the run completes successfully)"),

        envSection('PIPELINE'),
        envLine('Git commit', gitCommit ?: 'unknown (not a git checkout)'),
        envLine('Git branch', gitBranch),
        envLine('Git remote', gitRemote),
        envLine('Uncommitted changes', gitChanges ? 'YES - see "Uncommitted pipeline changes" at the end of this record' : 'NONE'),
        envLine('Pipeline folder', gitDir),

        envSection('CONTAINER'),
        envLine('Image', imageName),
        envLine('Image run by ENV_LOG', '@@TASK_CONTAINER@@'),
        envLine('Image ID', imageId ?: 'unknown (image not found locally)'),
        envLine('Image digest', imageDigests ?: 'none (image built locally, not pulled from a registry)'),
        envLine('Created', imageCreated),
        envLine('Image platform', imagePlatform),

        envSection('COMPUTE'),
        envLine('Docker', dockerVersion),
        envLine('Docker system', dockerSystem),
        envLine('Host system', "${System.getProperty('os.name')} ${System.getProperty('os.version')} (${System.getProperty('os.arch')})"),
        envLine('Java', System.getProperty('java.version')),
        envLine('Pipeline CPUs', params.cpus),
        envLine('Pipeline start memory', params.memory),
        envLine('Pipeline max memory', params.max_memory ?: 'auto (Docker memory minus --memory_margin)'),
        envLine('Memory margin', params.memory_margin),
        envLine('Note', 'requested settings; each task can get more memory on retries (see nextflow_trace.txt)')
    ]
    def runHeader = headerLines.join('\n')

    // env_log.txt / env_log.pdf in the working folder always describe the
    // LATEST run. The old ones are removed first, so a stale env_log.pdf can
    // never be mistaken for (or copied into run_history/ as) this run's.
    // Successful runs keep their own copy in run_history/.
    ['env_log.txt', 'env_log.pdf'].each { n -> new File(WORKDIR_ABS, n).delete() }

    println "Docker image      : ${imageName} (${imageId ?: 'ID unknown'})"

    env_ch = ENV_LOG(
        Channel.of(
            tuple(
                WORKDIR_ABS,
                runHeader,
                runChanges
            )
        )
    )

    // INIT_PARAMS starts only after env_log.txt has been written.
    init_ch = INIT_PARAMS(
        env_ch.map { env_log ->
            tuple(
                WORKDIR_ABS,
                PIPELINE_LOG_ABS,
                file("${projectDir}/Initialization/Initialize.py"),
                file("${projectDir}/Initialization/dict_flags.json"),
                file("${projectDir}/Initialization/dict_parameters.json"),
                file("${WORKDIR_ABS}/Parameters.txt"),
                file(original_data.absolutePath),
                file(metaFile.absolutePath),
                CONTAINER_MOUNTS
            )
        }
    )

    run_input_ch = init_ch.map { paramdir, setup_log, init_done ->
        tuple(
            WORKDIR_ABS,
            steps_tsv,
            file(codePyDir.absolutePath),
            file(codeRDir.absolutePath),
            file(reportPy.absolutePath),
            file(qaScript.absolutePath),
            file(evalConfig.absolutePath),
            file(figscptn.absolutePath),
            file(seed_data.absolutePath),
            file(metaFile.absolutePath),
            file(seed_idmap.absolutePath),
            file(seed_mask.absolutePath),
            file(seed_maskdp.absolutePath),
            PIPELINE_LOG_ABS,
            file(reportSeed.absolutePath),
            file(methodDoc.absolutePath),
            paramdir,
            IMPUTATION_PROCEDURES_ABS ?: '',
            CONTAINER_MOUNTS
        )
    }

    RUN_ALL_MODULES(run_input_ch)

    workflow.onComplete = {
        // Always run cleanup, even when the pipeline failed.
        println ''
        if( !workflow.success ) {
            println '[CLEANUP] Running cleanup after pipeline failure.'
        }
        else {
            println '[CLEANUP] Starting post-run cleanup.'
        }

        // Run history: only for a successful run. A failed or interrupted run
        // writes nothing; the settings stay in the working folder as they are.
        if( workflow.success ) {
            def historyDir = writeRunHistory(
                WORKDIR_ABS,
                historyFolderName,
                capturedSettings,
                historyInfo + [
                    envLine('Run finished', workflow.complete),
                    envLine('Duration', workflow.duration),
                    envLine('Status', 'COMPLETED')
                ],
                [
                    new File(WORKDIR_ABS, 'env_log.txt'),
                    new File(WORKDIR_ABS, 'env_log.pdf'),
                    new File(WORKDIR_ABS, 'Parameter_Files'),
                    new File(WORKDIR_ABS, 'logfile_setup.txt'),
                    new File(WORKDIR_ABS, 'report.pdf'),
                    new File(params.outdir.toString(), 'nextflow_trace.txt')
                ]
            )
            if( historyDir ) {
                // Only the part of pipeline_log.txt that this run wrote.
                try {
                    // (No Java array syntax here: not allowed in strict Nextflow syntax.)
                    def logBytes = new File(PIPELINE_LOG_ABS).bytes
                    def startAt = (pipelineLogStart as long) <= logBytes.length ? (pipelineLogStart as int) : 0
                    new File(historyDir, 'pipeline_log.txt').bytes =
                        java.util.Arrays.copyOfRange(logBytes, startAt, logBytes.length)
                }
                catch( Exception e ) {
                    println "[RUN HISTORY WARNING] Could not save this run's pipeline log: ${e.message}"
                }
                println "[RUN HISTORY] Settings of this run saved to: ${historyDir}"
            }
        }
        else {
            println '[RUN HISTORY] Run failed: nothing saved to run_history/ (settings stay in the working folder).'
        }

        def cfgDataPath2 = resolveLaunch(config_data as String)
        def cfgData2 = parseKeyValueFile(cfgDataPath2)
        def workingPath2 = cfgData2.get('WORKING_PATH') ?: '.'
        def workdirAbs2 = resolveFrom(launchDir.toString(), workingPath2)

        // Nextflow's own folders are never touched by the file cleanup below.
        def nfWorkDir = workDir.toString()
        def nfMetaDir = new File(launchDir.toString(), '.nextflow').absolutePath

        println '[CLEANUP] Remove 0-byte files from output folder'
        runCommand([
            'find', workdirAbs2,
            '-type', 'f', '-size', '0',
            '-not', '-path', "${nfWorkDir}/*",
            '-not', '-path', "${nfMetaDir}/*",
            // Keep the run history complete, even an empty settings file.
            '-not', '-path', "${new File(workdirAbs2, 'run_history').absolutePath}/*",
            '-print', '-delete'
        ])

        println '[CLEANUP] Remove common hidden junk files'
        runCommand([
            'find', workdirAbs2,
            '-type', 'f',
            '(', '-name', '.DS_Store', '-o', '-name', '._*', '-o', '-name', 'Thumbs.db', ')',
            '-not', '-path', "${nfWorkDir}/*",
            '-not', '-path', "${nfMetaDir}/*",
            '-print', '-delete'
        ])

        if( cleanup_aggressive ) {
            println '[CLEANUP] Aggressive cleanup enabled'
            // rm -rf removes symlinks themselves and never follows them,
            // so linked input folders (e.g. modules_python) stay untouched.
            runCommand(['rm', '-rf', nfWorkDir])
            runCommand(['rm', '-rf', nfMetaDir, new File(launchDir.toString(), '.nextflow.log').absolutePath])
        }
        else {
            println '[CLEANUP] Aggressive cleanup disabled'
        }

        println '[CLEANUP] Done.'

        if( !workflow.success ) {
            println ''
            println '[PIPELINE FAILED]'
            println 'Cleanup completed, but the pipeline failed.'
            println ''
        }
    }


    workflow.onError = {
        println ''
        println '============================================================'
        println '[PIPELINE FAILED]'
        println 'The pipeline stopped before completion.'
        println 'Please check pipeline_log.txt and nextflow_trace.txt for details.'
        println 'Nextflow task diagnostics remain available in the work directory.'
        println '============================================================'
        println ''

        try {
            def failureLog = new File(PIPELINE_LOG_ABS)
            failureLog.parentFile?.mkdirs()

            failureLog << System.lineSeparator()
            failureLog << '============================================================' + System.lineSeparator()
            failureLog << '[PIPELINE STATUS] FAILED' + System.lineSeparator()
            failureLog << "[NEXTFLOW ERROR] ${workflow.errorMessage ?: 'Unknown error'}" + System.lineSeparator()
            failureLog << "[NEXTFLOW EXIT STATUS] ${workflow.exitStatus}" + System.lineSeparator()
            failureLog << '============================================================' + System.lineSeparator()

            println "[PIPELINE LOG] Failure written to ${PIPELINE_LOG_ABS}"
        }
        catch( Throwable logError ) {
            println "[PIPELINE LOG ERROR] Could not update pipeline_log.txt: ${logError.message}"
        }
    }

}
