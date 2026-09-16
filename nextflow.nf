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

String bashQuote(String s) {
    return s.replace("'", "'\\''")
}

process INIT_PARAMS {

    tag "INIT_PARAMS | attempt=${task.attempt} | mem=${task.memory} | cpus=${task.cpus}"

    maxForks 1

    publishDir { publish_root }, mode: 'copy', overwrite: true

    input:
    tuple val(publish_root),
          val(pipeline_log_path),
          path(init_script),
          path(flags_pkl),
          path(params_pkl),
          path(params_txt),
          path(original_data),
          path(meta_file)

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

    python "${init_script}" \
        -w "." \
        -m "${meta_file}" \
        -d "${original_data}" 2>&1 |
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
          path(maxlfqdll),
          path(maxlfqdylib),
          val(maxlfqso),
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
          val(imputation_procedures_input)

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

    current_data="\$BASE/${seed_data}"
    current_id="\$BASE/${id_file}"
    current_mask="\$BASE/${mask_file}"
    current_mask_dp="\$BASE/${mask_dp_file}"
    current_report="\$BASE/${prev_report}"


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

        if [ -f Normalization_factors.tsv ]; then
            cp -f Normalization_factors.tsv "\$publish_step/" || true
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

        # Prepare the MaxLFQ shared library.
        set +e
        cp -f "${maxlfqso}" ./maxlfq_docker.so
        maxlfq_cp_rc=\$?
        set -e

        if [ "\$maxlfq_cp_rc" -ne 0 ] ||
           [ ! -f ./maxlfq_docker.so ]; then

            log_both "###MOD_ERROR### maxlfq_docker.so could not be prepared: ${maxlfqso}"
            log_both "[MODULE STATUS] FAILED"
            log_both "[QA STATUS] SKIPPED because the module could not start."
            log_both "----- MODULE END -----"

            touch MOD.fail
            touch QA.skip

            printf '%s\\t%s\\t%s\\n' \
                "\$module_name" "SKIPPED" "NA" >> qa_status.tsv

            cp -f step_log.txt log.txt || true
            publish_current_step "\$publish_step" "\$counter"
            printf '%s\\n' "[FATAL] \${module_name}: maxlfq_docker.so is missing" >> "\$persistent_log"
            exit 100
        fi

        # Resolve module runner and script.
        if [ "\$engine" = "py" ]; then
            script_path="\$BASE/${code_py_dir}/\${script_name}"
            runner="python"

        elif [ "\$engine" = "r" ] ||
             [ "\$engine" = "R" ]; then

            script_path="\$BASE/${code_r_dir}/\${script_name}"
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

        param_file="\$BASE/${paramdir}/\${param_name}"

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
            -m "\$BASE/${meta_file}" \
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

        python "\$BASE/${qa_script}" \
            -w . \
            -m "\$BASE/${meta_file}" \
            -s "\${module_name}" \
            -e "\$BASE/${eval_cfg}" \
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

            python "\$BASE/${report_script}" \
                --log step_log.txt \
                --in_pdf "\$current_report" \
                --out report_new.pdf \
                --doc_tsv "\$BASE/${doc_tsv}" 2>&1 |
                tee -a \
                    step_log.txt \
                    "\$persistent_log"

            rep_rc=\${PIPESTATUS[0]}

        else
            python "\$BASE/${report_script}" \
                --log step_log.txt \
                --out report_new.pdf \
                --doc_tsv "\$BASE/${doc_tsv}" 2>&1 |
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
        current_norm=""

        if [ -f "\$BASE/data_\${counter}_mask.tsv" ]; then
            current_mask="\$BASE/data_\${counter}_mask.tsv"
        fi

        if [ -f "\$BASE/data_\${counter}_mask_description.tsv" ]; then
            current_mask_dp="\$BASE/data_\${counter}_mask_description.tsv"
        fi

        if [ -f "\$BASE/Normalization_factors.tsv" ]; then
            current_norm="\$BASE/Normalization_factors.tsv"
        fi

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
    def maxlfq_dll_param = pget('maxlfq_dll', 'modules_python/maxlfq.dll')
    def maxlfq_dylib_param = pget('maxlfq_dylib', 'modules_python/maxlfq.dylib')
    def maxlfq_so_param = pget('maxlfq_so', '/workflow/modules_python/maxlfq_docker.so')
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

    run_value = lines.find {
        it.trim().toLowerCase().startsWith('run')
    }?.split('=', 2)?.getAt(1)?.trim()

    if( !run_value ) {
        die("Config file missing 'Run =' entry")
    }

    order_line = lines.find {
        it.trim().toLowerCase().startsWith('order')
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
    maxlfqdll = mustExist(resolveCode(maxlfq_dll_param as String), 'maxlfq dll')
    maxlfqdylib = mustExist(resolveCode(maxlfq_dylib_param as String), 'maxlfq dylib')

    figscptn = mustExist(resolveCode(figs_cptn_param as String), 'Figures Captions')
    codePyDir = mustExist(resolveCode(code_py as String), 'modules_python')
    codeRDir = mustExist(resolveCode(code_r as String), 'modules_R')


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

    init_ch = INIT_PARAMS(
        Channel.of(
            tuple(
                WORKDIR_ABS,
                PIPELINE_LOG_ABS,
                file("${projectDir}/Initialization/Initialize.py"),
                file("${projectDir}/Initialization/dict_flags.json"),
                file("${projectDir}/Initialization/dict_parameters.json"),
                file("${WORKDIR_ABS}/Parameters.txt"),
                file(original_data.absolutePath),
                file(metaFile.absolutePath)
            )
        )
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
            file(maxlfqdll.absolutePath),
            file(maxlfqdylib.absolutePath),
            maxlfq_so_param as String,
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
            IMPUTATION_PROCEDURES_ABS ?: ''
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

        def cfgDataPath2 = resolveLaunch(config_data as String)
        def cfgData2 = parseKeyValueFile(cfgDataPath2)
        def workingPath2 = cfgData2.get('WORKING_PATH') ?: '.'
        def workdirAbs2 = resolveFrom(launchDir.toString(), workingPath2)

        println '[CLEANUP] Remove 0-byte files from output folder'
        "bash -lc \"find '${bashQuote(workdirAbs2)}' -type f -size 0 -print -delete || true\"".execute().waitFor()

        println '[CLEANUP] Remove common hidden junk files'
        "bash -lc \"find '${bashQuote(workdirAbs2)}' -type f \\( -name '.DS_Store' -o -name '._*' -o -name 'Thumbs.db' \\) -print -delete || true\"".execute().waitFor()

        if( cleanup_aggressive ) {
            println '[CLEANUP] Aggressive cleanup enabled'
            "bash -lc \"rm -rf '${bashQuote(workDir.toString())}' || true\"".execute().waitFor()
            "bash -lc \"rm -rf '${bashQuote(launchDir.toString())}/.nextflow' '${bashQuote(launchDir.toString())}/.nextflow.log' || true\"".execute().waitFor()
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
