nextflow.enable.dsl=2
nextflow.preview.recursion=true

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

    set +e
    python "${init_script}" \\
      -w "." \\
      -m "${meta_file}" \\
      -d "${original_data}"
    py_exit=\$?
    set -e

    if [ "\$py_exit" -ne 0 ]; then
        echo "[FATAL] Initialize.py failed with exit code \$py_exit" >&2
        exit 100
    fi

    if [ ! -f logfile_setup.txt ]; then
        echo "[FATAL] logfile_setup.txt was not created by Initialize.py" >&2
        exit 100
    fi

    if grep -Fq '###ERROR###' logfile_setup.txt; then
        echo "[FATAL] Initialize.py produced ###ERROR### in logfile_setup.txt" >&2
        exit 100
    fi

    if [ ! -d Parameter_Files ]; then
        echo "[FATAL] Parameter_Files was not created by Initialize.py" >&2
        exit 100
    fi

    touch init.done
    """
}


/*
 * RUN_STEP — runs exactly ONE pipeline step, then emits the SAME tuple shape it
 * received (mandatory for .recurse()). The step to run is steps[idx]; idx is
 * advanced by 1 on output so the next recursion runs the next step.
 *
 * Sequential guarantee: .recurse() feeds this iteration's OUTPUT as the next
 * iteration's INPUT. The data file emitted here (out_data.tsv) literally becomes
 * the next step's in_data. Steps therefore run strictly in order, each consuming
 * the previous step's result. Nothing runs independently (maxForks 1, local).
 *
 * Display: each iteration is a separate task, so the executor shows
 *   [xx/aaaa] process > RUN_STEP (1)  ✔
 *   [xx/bbbb] process > RUN_STEP (2)  ✔
 *
 * Failure policy: module / QA / report failures are NON-FATAL. The step never
 * exits non-zero for those and ALWAYS emits a complete tuple; if the module
 * produces no new data, previous data/id/mask/report are CARRIED FORWARD. The
 * ONLY hard failure is an OOM kill (137/143), which cannot carry itself forward.
 */
process RUN_STEP {

    tag {
        def f = steps_tsv.split('\n')[(idx as int)].split('\t')
        "${f[0]}_${f[1]} | attempt=${task.attempt} | mem=${task.memory} | cpus=${task.cpus}"
    }

    maxForks 1

    beforeScript """
    echo "[RESOURCE INFO] Iteration ${(idx as int) + 1} | Attempt ${task.attempt} | CPUs=${task.cpus} | Memory=${task.memory}"
    """

    input:
    tuple val(publish_root),
          val(idx),
          val(steps_tsv),
          path(in_data, stageAs: 'in_data.tsv'),
          path(in_id, stageAs: 'in_id.tsv'),
          path(in_mask, stageAs: 'in_mask.tsv'),
          path(in_mask_dp, stageAs: 'in_mask_dp.tsv'),
          path(in_report, stageAs: 'in_report.pdf'),
          path(code_py_dir),
          path(code_r_dir),
          path(report_script),
          path(qa_script),
          path(eval_cfg),
          path(meta_file),
          path(doc_tsv),
          path(paramdir),
          val(maxlfqso)

    output:
    tuple val(publish_root),
          val(idx + 1),
          val(steps_tsv),
          path("out_data.tsv"),
          path("out_id_mapping.tsv"),
          path("out_mask.tsv"),
          path("out_mask_description.tsv"),
          path("out_report.pdf"),
          path(code_py_dir),
          path(code_r_dir),
          path(report_script),
          path(qa_script),
          path(eval_cfg),
          path(meta_file),
          path(doc_tsv),
          path(paramdir),
          val(maxlfqso)

    script:
    """
    set -uo pipefail

    BASE="\$PWD"
    IDX=${idx}

    cat > steps.tsv <<'EOF_STEPS'
${steps_tsv}
EOF_STEPS

    line=\$(sed -n "\$((IDX + 1))p" steps.tsv)
    if [ -z "\$line" ]; then
        echo "[FATAL] No step at index \$IDX" >&2
        exit 100
    fi

    counter=\$(printf '%s' "\$line" | cut -f1)
    base=\$(printf '%s' "\$line" | cut -f2)
    engine=\$(printf '%s' "\$line" | cut -f3)
    script_name=\$(printf '%s' "\$line" | cut -f4)
    param_name=\$(printf '%s' "\$line" | cut -f5)

    step_name="\${base}"

    echo "[MODULE INFO] Running \${counter}_\${step_name}"

    mkdir -p report_folder
    : > step_log.txt
    echo "----- MODULE START (\${counter}_\${step_name}) -----" >> step_log.txt

    if cp -f "${maxlfqso}" ./maxlfq_docker.so 2>>step_log.txt; then
        :
    else
        echo "[WARN] could not copy maxlfq_docker.so from ${maxlfqso}" >> step_log.txt
    fi

    if [ "\$engine" = "py" ]; then
        script_path="\$BASE/${code_py_dir}/\${script_name}"
        runner="python"
    else
        script_path="\$BASE/${code_r_dir}/\${script_name}"
        runner="Rscript"
    fi

    param_file="\$BASE/${paramdir}/\${param_name}"

    rm -f ./*_Parameters.json
    if [ -f "\$param_file" ]; then
        cp -f "\$param_file" .
    else
        echo "[WARN] Missing parameter file: \$param_file" >> step_log.txt
    fi

    mod_rc=0
    if [ ! -f "\$script_path" ]; then
        echo "[WARN] Missing module script: \$script_path" >> step_log.txt
        echo "###MOD_ERROR### module script missing for \${counter}_\${step_name}" >> step_log.txt
        mod_rc=127
    else
        set +e
        "\$runner" "\$script_path" \\
          -d "\$BASE/in_data.tsv" \\
          -p "\$param_file" \\
          -o . \\
          -m "\$BASE/${meta_file}" \\
          -l step_log.txt \\
          -i "\$BASE/in_id.tsv" \\
          -c "\$counter" >> step_log.txt 2>&1
        mod_rc=\$?
        set -e
    fi

    if [ \$mod_rc -ne 0 ]; then
        echo "###MOD_ERROR### module failed exit=\$mod_rc for \${counter}_\${step_name}" >> step_log.txt
        touch MOD.fail

        if [ \$mod_rc -eq 137 ] || [ \$mod_rc -eq 143 ]; then
            echo "[MEMORY ERROR] \${counter}_\${step_name} killed with exit \$mod_rc (out of memory)." >> step_log.txt
            echo "[MEMORY ERROR] Nextflow will retry with more memory if retries remain." >> step_log.txt
            exit \$mod_rc
        fi

        echo "[INFO] module failed but continuing (non-fatal); previous data carried forward." >> step_log.txt
    else
        touch MOD.ok
    fi

    echo "----- MODULE END -----" >> step_log.txt

    echo "" >> step_log.txt
    echo "----- QA START (\${counter}_\${step_name}) -----" >> step_log.txt

    if [ \$mod_rc -eq 0 ]; then
        set +e
        python "\$BASE/${qa_script}" \\
          -w . \\
          -m "\$BASE/${meta_file}" \\
          -s "\${counter}_\${step_name}" \\
          -e "\$BASE/${eval_cfg}" \\
          --docker >> step_log.txt 2>&1
        qa_rc=\$?
        set -e

        if [ \$qa_rc -ne 0 ]; then
            echo "###QA_ERROR### QA failed exit=\$qa_rc for \${counter}_\${step_name}" >> step_log.txt
            echo "[WARN] QA failed but continuing" >> step_log.txt
            touch QA.fail
        else
            touch QA.ok
        fi
    else
        echo "[SKIP] QA skipped because module failed." >> step_log.txt
        touch QA.skip
    fi

    echo "----- QA END -----" >> step_log.txt

    echo "" >> step_log.txt
    echo "----- REPORT START (\${counter}_\${step_name}) -----" >> step_log.txt

    set +e
    if [ -f "\$BASE/in_report.pdf" ] && [ -s "\$BASE/in_report.pdf" ]; then
        python "\$BASE/${report_script}" \\
          --log step_log.txt \\
          --in_pdf "\$BASE/in_report.pdf" \\
          --out report_new.pdf \\
          --doc_tsv "\$BASE/${doc_tsv}" >> step_log.txt 2>&1
    else
        python "\$BASE/${report_script}" \\
          --log step_log.txt \\
          --out report_new.pdf \\
          --doc_tsv "\$BASE/${doc_tsv}" >> step_log.txt 2>&1
    fi
    rep_rc=\$?
    set -e

    if [ \$rep_rc -ne 0 ]; then
        echo "###REP_ERROR### report failed exit=\$rep_rc for \${counter}_\${step_name}" >> step_log.txt
        echo "[WARN] report failed but continuing" >> step_log.txt
        touch REP.fail
    else
        mv -f report_new.pdf report.pdf
        touch REP.ok
    fi

    echo "----- REPORT END -----" >> step_log.txt

    cp -f step_log.txt log.txt || true

    if [ -f "data_\${counter}.tsv" ]; then
        cp -f "data_\${counter}.tsv" out_data.tsv
    else
        echo "[INFO] carrying forward previous data (no data_\${counter}.tsv produced)" >> step_log.txt
        cp -f "\$BASE/in_data.tsv" out_data.tsv
    fi

    if [ -f "id_mapping.tsv" ]; then
        cp -f "id_mapping.tsv" out_id_mapping.tsv
    else
        cp -f "\$BASE/in_id.tsv" out_id_mapping.tsv
    fi

    if [ -f "data_\${counter}_mask.tsv" ]; then
        cp -f "data_\${counter}_mask.tsv" out_mask.tsv
    else
        cp -f "\$BASE/in_mask.tsv" out_mask.tsv
    fi

    if [ -f "data_\${counter}_mask_description.tsv" ]; then
        cp -f "data_\${counter}_mask_description.tsv" out_mask_description.tsv
    else
        cp -f "\$BASE/in_mask_dp.tsv" out_mask_description.tsv
    fi

    if [ -f "report.pdf" ]; then
        cp -f "report.pdf" out_report.pdf
    elif [ -f "\$BASE/in_report.pdf" ]; then
        cp -f "\$BASE/in_report.pdf" out_report.pdf
    else
        : > out_report.pdf
    fi

    publish_step="${publish_root}/\${counter}_\${step_name}"
    mkdir -p "\$publish_step"
    cp -f step_log.txt "\$publish_step/" || true
    cp -f log.txt "\$publish_step/" || true
    [ -f report.pdf ] && cp -f report.pdf "\$publish_step/" || true
    cp -rf report_folder "\$publish_step/" || true
    [ -f "data_\${counter}.tsv" ] && cp -f "data_\${counter}.tsv" "\$publish_step/" || true
    [ -f id_mapping.tsv ] && cp -f id_mapping.tsv "\$publish_step/" || true
    [ -f "data_\${counter}_mask.tsv" ] && cp -f "data_\${counter}_mask.tsv" "\$publish_step/" || true
    [ -f "data_\${counter}_mask_description.tsv" ] && cp -f "data_\${counter}_mask_description.tsv" "\$publish_step/" || true
    [ -f Normalization_factors.tsv ] && cp -f Normalization_factors.tsv "\$publish_step/" || true

    echo "[MODULE INFO] Finished \${counter}_\${step_name}"
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

    CONFIG_TXT_ABS = new File(WORKDIR_ABS, 'Config.txt').absolutePath

    ORIGINALDATA_ABS = cfgData.get('ORIGINALDATA') ? resolveFrom(WORKDIR_ABS, cfgData['ORIGINALDATA']) : null
    METADATA_ABS = cfgData.get('METADATA') ? resolveFrom(WORKDIR_ABS, cfgData['METADATA']) : null

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

    logSeed = ensureExistsOrCreateEmpty(new File(WORKDIR_ABS, 'pipeline_log.txt').absolutePath, 'pipeline_log.txt')
    reportSeed = ensureExistsOrCreateEmpty(new File(start_dir, 'report.pdf').absolutePath, 'report seed')

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

    def steps_tsv = stepsResolved.collect {
        "${it.counter}\t${it.base}\t${it.engine}\t${it.scriptName}\t${it.paramName}"
    }.join('\n')

    def n_steps = stepsResolved.size()

    println "Working path      : ${WORKDIR_ABS}"
    println "Original data     : ${ORIGINALDATA_ABS}"
    println "Metadata          : ${METADATA_ABS}"
    println "Steps             : ${n_steps}"

    if( n_steps == 0 ) {
        die("No steps to run (Order = is empty)")
    }

    init_ch = INIT_PARAMS(
        Channel.of(
            tuple(
                WORKDIR_ABS,
                file("${projectDir}/Initialization/Initialize.py"),
                file("${projectDir}/Initialization/dict_flags.pkl"),
                file("${projectDir}/Initialization/dict_parameters.pkl"),
                file("${WORKDIR_ABS}/Parameters.txt"),
                file(original_data.absolutePath),
                file(metaFile.absolutePath)
            )
        )
    )

    // Seed tuple (iteration 0 input). MUST match RUN_STEP input shape exactly.
    // .first() converts the queue channel into a value channel, required by
    // .recurse(). Mapping INIT_PARAMS output also gates start on init completion.
    seed_ch = init_ch
        .map { paramdir, setup_log, init_done ->
            tuple(
                WORKDIR_ABS,
                0,
                steps_tsv,
                file(seed_data.absolutePath),
                file(seed_idmap.absolutePath),
                file(seed_mask.absolutePath),
                file(seed_maskdp.absolutePath),
                file(reportSeed.absolutePath),
                file(codePyDir.absolutePath),
                file(codeRDir.absolutePath),
                file(reportPy.absolutePath),
                file(qaScript.absolutePath),
                file(evalConfig.absolutePath),
                file(metaFile.absolutePath),
                file(methodDoc.absolutePath),
                paramdir,
                maxlfq_so_param as String
            )
        }
        .first()

    // Run all steps sequentially. .times(n_steps) => one iteration per step;
    // each iteration's OUTPUT feeds the next iteration's INPUT (strict order).
    RUN_STEP
        .recurse(seed_ch)
        .times(n_steps)


    workflow.onComplete {
        if( !workflow.success ) {
            println ''
            println '[CLEANUP] Skipped because pipeline failed.'
            println ''
        }
        else {
            try {
                def out = new File(WORKDIR_ABS, 'pipeline_log.txt')
                out.text = ''
                stepsResolved.each { st ->
                    def slog = new File("${WORKDIR_ABS}/${st.counter}_${st.base}", 'step_log.txt')
                    out << "\n===== ${st.counter}_${st.base} (${st.engine}) =====\n"
                    out << (slog.exists() ? slog.text : "[WARN] step_log.txt not found for ${st.counter}_${st.base}\n")
                }
                println "[LOG] Rebuilt aggregated log: ${out.absolutePath}"
            } catch( Throwable t ) {
                println "[WARN] Could not rebuild pipeline_log.txt: ${t.message}"
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
        }
    }

    workflow.onError {
        println ''
        println '[PIPELINE FAILED]'
        println 'The pipeline stopped before completion.'
        println 'If the failed exit status was 137 or 143, the likely reason is insufficient memory.'
        println 'Nextflow already retried with larger memory if retries were available.'
        println 'Please rerun with more memory, for example:'
        println '  nextflow run nextflow.nf --memory "32 GB" --cpus 8'
        println 'For initialization memory:'
        println '  nextflow run nextflow.nf --init_memory "16 GB"'
        println 'If using Docker Desktop, increase Docker Desktop memory in Settings > Resources.'
        println ''
    }
}

