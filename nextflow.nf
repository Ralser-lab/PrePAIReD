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

    // tag "INIT_PARAMS"
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


process MODULE {

    tag "${step_id}_${step_name} | attempt=${task.attempt} | mem=${task.memory} | cpus=${task.cpus}"
    maxForks 1

    input:
    tuple val(publish_root),
          val(step_id),
          val(step_name),
          val(engine),
          path(script_file),
          path(param_file),
          path(code_py_dir),
          path(code_r_dir),
          path(report_script),
          path(qa_script),
          path(eval_cfg),
          val(maxlfqso),
          path(meta_file),
          path(input_data),
          path(input_id),
          path(input_mask),
          path(input_mask_dp),
          path(input_report),
          path(doc_tsv)

    output:
    tuple val(step_id),
          val(step_name),
          path("data_${step_id}.tsv"),
          path("id_mapping.tsv"),
          path("data_${step_id}_mask.tsv"),
          path("data_${step_id}_mask_description.tsv"),
          path("report.pdf"),
          path("step_log.txt"),
          emit: step_out

    script:
    """
    set -euo pipefail

    mkdir -p report_folder

    echo "----- MODULE START (${step_id}_${step_name}) -----" > step_log.txt
    echo "[RESOURCE INFO] attempt=${task.attempt} cpus=${task.cpus} memory=${task.memory}" >> step_log.txt
    echo "[MODULE INFO] Running ${step_id}_${step_name}" >> step_log.txt

    cp -f "${maxlfqso}" ./maxlfq_docker.so || true

    if [ "${engine}" = "py" ]; then
        runner="python"
    else
        runner="Rscript"
    fi

    set +e
    "\$runner" "${script_file}" \\
      -d "${input_data}" \\
      -p "${param_file}" \\
      -o . \\
      -m "${meta_file}" \\
      -l step_log.txt \\
      -i "${input_id}" \\
      -c "${step_id}" >> step_log.txt 2>&1
    mod_rc=\$?
    set -e

    if [ \$mod_rc -ne 0 ]; then
        echo "###MOD_ERROR### module failed exit=\$mod_rc for ${step_id}_${step_name}" >> step_log.txt
        echo "[SKIP] QA skipped because module failed." >> step_log.txt
        touch MOD.fail
        touch QA.skip
        if [ \$mod_rc -eq 137 ] || [ \$mod_rc -eq 143 ]; then
            echo "[MEMORY ERROR] ${step_id}_${step_name} was killed with exit \$mod_rc." >> step_log.txt
            exit \$mod_rc
        fi

        exit 100
    else
        touch MOD.ok
    fi
        

    echo "----- MODULE END -----" >> step_log.txt

    echo "----- QA START (${step_id}_${step_name}) -----" >> step_log.txt
    set +e
    python "${qa_script}" \\
      -w . \\
      -m "${meta_file}" \\
      -s "${step_id}_${step_name}" \\
      -e "${eval_cfg}" \\
      --docker >> step_log.txt 2>&1
    qa_rc=\$?
    set -e

    if [ \$qa_rc -ne 0 ]; then
        echo "###QA_ERROR### QA failed exit=\$qa_rc for ${step_id}_${step_name}" >> step_log.txt
        echo "[WARN] QA failed but continuing" >> step_log.txt
        touch QA.fail
    else
        touch QA.ok
    fi
    
    echo "----- QA END -----" >> step_log.txt

    echo "----- REPORT START (${step_id}_${step_name}) -----" >> step_log.txt
    set +e
    if [ -f "${input_report}" ] && [ -s "${input_report}" ]; then
        python "${report_script}" \\
          --log step_log.txt \\
          --in_pdf "${input_report}" \\
          --out report_new.pdf \\
          --doc_tsv "${doc_tsv}" >> step_log.txt 2>&1
    else
        python "${report_script}" \\
          --log step_log.txt \\
          --out report_new.pdf \\
          --doc_tsv "${doc_tsv}" >> step_log.txt 2>&1
    fi
    rep_rc=\$?
    set -e

    if [ \$rep_rc -ne 0 ]; then
        echo "###REP_ERROR### report failed exit=\$rep_rc for ${step_id}_${step_name}" >> step_log.txt
        touch REP.fail
    else
        mv report_new.pdf report.pdf
        touch REP.ok
    fi

    
    echo "----- REPORT END -----" >> step_log.txt

    touch pipeline_log.txt
    {
        echo ""
        echo "===== ${step_id}_${step_name} (${engine}) ====="
        cat step_log.txt
    } >> pipeline_log.txt

    cp -f step_log.txt log.txt

    publish_step="${publish_root}/${step_id}_${step_name}"
    mkdir -p "\$publish_step"

    cp -f step_log.txt "\$publish_step/"
    cp -f log.txt "\$publish_step/" || true
    cp -f pipeline_log.txt "${publish_root}/" || true
    cp -f report.pdf "\$publish_step/" || true
    cp -f report.pdf "\$publish_step/" || true
    cp -f data_${step_id}.tsv "\$publish_step/" || true
    cp -f id_mapping.tsv "\$publish_step/" || true
    cp -f data_${step_id}_mask.tsv "\$publish_step/" || true
    cp -f data_${step_id}_mask_description.tsv "\$publish_step/" || true
    [ -f Normalization_factors.tsv ] && cp -f Normalization_factors.tsv "\$publish_step/" || true

    test -f data_${step_id}.tsv || {
        echo "[FATAL] Expected output data_${step_id}.tsv was not created." >> step_log.txt
        exit 100
    }

    test -f id_mapping.tsv || {
        echo "[FATAL] Expected output id_mapping.tsv was not created." >> step_log.txt
        exit 100
    }

    touch step.done
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
                file("${projectDir}/Initialization/Initialize.py"),
                file("${projectDir}/Initialization/dict_flags.pkl"),
                file("${projectDir}/Initialization/dict_parameters.pkl"),
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
            file(logSeed.absolutePath),
            file(reportSeed.absolutePath),
            file(methodDoc.absolutePath),
            paramdir
        )
    }

    step_tuples = []

    def prev_folder = start_dir
    def prev_counter = precounter

    stepsResolved.each { s ->

        def step_folder = new File(WORKDIR_ABS, "${s.counter}_${s.base}").absolutePath

        step_tuples << tuple(
            WORKDIR_ABS,
            s.counter,
            s.base,
            s.engine,
            file(new File(s.engine == 'py' ? codePyDir : codeRDir, s.scriptName).absolutePath),
            file(new File(WORKDIR_ABS, "Parameter_Files/${s.paramName}").absolutePath),
            file(codePyDir.absolutePath),
            file(codeRDir.absolutePath),
            file(reportPy.absolutePath),
            file(qaScript.absolutePath),
            file(evalConfig.absolutePath),
            maxlfq_so_param as String,
            file(metaFile.absolutePath),
            new File(prev_folder, "data_${prev_counter}.tsv").absolutePath,
            new File(prev_folder, "id_mapping.tsv").absolutePath,
            file(new File(prev_folder, "data_${prev_counter}_mask.tsv").absolutePath),
            file(new File(prev_folder, "data_${prev_counter}_mask_description.tsv").absolutePath),
            new File(prev_folder, "report.pdf").absolutePath,
            file(methodDoc.absolutePath)
        )

        prev_folder = step_folder
        prev_counter = s.counter
    }

    MODULE(Channel.fromList(step_tuples))

    onComplete:
    
    if( !workflow.success ) {
        println ''
        println '[CLEANUP] Skipped because pipeline failed.'
        println ''
        return
    }
    else {
        cfgDataPath2 = resolveLaunch(config_data as String)
        cfgData2 = parseKeyValueFile(cfgDataPath2)
        workingPath2 = cfgData2.get('WORKING_PATH') ?: '.'
        workdirAbs2 = resolveFrom(launchDir.toString(), workingPath2)
        
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

    onError:
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


