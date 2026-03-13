nextflow.enable.dsl=2

def pget = { k, d -> params.keySet().contains(k) ? params[k] : d }
params.config_data    = pget('config_data','Config_data.txt')
// params.config_file    = pget('config_file','Config.txt')
// params.param_dir      = pget('param_dir','Parameter_Files')
// params.initial_dir    = pget('initial_dir','Data')
// params.meta_path      = pget('meta_path', params.initial_dir)
params.code_py        = pget('code_py','modules_python')
params.code_r         = pget('code_r','modules_R')
params.method_doc     = pget('method_doc','docs/Method_description.tsv')
params.keep_generated = pget('keep_generated', false)
params.qa_script   = pget('qa_script', 'modules_python/module_quality_assessment.py')

params.cleanup_aggressive = pget('cleanup_aggressive', false)

params.eval_config = pget('eval_config', 'Config_evaluation.txt')

params.maxlfq_dll = pget('maxlfq_dll', 'modules_python/maxlfq.dll')
params.maxlfq_dylib = pget('maxlfq_dylilb', 'modules_python/maxlfq.dylib')
params.maxlfq_so = pget('maxlfq_so', '/workflow/modules_python/maxlfq_docker.so')

params.figs_cptn = pget('figs_cptn', 'modules_python/Figures_Captions.tsv')


def PARAM_DIR_NAME = "Parameter_Files"


def die = { msg -> throw new RuntimeException(msg as String) }

def dieSoft = { String msg ->
  log.error(msg)     // prints nicely in Nextflow logs
  System.err.println("\n[FATAL] " + msg + "\n")
  System.exit(1)     // clean stop before running tasks
}
// --- path helpers

def parseKeyValueFile = { String path ->
  def m = [:]
  def f = new File(path)
  if( !f.exists() ) die("Missing config_data file: ${path}")
  f.eachLine { line ->
    def s = line.trim()
    if( !s || s.startsWith('#') ) return
    def parts = s.split('=',2)
    if( parts.size() == 2 )
      m[parts[0].trim()] = parts[1].trim().replaceAll(/^["']|["']$/, '')
  }
  return m
}

def resolveLaunch = { String p ->
  def f = new File(p)
  f.isAbsolute() ? f.getAbsolutePath() : new File(launchDir.toString(), p).getAbsolutePath()
}
def resolveFrom = { String base, String p ->
  def f = new File(p)
  f.isAbsolute() ? f.getAbsolutePath() : new File(base, p).getAbsolutePath()
}
def resolveCode = { String p ->
  def f = new File(p)
  f.isAbsolute() ? f.getAbsolutePath() : new File(projectDir.toString(), p).getAbsolutePath()
}

def cfgDataPath = resolveLaunch(params.config_data as String)
def cfgData     = parseKeyValueFile(cfgDataPath)

def WORKING_PATH = cfgData.get('WORKING_PATH') ?: '.'
def WORKDIR_ABS  = resolveFrom(launchDir.toString(), WORKING_PATH)


// everything instance-relative goes here
def resolveInstance = { String p ->
  def f = new File(p)
  f.isAbsolute() ? f.getAbsolutePath() : new File(WORKDIR_ABS, p).getAbsolutePath()
}

// --- file helpers
def mustExist = { String path, String label ->
  def f = new File(path)
  if( !f.exists() ) die("Missing ${label}: ${path}")
  return f
}

def ensureExistsOrCreateEmpty = { String path, String label ->
  def f = new File(path)
  if( !f.exists() ) {
    f.parentFile?.mkdirs()
    f.text = ""
    println "Created empty ${label}: ${path}"
  }
  return f
}

// Always inside working dir:
// def PARAM_DIR_ABS   = new File(WORKDIR_ABS, "Parameter_Files").getAbsolutePath()
def CONFIG_TXT_ABS  = new File(WORKDIR_ABS, "Config.txt").getAbsolutePath()

// Data/metadata can be absolute or relative-to-workingdir
def ORIGINALDATA_ABS = cfgData.get('ORIGINALDATA') ? resolveFrom(WORKDIR_ABS, cfgData['ORIGINALDATA']) : null
def METADATA_ABS     = cfgData.get('METADATA')     ? resolveFrom(WORKDIR_ABS, cfgData['METADATA'])     : null

params.config_file = CONFIG_TXT_ABS
// params.param_dir   = PARAM_DIR_ABS
params.initial_dir = new File(WORKDIR_ABS, "Data").getAbsolutePath()
params.meta_path   = params.initial_dir


// -------------------- parse config --------------------
def cfgFile = mustExist(params.config_file as String, "Config.txt")

def lines = cfgFile.readLines()
def run_value = lines.find { it.trim().toLowerCase().startsWith('run') }?.split('=',2)?.getAt(1)?.trim()
if( !run_value ) die("Config file missing 'Run =' entry")

def order_line = lines.find { it.trim().toLowerCase().startsWith('order') }
if( !order_line ) die("Config file missing 'Order =' entry")

def ordered = order_line.split('=',2)[1].replaceAll(/\s+/,'').split(',').findAll{it}.toList()

println "Run mode          : ${run_value}"
println "Ordered processes : ${ordered}"

// -------------------- rerun handling --------------------
def start_dir
def precounter
if( run_value.toLowerCase() == 'new' ) {
  start_dir  = WORKDIR_ABS
  precounter = "00"
} else {
  // remove trailing _py or _r from Run value for folder lookup
  def rerun_name = run_value.replaceFirst(/_(?i)(py|r)$/, '')
  
  def rerun_folder = new File(WORKDIR_ABS, rerun_name)
  if( !rerun_folder.isDirectory() )
      die("Rerun folder '${rerun_name}' not found (resolved to ${rerun_folder})")
  start_dir  = rerun_folder.getAbsolutePath()
  precounter = rerun_name.split('_',2)[0]
  
  // --- PRE-FLIGHT CHECK: block rerun if the chosen step contains ###ERROR###
  def stepLog = new File(rerun_folder, "step_log.txt")
  if( stepLog.exists() ) {
    def bad = stepLog.readLines().any { it.contains('###ERROR###') }
    if( bad ) {
      def hits = stepLog.readLines()
                        .withIndex()
                        .findAll { line, idx -> line.contains('###ERROR###') }
                        .takeRight(20)
                        .collect { line, idx -> "${idx+1}: ${line}" }
                        .join('\n')

      dieSoft(
        "Cannot rerun from '${rerun_name}' because it contains error markers in step_log.txt.\n" +
        "Fix that step (or choose an earlier step), then rerun.\n\n" +
        "---- ###ERROR### lines (last 20) ----\n${hits}\n" +
        "-------------------------------------"
      )
    }
  } else {
    // optional: if missing step_log should also block rerun
    // dieSoft("Cannot rerun from '${rerun_name}' because step_log.txt is missing in ${rerun_folder}")
    println "[WARN] step_log.txt not found in ${rerun_folder} (skipping pre-flight error check)"
  }
}

println "Start dir         : ${start_dir}"
println "Precounter        : ${precounter}"

// -------------------- validation --------------------
// def seed_data   = mustExist(new File(start_dir, "data_${precounter}.tsv").getAbsolutePath(), "seed data")
def seed_data_path = run_value.equalsIgnoreCase('new')
  ? ORIGINALDATA_ABS
  : new File(start_dir, "data_${precounter}.tsv").getAbsolutePath()

def seed_data = mustExist(seed_data_path, "seed data")

// always original
def original_data = mustExist(ORIGINALDATA_ABS, "ORIGINALDATA")

// def seed_mask   = mustExist(new File(start_dir, "data_${precounter}_mask.tsv").getAbsolutePath(), "seed mask")
// def seed_maskdp = mustExist(new File(start_dir, "data_${precounter}_mask_description.tsv").getAbsolutePath(), "seed mask description")

def seed_idmap_path  = new File(start_dir, "id_mapping.tsv").getAbsolutePath()

def seed_mask_path   = new File(start_dir, "data_${precounter}_mask.tsv").getAbsolutePath()
def seed_maskdp_path = new File(start_dir, "data_${precounter}_mask_description.tsv").getAbsolutePath()

def seed_mask   = ensureExistsOrCreateEmpty(seed_mask_path,   "seed mask (auto-created)")
def seed_maskdp = ensureExistsOrCreateEmpty(seed_maskdp_path, "seed mask description (auto-created)")
// def seed_idmap = ensureExistsOrCreateEmpty(seed_idmap_path, "seed idmap (auto-created)")
def seed_idmap = (run_value.toLowerCase() == 'new')
  ? ensureExistsOrCreateEmpty(seed_idmap_path, "seed idmap (auto-created)")
  : mustExist(seed_idmap_path, "seed idmap (id_mapping.tsv)")
  
  
def logSeed = ensureExistsOrCreateEmpty(resolveInstance("pipeline_log.txt"), "pipeline_log.txt")

def reportSeed = ensureExistsOrCreateEmpty(new File(start_dir, "report.pdf").getAbsolutePath(), "report seed")
println "Using report seed : ${reportSeed}"

// def metaFile  = mustExist(resolveInstance("${params.meta_path}/Metadata.tsv"), "Metadata.tsv")

def meta_path = (run_value.equalsIgnoreCase('new') && METADATA_ABS)
  ? METADATA_ABS
  : resolveInstance("${params.meta_path}/Metadata.tsv")

// def metaFile = mustExist(meta_path, "Metadata.tsv")

def metaFile = mustExist(
    run_value.equalsIgnoreCase('new') ? METADATA_ABS : METADATA_ABS,
    "Metadata.tsv"
)

def methodDoc = mustExist(resolveCode(params.method_doc as String), "method doc TSV")
def reportPy  = mustExist(resolveCode("${params.code_py}/report.py"), "modules_python/report.py")

def qaScript   = mustExist(resolveCode(params.qa_script as String), "QA script")
def evalConfig = mustExist(resolveInstance(params.eval_config as String), "Config_evaluation.txt")

def maxlfqdll    = mustExist(resolveCode(params.maxlfq_dll as String), "maxlfq dll")
def maxlfqdylib  = mustExist(resolveCode(params.maxlfq_dylib as String), "maxlfq dylib")
def maxlfqso     = params.maxlfq_so as String

def figscptn = mustExist(resolveCode(params.figs_cptn as String), "Figures Captions")

// validate steps
def stepsResolved = []
ordered.each { item ->
  def (counter, rest) = item.split('_',2)
  def suffix = rest.tokenize('_')[-1].toLowerCase()
  def engine = (suffix == 'r') ? 'r' : 'py'
  def base   = rest.replaceAll(/_(?i)(py|r)$/,'')

  def codeDirRel = (engine == 'py') ? (params.code_py as String) : (params.code_r as String)
  def ext        = (engine == 'py') ? 'py' : 'R'

  def scriptPath = resolveCode("${codeDirRel}/module_${base}.${ext}")
  // def paramPath  = resolveInstance("${params.param_dir}/${counter}_${base}_parameters.json")
  
  def paramRel = "${counter}_${base}_parameters.json"
  
  mustExist(scriptPath, "step script for ${item}")
  // mustExist(paramPath,  "param JSON for ${item}")

  stepsResolved << [counter:counter, base:base, engine:engine, script:scriptPath, paramRel:paramRel] //, param:paramPath]
}

println "Working path      : ${WORKDIR_ABS}"
println "Original data     : ${ORIGINALDATA_ABS}"
println "Metadata          : ${METADATA_ABS}"
def param_dir = new File(WORKDIR_ABS, "Parameter_Files").getAbsolutePath()

// -------------------- temp generation dir --------------------
def genDir = new File("${workDir}/_generated_${System.currentTimeMillis()}")
genDir.mkdirs()

workflow.onComplete {
  if( params.keep_generated ) {
    println "Generated files kept at: ${genDir}"
  } else {
    try { genDir.deleteDir() } catch(e) { }
  }
}

def modulePath = new File(genDir, "RUN_MODULE.nf").getAbsolutePath()
def wfPath     = new File(genDir, "workflow_generated.nf").getAbsolutePath()

// -------------------- module file (container comes from the config/profile) --------------------
new File(modulePath).text = '''\
nextflow.enable.dsl=2

process INIT_PARAMS {

  tag "INIT_PARAMS"
  maxForks 1

  publishDir { "${publish_root}" }, mode: 'copy'

  input:
    tuple val(publish_root),
          path(init_script),
          path(flags_pkl),
          path(params_pkl),
          path(params_txt),
          path(original_data),
          path(meta_file)

  output:
    tuple path("Parameter_Files"), path("init.done"), path("logfile_setup.txt")

  script:
  """
  set -euo pipefail

  python "${init_script}" \\
    -w "." \\
    -m "${meta_file}" \\
    -d "${original_data}"

  test -d Parameter_Files
  touch init.done
  test -f logfile_setup.txt

  if grep -q '###ERROR###' logfile_setup.txt; then
      mkdir -p "${publish_root}"
      cp logfile_setup.txt "${publish_root}/"
      cp -r Parameter_Files "${publish_root}/"
      echo "Found ###ERROR### in logfile_setup.txt"
      exit 100
  fi
  """
}


process RUN_MODULE {

  tag { "${counter}_${step_name}" }
  maxForks 1

  publishDir { "${publish_root}/${counter}_${step_name}" }, mode:'copy', saveAs: { fname ->
      if( fname in ['log.txt','pipeline_log.txt'] ) return null
      return fname
    }

  input:
    tuple val(publish_root),
          val(engine),
          val(counter),
          val(step_name),
          val(param_dir),
          path(paramdir),
          path(script),
          path(report_script),
          path(qa_script),
          path(eval_cfg),
          path(maxlfqdll),
          path(maxlfqdylib),
          val(maxlfqso),
          path(figscptn),
          path(data_file),
          path(param_file),
          path(meta_file),
          path(id_file),
          path(mask_file),
          path(mask_dp_file),
          path(log_file),
          path(prev_report),
          path(doc_tsv)

  output:
    tuple path("data_${counter}.tsv"),
          path("id_mapping.tsv"),
          path("data_${counter}_mask.tsv", optional: true),
          path("data_${counter}_mask_description.tsv", optional: true),
          path(log_file),
          path("step_log.txt"),
          path("report_folder"),
          path("report.pdf"),
          val(param_dir)

  script:
  """
  set -euo pipefail
  mkdir -p report_folder
  : > step_log.txt

  echo "----- MODULE START (${counter}_${step_name}) -----" >> step_log.txt
  
  cp -f /workflow/modules_python/maxlfq_docker.so ./maxlfq_docker.so
  test -f ./maxlfq_docker.so || { echo "maxlfq_docker.so not found in container"; exit 1; }

  set +e
  if [ "${engine}" = "py" ]; then
    python "${script}" \\
      -d "${data_file}" \\
      -p "${param_file}" \\
      -o . \\
      -m "${meta_file}" \\
      -l step_log.txt \\
      -i "${id_file}" \\
      -c "${counter}" >> step_log.txt 2>&1
  else
    Rscript "${script}" \\
      -d "${data_file}" \\
      -p "${param_file}" \\
      -o . \\
      -m "${meta_file}" \\
      -l step_log.txt \\
      -i "${id_file}" \\
      -c "${counter}" >> step_log.txt 2>&1
  fi
  mod_rc=\\$?
  set -e

  if [ \\$mod_rc -ne 0 ]; then
    echo "###MOD_ERROR### module failed (exit=\\$mod_rc) for ${counter}_${step_name}" >> step_log.txt
    touch MOD.fail
  else
    touch MOD.ok
  fi

  echo "----- MODULE END -----" >> step_log.txt


  echo "" >> step_log.txt
  echo "----- QA START (${counter}_${step_name}) -----" >> step_log.txt

  set +e
  python "${qa_script}" \\
    -w "." \\
    -m "${meta_file}" \\
    -s "${counter}_${step_name}" \\
    -e "${eval_cfg}" \\
    --docker >> step_log.txt 2>&1
  qa_rc=\\$?
  set -e

  if [ \\$qa_rc -ne 0 ]; then
    echo "###QA_ERROR### QA failed (exit=\\$qa_rc) for ${counter}_${step_name}" >> step_log.txt
    echo "[WARN] QA failed but continuing" >> step_log.txt
    touch QA.fail
  else
    touch QA.ok
  fi

  echo "----- QA END -----" >> step_log.txt


  echo "" >> step_log.txt
  echo "----- REPORT START (${counter}_${step_name}) -----" >> step_log.txt

  set +e
  if [ -f "${prev_report}" ] && [ -s "${prev_report}" ]; then
    python "${report_script}" \\
      --log step_log.txt \\
      --in_pdf "${prev_report}" \\
      --out report_new.pdf \\
      --doc_tsv "${doc_tsv}" >> step_log.txt 2>&1
  else
    python "${report_script}" \\
      --log step_log.txt \\
      --out report_new.pdf \\
      --doc_tsv "${doc_tsv}" >> step_log.txt 2>&1
  fi
  rep_rc=\\$?
  set -e

  if [ \\$rep_rc -ne 0 ]; then
    echo "###REP_ERROR### report failed (exit=\\$rep_rc) for ${counter}_${step_name}" >> step_log.txt
    touch REP.fail
  else
    mv report_new.pdf report.pdf
    touch REP.ok
  fi

  echo "----- REPORT END -----" >> step_log.txt


  touch "${log_file}"
  {
    echo ""
    echo "===== ${counter}_${step_name} (${engine}) ====="
    cat step_log.txt
  } >> "${log_file}"

  cp -f step_log.txt log.txt


  # STOP only on module/report errors
  if grep -Eq "###MOD_ERROR###|###REP_ERROR###" step_log.txt; then
    echo "[FATAL] ${counter}_${step_name}: module/report failed; see step_log.txt" >&2
    exit 100
  fi

  # QA errors: warn only
  if grep -q "###QA_ERROR###" step_log.txt; then
    echo "[WARN] ${counter}_${step_name}: QA failed; continuing" >&2
  fi
  """
}
'''
// -------------------- workflow file with aliases + named workflow component + main: --------------------
def w = new File(wfPath).newWriter()
w << "nextflow.enable.dsl=2\n\n"

w << "include { INIT_PARAMS } from '${modulePath}'\n"
stepsResolved.eachWithIndex { _, i ->
  def idx = String.format('%02d', i+1)
  w << "include { RUN_MODULE as RUN_${idx} } from '${modulePath}'\n"
}

w << """

/*
 * === AUTO-GENERATED (in workDir) ===
 * Source config : ${cfgFile}
 * Start dir     : ${start_dir}
 * Precounter    : ${precounter}
 * Steps         : ${stepsResolved.size()}
 */


workflow GENERATED {

  main:

    def meta       = file("${metaFile.getAbsolutePath()}")
    def log        = file("${logSeed.getAbsolutePath()}")
    def report0    = file("${reportSeed.getAbsolutePath()}")
    def method_doc = file("${methodDoc.getAbsolutePath()}")
    
    def init_ch = INIT_PARAMS(
      Channel.of( tuple(
        "${WORKDIR_ABS}",
        file("${projectDir}/Initialization/Initialize.py"),
        file("${projectDir}/Initialization/dict_flags.pkl"),
        file("${projectDir}/Initialization/dict_parameters.pkl"),
        file("${WORKDIR_ABS}/Parameters.txt"),
        file("${original_data.getAbsolutePath()}"),
        meta
      ))
    )
    
    def ch0 = init_ch.map { param_dir, init_done, setup_log ->
      [
        file("${seed_data.getAbsolutePath()}"),
        file("${seed_idmap.getAbsolutePath()}"),
        file("${seed_mask.getAbsolutePath()}"),
        file("${seed_maskdp.getAbsolutePath()}"),
        file("${logSeed.getAbsolutePath()}"),
        null,
        null,
        report0,
        param_dir
      ]
    }
    def ch00 = ch0//init_ch.map { param_dir, init_done -> ch0 + [param_dir] }
    
    //def ch00 = ch0.combine(init_ch).map { base_tuple, init_tuple ->
     // def (data, id, mask, mask_dp, log_fp, step_log, report_dir, report_pdf) = base_tuple
      // def (param_dir, init_done) = init_tuple
      // pass param_dir forward by appending it
      //[data, id, mask, mask_dp, log_fp, step_log, report_dir, report_pdf, param_dir]
    //}
    

"""

stepsResolved.eachWithIndex { s, i ->
  def idx  = String.format('%02d', i+1)
  def prev = (i == 0) ? '00' : String.format('%02d', i)

  w << """\
    def ch${idx} = RUN_${idx}(
      ch${prev}.map { data, id, mask, mask_dp, log_fp, step_log, report_dir, report_pdf, param_dir ->
        [
          '${WORKDIR_ABS}',
          '${s.engine}',
          '${s.counter}',
          '${s.base}',
          param_dir,
          param_dir,
          file("${s.script}"),
          file("${reportPy.getAbsolutePath()}"),
          file("${qaScript.getAbsolutePath()}"),
          file("${evalConfig.getAbsolutePath()}"),
          file("${maxlfqdll.getAbsolutePath()}"),
          file("${maxlfqdylib.getAbsolutePath()}"),
          "${maxlfqso}",
          file("${figscptn.getAbsolutePath()}"),
          data,
          file("${param_dir}/${s.counter}_${s.base}_parameters.json"),
          meta,
          id,
          mask,
          mask_dp,
          log_fp,
          report_pdf,
          method_doc
        ]
      }
    )

"""
}

w << "}\n"
w.close()

// -------------------- include and run generated workflow --------------------
include { GENERATED } from "${wfPath}"

workflow {
  GENERATED()
}


workflow.onError {
    println ""
    println "[PIPELINE FAILED]"
    println ""
}

workflow.onComplete {

  def runCmd = { String label, String cmd ->
    try {
      println "[CLEANUP] ${label}"
      def proc = ["bash","-lc", cmd].execute()
      proc.waitFor()
      if( proc.exitValue() != 0 ) {
        println "[CLEANUP] WARN: ${label} failed (exit=${proc.exitValue()})"
        def err = proc.err?.text
        if( err ) println err
      }
    } catch(e) {
      println "[CLEANUP] WARN: ${label} exception: ${e}"
    }
  }

  def q = { String s -> s.replace("'", "'\\''") }  // simple single-quote escape for bash strings

  // --------------------
  // SAFE CLEANUP (default)
  // --------------------

  // A) Remove your generated NF modules (safe; they are recreated every run)
  runCmd("Remove _generated_* folders in workDir",
    "find '${q(workDir.toString())}' -maxdepth 1 -type d -name '_generated_*' -print -exec rm -rf {} + || true"
  )

  // B) Remove 0-byte files ONLY from the user-facing output tree (NOT from workDir)
  runCmd("Remove 0-byte files from WORKDIR_ABS",
    "find '${q(WORKDIR_ABS.toString())}' -type f -size 0 -print -delete || true"
  )

  // C) Remove common OS “hidden junk” files from user-facing output tree
  // (Keeps Nextflow internals and logs intact.)
  runCmd("Remove common hidden junk files from WORKDIR_ABS",
    """
      set -euo pipefail
      find '${q(WORKDIR_ABS.toString())}' -type f \\( \
        -name '.DS_Store' -o \
        -name '._*' -o \
        -name 'Thumbs.db' \
      \\) -print -delete || true
    """.stripIndent().trim()
  )

  // --------------------
  // OPTIONAL AGGRESSIVE CLEANUP (breaks -resume)
  // --------------------
  if( params.cleanup_aggressive ) {
    println "[CLEANUP] Aggressive cleanup enabled (this will break -resume)."

    // D) Remove Nextflow work cache completely
    runCmd("Remove workDir (Nextflow cache)",
      "rm -rf '${q(workDir.toString())}' || true"
    )

    // E) Remove Nextflow metadata + common reports in launchDir
    runCmd("Remove Nextflow metadata and reports in launchDir",
      """
        rm -rf '${q(launchDir.toString())}/.nextflow' || true
        rm -f  '${q(launchDir.toString())}/.nextflow.log' || true
        rm -f  '${q(launchDir.toString())}/trace.txt' || true
        rm -f  '${q(launchDir.toString())}/report.html' || true
        rm -f  '${q(launchDir.toString())}/timeline.html' || true
        rm -f  '${q(launchDir.toString())}/dag.html' || true
      """.stripIndent().trim()
    )
  } else {
    println "[CLEANUP] Aggressive cleanup disabled (keeps workDir for -resume)."
  }

  println "[CLEANUP] Done."
}

// nextflow run nextflow.nf --cleanup_aggressive true
