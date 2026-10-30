nextflow.enable.dsl=2

params.fixture = null
params.pipeline_bin = "${projectDir}/../bin"
params.rscript_bin = 'Rscript'
params.qc_project = 'scLRA-LLM'
params.qc_llm_strict = false
params.qc_read_source = 'auto'
params.llm_enabled = false
params.llm_base_url = ''
params.llm_model = 'test-model'
params.use_read_summary = false
params.outdir = 'qc-smoke-results'

include { EXTRACT_SC_QC; LLM_SC_QC_REPORT; RENDER_SC_QC_REPORT } from '../modules/PreliminarySingleCellAnalysis.nf'

workflow {
  fixture = file(params.fixture, checkIfExists: true)
  pipeline_bin = file(params.pipeline_bin, checkIfExists: true)
  summary = params.use_read_summary ? fixture.resolve('bam_summary.tsv') : []
  extracted = EXTRACT_SC_QC(fixture.resolve('gene.rds'), fixture.resolve('transcript.rds'),
    fixture.resolve('isoquant'), fixture.resolve('read_groups.tsv'), fixture.resolve('reference.gtf'),
    summary, [], pipeline_bin)
  interpretation = []
  if (params.llm_enabled) {
    interpreted = LLM_SC_QC_REPORT(extracted.qc_data, pipeline_bin)
    interpretation = interpreted.interpretation
  }
  RENDER_SC_QC_REPORT(extracted.qc_data, interpretation, pipeline_bin)
}
