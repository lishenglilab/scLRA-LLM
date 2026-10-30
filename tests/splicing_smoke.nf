nextflow.enable.dsl=2

params.fixture = null
params.pipeline_bin = "${projectDir}/../bin"
params.rscript_bin = 'Rscript'
params.outdir = 'splicing-smoke-results'
params.split_size = 2
params.ase_fdr = 0.05
params.ase_fold_change = 1.5
params.ase_pseudocount = 1e-6
params.sample_stage_mapping = [['alpha', 'case'], ['beta', 'control']]

include { SPLIT_EXPRESSION; CALCULATE_PSI; MERGE_AS_RESULTS; ASE_DIFFERENTIAL } from '../modules/AlternativeSplicing.nf'

workflow {
  fixture = file(params.fixture, checkIfExists: true)
  pipeline_bin = file(params.pipeline_bin, checkIfExists: true)
  split = SPLIT_EXPRESSION(fixture.resolve('transcript.rds'), pipeline_bin)
  psi = CALCULATE_PSI(fixture.resolve('events.ioe'), split.split_files.flatten())
  merged = MERGE_AS_RESULTS(psi.psi_file.collect(), pipeline_bin)
  comparisons = channel.of(tuple('case', 'control', 'case_vs_control'), tuple('control', 'case', 'control_vs_case'))
  ASE_DIFFERENTIAL(merged.merged_psi, comparisons, pipeline_bin)
}
