nextflow.enable.dsl=2

/*
 * Prepare gene expression matrix from IsoQuant output
 */
process PREPARE_GENE_MATRIX {
  label 'sclong'
  tag "prepare_gene_matrix"
  
  publishDir "${params.outdir ?: 'Result'}/PreliminarySingleCellAnalysis/gene_matrix", mode: 'copy'
  
  input:
    path isoquant_dir
  
  output:
    path "gene_matrix", emit: gene_matrix_dir
  
  script:
    """
    # Create output directory
    mkdir -p gene_matrix
    
    # Check if required files exist
    if [[ ! -f "${isoquant_dir}/OUT/OUT.discovered_gene_grouped_counts.barcodes.tsv" ]]; then
      echo "ERROR: Gene barcodes file not found" >&2
      exit 1
    fi
    
    if [[ ! -f "${isoquant_dir}/OUT/OUT.discovered_gene_grouped_counts.features.tsv" ]]; then
      echo "ERROR: Gene features file not found" >&2
      exit 1
    fi
    
    if [[ ! -f "${isoquant_dir}/OUT/OUT.discovered_gene_grouped_counts.matrix.mtx" ]]; then
      echo "ERROR: Gene matrix file not found" >&2
      exit 1
    fi
    
    # Copy and rename files
    cp "${isoquant_dir}/OUT/OUT.discovered_gene_grouped_counts.barcodes.tsv" gene_matrix/barcodes.tsv
    cp "${isoquant_dir}/OUT/OUT.discovered_gene_grouped_counts.features.tsv" gene_matrix/genes.tsv
    cp "${isoquant_dir}/OUT/OUT.discovered_gene_grouped_counts.matrix.mtx" gene_matrix/matrix.mtx
    
    echo "Gene matrix prepared in: gene_matrix"
    echo "Files:"
    ls -la gene_matrix/
    """
}

/*
 * Prepare transcript expression matrix from IsoQuant output
 */
process PREPARE_TRANSCRIPT_MATRIX {
  label 'sclong'
  tag "prepare_transcript_matrix"
  
  publishDir "${params.outdir ?: 'Result'}/PreliminarySingleCellAnalysis/transcript_matrix", mode: 'copy'
  
  input:
    path isoquant_dir
  
  output:
    path "transcript_matrix", emit: transcript_matrix_dir
  
  script:
    """
    # Create output directory
    mkdir -p transcript_matrix
    
    # Check if required files exist
    if [[ ! -f "${isoquant_dir}/OUT/OUT.discovered_transcript_grouped_counts.barcodes.tsv" ]]; then
      echo "ERROR: Transcript barcodes file not found" >&2
      exit 1
    fi
    
    if [[ ! -f "${isoquant_dir}/OUT/OUT.discovered_transcript_grouped_counts.features.tsv" ]]; then
      echo "ERROR: Transcript features file not found" >&2
      exit 1
    fi
    
    if [[ ! -f "${isoquant_dir}/OUT/OUT.discovered_transcript_grouped_counts.matrix.mtx" ]]; then
      echo "ERROR: Transcript matrix file not found" >&2
      exit 1
    fi
    
    # Copy and rename files
    cp "${isoquant_dir}/OUT/OUT.discovered_transcript_grouped_counts.barcodes.tsv" transcript_matrix/barcodes.tsv
    cp "${isoquant_dir}/OUT/OUT.discovered_transcript_grouped_counts.features.tsv" transcript_matrix/genes.tsv
    cp "${isoquant_dir}/OUT/OUT.discovered_transcript_grouped_counts.matrix.mtx" transcript_matrix/matrix.mtx
    
    echo "Transcript matrix prepared in: transcript_matrix"
    echo "Files:"
    ls -la transcript_matrix/
    """
}

/*
 * Run single-cell analysis with R
 */
process RUN_SC_ANALYSIS {
  label 'sclong'
  tag "sc_analysis"
  
  publishDir "${params.outdir ?: 'Result'}/PreliminarySingleCellAnalysis", mode: 'copy'
  
  input:
    path gene_matrix_dir
    path transcript_matrix_dir
    path pipeline_bin
  
  output:
    path "sce.gene.rds", emit: seurat_gene
    path "sce.tr.rds", emit: seurat_tr
    path "*.pdf", emit: plots
    path "*.txt", emit: marker_tables
    path "qc_metrics.tsv", emit: qc_metrics
    path "*.csv", optional: true, emit: annotation_files
    path "*.requests.json", optional: true, emit: annotation_requests
  
  script:
    """
    # Run R analysis script
    "${params.rscript_bin ?: 'Rscript'}" "${pipeline_bin}/preliminary_single_cell_analysis.R" \
      --gene_matrix "${gene_matrix_dir}" \
      --transcript_matrix "${transcript_matrix_dir}" \
      --output_dir "." \
      --species "${params.species ?: 'Human'}" \
      --tissue "${params.tissue ?: 'Skin'}" \
      --project_dir "." \
      --min_cells ${params.min_cells ?: 10} \
      --min_features ${params.min_features ?: 200} \
      --llm_enabled "${params.llm_enabled}" \
      --llm_base_url "${params.llm_base_url}" \
      --llm_model "${params.llm_model}"
    
    echo "Single-cell analysis completed"
    echo "Generated files:"
    ls -la *.rds *.pdf *.txt *.tsv *.csv 2>/dev/null || true
    """
}

/*
 * Extract run-specific metrics from the aligned, post-filter cell cohort.
 * FASTQ read-group counts are explicitly distinguished from optional BAM summaries.
 */
def qcQuote(value) {
  return "'" + value.toString().replace("'", "'\"'\"'") + "'"
}

process EXTRACT_SC_QC {
  label 'sclong'
  tag "filtered_cell_qc"

  input:
    path seurat_gene
    path seurat_tr
    path isoquant_dir
    path read_groups
    path reference_gtf
    path read_summary
    path classification
    path pipeline_bin

  output:
    path "qc_data", emit: qc_data

  script:
    def summaryArg = read_summary ? "--read_summary ${qcQuote(read_summary)}" : "--read_groups ${qcQuote(read_groups)}"
    def classificationArg = classification ? "--classification ${qcQuote(classification)}" : ''
    def projectArg = qcQuote(params.qc_project ?: 'scLRA-LLM')
    def sourceArg = qcQuote(params.qc_read_source ?: 'auto')
    """
    "${params.rscript_bin ?: 'Rscript'}" "${pipeline_bin}/qc/extract_qc_metrics.R" \
      --gene "${seurat_gene}" --transcript "${seurat_tr}" \
      --isoquant_dir "${isoquant_dir}" --reference_gtf "${reference_gtf}" \
      --output_dir qc_data ${classificationArg}
    python3 "${pipeline_bin}/qc/prepare_qc_report.py" \
      --data_dir qc_data ${summaryArg} --project ${projectArg} \
      --isoquant_dir "${isoquant_dir}" --read_source ${sourceArg}
    "${params.rscript_bin ?: 'Rscript'}" "${pipeline_bin}/qc/plot_qc.R" \
      qc_data/cells.tsv.gz qc_data/qc_distributions.svg
    """
}

/* Only aggregate metrics and definitions are sent to the configured provider. */
process LLM_SC_QC_REPORT {
  label 'preprocess'
  tag "llm_sc_qc_report"

  // Retain the historical text-report path for existing consumers.
  publishDir "${params.outdir ?: 'Result'}/PreliminarySingleCellAnalysis/deepseek_qc", mode: 'copy',
    pattern: 'sc_qc_report.txt'

  input:
    path qc_data
    path pipeline_bin

  output:
    path "qc_interpretation", emit: interpretation
    path "sc_qc_report.txt", emit: sc_qc_report

  script:
    def strictArg = params.qc_llm_strict ? '--strict' : ''
    def baseUrlArg = qcQuote(params.llm_base_url)
    def modelArg = qcQuote(params.llm_model)
    """
    python3 "${pipeline_bin}/qc/interpret_qc.py" \
      --data "${qc_data}/qc_report.json" --output_dir qc_interpretation \
      --base_url ${baseUrlArg} --model ${modelArg} ${strictArg}
    cp qc_interpretation/sc_qc_report.txt sc_qc_report.txt
    """
}

/* The portable English HTML is generated with or without an API key. */
process RENDER_SC_QC_REPORT {
  label 'preprocess'
  tag "render_sc_qc_report"

  publishDir "${params.outdir ?: 'Result'}/PreliminarySingleCellAnalysis", mode: 'copy'

  input:
    path qc_data
    path interpretation
    path pipeline_bin

  output:
    path "qc_report", emit: report_dir

  script:
    def interpretedArg = interpretation ? "--interpreted ${qcQuote(interpretation.toString() + '/qc_report.json')}" : ''
    """
    python3 "${pipeline_bin}/qc/publish_qc_report.py" \
      --data_dir "${qc_data}" ${interpretedArg} --output_dir qc_report
    """
}

/*
 * Main workflow for single-cell analysis
 */
workflow PreliminarySingleCellAnalysis {
  take:
    isoquant_dir
    read_groups
    reference_gtf
  
  main:
    pipeline_bin = file("${projectDir}/bin", checkIfExists: true)

    // Prepare gene and transcript matrices
    gene_matrix = PREPARE_GENE_MATRIX(isoquant_dir)
    transcript_matrix = PREPARE_TRANSCRIPT_MATRIX(isoquant_dir)
    
    // Run single-cell analysis
    sc_results = RUN_SC_ANALYSIS(
      gene_matrix.gene_matrix_dir,
      transcript_matrix.transcript_matrix_dir,
      pipeline_bin
    )
    
    report_dirs = channel.empty()
    report_html = channel.empty()
    report_tsv = channel.empty()
    report_json = channel.empty()
    report_text = channel.empty()
    llm_report_text = channel.empty()
    if (params.qc_report_enabled) {
      read_summary = params.qc_read_summary ? file(params.qc_read_summary, checkIfExists: true) : []
      classification = params.qc_classification ? file(params.qc_classification, checkIfExists: true) : []
      qc_data = EXTRACT_SC_QC(sc_results.seurat_gene, sc_results.seurat_tr,
        isoquant_dir, read_groups, reference_gtf, read_summary, classification, pipeline_bin)
      interpretation = []
      if (params.llm_enabled) {
        interpreted = LLM_SC_QC_REPORT(qc_data.qc_data, pipeline_bin)
        interpretation = interpreted.interpretation
        llm_report_text = interpreted.sc_qc_report
      }
      rendered = RENDER_SC_QC_REPORT(qc_data.qc_data, interpretation, pipeline_bin)
      report_dirs = rendered.report_dir
      report_html = report_dirs.map { dir -> dir.resolve('scLRA_LLM_QC_report.html') }
      report_tsv = report_dirs.map { dir -> dir.resolve('scLRA_LLM_QC_report.tsv') }
      report_json = report_dirs.map { dir -> dir.resolve('data/qc_report.json') }
      report_text = report_dirs.map { dir -> dir.resolve('sc_qc_report.txt') }
    }
  
  emit:
    seurat_gene = sc_results.seurat_gene     // sce.gene.rds
    seurat_tr   = sc_results.seurat_tr       // sce.tr.rds
    plots = sc_results.plots
    marker_tables = sc_results.marker_tables
    qc_metrics = sc_results.qc_metrics
    annotation_files = sc_results.annotation_files  // optional
    qc_report_dir = report_dirs
    qc_report_html = report_html
    qc_report_tsv = report_tsv
    qc_report_json = report_json
    llm_qc_report = llm_report_text
    // Backward-compatible text output, also available without an LLM.
    sc_qc_report = report_text
}
