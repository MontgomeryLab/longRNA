cwlVersion: v1.0
class: Workflow

doc: |
  RNA-seq analysis workflow using RSEM + STAR
  
  This workflow performs:
  1. Quality control and trimming with fastp (single or paired-end)
  2. RSEM reference preparation with STAR aligner
  3. Quantification with RSEM + STAR
  4. Differential expression analysis with DESeq2
  
  Save as: cwl/workflows/rnaseq_rsem_workflow.cwl

requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: SubworkflowFeatureRequirement

inputs:
  # Input FASTQ files
  in_fastq_r1:
    type: File[]
    doc: "Input FASTQ R1 files (or single-end files)"
  
  in_fastq_r2:
    type: File[]?
    doc: "Input FASTQ R2 files (for paired-end, omit for single-end)"
  
  is_paired_end:
    type: boolean
    doc: "Whether data is paired-end"
  
  sample_basenames:
    type: string[]
    doc: "Sample names corresponding to input files"
  
  sample_groups:
    type: string[]
    doc: "Sample groups for differential expression"

  sample_csv_file:
    type: File?
    doc: "Optional CSV file with detailed sample metadata (group names, replicates, control designation)"

  run_name:
    type: string
    doc: "Name prefix for output files"
  
  # Reference files
  genome_fasta:
    type: File
    doc: "Genome FASTA file"
  
  gtf_annotation:
    type: File
    doc: "GTF annotation file"
  
  # RSEM parameters
  rsem_reference_name:
    type: string
    doc: "Name for RSEM reference"

  rsem_sjdboverhang:
    type: int
    doc: "STAR sjdbOverhang (read_length - 1)"

  rsem_genomeSAindexNbases:
    type: int?
    doc: "STAR genomeSAindexNbases (12 for C. elegans, 14 for human)"

  rsem_forward_prob:
    type: float
    doc: "Forward probability (0=unstranded, 1=forward stranded)"
  
  # Tool parameters
  fastp_threads:
    type: int?

  qualified_quality_phred:
    type: int?

  length_required:
    type: int?

  detect_adapter_for_pe:
    type: boolean?

  trim_poly_x:
    type: boolean?

  trim_poly_g:
    type: boolean?

  rsem_threads:
    type: int?

  perfect_match_only:
    type: boolean?
    doc: "Only allow perfectly matching reads (no mismatches)"

  unique_mappers_only:
    type: boolean?
    doc: "Only allow uniquely mapping reads"

  alpha:
    type: float?

  lfc_threshold:
    type: float?
  
  # Output directory names
  dir_name_fastp: string
  dir_name_rsem: string
  dir_name_rsem_bams: string
  dir_name_rsem_logs: string
  dir_name_analysis: string
  dir_name_gene_counts: string

steps:
  
  # Step 1: Quality control with fastp
  fastp_qc:
    run: ../tools/fastp.cwl
    scatter: [input_fastq_r1, input_fastq_r2, sample_name, output_name_r1, output_name_r2, html_report, json_report]
    scatterMethod: dotproduct
    in:
      input_fastq_r1: in_fastq_r1
      input_fastq_r2: in_fastq_r2
      sample_name: sample_basenames
      output_name_r1: sample_basenames
      output_name_r2: sample_basenames
      html_report: sample_basenames
      json_report: sample_basenames
      threads: fastp_threads
      qualified_quality_phred: qualified_quality_phred
      length_required: length_required
      detect_adapter_for_pe: detect_adapter_for_pe
      trim_poly_x: trim_poly_x
      trim_poly_g: trim_poly_g
    out: [cleaned_fastq_r1, cleaned_fastq_r2, html_report_file, json_report_file]
  
  # Step 2: Build RSEM reference with STAR
  rsem_prepare_ref:
    run: ../tools/build-rsem-star-reference.cwl
    in:
      genome_fasta: genome_fasta
      gtf_annotation: gtf_annotation
      reference_name: rsem_reference_name
      sjdboverhang: rsem_sjdboverhang
      genomeSAindexNbases: rsem_genomeSAindexNbases
      threads: fastp_threads
    out: [reference_dir]
  
  # Step 3: Quantify with RSEM + STAR
  rsem_quant:
    run: ../tools/rsem-calculate-expression.cwl
    scatter: [cleaned_fastq_r1, cleaned_fastq_r2, sample_name]
    scatterMethod: dotproduct
    in:
      reference_dir: rsem_prepare_ref/reference_dir
      reference_name: rsem_reference_name
      cleaned_fastq_r1: fastp_qc/cleaned_fastq_r1
      cleaned_fastq_r2: fastp_qc/cleaned_fastq_r2
      sample_name: sample_basenames
      is_paired_end: is_paired_end
      forward_prob: rsem_forward_prob
      threads: rsem_threads
      perfect_match_only: perfect_match_only
      unique_mappers_only: unique_mappers_only
    out: [genes_results, isoforms_results, transcript_bam, stat_dir, quant_dir,
          rsem_log, rsem_err, star_log_final, star_log_progress, star_log_full]
  
  # Step 4: Differential expression analysis with DESeq2
  deseq_analysis:
    run: ../tools/deseq-analysis-rsem.cwl
    in:
      r_script:
        default:
          class: File
          location: ../../scripts/deseq_analysis_rsem.R
      gtf_annotation: gtf_annotation
      sample_csv_file: sample_csv_file
      genes_results_files: rsem_quant/genes_results
      sample_names: sample_basenames
      sample_groups: sample_groups
      output_prefix: run_name
      alpha: alpha
      lfc_threshold: lfc_threshold
    out: [results_csv, normalized_counts, gene_counts_raw, gene_counts_normalized,
          ma_plot, volcano_plot, pca_plot, replicate_scatter, mean_scatter]
  
  # Step 5: Organize outputs into directories
  organize_fastp:
    run: ../tools/make-subdir.cwl
    in:
      dir_files:
        source: [fastp_qc/cleaned_fastq_r1, fastp_qc/cleaned_fastq_r2,
                 fastp_qc/html_report_file, fastp_qc/json_report_file]
        linkMerge: merge_flattened
      dir_name: dir_name_fastp
    out: [subdir]
  
  organize_rsem:
    run: ../tools/make-subdir.cwl
    in:
      dir_files: rsem_quant/genes_results
      dir_name: dir_name_rsem
    out: [subdir]

  # organize_bams step removed - not needed with --no-bam-output in strict filtering mode

  organize_logs:
    run: ../tools/make-subdir.cwl
    in:
      dir_files:
        source: [rsem_quant/rsem_log, rsem_quant/rsem_err,
                 rsem_quant/star_log_final, rsem_quant/star_log_progress,
                 rsem_quant/star_log_full]
        linkMerge: merge_flattened
      dir_name: dir_name_rsem_logs
    out: [subdir]

  organize_analysis:
    run: ../tools/make-subdir.cwl
    in:
      dir_files: 
        source: [deseq_analysis/results_csv, deseq_analysis/normalized_counts,
                 deseq_analysis/ma_plot, deseq_analysis/volcano_plot,
                 deseq_analysis/pca_plot, deseq_analysis/replicate_scatter,
                 deseq_analysis/mean_scatter]
        linkMerge: merge_flattened
      dir_name: dir_name_analysis
    out: [subdir]
  
  organize_gene_counts:
    run: ../tools/make-subdir.cwl
    in:
      dir_files:
        source: [deseq_analysis/gene_counts_raw, deseq_analysis/gene_counts_normalized]
        linkMerge: merge_flattened
      dir_name: dir_name_gene_counts
    out: [subdir]

outputs:
  
  # Organized directories
  fastp_dir:
    type: Directory
    outputSource: organize_fastp/subdir
    doc: "Directory containing fastp cleaned FASTQ files"
  
  rsem_dir:
    type: Directory
    outputSource: organize_rsem/subdir
    doc: "Directory containing RSEM gene quantification results"

  # bam_dir removed - not generated with --no-bam-output in strict filtering mode

  logs_dir:
    type: Directory
    outputSource: organize_logs/subdir
    doc: "Directory containing RSEM and STAR log files"

  analysis_dir:
    type: Directory
    outputSource: organize_analysis/subdir
    doc: "Directory containing DESeq2 analysis results"
  
  gene_counts_dir:
    type: Directory
    outputSource: organize_gene_counts/subdir
    doc: "Directory containing gene count matrices (raw and normalized)"
