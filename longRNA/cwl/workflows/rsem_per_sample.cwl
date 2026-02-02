cwlVersion: v1.0
class: Workflow

doc: |
  Per-sample RSEM workflow - processes ONE sample through fastp and RSEM
  Used for sequential processing with cleanup between samples

requirements:
  - class: StepInputExpressionRequirement

inputs:
  # Input FASTQ files for ONE sample
  fastq_r1:
    type: File
    doc: "Input FASTQ R1 file"

  fastq_r2:
    type: File?
    doc: "Input FASTQ R2 file (for paired-end)"

  sample_name:
    type: string
    doc: "Sample name"

  is_paired_end:
    type: boolean
    doc: "Whether data is paired-end"

  # Reference directory (pre-built)
  reference_dir:
    type: Directory
    doc: "Pre-built RSEM reference directory"

  reference_name:
    type: string
    doc: "RSEM reference name"

  # Parameters
  forward_prob:
    type: float
    doc: "Forward probability for strand-specificity"

  threads:
    type: int?
    default: 8

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

  perfect_match_only:
    type: boolean?

  unique_mappers_only:
    type: boolean?

steps:

  # Step 1: Quality control with fastp
  fastp_qc:
    run: ../tools/fastp.cwl
    in:
      input_fastq_r1: fastq_r1
      input_fastq_r2: fastq_r2
      sample_name: sample_name
      output_name_r1: sample_name
      output_name_r2: sample_name
      html_report: sample_name
      json_report: sample_name
      threads: threads
      qualified_quality_phred: qualified_quality_phred
      length_required: length_required
      detect_adapter_for_pe: detect_adapter_for_pe
      trim_poly_x: trim_poly_x
      trim_poly_g: trim_poly_g
    out: [cleaned_fastq_r1, cleaned_fastq_r2, html_report_file, json_report_file]

  # Step 2: Quantify with RSEM + STAR
  rsem_quant:
    run: ../tools/rsem-calculate-expression.cwl
    in:
      reference_dir: reference_dir
      reference_name: reference_name
      cleaned_fastq_r1: fastp_qc/cleaned_fastq_r1
      cleaned_fastq_r2: fastp_qc/cleaned_fastq_r2
      sample_name: sample_name
      is_paired_end: is_paired_end
      forward_prob: forward_prob
      threads: threads
      perfect_match_only: perfect_match_only
      unique_mappers_only: unique_mappers_only
      star_gzipped:
        default: true
    out: [genes_results, isoforms_results, sorted_bam, stat_dir, quant_dir,
          rsem_log, rsem_err, star_log_final, star_log_progress, star_log_full]

outputs:
  # Fastp outputs
  cleaned_fastq_r1:
    type: File
    outputSource: fastp_qc/cleaned_fastq_r1
  cleaned_fastq_r2:
    type: File?
    outputSource: fastp_qc/cleaned_fastq_r2
  fastp_html:
    type: File
    outputSource: fastp_qc/html_report_file
  fastp_json:
    type: File
    outputSource: fastp_qc/json_report_file

  # RSEM outputs
  genes_results:
    type: File
    outputSource: rsem_quant/genes_results
  isoforms_results:
    type: File
    outputSource: rsem_quant/isoforms_results
  sorted_bam:
    type: File?
    outputSource: rsem_quant/sorted_bam
  rsem_log:
    type: File
    outputSource: rsem_quant/rsem_log
  rsem_err:
    type: File
    outputSource: rsem_quant/rsem_err
  star_log_final:
    type: File?
    outputSource: rsem_quant/star_log_final
  star_log_progress:
    type: File?
    outputSource: rsem_quant/star_log_progress
  star_log_full:
    type: File?
    outputSource: rsem_quant/star_log_full
