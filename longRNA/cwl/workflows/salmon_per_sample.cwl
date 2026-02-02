cwlVersion: v1.0
class: Workflow

doc: |
  Per-sample Salmon workflow - processes ONE sample through fastp and Salmon
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

  # Salmon index directory (pre-built)
  index_dir:
    type: Directory
    doc: "Pre-built Salmon index directory"

  # Parameters
  library_type:
    type: string?
    default: "A"
    doc: "Library type (A=auto-detect)"

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

  # Step 2: Quantify with Salmon
  salmon_quant:
    run: ../tools/salmon-quant.cwl
    in:
      index_dir: index_dir
      cleaned_fastq_r1: fastp_qc/cleaned_fastq_r1
      cleaned_fastq_r2: fastp_qc/cleaned_fastq_r2
      sample_name: sample_name
      output_name: sample_name
      library_type: library_type
      threads: threads
    out: [quant_dir, quant_sf]

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

  # Salmon outputs
  quant_dir:
    type: Directory
    outputSource: salmon_quant/quant_dir
  quant_sf:
    type: File
    outputSource: salmon_quant/quant_sf
