#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

doc: |
  Quality control and adapter trimming with fastp
  Supports both single-end and paired-end reads

baseCommand: fastp

inputs:
  input_fastq_r1:
    type: File
    inputBinding:
      prefix: -i
    doc: "Input FASTQ file (R1 for paired-end, or single-end)"
  
  input_fastq_r2:
    type: File?
    inputBinding:
      prefix: -I
    doc: "Input FASTQ file R2 (for paired-end only)"
  
  output_name_r1:
    type: string
    inputBinding:
      prefix: -o
      valueFrom: $(inputs.sample_name)_R1_clean.fastq
    doc: "Output filename for cleaned R1 FASTQ"
  
  output_name_r2:
    type: string?
    inputBinding:
      prefix: -O
      valueFrom: $(inputs.sample_name)_R2_clean.fastq
    doc: "Output filename for cleaned R2 FASTQ (paired-end only)"
  
  sample_name:
    type: string
    doc: "Sample identifier for output naming"
  
  html_report:
    type: string
    inputBinding:
      prefix: -h
      valueFrom: $(inputs.sample_name)_fastp.html
    doc: "HTML report filename"
  
  json_report:
    type: string
    inputBinding:
      prefix: -j
      valueFrom: $(inputs.sample_name)_fastp.json
    doc: "JSON report filename"
  
  threads:
    type: int?
    inputBinding:
      prefix: --thread
    doc: "Number of threads to use"
  
  qualified_quality_phred:
    type: int?
    inputBinding:
      prefix: --qualified_quality_phred
    doc: "Quality value that a base is qualified"

  length_required:
    type: int?
    inputBinding:
      prefix: --length_required
    doc: "Minimum read length required"

  detect_adapter_for_pe:
    type: boolean?
    inputBinding:
      prefix: --detect_adapter_for_pe
    doc: "Enable adapter auto-detection for paired-end reads"

  trim_poly_x:
    type: boolean?
    inputBinding:
      prefix: --trim_poly_x
    doc: "Enable polyX (polyA/polyT/polyG/polyC) tail trimming"

  trim_poly_g:
    type: boolean?
    inputBinding:
      prefix: --trim_poly_g
    doc: "Enable polyG tail trimming (for NovaSeq/NextSeq)"

requirements:
  - class: ResourceRequirement
    coresMin: $(inputs.threads)
    coresMax: $(inputs.threads)

outputs:
  cleaned_fastq_r1:
    type: File
    outputBinding:
      glob: $(inputs.sample_name)_R1_clean.fastq
    doc: "Cleaned R1 FASTQ file"
  
  cleaned_fastq_r2:
    type: File?
    outputBinding:
      glob: $(inputs.sample_name)_R2_clean.fastq
    doc: "Cleaned R2 FASTQ file (paired-end only)"
  
  html_report_file:
    type: File
    outputBinding:
      glob: $(inputs.sample_name)_fastp.html
    doc: "HTML QC report"
  
  json_report_file:
    type: File
    outputBinding:
      glob: $(inputs.sample_name)_fastp.json
    doc: "JSON QC report"
