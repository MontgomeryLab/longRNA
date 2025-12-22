#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

doc: |
  Quantify transcript expression with Salmon
  Supports both single-end and paired-end reads

baseCommand: [salmon, quant]

inputs:
  index_dir:
    type: Directory
    inputBinding:
      prefix: -i
      position: 1
    doc: "Salmon index directory"
  
  library_type:
    type: string
    inputBinding:
      prefix: -l
      position: 2
    doc: "Library type (A = automatic detection)"
  
  cleaned_fastq_r1:
    type: File
    doc: "Input FASTQ file R1 (or single-end)"
  
  cleaned_fastq_r2:
    type: File?
    doc: "Input FASTQ file R2 (for paired-end, leave null for single-end)"
  
  sample_name:
    type: string
    doc: "Sample identifier for output naming"
  
  output_name:
    type: string
    inputBinding:
      prefix: -o
      position: 10
      valueFrom: $(inputs.sample_name)_quant
    doc: "Output directory name"
  
  threads:
    type: int?
    inputBinding:
      prefix: --threads
      position: 11
    doc: "Number of threads to use"
  
  validate_mappings:
    type: boolean?
    default: true
    inputBinding:
      prefix: --validateMappings
      position: 12
    doc: "Validate mappings using alignment-based verification"
  
  gc_bias:
    type: boolean?
    default: true
    inputBinding:
      prefix: --gcBias
      position: 13
    doc: "Perform GC bias correction"

outputs:
  quant_dir:
    type: Directory
    outputBinding:
      glob: $(inputs.sample_name)_quant
    doc: "Directory containing quantification results"
  
  quant_sf:
    type: File
    outputBinding:
      glob: $(inputs.sample_name)_quant/quant.sf
    doc: "Quantification file in Salmon format"

requirements:
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.threads)
    coresMax: $(inputs.threads)

arguments:
  - position: 3
    prefix: "-1"
    valueFrom: |
      ${
        if (inputs.cleaned_fastq_r2) {
          return inputs.cleaned_fastq_r1.path;
        } else {
          return null;
        }
      }
  - position: 4
    prefix: "-2"
    valueFrom: |
      ${
        if (inputs.cleaned_fastq_r2) {
          return inputs.cleaned_fastq_r2.path;
        } else {
          return null;
        }
      }
  - position: 3
    prefix: "-r"
    valueFrom: |
      ${
        if (inputs.cleaned_fastq_r2) {
          return null;
        } else {
          return inputs.cleaned_fastq_r1.path;
        }
      }
