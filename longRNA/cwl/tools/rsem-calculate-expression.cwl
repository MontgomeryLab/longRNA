#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

doc: |
  Quantify gene and transcript expression using RSEM with STAR aligner
  Supports both single-end and paired-end reads

baseCommand: [bash]

arguments:
  - $(inputs.script_file.path)

inputs:
  script_file:
    type: File
    default:
      class: File
      location: rsem_wrapper.sh

  reference_dir:
    type: Directory

  reference_name:
    type: string

  cleaned_fastq_r1:
    type: File

  cleaned_fastq_r2:
    type: File?

  sample_name:
    type: string

  is_paired_end:
    type: boolean

  use_star:
    type: boolean
    default: true

  star_gzipped:
    type: boolean
    default: false

  forward_prob:
    type: float?

  estimate_rspd:
    type: boolean?
    default: true

  calc_ci:
    type: boolean?
    default: false

  seed:
    type: int?
    default: 12345

  threads:
    type: int?

  perfect_match_only:
    type: boolean?
    default: false

  unique_mappers_only:
    type: boolean?
    default: false

outputs:
  genes_results:
    type: File
    outputBinding:
      glob: $(inputs.sample_name).genes.results

  isoforms_results:
    type: File
    outputBinding:
      glob: $(inputs.sample_name).isoforms.results

  transcript_bam:
    type: File?
    outputBinding:
      glob: $(inputs.sample_name).transcript.bam

  stat_dir:
    type: Directory
    outputBinding:
      glob: $(inputs.sample_name).stat

  quant_dir:
    type: Directory
    outputBinding:
      glob: "."

  rsem_log:
    type: File
    outputBinding:
      glob: $(inputs.sample_name).rsem.log

  rsem_err:
    type: File
    outputBinding:
      glob: $(inputs.sample_name).rsem.err

  star_log_final:
    type: File?
    outputBinding:
      glob: $(inputs.sample_name).temp/Log.final.out

  star_log_progress:
    type: File?
    outputBinding:
      glob: $(inputs.sample_name).temp/Log.progress.out

  star_log_full:
    type: File?
    outputBinding:
      glob: $(inputs.sample_name).temp/Log.out

requirements:
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.threads)
    coresMax: $(inputs.threads)
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.reference_dir)
  - class: EnvVarRequirement
    envDef:
      USE_STRICT_FILTERING: "$(inputs.perfect_match_only || inputs.unique_mappers_only ? 'true' : 'false')"
      PERFECT_MATCH_ONLY: "$(inputs.perfect_match_only ? 'true' : 'false')"
      UNIQUE_MAPPERS_ONLY: "$(inputs.unique_mappers_only ? 'true' : 'false')"
      USE_STAR: "$(inputs.use_star ? 'true' : 'false')"
      STAR_GZIPPED: "$(inputs.star_gzipped ? 'true' : 'false')"
      ESTIMATE_RSPD: "$(inputs.estimate_rspd ? 'true' : 'false')"
      CALC_CI: "$(inputs.calc_ci ? 'true' : 'false')"
      REF_DIR: "$(inputs.reference_dir.path)"
      REF_NAME: "$(inputs.reference_name)"
      R1: "$(inputs.cleaned_fastq_r1.path)"
      R2: "$(inputs.is_paired_end && inputs.cleaned_fastq_r2 ? inputs.cleaned_fastq_r2.path : '')"
      SAMPLE_NAME: "$(inputs.sample_name)"
      THREADS: "$(inputs.threads.toString())"
      FORWARD_PROB: "$(inputs.forward_prob !== null && inputs.forward_prob !== undefined ? inputs.forward_prob.toString() : '')"
      SEED: "$(inputs.seed !== null && inputs.seed !== undefined ? inputs.seed.toString() : '')"

stdout: $(inputs.sample_name).wrapper.log
stderr: $(inputs.sample_name).wrapper.err
