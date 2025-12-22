#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

doc: |
  Build Salmon index from transcriptome FASTA

baseCommand: [salmon, index]

inputs:
  transcriptome:
    type: File
    inputBinding:
      prefix: -t
    doc: "Transcriptome FASTA file (can be gzipped)"

  index_name:
    type: string
    default: salmon_index
    inputBinding:
      prefix: -i
    doc: "Name for the index directory"

  kmer_size:
    type: int?
    default: 31
    inputBinding:
      prefix: -k
    doc: "K-mer size for indexing"

  threads:
    type: int?
    inputBinding:
      prefix: --threads
    doc: "Number of threads to use"

outputs:
  index_dir:
    type: Directory
    outputBinding:
      glob: $(inputs.index_name)
    doc: "Salmon index directory"

requirements:
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.threads)
    coresMax: $(inputs.threads)
