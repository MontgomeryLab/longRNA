#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

doc: |
  Standalone DESeq2 analysis for RSEM output
  Reads .genes.results files from a directory (not from cache)

baseCommand: Rscript

inputs:
  r_script:
    type: File
    inputBinding:
      position: 1
    doc: "R script for DESeq2 analysis with RSEM data"

  gtf_annotation:
    type: File
    inputBinding:
      position: 2
    doc: "GTF annotation file for gene name mapping"

  sample_csv_file:
    type: File?
    inputBinding:
      position: 3
    doc: "Optional CSV file with detailed sample metadata"

  geneid2name_file:
    type: File?
    inputBinding:
      position: 4
    doc: "Optional Excel file mapping gene IDs to gene names"

  genes_results_files:
    type: File[]
    doc: "Array of RSEM .genes.results files"

  sample_names:
    type: string[]
    doc: "Sample names corresponding to results files"

  sample_groups:
    type: string[]
    doc: "Sample groups for differential expression"

  output_prefix:
    type: string
    default: deseq2_rsem
    doc: "Prefix for output files"

outputs:
  results_csv:
    type: File[]
    outputBinding:
      glob: "*_results*.csv"

  normalized_counts:
    type: File
    outputBinding:
      glob: "*_normalized_counts.csv"

  gene_counts_raw:
    type: File
    outputBinding:
      glob: "gene_counts_raw.tsv"

  gene_counts_normalized:
    type: File
    outputBinding:
      glob: "gene_counts_normalized.tsv"

  genename_counts_raw:
    type: File?
    outputBinding:
      glob: "genename_counts_raw.tsv"

  genename_counts_normalized:
    type: File?
    outputBinding:
      glob: "genename_counts_normalized.tsv"

  ma_plot:
    type: File[]
    outputBinding:
      glob: "*MA_plot*.pdf"

  volcano_plot:
    type: File[]
    outputBinding:
      glob: "*volcano_plot*.pdf"

  pca_plot:
    type: File
    outputBinding:
      glob: "*PCA_plot*.pdf"

  replicate_scatter:
    type: File
    outputBinding:
      glob: "*replicate_scatter*.pdf"

  mean_scatter:
    type: File
    outputBinding:
      glob: "*mean_scatter*.pdf"

requirements:
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entryname: sample_info.csv
        entry: |
          sample,condition
          ${
            var lines = [];
            for (var i = 0; i < inputs.sample_names.length; i++) {
              lines.push(inputs.sample_names[i] + "," + inputs.sample_groups[i]);
            }
            return lines.join("\n");
          }
      - entryname: genes_files.txt
        entry: |
          ${
            var lines = [];
            for (var i = 0; i < inputs.genes_results_files.length; i++) {
              lines.push(inputs.genes_results_files[i].path);
            }
            return lines.join("\n");
          }
