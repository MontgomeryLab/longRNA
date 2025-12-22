#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

doc: |
  Perform differential expression analysis using DESeq2 with RSEM output
  
  Save as: cwl/tools/deseq-analysis-rsem.cwl

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
    doc: "Optional CSV file with detailed sample metadata (group names, replicates, control designation)"

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
  
  alpha:
    type: float?
    default: 0.05
    doc: "Significance threshold (FDR)"
  
  lfc_threshold:
    type: float?
    default: 1.0
    doc: "Log2 fold change threshold"

outputs:
  results_csv:
    type: File[]
    outputBinding:
      glob: "*_results*.csv"
    doc: "DESeq2 results tables (one per comparison)"

  normalized_counts:
    type: File
    outputBinding:
      glob: "*_normalized_counts.csv"
    doc: "Normalized count matrix"

  gene_counts_raw:
    type: File
    outputBinding:
      glob: "gene_counts_raw.tsv"
    doc: "Raw gene counts in TSV format"

  gene_counts_normalized:
    type: File
    outputBinding:
      glob: "gene_counts_normalized.tsv"
    doc: "Normalized gene counts in TSV format"

  ma_plot:
    type: File[]
    outputBinding:
      glob: "*MA_plot*.pdf"
    doc: "MA plots (one per comparison)"

  volcano_plot:
    type: File[]
    outputBinding:
      glob: "*volcano_plot*.pdf"
    doc: "Volcano plots (one per comparison)"
  
  pca_plot:
    type: File
    outputBinding:
      glob: "*PCA_plot*.pdf"
    doc: "PCA plot"

  replicate_scatter:
    type: File
    outputBinding:
      glob: "*replicate_scatter*.pdf"
    doc: "Replicate scatter plots"
  
  mean_scatter:
    type: File
    outputBinding:
      glob: "*mean_scatter*.pdf"
    doc: "Mean expression scatter plot"

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
