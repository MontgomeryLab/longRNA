cwlVersion: v1.0
class: Workflow

doc: |
  RNA-seq analysis workflow
  
  This workflow performs:
  1. Quality control and trimming with fastp (single or paired-end)
  2. Transcriptome indexing with Salmon (if needed)
  3. Quantification with Salmon  
  4. Differential expression analysis with DESeq2

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
  transcriptome:
    type: File
    doc: "Transcriptome FASTA file"

  gtf_annotation:
    type: File
    doc: "GTF annotation file"
  
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

  salmon_threads:
    type: int?
  
  library_type:
    type: string?
  
  alpha:
    type: float?

  lfc_threshold:
    type: float?
  
  # Output directory names
  dir_name_fastp: string
  dir_name_salmon: string
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
  
  # Step 2: Build Salmon index
  salmon_index:
    run: ../tools/salmon-index.cwl
    in:
      transcriptome: transcriptome
      threads: salmon_threads
    out: [index_dir]
  
  # Step 3: Quantify with Salmon
  salmon_quant:
    run: ../tools/salmon-quant.cwl
    scatter: [cleaned_fastq_r1, cleaned_fastq_r2, sample_name, output_name]
    scatterMethod: dotproduct
    in:
      index_dir: salmon_index/index_dir
      cleaned_fastq_r1: fastp_qc/cleaned_fastq_r1
      cleaned_fastq_r2: fastp_qc/cleaned_fastq_r2
      sample_name: sample_basenames
      output_name: sample_basenames
      library_type: library_type
      threads: salmon_threads
    out: [quant_dir, quant_sf]
  
  # Step 4: Differential expression analysis
  deseq_analysis:
    run: ../tools/deseq-analysis.cwl
    in:
      r_script:
        default:
          class: File
          location: ../../scripts/deseq_analysis_salmon.R
      gtf_annotation: gtf_annotation
      sample_csv_file: sample_csv_file
      quant_dirs: salmon_quant/quant_dir
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
  
  organize_salmon:
    run: ../tools/make-subdir.cwl
    in:
      dir_files: salmon_quant/quant_dir
      dir_name: dir_name_salmon
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
  
  salmon_dir:
    type: Directory
    outputSource: organize_salmon/subdir
    doc: "Directory containing Salmon quantification results"
  
  analysis_dir:
    type: Directory
    outputSource: organize_analysis/subdir
    doc: "Directory containing DESeq2 analysis results"
  
  gene_counts_dir:
    type: Directory
    outputSource: organize_gene_counts/subdir
    doc: "Directory containing gene count matrices (raw and normalized)"
