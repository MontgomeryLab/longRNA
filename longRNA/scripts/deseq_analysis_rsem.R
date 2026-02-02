#!/usr/bin/env Rscript
# DESeq2 Differential Expression Analysis with RSEM Gene-Level Counts
# This script imports RSEM quantifications and performs DE analysis
#
# Save as: scripts/deseq_analysis_rsem.R

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(tximport)
  library(ggplot2)
  library(pheatmap)
  library(readxl)
})

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
gtf_file <- if (length(args) > 0) args[1] else "annotation.gtf"
sample_csv_file <- if (length(args) > 1) args[2] else NULL
geneid2name_file <- if (length(args) > 2) args[3] else NULL

cat("=======================================================\n")
cat("RNA-seq Gene-Level Analysis Pipeline (RSEM)\n")
cat("=======================================================\n\n")

# Read sample information
cat("Reading sample information...\n")

# Read basic sample_info.csv
sample_info <- read.csv("sample_info.csv", stringsAsFactors = FALSE)

# If sample_csv_file is provided, read it for detailed metadata
experimental_design <- NULL
if (!is.null(sample_csv_file) && file.exists(sample_csv_file)) {
  cat(sprintf("Reading detailed sample metadata from: %s\n", sample_csv_file))
  sample_metadata <- read.csv(sample_csv_file, stringsAsFactors = FALSE)

  # Extract sample names from forward reads filenames (remove _1.fq.gz or _R1.fq.gz suffix)
  sample_metadata$sample_name <- sub("_R?1\\.fq\\.gz$", "", sample_metadata$Forward.Reads)

  # Update sample_info with condition information from metadata
  sample_info$condition <- sapply(sample_info$sample, function(s) {
    idx <- which(sample_metadata$sample_name == s)
    if (length(idx) > 0) {
      return(sample_metadata$Sample.Group.Name[idx[1]])
    } else {
      return(NA)
    }
  })

  # Detect experimental design
  if ("Control" %in% colnames(sample_metadata)) {
    control_rows <- sample_metadata[!is.na(sample_metadata$Control) &
                                    sample_metadata$Control != "", ]
    if (nrow(control_rows) > 0) {
      control_condition <- unique(control_rows$Sample.Group.Name)
      if (length(control_condition) == 1) {
        experimental_conditions <- unique(sample_metadata$Sample.Group.Name)
        experimental_conditions <- experimental_conditions[experimental_conditions != control_condition]

        experimental_design <- list(
          control = control_condition,
          experimental = experimental_conditions
        )

        cat(sprintf("\n✓ Detected experimental design:\n"))
        cat(sprintf("  Control: %s\n", control_condition))
        cat(sprintf("  Experimental conditions: %s\n", paste(experimental_conditions, collapse=", ")))
        cat(sprintf("  Total comparisons to perform: %d\n\n", length(experimental_conditions)))
      }
    }
  }
}

# Read RSEM .genes.results file paths
genes_files <- readLines("genes_files.txt")
names(genes_files) <- sample_info$sample

# Check that all files exist
if (!all(file.exists(genes_files))) {
  stop("Some RSEM genes.results files not found!")
}

cat("Samples found:\n")
print(sample_info)

# Import RSEM quantifications at GENE level
cat("\nImporting RSEM quantifications (gene-level)...\n")
cat("This may take a few minutes...\n")

txi.rsem <- tximport(genes_files,
                     type = "rsem",
                     txIn = FALSE,    # Gene-level import
                     txOut = FALSE)

cat(sprintf("Imported %d genes\n", nrow(txi.rsem$counts)))

# Handle genes with zero effective length to ensure DESeq2 compatibility
# DESeq2 requires all lengths > 0, so we replace zeros with small positive values
genes_with_zero_length <- apply(txi.rsem$length, 1, function(x) any(x <= 0))
if (sum(genes_with_zero_length) > 0) {
  cat(sprintf("Found %d genes with zero effective length in at least one sample\n", sum(genes_with_zero_length)))
  cat("Replacing zero lengths with minimum positive length for each gene\n")

  # For each gene, replace zero/negative lengths with the minimum positive length for that gene
  # If all lengths are zero, use 1 as default
  for (i in which(genes_with_zero_length)) {
    gene_lengths <- txi.rsem$length[i, ]
    positive_lengths <- gene_lengths[gene_lengths > 0]
    replacement_value <- if (length(positive_lengths) > 0) min(positive_lengths) else 1
    txi.rsem$length[i, gene_lengths <= 0] <- replacement_value
  }

  cat(sprintf("All %d genes retained with adjusted lengths\n", nrow(txi.rsem$counts)))
}

# Read GTF for gene ID to gene name mapping
cat("\nCreating gene ID to gene name mapping from GTF...\n")
cat(sprintf("GTF file: %s\n", gtf_file))

if (!file.exists(gtf_file)) {
  stop(sprintf("GTF file not found: %s", gtf_file))
}

# Read GTF file
if (grepl("\\.gz$", gtf_file)) {
  gtf_data <- read.table(gzfile(gtf_file), sep = "\t", comment.char = "#", 
                         stringsAsFactors = FALSE, quote = "")
} else {
  gtf_data <- read.table(gtf_file, sep = "\t", comment.char = "#", 
                         stringsAsFactors = FALSE, quote = "")
}

colnames(gtf_data) <- c("chr", "source", "feature", "start", "end", 
                        "score", "strand", "frame", "attributes")

# Function to extract attribute values from GTF
extract_gtf_attr <- function(attr_string, attr_name) {
  pattern <- paste0(attr_name, ' "([^"]+)"')
  m <- regexpr(pattern, attr_string, perl = TRUE)
  if (m > 0) {
    matched_str <- regmatches(attr_string, m)
    value_pattern <- '"([^"]+)"'
    value_match <- regexpr(value_pattern, matched_str, perl = TRUE)
    value_str <- regmatches(matched_str, value_match)
    return(gsub('"', '', value_str))
  } else {
    return(NA_character_)
  }
}

# Extract gene_id to gene_name mapping from gene features
genes <- gtf_data[gtf_data$feature == "gene", ]
gene_ids_all <- sapply(genes$attributes, function(x) extract_gtf_attr(x, "gene_id"))
gene_names <- sapply(genes$attributes, function(x) extract_gtf_attr(x, "gene_name"))

# Remove version numbers from gene IDs
gene_ids_all_clean <- sub("\\.\\d+$", "", gene_ids_all)

# Create gene_id to gene_name mapping
gene_id_to_name <- data.frame(
  gene_id = gene_ids_all_clean,
  gene_name = gene_names,
  stringsAsFactors = FALSE
)

# Remove rows with NA
gene_id_to_name <- gene_id_to_name[!is.na(gene_id_to_name$gene_id) & 
                                   !is.na(gene_id_to_name$gene_name), ]

# Remove duplicates
gene_id_to_name <- gene_id_to_name[!duplicated(gene_id_to_name$gene_id), ]

cat(sprintf("Created mapping for %d genes to gene symbols\n\n", nrow(gene_id_to_name)))

# Convert Ensembl gene IDs to gene symbols
cat("Converting gene IDs to gene symbols...\n")
gene_ids_in_counts <- rownames(txi.rsem$counts)
gene_ids_clean <- sub("\\.\\d+$", "", gene_ids_in_counts)

# Match gene IDs to gene names
gene_names_mapped <- gene_id_to_name$gene_name[match(gene_ids_clean, gene_id_to_name$gene_id)]

# For genes without a name, keep the gene ID
gene_names_mapped[is.na(gene_names_mapped)] <- gene_ids_clean[is.na(gene_names_mapped)]

# Handle duplicate gene names by appending gene ID
if (any(duplicated(gene_names_mapped))) {
  dup_names <- unique(gene_names_mapped[duplicated(gene_names_mapped)])
  cat(sprintf("Found %d duplicate gene names, adding gene IDs to make unique\n", length(dup_names)))
  for (dup in dup_names) {
    idx <- which(gene_names_mapped == dup)
    gene_names_mapped[idx] <- paste0(gene_names_mapped[idx], "_", gene_ids_clean[idx])
  }
}

# Update rownames in tximport object
rownames(txi.rsem$counts) <- gene_names_mapped
rownames(txi.rsem$abundance) <- gene_names_mapped
rownames(txi.rsem$length) <- gene_names_mapped

cat(sprintf("Successfully converted to gene symbols\n"))
cat(sprintf("Genes with symbols: %d\n", sum(gene_names_mapped %in% gene_id_to_name$gene_name)))
cat(sprintf("Genes keeping gene IDs: %d\n", sum(!gene_names_mapped %in% gene_id_to_name$gene_name)))

# Save RAW gene counts as TSV
cat("\nSaving raw gene counts...\n")
gene_counts <- as.data.frame(round(txi.rsem$counts))
gene_counts <- cbind(gene = rownames(gene_counts), gene_counts)
write.table(gene_counts, "gene_counts_raw.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  - gene_counts_raw.tsv (%d genes, raw integer counts)\n", nrow(gene_counts) - 1))

# Prepare sample data
sample_data <- data.frame(
  condition = factor(sample_info$condition),
  row.names = sample_info$sample
)

# Create DESeq2 dataset
cat("\nCreating DESeq2 dataset...\n")
dds <- DESeqDataSetFromTximport(txi.rsem, colData = sample_data, 
                                design = ~ condition)

# Filter low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

cat(sprintf("Retained %d genes after filtering (rowSums >= 10)\n", sum(keep)))

# Run DESeq2 analysis
cat("\nRunning DESeq2 analysis...\n")
dds <- DESeq(dds)

# Perform pairwise comparisons
all_comparisons <- list()

if (!is.null(experimental_design)) {
  # Multi-condition experimental design: compare each condition vs control
  cat(sprintf("\nPerforming pairwise comparisons (%s as control):\n", experimental_design$control))

  for (exp_cond in experimental_design$experimental) {
    cat(sprintf("  - %s vs %s\n", exp_cond, experimental_design$control))

    res <- results(dds,
                   contrast = c("condition", exp_cond, experimental_design$control),
                   alpha = 0.05)
    res_ordered <- res[order(res$padj),]

    # Store results
    all_comparisons[[exp_cond]] <- list(
      name = exp_cond,
      results = res,
      results_ordered = res_ordered
    )

    # Save results
    output_file <- sprintf("deseq2_rsem_results_%s_vs_%s.csv", exp_cond, experimental_design$control)
    write.csv(as.data.frame(res_ordered), output_file)
    cat(sprintf("    Saved: %s\n", output_file))

    # Summary
    cat(sprintf("    Summary for %s vs %s:\n", exp_cond, experimental_design$control))
    print(summary(res))
    cat("\n")
  }

  # Use first comparison as default for backward compatibility
  res <- all_comparisons[[1]]$results
  res_ordered <- all_comparisons[[1]]$results_ordered

} else {
  # Simple two-condition comparison
  cat("\nSimple two-condition comparison\n")
  res <- results(dds, alpha = 0.05)
  res_ordered <- res[order(res$padj),]

  # Summary
  cat("\nDESeq2 Results Summary:\n")
  summary(res)

  # Save results
  cat("\nSaving results...\n")
  write.csv(as.data.frame(res_ordered), "deseq2_rsem_results.csv")
}

# Get normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, "deseq2_rsem_normalized_counts.csv")

# Save normalized counts as TSV with ALL genes
cat("\nCalculating normalized counts for all genes...\n")

# Create a full DESeq2 dataset with ALL genes for normalization
dds_full <- DESeqDataSetFromTximport(txi.rsem, colData = sample_data, design = ~ condition)
dds_full <- estimateSizeFactors(dds_full)

# Get normalized counts for all genes
normalized_counts_full <- counts(dds_full, normalized = TRUE)
normalized_counts_tsv <- as.data.frame(round(normalized_counts_full, 3))
normalized_counts_tsv <- cbind(gene = rownames(normalized_counts_tsv), normalized_counts_tsv)
write.table(normalized_counts_tsv, "gene_counts_normalized.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  - gene_counts_normalized.tsv (%d genes, 3 decimal places)\n",
            nrow(normalized_counts_tsv) - 1))

# Create gene name versions if geneid2name_file is provided
if (!is.null(geneid2name_file) && file.exists(geneid2name_file)) {
  cat(sprintf("\nReading gene ID to name mapping from: %s\n", geneid2name_file))

  # Read the xlsx file
  geneid2name <- read_excel(geneid2name_file)
  colnames(geneid2name) <- c("gene_id", "gene_name")

  cat(sprintf("Loaded %d gene ID to name mappings\n", nrow(geneid2name)))

  # Function to convert gene IDs to gene names
  convert_ids_to_names <- function(gene_ids, mapping) {
    gene_names <- sapply(gene_ids, function(gid) {
      idx <- which(mapping$gene_id == gid)
      if (length(idx) > 0) {
        return(mapping$gene_name[idx[1]])
      } else {
        return(gid)  # Keep original ID if no mapping found
      }
    })
    return(as.character(gene_names))
  }

  # Create gene name versions
  cat("Creating gene name count files...\n")

  # Raw counts with gene names
  genename_counts_raw <- gene_counts
  genename_counts_raw$gene <- convert_ids_to_names(gene_counts$gene, geneid2name)
  write.table(genename_counts_raw, "genename_counts_raw.tsv",
              sep = "\t", quote = FALSE, row.names = FALSE)

  # Normalized counts with gene names
  genename_counts_normalized <- normalized_counts_tsv
  genename_counts_normalized$gene <- convert_ids_to_names(normalized_counts_tsv$gene, geneid2name)
  write.table(genename_counts_normalized, "genename_counts_normalized.tsv",
              sep = "\t", quote = FALSE, row.names = FALSE)

  # Count how many were converted
  converted_count <- sum(genename_counts_raw$gene != gene_counts$gene)
  cat(sprintf("  - genename_counts_raw.tsv (%d genes converted to names)\n", converted_count))
  cat(sprintf("  - genename_counts_normalized.tsv (%d genes converted to names)\n", converted_count))
} else {
  if (!is.null(geneid2name_file)) {
    cat(sprintf("\nWarning: geneid2name file not found: %s\n", geneid2name_file))
  }
}

# Convert gene IDs to transcript IDs using tx2gene mapping
if (file.exists("tx2gene_mapping.tsv")) {
  cat("\nConverting gene IDs to transcript IDs...\n")

  # Read tx2gene mapping
  tx2gene <- read.table("tx2gene_mapping.tsv", header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE)

  # Create gene to transcript mapping (use first transcript for each gene)
  gene2tx <- tx2gene[!duplicated(tx2gene$GENEID), ]
  rownames(gene2tx) <- gene2tx$GENEID

  # Function to convert gene IDs to transcript IDs
  convert_to_transcript_ids <- function(gene_ids) {
    transcript_ids <- sapply(gene_ids, function(gid) {
      if (gid %in% gene2tx$GENEID) {
        return(gene2tx[gid, "TXNAME"])
      } else {
        return(gid)  # Keep gene ID if no transcript found
      }
    })
    return(as.character(transcript_ids))
  }

  # Convert raw counts
  gene_counts_tx <- gene_counts
  gene_counts_tx$gene <- convert_to_transcript_ids(gene_counts$gene)
  write.csv(gene_counts_tx, "gene_counts_raw_names.csv", row.names = FALSE)

  # Convert normalized counts
  normalized_counts_tx <- normalized_counts_tsv
  normalized_counts_tx$gene <- convert_to_transcript_ids(normalized_counts_tsv$gene)
  write.csv(normalized_counts_tx, "gene_counts_normalized_names.csv", row.names = FALSE)

  cat(sprintf("  - gene_counts_raw_names.csv (%d genes with transcript IDs)\n",
              nrow(gene_counts_tx) - 1))
  cat(sprintf("  - gene_counts_normalized_names.csv (%d genes with transcript IDs)\n",
              nrow(normalized_counts_tx) - 1))
  cat(sprintf("  Genes converted to transcript IDs: %d\n",
              sum(gene_counts_tx$gene != gene_counts$gene)))
}

cat("Results saved to:\n")
cat("  - deseq2_rsem_results.csv (differential expression results)\n")
cat("  - deseq2_rsem_normalized_counts.csv (normalized counts for filtered genes)\n")
cat("  - gene_counts_normalized.tsv (normalized counts for ALL genes, 3 decimals)\n")

# Generate plots
cat("\nGenerating plots...\n")

# If we have multiple comparisons, generate plots for each
if (!is.null(experimental_design) && length(all_comparisons) > 0) {
  cat("Generating condition-specific plots...\n")

  for (comp_name in names(all_comparisons)) {
    comp <- all_comparisons[[comp_name]]
    res <- comp$results
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", comp_name)  # Replace special chars

    cat(sprintf("\n  Creating plots for: %s vs %s\n", comp_name, experimental_design$control))

    # 1. MA plot with LFC shrinkage
    pdf_name <- sprintf("deseq2_rsem_MA_plot_%s_vs_%s.pdf", safe_name, experimental_design$control)
    pdf(pdf_name, width = 10, height = 6)
    par(mfrow=c(1,2))
    plotMA(res, main = sprintf("MA Plot (Original)\n%s vs %s", comp_name, experimental_design$control),
           ylim = c(-5, 5))
    # Find the coefficient name for this contrast
    coef_idx <- which(grepl(paste0("condition.*", gsub("[^A-Za-z0-9]", ".", comp_name)), resultsNames(dds)))
    if (length(coef_idx) > 0) {
      resLFC <- lfcShrink(dds, coef=coef_idx[1], type="normal", res=res)
      plotMA(resLFC, main = sprintf("MA Plot (LFC Shrinkage)\n%s vs %s", comp_name, experimental_design$control),
             ylim = c(-5, 5))
    } else {
      plot.new()
      text(0.5, 0.5, "LFC shrinkage not available", cex = 1.5)
    }
    dev.off()

    # 2. Volcano plot
    pdf_name <- sprintf("deseq2_rsem_volcano_plot_%s_vs_%s.pdf", safe_name, experimental_design$control)
    pdf(pdf_name, width = 8, height = 6)
    res_df <- as.data.frame(res)
    res_df$significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1

    print(ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
      geom_point(alpha = 0.5, size = 1) +
      scale_color_manual(values = c("gray", "red")) +
      theme_bw() +
      labs(title = sprintf("Volcano Plot: %s vs %s (RSEM)", comp_name, experimental_design$control),
           x = "Log2 Fold Change",
           y = "-Log10 Adjusted P-value") +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed"))
    dev.off()
  }

} else {
  # Original single comparison plots
  # 1. MA plot with LFC shrinkage
  pdf("deseq2_rsem_MA_plot.pdf", width = 10, height = 6)
  par(mfrow=c(1,2))
  plotMA(res, main = "MA Plot (Original)", ylim = c(-5, 5))
  # LFC shrinkage for better visualization
  resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="normal", res=res)
  plotMA(resLFC, main = "MA Plot (LFC Shrinkage)", ylim = c(-5, 5))
  dev.off()

  # 2. Volcano plot
  pdf("deseq2_rsem_volcano_plot.pdf", width = 8, height = 6)
  res_df <- as.data.frame(res)
  res_df$significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1

  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c("gray", "red")) +
    theme_bw() +
    labs(title = "Volcano Plot (RSEM)",
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-value") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed")
  dev.off()
}

# Create res_df for any subsequent use (use first comparison if experimental_design exists)
if (!is.null(experimental_design) && length(all_comparisons) > 0) {
  # Use the first experimental comparison
  first_comp_name <- names(all_comparisons)[1]
  first_comp <- all_comparisons[[first_comp_name]]
  res <- first_comp$results
  res_df <- as.data.frame(res)
  res_df$significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1
  cat(sprintf("\nUsing %s vs %s comparison for any subsequent plots\n",
              first_comp_name, experimental_design$control))
} else if (!exists("res_df")) {
  # Fallback: ensure res_df exists
  cat("\nCreating res_df from default results\n")
  res_df <- as.data.frame(res)
  res_df$significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1
}

# Verify res_df exists
if (!exists("res_df")) {
  stop("ERROR: res_df not defined!")
}
cat(sprintf("res_df contains %d genes\n", nrow(res_df)))

# 3. PCA plot
pdf("deseq2_rsem_PCA_plot.pdf", width = 8, height = 6)
# Use varianceStabilizingTransformation for small datasets
vsd <- tryCatch({
  vst(dds, blind = FALSE)
}, error = function(e) {
  if (grepl("less than 'nsub' rows", e$message)) {
    cat("Using varianceStabilizingTransformation() for small dataset\n")
    varianceStabilizingTransformation(dds, blind = FALSE)
  } else {
    stop(e)
  }
})
plotPCA(vsd, intgroup = "condition") +
  theme_bw() +
  labs(title = "PCA Plot of Samples (RSEM)")
dev.off()

# 4. Replicate scatter plots
cat("Generating replicate scatter plots...\n")
pdf("deseq2_rsem_replicate_scatter.pdf", width = 12, height = 10)

norm_counts <- counts(dds, normalized=TRUE)
conditions <- unique(sample_info$condition)

# Calculate total number of plots needed (one plot per condition with >=2 replicates)
total_plots <- 0
for (cond in conditions) {
  idx <- which(sample_info$condition == cond)
  if (length(idx) >= 2) {
    total_plots <- total_plots + 1
  }
}

# Set up plot layout based on number of plots
if (total_plots <= 4) {
  par(mfrow=c(2,2), mar = c(4, 4, 2, 1), mgp = c(2, 0.5, 0), pty="s")
} else if (total_plots <= 6) {
  par(mfrow=c(2,3), mar = c(4, 4, 2, 1), mgp = c(2, 0.5, 0), pty="s")
} else if (total_plots <= 9) {
  par(mfrow=c(3,3), mar = c(4, 4, 2, 1), mgp = c(2, 0.5, 0), pty="s")
} else {
  par(mfrow=c(3,4), mar = c(3, 3, 2, 1), mgp = c(1.5, 0.5, 0), pty="s")
}

# Plot replicates for each condition
for (cond in conditions) {
  idx <- which(sample_info$condition == cond)

  if (length(idx) >= 2) {
    # For conditions with exactly 2 replicates, plot them directly
    if (length(idx) == 2) {
      plot(log2(1 + norm_counts[, idx[1]]),
           log2(1 + norm_counts[, idx[2]]),
           col="#00000020", pch=19, cex=0.3,
           xlab=colnames(norm_counts)[idx[1]],
           ylab=colnames(norm_counts)[idx[2]],
           main=paste(cond, "replicates"))
      abline(0, 1, col="red", lwd=2)
    } else {
      # For conditions with >2 replicates, plot all pairwise combinations
      for (i in 1:(length(idx)-1)) {
        for (j in (i+1):length(idx)) {
          plot(log2(1 + norm_counts[, idx[i]]),
               log2(1 + norm_counts[, idx[j]]),
               col="#00000020", pch=19, cex=0.3,
               xlab=colnames(norm_counts)[idx[i]],
               ylab=colnames(norm_counts)[idx[j]],
               main=paste(cond, "replicates"))
          abline(0, 1, col="red", lwd=2)
        }
      }
    }
  }
}

dev.off()

# 6. Mean scatter plots
cat("Generating mean scatter plots...\n")

# Helper function for log2 tick marks
plot_log2_ticks <- function(min, max) {
  tick_sep <- seq(from=min, to=max, by=4)
  tick_sep <- 2^tick_sep

  for (i in 1:length(tick_sep)) {
    axis(1, at=log2(tick_sep[i]),
         labels=prettyNum(tick_sep[i], digits=2, format="f"),
         las=1, cex.axis=1.0)
    axis(2, at=log2(tick_sep[i]),
         labels=prettyNum(tick_sep[i], digits=2, format="f"),
         las=1, cex.axis=1.0)
  }

  for (i in 2:length(tick_sep)) {
    minor_sep <- (tick_sep[i] - tick_sep[i-1]) / 8
    minor_ticks <- seq(from=tick_sep[i-1], to=tick_sep[i], by=minor_sep)
    axis(1, at=log2(minor_ticks[1:8]), labels=FALSE, tck=-0.005)
    axis(2, at=log2(minor_ticks[1:8]), labels=FALSE, tck=-0.005)
    axis(1, at=log2(minor_ticks[5]), tck=-0.009, labels=FALSE)
    axis(2, at=log2(minor_ticks[5]), tck=-0.009, labels=FALSE)
  }
}

cnts <- counts(dds, normalized=TRUE)

if (!is.null(experimental_design) && length(all_comparisons) > 0) {
  # Generate mean scatter plots for each pairwise comparison
  pdf("deseq2_rsem_mean_scatter.pdf", width = 10, height = 10)

  for (exp_cond in names(all_comparisons)) {
    cat(sprintf("  - Creating mean scatter: %s vs %s\n", exp_cond, experimental_design$control))

    # Get indices for this comparison
    idx_control <- which(sample_info$condition == experimental_design$control)
    idx_exp <- which(sample_info$condition == exp_cond)

    # Calculate average expression
    cnts_avg <- data.frame(
      control = log2(rowMeans(cnts[, idx_control])),
      experimental = log2(rowMeans(cnts[, idx_exp]))
    )
    colnames(cnts_avg) <- c(experimental_design$control, exp_cond)
    cnts_avg <- replace(cnts_avg, cnts_avg == -Inf, -4)

    # Get significant genes for this comparison
    res_comp <- all_comparisons[[exp_cond]]$results
    res_comp_df <- as.data.frame(res_comp)
    res_comp_df$significant <- res_comp_df$padj < 0.05 & abs(res_comp_df$log2FoldChange) > 1
    sig_genes_comp <- subset(res_comp_df, significant == TRUE)

    # Create plot
    par(mar=c(5,5,5,2))
    plot(-4:20, -4:20,
         tck = 0, xaxt='n', yaxt='n', type="n",
         main=paste(exp_cond, "vs", experimental_design$control, "(RSEM)"),
         xlab=paste("Average reads in", experimental_design$control),
         ylab=paste("Average reads in", exp_cond),
         font.lab=2, las=1, tcl=-0.3,
         cex.lab=1.5, cex.main=1.8, cex.axis=1.2)

    plot_log2_ticks(-4, 20)

    points(cnts_avg[,1], cnts_avg[,2], pch=20, col="grey80", cex=0.5)

    abline(h=log2(10), v=log2(10), lty=2, col="grey60")
    abline(0, 1, col="grey60")
    abline(log2(0.5), 1, col="grey60")
    abline(log2(2), 1, col="grey60")

    if (nrow(sig_genes_comp) > 0) {
      points(cnts_avg[rownames(sig_genes_comp), 1],
             cnts_avg[rownames(sig_genes_comp), 2],
             pch=19, col="deepskyblue", cex=0.6)
    }

    legend_labels <- c(
      paste0("p < 0.05 & |log2FC| > 1"),
      paste0("p >= 0.05 or |log2FC| <= 1")
    )
    legend('topleft',
           legend = legend_labels,
           fill = c("deepskyblue", "grey80"),
           bty = 'n',
           cex = 1.2)
  }

  dev.off()

} else {
  # Simple two-condition comparison - single mean scatter plot
  pdf("deseq2_rsem_mean_scatter.pdf", width = 10, height = 10)
  par(mar=c(5,5,5,2))

  conditions <- unique(sample_info$condition)
  condition1 <- conditions[1]
  condition2 <- conditions[2]

  idx1 <- which(sample_info$condition == condition1)
  idx2 <- which(sample_info$condition == condition2)

  cnts_avg <- data.frame(
    cond1 = log2(rowMeans(cnts[, idx1])),
    cond2 = log2(rowMeans(cnts[, idx2]))
  )
  colnames(cnts_avg) <- c(condition1, condition2)
  cnts_avg <- replace(cnts_avg, cnts_avg == -Inf, -4)

  res_df <- as.data.frame(res)
  res_df$significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1
  sig_genes <- subset(res_df, significant == TRUE)

  plot(-4:20, -4:20,
       tck = 0, xaxt='n', yaxt='n', type="n",
       main=paste(condition2, "vs", condition1, "(RSEM)"),
       xlab=paste("Average reads in", condition1),
       ylab=paste("Average reads in", condition2),
       font.lab=2, las=1, tcl=-0.3,
       cex.lab=1.5, cex.main=1.8, cex.axis=1.2)

  plot_log2_ticks(-4, 20)

  points(cnts_avg[,1], cnts_avg[,2], pch=20, col="grey80", cex=0.5)

  abline(h=log2(10), v=log2(10), lty=2, col="grey60")
  abline(0, 1, col="grey60")
  abline(log2(0.5), 1, col="grey60")
  abline(log2(2), 1, col="grey60")

  if (nrow(sig_genes) > 0) {
    points(cnts_avg[rownames(sig_genes), 1],
           cnts_avg[rownames(sig_genes), 2],
           pch=19, col="deepskyblue", cex=0.6)
  }

  legend_labels <- c(
    paste0("p < 0.05 & |log2FC| > 1"),
    paste0("p >= 0.05 or |log2FC| <= 1")
  )
  legend('topleft',
         legend = legend_labels,
         fill = c("deepskyblue", "grey80"),
         bty = 'n',
         cex = 1.2)

  dev.off()
}

cat("Plots saved to:\n")
if (!is.null(experimental_design) && length(all_comparisons) > 0) {
  cat("  MA plots and Volcano plots for each comparison:\n")
  for (comp_name in names(all_comparisons)) {
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", comp_name)
    cat(sprintf("    - deseq2_rsem_MA_plot_%s_vs_%s.pdf\n", safe_name, experimental_design$control))
    cat(sprintf("    - deseq2_rsem_volcano_plot_%s_vs_%s.pdf\n", safe_name, experimental_design$control))
  }
} else {
  cat("  - deseq2_rsem_MA_plot.pdf (original + LFC shrinkage)\n")
  cat("  - deseq2_rsem_volcano_plot.pdf\n")
}
cat("  - deseq2_rsem_PCA_plot.pdf\n")
cat("  - deseq2_rsem_replicate_scatter.pdf (all conditions)\n")
cat("  - deseq2_rsem_mean_scatter.pdf (all comparisons)\n")

cat("\n✅ Analysis complete!\n")

# Print summary for all comparisons
if (!is.null(experimental_design) && length(all_comparisons) > 0) {
  cat("\n=======================================================\n")
  cat("Summary of All Comparisons\n")
  cat("=======================================================\n")

  for (comp_name in names(all_comparisons)) {
    res <- all_comparisons[[comp_name]]$results
    res_df <- as.data.frame(res)
    sig_genes <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)

    cat(sprintf("\n--- %s vs %s ---\n", comp_name, experimental_design$control))
    cat(sprintf("Total significant genes (padj < 0.05, |log2FC| > 1): %d\n", nrow(sig_genes)))

    if (nrow(sig_genes) > 0) {
      cat("Top 5 differentially expressed genes:\n")
      top_genes <- head(sig_genes[order(sig_genes$padj),], 5)
      print(top_genes[, c("log2FoldChange", "padj")])
    } else {
      cat("No significantly differentially expressed genes found.\n")
    }
  }
  cat("\n=======================================================\n")
} else {
  # Print top 10 significant genes (original behavior)
  cat("\nTop 10 differentially expressed genes:\n")
  if (exists("sig_genes") && nrow(sig_genes) > 0) {
    top_genes <- head(sig_genes[order(sig_genes$padj),], 10)
    print(top_genes[, c("log2FoldChange", "padj")])
    cat(sprintf("\nTotal significant genes: %d\n", nrow(sig_genes)))
  } else {
    cat("No significantly differentially expressed genes found.\n")
  }
}

cat("\n=======================================================\n")
cat("Gene count matrix info:\n")
cat(sprintf("  Total genes in TSV files: %d\n", nrow(txi.rsem$counts)))
cat(sprintf("  Genes used in DESeq2 analysis: %d (after filtering)\n", nrow(dds)))
cat(sprintf("  Samples: %d\n", ncol(txi.rsem$counts)))
cat("=======================================================\n")
