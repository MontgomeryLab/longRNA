#!/usr/bin/env Rscript
# DESeq2 Differential Expression Analysis with Gene-Level Counts
# Modified for C. elegans WormBase annotations

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
gtf_file <- if (length(args) > 0) args[1] else "annotation.gtf.gz"
sample_csv_file <- if (length(args) > 1) args[2] else NULL
geneid2name_file <- if (length(args) > 2) args[3] else NULL

cat("=======================================================\n")
cat("RNA-seq Gene-Level Analysis Pipeline\n")
cat("=======================================================\n\n")

# Read sample information from CSV
cat("Reading sample information...\n")

# First read the basic sample_info.csv for sample names and quant files
sample_info <- read.csv("sample_info.csv", stringsAsFactors = FALSE)

# If sample_csv_file is provided, read it to get detailed experimental design
if (!is.null(sample_csv_file) && file.exists(sample_csv_file)) {
  cat(sprintf("Reading detailed sample metadata from: %s\n", sample_csv_file))
  sample_metadata <- read.csv(sample_csv_file, stringsAsFactors = FALSE)

  # Extract sample names from forward reads filenames (remove _1.fq.gz or _R1.fq.gz suffix)
  sample_metadata$sample_name <- sub("_R?1\\.fq\\.gz$", "", sample_metadata$Forward.Reads)

  # Identify control condition
  control_samples <- sample_metadata[sample_metadata$Control == "TRUE", ]
  if (nrow(control_samples) > 0) {
    control_condition <- unique(control_samples$Sample.Group.Name)[1]
    cat(sprintf("Control condition identified: %s\n", control_condition))
  } else {
    stop("No control samples found in CSV file! Please mark control samples with Control=TRUE")
  }

  # Get all experimental conditions (non-control)
  experimental_conditions <- setdiff(unique(sample_metadata$Sample.Group.Name), control_condition)
  cat(sprintf("Experimental conditions: %s\n", paste(experimental_conditions, collapse=", ")))

  # Update sample_info with condition information from metadata
  # Match by sample name
  sample_info$condition <- sapply(sample_info$sample, function(s) {
    idx <- which(sample_metadata$sample_name == s)
    if (length(idx) > 0) {
      return(sample_metadata$Sample.Group.Name[idx[1]])
    } else {
      return(NA)
    }
  })

  # Store metadata for later use
  experimental_design <- list(
    control = control_condition,
    experimental = experimental_conditions,
    metadata = sample_metadata
  )
} else {
  cat("No sample CSV file provided, using condition info from sample_info.csv\n")
  experimental_design <- NULL
}

# Read quant file paths
quant_files <- readLines("quant_files.txt")
names(quant_files) <- sample_info$sample

# Check that all files exist
if (!all(file.exists(quant_files))) {
  stop("Some quantification files not found!")
}

cat("Samples found:\n")
print(sample_info)

# Create tx2gene mapping from GTF
cat("\nCreating transcript-to-gene mapping from GTF...\n")
cat(sprintf("GTF file: %s\n", gtf_file))

# Check if GTF file exists
if (!file.exists(gtf_file)) {
  stop(sprintf("GTF file not found: %s", gtf_file))
}

# Read GTF file to extract transcript-to-gene and gene_id-to-gene_name mappings
cat("Parsing GTF file for transcript and gene name mappings...\n")

# Read GTF (handle gzipped files)
if (grepl("\\.gz$", gtf_file)) {
  gtf_data <- read.table(gzfile(gtf_file), sep = "\t", comment.char = "#", 
                         stringsAsFactors = FALSE, quote = "")
} else {
  gtf_data <- read.table(gtf_file, sep = "\t", comment.char = "#", 
                         stringsAsFactors = FALSE, quote = "")
}

colnames(gtf_data) <- c("chr", "source", "feature", "start", "end", 
                        "score", "strand", "frame", "attributes")

# Function to extract attribute values from GTF (handles multiple formats)
extract_gtf_attr <- function(attr_string, attr_name) {
  # Try quoted format first: attr_name "value"
  pattern <- paste0(attr_name, ' "([^"]+)"')
  m <- regexpr(pattern, attr_string, perl = TRUE)
  if (m > 0) {
    matched_str <- regmatches(attr_string, m)
    value_pattern <- '"([^"]+)"'
    value_match <- regexpr(value_pattern, matched_str, perl = TRUE)
    value_str <- regmatches(matched_str, value_match)
    return(gsub('"', '', value_str))
  }
  
  # Try unquoted format: attr_name=value or attr_name value
  pattern2 <- paste0(attr_name, '[= ]([^;]+)')
  m2 <- regexpr(pattern2, attr_string, perl = TRUE)
  if (m2 > 0) {
    matched_str <- regmatches(attr_string, m2)
    value <- sub(paste0(attr_name, '[= ]'), '', matched_str)
    value <- gsub('["; ]', '', value)
    return(value)
  }
  
  return(NA_character_)
}

# Extract transcript-to-gene mapping from transcript features
transcripts <- gtf_data[gtf_data$feature == "transcript", ]

# Extract transcript_id and gene_id
transcript_ids <- sapply(transcripts$attributes, function(x) extract_gtf_attr(x, "transcript_id"))
gene_ids <- sapply(transcripts$attributes, function(x) extract_gtf_attr(x, "gene_id"))

# DO NOT strip versions manually - let tximport handle it with ignoreTxVersion=TRUE
# This prevents issues with C. elegans IDs like Y38C1AB.9 being incorrectly stripped to Y38C1AB

# Create tx2gene mapping (must have columns TXNAME and GENEID)
tx2gene <- data.frame(
  TXNAME = transcript_ids,
  GENEID = gene_ids,
  stringsAsFactors = FALSE
)

# Remove rows with NA
tx2gene <- tx2gene[!is.na(tx2gene$TXNAME) & !is.na(tx2gene$GENEID), ]

# Remove duplicates (keep first occurrence for each transcript)
tx2gene <- tx2gene[!duplicated(tx2gene$TXNAME), ]

cat(sprintf("Mapped %d unique transcripts to %d genes\n",
            nrow(tx2gene), length(unique(tx2gene$GENEID))))

# Show examples for debugging
cat("Example transcript IDs from tx2gene (should match Salmon output):\n")
cat(paste(head(tx2gene$TXNAME, 5), collapse = ", "), "\n")
cat("Example gene IDs from tx2gene:\n")
cat(paste(head(tx2gene$GENEID, 5), collapse = ", "), "\n\n")

# Save tx2gene mapping for debugging
write.table(tx2gene, "tx2gene_mapping.tsv", sep="\t", quote=FALSE, row.names=FALSE)
cat("Saved tx2gene mapping to tx2gene_mapping.tsv for inspection\n")

# Import Salmon quantifications at GENE level
cat("Importing Salmon quantifications (gene-level)...\n")
cat("This may take a few minutes...\n")

# Wrap in tryCatch for better error messages
txi <- tryCatch({
  tximport(quant_files,
           type = "salmon",
           tx2gene = tx2gene,
           ignoreTxVersion = FALSE,
           ignoreAfterBar = FALSE)
}, error = function(e) {
  cat("\nError during tximport:\n")
  cat(conditionMessage(e), "\n")
  
  # Debug information
  cat("\nDebug information:\n")
  cat("Number of quant files:", length(quant_files), "\n")
  cat("Number of tx2gene rows:", nrow(tx2gene), "\n")
  cat("Number of unique transcripts in tx2gene:", length(unique(tx2gene$TXNAME)), "\n")
  cat("Number of unique genes in tx2gene:", length(unique(tx2gene$GENEID)), "\n")
  
  # Check first quant file
  if (length(quant_files) > 0 && file.exists(quant_files[1])) {
    cat("\nChecking first quant file:", quant_files[1], "\n")
    quant_sample <- read.table(quant_files[1], header = TRUE, sep = "\t", nrows = 5)
    cat("First few transcript Names from Salmon:\n")
    print(quant_sample$Name)
    
    # Extract transcript IDs from Salmon (before first |)
    quant_tx_full <- read.table(quant_files[1], header = TRUE, sep = "\t")$Name
    quant_tx_first <- sub("\\|.*", "", as.character(quant_tx_full))
    quant_tx_clean <- sub("\\.\\d+$", "", quant_tx_first)
    
    cat("\nFirst 5 Salmon transcript IDs (cleaned):\n")
    cat(paste(head(quant_tx_clean, 5), collapse = ", "), "\n")
    
    cat("\nFirst 5 tx2gene transcript IDs:\n")
    cat(paste(head(tx2gene$TXNAME, 5), collapse = ", "), "\n")
    
    n_match <- sum(quant_tx_clean %in% tx2gene$TXNAME)
    cat(sprintf("\nTranscripts in quant file matching tx2gene: %d out of %d\n", 
                n_match, length(quant_tx_clean)))
  }
  
  stop("tximport failed. See debug information above.")
})

cat(sprintf("Imported %d genes\n", nrow(txi$counts)))

# DEBUG: Check what gene IDs we actually got
cat("\nDEBUG: First 10 gene IDs from tximport:\n")
cat(paste(head(rownames(txi$counts), 10), collapse = ", "), "\n")
cat("DEBUG: Last 10 gene IDs from tximport:\n")
cat(paste(tail(rownames(txi$counts), 10), collapse = ", "), "\n\n")

# SKIP gene symbol conversion - keep WBGene IDs for accurate counts
cat("Keeping gene IDs as-is (not converting to symbols)...\n")
cat(sprintf("Total genes: %d\n", nrow(txi$counts)))

# Save RAW gene counts as TSV (matching your example format)
cat("\nSaving raw gene counts...\n")
gene_counts <- as.data.frame(txi$counts)
# Round to integers (raw counts should be whole numbers)
gene_counts <- round(gene_counts)
gene_counts <- cbind(gene = rownames(gene_counts), gene_counts)
write.table(gene_counts, "gene_counts_raw.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  - gene_counts_raw.tsv (%d genes, raw integer counts)\n", nrow(gene_counts) - 1))
cat(sprintf("  Sample order: %s, %s, ... %s\n", 
            colnames(gene_counts)[2], colnames(gene_counts)[3], tail(colnames(gene_counts), 1)))

# Prepare sample data
sample_data <- data.frame(
  condition = factor(sample_info$condition),
  row.names = sample_info$sample
)

# Create DESeq2 dataset with ALL samples
cat("\nCreating DESeq2 dataset...\n")
dds_full <- DESeqDataSetFromTximport(txi, colData = sample_data,
                                      design = ~ condition)

# Filter low count genes
keep <- rowSums(counts(dds_full)) >= 10
dds_full <- dds_full[keep,]

cat(sprintf("Retained %d genes after filtering (rowSums >= 10)\n", sum(keep)))

# If experimental design is available, perform pairwise comparisons
if (!is.null(experimental_design)) {
  cat("\n=======================================================\n")
  cat("Performing pairwise comparisons: Each condition vs Control\n")
  cat("=======================================================\n")

  # Run DESeq2 once on the full dataset
  cat("\nRunning DESeq2 on full dataset...\n")
  dds_full <- DESeq(dds_full)

  # For each experimental condition, extract results vs control
  all_comparisons <- list()

  for (exp_cond in experimental_design$experimental) {
    cat(sprintf("\n--- Analyzing: %s vs %s ---\n", exp_cond, experimental_design$control))

    # Get results for this specific comparison
    res <- results(dds_full,
                   contrast = c("condition", exp_cond, experimental_design$control),
                   alpha = 0.05)
    res_ordered <- res[order(res$padj),]

    # Summary
    cat(sprintf("\n%s vs %s Results Summary:\n", exp_cond, experimental_design$control))
    summary(res)

    # Save results with condition-specific filename
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", exp_cond)  # Replace special chars
    results_file <- sprintf("deseq2_results_%s_vs_%s.csv",
                            safe_name, experimental_design$control)
    cat(sprintf("Saving results to: %s\n", results_file))
    write.csv(as.data.frame(res_ordered), results_file)

    # Store for later plotting
    all_comparisons[[exp_cond]] <- list(
      results = res,
      name = exp_cond,
      safe_name = safe_name
    )
  }

  # Also save the full DESeq2 object for reference
  dds <- dds_full

} else {
  # Fallback to original behavior if no experimental design
  cat("\nRunning DESeq2 analysis (all conditions)...\n")
  dds <- DESeq(dds_full)

  # Get results
  res <- results(dds, alpha = 0.05)
  res_ordered <- res[order(res$padj),]

  # Summary
  cat("\nDESeq2 Results Summary:\n")
  summary(res)

  # Save results
  cat("\nSaving results...\n")
  write.csv(as.data.frame(res_ordered), "deseq2_results.csv")
}

# Get normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, "deseq2_normalized_counts.csv")

# Save normalized counts as TSV with ALL genes (not just filtered)
cat("\nCalculating normalized counts for all genes...\n")

# Create a full DESeq2 dataset with ALL genes for normalization
dds_full <- DESeqDataSetFromTximport(txi, colData = sample_data, design = ~ condition)

# Estimate size factors for the full dataset
dds_full <- estimateSizeFactors(dds_full)

# Get normalized counts for all genes
normalized_counts_full <- counts(dds_full, normalized = TRUE)

# Format to 3 decimal places
normalized_counts_tsv <- as.data.frame(round(normalized_counts_full, 3))
normalized_counts_tsv <- cbind(gene = rownames(normalized_counts_tsv), normalized_counts_tsv)
write.table(normalized_counts_tsv, "gene_counts_normalized.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  - gene_counts_normalized.tsv (%d genes, 3 decimal places)\n", nrow(normalized_counts_tsv) - 1))

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

# Convert gene IDs to transcript IDs if tx2gene_mapping.tsv exists
if (file.exists("tx2gene_mapping.tsv")) {
  cat("\nConverting gene IDs to transcript IDs...\n")

  # Read tx2gene mapping
  tx2gene <- read.table("tx2gene_mapping.tsv", header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE)

  # Create gene-to-transcript mapping (use first transcript for each gene)
  gene2tx <- tx2gene[!duplicated(tx2gene$GENEID), ]
  rownames(gene2tx) <- gene2tx$GENEID

  # Function to convert gene IDs to transcript IDs
  convert_to_transcript_ids <- function(gene_ids) {
    transcript_ids <- sapply(gene_ids, function(gid) {
      if (gid %in% gene2tx$GENEID) {
        return(gene2tx[gid, "TXNAME"])
      } else {
        return(gid)
      }
    })
    return(as.character(transcript_ids))
  }

  # Convert raw counts
  gene_counts_tx <- gene_counts
  gene_counts_tx$gene <- convert_to_transcript_ids(gene_counts$gene)
  write.csv(gene_counts_tx, "TransID_counts_raw.csv", row.names = FALSE)
  cat(sprintf("  - TransID_counts_raw.csv (%d genes with transcript IDs)\n", nrow(gene_counts_tx) - 1))

  # Convert normalized counts
  normalized_counts_tx <- normalized_counts_tsv
  normalized_counts_tx$gene <- convert_to_transcript_ids(normalized_counts_tsv$gene)
  write.csv(normalized_counts_tx, "TransID_counts_normalized.csv", row.names = FALSE)
  cat(sprintf("  - TransID_counts_normalized.csv (%d genes with transcript IDs)\n", nrow(normalized_counts_tx) - 1))

  cat(sprintf("  Genes converted to transcript IDs: %d\n", nrow(gene2tx)))
}

cat("Results saved to:\n")
cat("  - deseq2_results.csv (differential expression results)\n")
cat("  - deseq2_normalized_counts.csv (normalized counts for filtered genes)\n")
cat("  - gene_counts_normalized.tsv (normalized counts for ALL genes, 3 decimals)\n")
cat(sprintf("  Sample order in TSV files: %s ... %s\n", 
            tail(colnames(gene_counts), 1), colnames(gene_counts)[2]))

# Generate plots
cat("\nGenerating plots...\n")

# If we have multiple comparisons, generate plots for each
if (!is.null(experimental_design) && length(all_comparisons) > 0) {
  cat("Generating condition-specific plots...\n")

  for (comp_name in names(all_comparisons)) {
    comp <- all_comparisons[[comp_name]]
    res <- comp$results
    safe_name <- comp$safe_name

    cat(sprintf("\n  Creating plots for: %s vs %s\n", comp_name, experimental_design$control))

    # 1. MA plot with LFC shrinkage
    pdf_name <- sprintf("deseq2_MA_plot_%s_vs_%s.pdf", safe_name, experimental_design$control)
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
    pdf_name <- sprintf("deseq2_volcano_plot_%s_vs_%s.pdf", safe_name, experimental_design$control)
    pdf(pdf_name, width = 8, height = 6)
    res_df <- as.data.frame(res)
    res_df$significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1

    print(ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
      geom_point(alpha = 0.5, size = 1) +
      scale_color_manual(values = c("gray", "red")) +
      theme_bw() +
      labs(title = sprintf("Volcano Plot: %s vs %s", comp_name, experimental_design$control),
           x = "Log2 Fold Change",
           y = "-Log10 Adjusted P-value") +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed"))
    dev.off()
  }

} else {
  # Original single comparison plots
  # 1. MA plot with LFC shrinkage
  pdf("deseq2_MA_plot.pdf", width = 10, height = 6)
  par(mfrow=c(1,2))
  plotMA(res, main = "MA Plot (Original)", ylim = c(-5, 5))
  # LFC shrinkage for better visualization
  resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="normal", res=res)
  plotMA(resLFC, main = "MA Plot (LFC Shrinkage)", ylim = c(-5, 5))
  dev.off()

  # 2. Volcano plot
  pdf("deseq2_volcano_plot.pdf", width = 8, height = 6)
  res_df <- as.data.frame(res)
  res_df$significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1

  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c("gray", "red")) +
    theme_bw() +
    labs(title = "Volcano Plot",
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-value") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed")
  dev.off()
}

# Create res_df for scatter plots (use first comparison if experimental_design exists)
if (!is.null(experimental_design) && length(all_comparisons) > 0) {
  # Use the first experimental comparison for scatter plots
  first_comp_name <- names(all_comparisons)[1]
  first_comp <- all_comparisons[[first_comp_name]]
  res <- first_comp$results
  res_df <- as.data.frame(res)
  res_df$significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1
  cat(sprintf("\nUsing %s vs %s comparison for scatter plots\n",
              first_comp$name, experimental_design$control))
} else if (!exists("res_df")) {
  # Fallback: ensure res_df exists
  cat("\nCreating res_df from default results for scatter plots\n")
  res_df <- as.data.frame(res)
  res_df$significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1
}

# Verify res_df exists before scatter plots
if (!exists("res_df")) {
  stop("ERROR: res_df not defined before scatter plots!")
}
cat(sprintf("res_df contains %d genes for scatter plots\n", nrow(res_df)))

# 3. PCA plot
pdf("deseq2_PCA_plot.pdf", width = 8, height = 6)
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
  labs(title = "PCA Plot of Samples")
dev.off()

# 4. Replicate scatter plots
cat("Generating replicate scatter plots...\n")
pdf("deseq2_salmon_replicate_scatter.pdf", width = 12, height = 10)

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
cat("\nGenerating mean scatter plots...\n")

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
  pdf("deseq2_salmon_mean_scatter.pdf", width = 10, height = 10)

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
         main=paste(exp_cond, "vs", experimental_design$control, "(Salmon)"),
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
  pdf("deseq2_salmon_mean_scatter.pdf", width = 10, height = 10)
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
       main=paste(condition2, "vs", condition1, "(Salmon)"),
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
    safe_name <- all_comparisons[[comp_name]]$safe_name
    cat(sprintf("    - deseq2_MA_plot_%s_vs_%s.pdf\n", safe_name, experimental_design$control))
    cat(sprintf("    - deseq2_volcano_plot_%s_vs_%s.pdf\n", safe_name, experimental_design$control))
  }
} else {
  cat("  - deseq2_MA_plot.pdf (original + LFC shrinkage)\n")
  cat("  - deseq2_volcano_plot.pdf\n")
}
cat("  - deseq2_PCA_plot.pdf\n")
cat("  - deseq2_salmon_replicate_scatter.pdf\n")
cat("  - deseq2_salmon_mean_scatter.pdf\n")

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
cat(sprintf("  Total genes in TSV files: %d\n", nrow(txi$counts)))
cat(sprintf("  Genes used in DESeq2 analysis: %d (after filtering low counts)\n", nrow(dds)))
cat(sprintf("  Samples: %d\n", ncol(txi$counts)))
cat("\nTSV file formats:\n")
cat("  - gene_counts_raw.tsv: Integer counts (no decimals)\n")
cat("  - gene_counts_normalized.tsv: Normalized counts (3 decimal places)\n")
cat("=======================================================\n")
