#!/usr/bin/env Rscript
# run_deseq2_from_counts.R
# Usage: Rscript run_deseq2_from_counts.R
# Expects a folder "counts/" with per-sample count files (csv or two-column).
# Example filenames: pvy_control_Rep1.csv, pvy_infected_Rep2.csv, etc.

library(utils)

# ---------------------------
# Params
# ---------------------------
counts_dir <- "counts"
outdir <- "deseq2_results"
plots_dir <- file.path(outdir, "plots")
top_n_heatmap <- 50
alpha <- 0.05
lfc_cutoff <- 1

# Required packages
required_pkgs <- c("DESeq2", "ggplot2", "pheatmap", "RColorBrewer", "data.table", "apeglm")

# ---------------------------
# Helper: check packages and instruct if missing
# ---------------------------
missing_pkgs <- required_pkgs[!sapply(required_pkgs, function(p) requireNamespace(p, quietly = TRUE))]
if (length(missing_pkgs) > 0) {
  cat("ERROR: The following required R packages are missing:\n")
  cat(paste0("  - ", missing_pkgs, collapse = "\n"), "\n\n")
  cat("Please install them before running this script. Two recommended ways:\n\n")
  cat("1) Using conda (recommended):\n")
  cat("   conda install -n <your_env> -c bioconda -c conda-forge \\\n")
  cat("       r-deseq2 r-apeglm r-pheatmap r-ggplot2 r-data.table r-rcolorbrewer\n\n")
  cat("2) From within R (Bioconductor) — run these in an R session:\n")
  cat("   install.packages('BiocManager', repos='https://cloud.r-project.org')\n")
  cat("   BiocManager::install(c('DESeq2','apeglm'))\n")
  cat("   install.packages(c('ggplot2','pheatmap','RColorBrewer','data.table'), repos='https://cloud.r-project.org')\n\n")
  stop("Missing packages. Install them and re-run the script.")
}

# load packages for the data analyses!!
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(data.table)
library(apeglm)

# ---------------------------
# Create output directories
# ---------------------------
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------
# Read count files
# ---------------------------
files <- list.files(counts_dir, pattern = "\\.csv$|\\.counts$|\\.txt$|\\.tsv$", full.names = TRUE)
if (length(files) == 0) stop("No count files found in 'counts/' directory. Place per-sample count files there.")

# Read each file robustly
read_count_file <- function(path) {
  # try to detect separator and header
  # common formats:
  #   gene_id,count   (CSV)
  #   gene_id <tab> count  (two column)
  # handle both
  df <- tryCatch({
    # try read.csv first
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }, error = function(e) {
    tryCatch(read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE),
             error = function(e2) stop("Cannot read file: ", path))
  })
  # Normalize to two columns: gene, count
  if (ncol(df) >= 2) {
    # ensure numeric counts
    counts <- as.numeric(df[[2]])
    gene_ids <- as.character(df[[1]])
    return(data.frame(gene_id = gene_ids, count = counts, stringsAsFactors = FALSE))
  } else {
    stop("File has fewer than 2 columns: ", path)
  }
}

# build merged matrix
message("Reading count files...")
all_genes <- NULL
count_list <- list()
sample_names <- character(length(files))
for (i in seq_along(files)) {
  p <- files[i]
  df <- read_count_file(p)
  if (is.null(all_genes)) {
    all_genes <- df$gene_id
  } else {
    all_genes <- union(all_genes, df$gene_id)
  }
  count_list[[i]] <- df
  nm <- basename(p)
  nm <- sub("\\.(csv|counts|txt|tsv)$", "", nm)
  sample_names[i] <- nm
}

# create count matrix with union of genes (fill missing with 0)
count_mat <- matrix(0L, nrow = length(all_genes), ncol = length(files))
rownames(count_mat) <- all_genes
colnames(count_mat) <- sample_names

for (i in seq_along(count_list)) {
  df <- count_list[[i]]
  # match
  idx <- match(df$gene_id, all_genes)
  count_mat[idx, i] <- as.integer(round(df$count))
}

cat(sprintf("Built count matrix: %d genes x %d samples\n", nrow(count_mat), ncol(count_mat)))

# --- Handle missing values if any (NA) ---
na_total <- sum(is.na(count_mat))
if (na_total > 0) {
  message(sprintf("Warning: found %d NA entries in count matrix — replacing with 0 (assumes missing gene = 0 counts).", na_total))
  # optional: show per-sample NA counts
  print(data.frame(sample = colnames(count_mat), NA_counts = colSums(is.na(count_mat))))
  # Replace NA with 0
  count_mat[is.na(count_mat)] <- 0
}
# coerce to integer (DESeq2 expects integer counts)
mode(count_mat) <- "integer"
# confirm
cat(sprintf("After NA fix: %d genes x %d samples (integer matrix)\n", nrow(count_mat), ncol(count_mat)))


# -----------------------------------------------------------------
# Build sample table from filenames
# ---------------------------

infer_condition <- function(name) {
  if (grepl("control", name, ignore.case = TRUE)) return("control")
  if (grepl("infected|treated|pvY|pvy_infected|infect", name, ignore.case = TRUE)) return("infected")
  # fallback NA
  return(NA)
}
conditions <- sapply(sample_names, infer_condition)
if (any(is.na(conditions))) {
  cat("Warning: some samples have unknown condition inferred. Please ensure filenames contain 'control' or 'infected'.\n")
}

coldata <- data.frame(row.names = sample_names,
                      condition = factor(conditions))
print(coldata)

# --------------------------------------------------------------
# Create DESeq2 dataset and run DE
# ---------------------------
dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = coldata, design = ~ condition)

# filter low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
cat(sprintf("After filter: %d genes remain\n", nrow(dds)))

# set reference level if present
if ("control" %in% levels(dds$condition)) {
  dds$condition <- relevel(dds$condition, ref = "control")
}

dds <- DESeq(dds)

# get results (infected vs control)
res <- tryCatch({
  results(dds, contrast = c("condition", "infected", "control"))
}, error = function(e) {
  # if contrast fails because factor levels unknown, try reversed or default
  message("results() failed: ", e$message)
  results(dds)
})
# shrink LFC using apeglm if available
if ("apeglm" %in% rownames(installed.packages())) {
  # try to find the coef name safely
  coefs <- resultsNames(dds)
  # prefer "condition_infected_vs_control" style or fallback to coef index 2
  if ("condition_infected_vs_control" %in% coefs) {
    res <- lfcShrink(dds, coef="condition_infected_vs_control", type="apeglm")
  } else if (length(coefs) >= 2) {
    res <- lfcShrink(dds, coef=2, type="apeglm")
  } else {
    res <- as.data.frame(res)
  }
} else {
  # try built-in shrink if lfcShrink supports default
  try({ res <- lfcShrink(dds, coef=2, type="normal") }, silent=TRUE)
}

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]

# Reorder so gene_id is first column
if ("gene_id" %in% names(res_df)) {
  res_df <- res_df[, c("gene_id", setdiff(names(res_df), "gene_id"))]
}

# --------------------------------------------------------
# Save outputs (with gene_id as first column)
# ---------------------------
dir.create(outdir, showWarnings = FALSE)

# 1) full results (gene_id first)
write.csv(res_df, file = file.path(outdir, "deseq2_results_full.csv"), row.names = FALSE, quote = FALSE)

# 2) normalized counts: move rownames into first column 'gene_id'
norm_counts <- counts(dds, normalized=TRUE)
norm_df <- as.data.frame(norm_counts, check.names = FALSE)
norm_df <- cbind(gene_id = rownames(norm_df), norm_df)
write.csv(norm_df, file = file.path(outdir, "normalized_counts.csv"), row.names = FALSE, quote = FALSE)

# 3) transforms: vst and rlog; save with gene_id first
vsd <- vst(dds, blind = FALSE)
rld <- rlog(dds, blind = FALSE)

vsd_df <- as.data.frame(assay(vsd), check.names = FALSE)
vsd_df <- cbind(gene_id = rownames(vsd_df), vsd_df)
write.csv(vsd_df, file = file.path(outdir, "vsd_matrix.csv"), row.names = FALSE, quote = FALSE)

rld_df <- as.data.frame(assay(rld), check.names = FALSE)
rld_df <- cbind(gene_id = rownames(rld_df), rld_df)
write.csv(rld_df, file = file.path(outdir, "rld_matrix.csv"), row.names = FALSE, quote = FALSE)

# ===================================
# Plots
# ===========================
dir.create(plots_dir, showWarnings = FALSE)

# PCA
pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))
p <- ggplot(pca, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 3) +
  geom_text(vjust = -1, hjust = 0.5, size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  ggtitle("PCA (VST)") + theme_minimal()
ggsave(filename = file.path(plots_dir, "PCA_vst.png"), plot = p, width = 7, height = 6, dpi = 200)

# MA plot
png(file.path(plots_dir, "MA_plot.png"), width = 1200, height = 900, res = 150)
plotMA(res, main="DESeq2 MA-plot", ylim=c(-5,5))
dev.off()

# Volcano plot
res_plot <- res_df
# ensure padj and log2FoldChange columns exist before computing
if (!("padj" %in% names(res_plot))) res_plot$padj <- NA
if (!("log2FoldChange" %in% names(res_plot))) res_plot$log2FoldChange <- NA
res_plot$logp <- -log10(res_plot$padj + 1e-300)
res_plot$significant <- ifelse(!is.na(res_plot$padj) & res_plot$padj < alpha & abs(res_plot$log2FoldChange) >= lfc_cutoff, "yes", "no")
# remove NA/Inf rows for plotting
res_plot_plotable <- subset(res_plot, !is.na(logp) & is.finite(log2FoldChange))
pvol <- ggplot(res_plot_plotable, aes(x = log2FoldChange, y = logp, color = significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("no" = "#B0B0B0", "yes" = "#D62728")) +
  theme_minimal() +
  xlab("log2 Fold Change (infected vs control)") +
  ylab("-log10 adjusted p-value") +
  ggtitle("Volcano plot")
ggsave(filename = file.path(plots_dir, "volcano.png"), plot = pvol, width = 7, height = 6, dpi = 200)

# Heatmap of top genes
topgenes <- head(na.omit(res_df)$gene_id, n = top_n_heatmap)
if (length(topgenes) > 1) {
  mat <- assay(vsd)[topgenes, , drop = FALSE]
  mat_s <- t(scale(t(mat)))
  anno <- as.data.frame(colData(dds)[, "condition", drop = FALSE])
  rownames(anno) <- colnames(mat_s)
  pal <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  pheatmap(mat_s, cluster_rows = TRUE, cluster_cols = TRUE,
           annotation_col = anno, show_rownames = TRUE,
           filename = file.path(plots_dir, "heatmap_top50.png"),
           color = pal, fontsize_row = 6, main = "Top DE genes (VST)")
}

# significant genes table (gene_id first)
sig <- subset(res_df, !is.na(padj) & padj < alpha)
# if sig has gene_id as first column already, keep order; ensure gene_id column exists
if (!("gene_id" %in% names(sig))) {
  sig$gene_id <- rownames(sig)
  sig <- sig[, c("gene_id", setdiff(names(sig), "gene_id"))]
}
write.csv(sig, file = file.path(outdir, "deseq2_significant_genes.csv"), row.names = FALSE, quote = FALSE)

# session info
writeLines(capture.output(sessionInfo()), con = file.path(outdir, "sessionInfo.txt"))

cat("DESeq2 pipeline finished. Results and plots are in", outdir, "\n")
