# Sample-size analysis (full draws): dataset draw -> ComBat -> Pearson -> mutual rank ->
# rank-standardised MR network -> EGAD + complex AUPRC
# Builds a ranked MR network from every sample of the datasets in one P. falciparum draw and scores it
# with EGAD neighbour voting (GO AUROC, degree null) and complex-pair AUPRC.
# Run once per draw, e.g. as an array job on a cluster.
# Input:
#   - Plasmo_combined.tsv              combined RNA-seq TPM table, first column "Gene ID"
#   - draws_plas.tsv                   dataset draws (from dataset_draws.py)
#   - plasmo_computed_propogated.tsv   propagated computed GO annotations (see propagate_go.py)
#   - plasmo_complex_pairs.tsv         complex gene pairs (see complex_pairs.py)
#
# Output:
#   - results_mrrank_full_plas/draw_<draw_id>.tsv   one row: draw_id, n_datasets, rep, n_samples, nv_auroc,
#                                                   degree_null_auroc, auprc, prior, n_go_terms, n_na_corr

suppressPackageStartupMessages({
  library(sva)
  library(EGAD)
  library(ROCR)
})

# draw number from the command line, e.g. Rscript subsample_egad_mrrank_full.R 12
draw_id <- as.integer(commandArgs(trailingOnly = TRUE)[1])
cat("Draw:", draw_id, "\n")

# load expression data
data <- read.delim("Plasmo_combined.tsv", check.names = FALSE)
data <- data[!duplicated(data[["Gene ID"]]), ]
data <- data[, -c(2, 3)]          # remove non-sample columns
data <- na.omit(data)
gene_ids    <- data[[1]]
expr_matrix <- as.matrix(data[, -1])
rownames(expr_matrix) <- gene_ids

# filter lowly expressed genes on the FULL dataset, so every draw has the same genes
keep <- rowSums(expr_matrix >= 1) >= 0.1 * ncol(expr_matrix)  # >= 1 TPM in at least 10% of samples
expr_matrix <- expr_matrix[keep, ]

# dataset label for every sample, in column order
# (same order as names/times in dataset_draws.py, so column positions line up)
batch_vector <- rep(c("Plas_1", "Plas_10", "Plas_11", "Plas_12", "Plas_13", "Plas_14", "Plas_15", "Plas_16", "Plas_17",
                      "Plas_18", "Plas_19", "Plas_2", "Plas_20", "Plas_21", "Plas_22", "Plas_23", "Plas_24", "Plas_25",
                      "Plas_26", "Plas_27", "Plas_3", "Plas_4", "Plas_5", "Plas_6", "Plas_7", "Plas_8", "Plas_9"),
                    times = c(3, 2, 7, 43, 3, 8, 8, 6, 44, 25, 27, 8, 3, 4, 7, 5, 3, 2, 7, 24, 8, 21, 2, 2, 4, 6, 3))
stopifnot(length(batch_vector) == ncol(expr_matrix))

# keep every sample from the datasets in this draw
draws <- read.delim("draws_plas.tsv")
this_draw <- draws[draws$draw_id == draw_id, ]
chosen <- strsplit(this_draw$datasets, ",")[[1]]

cols <- batch_vector %in% chosen
expr_sub  <- expr_matrix[, cols]
batch_sub <- batch_vector[cols]
cat("Datasets:", length(chosen), " Samples:", ncol(expr_sub), "\n")
stopifnot(ncol(expr_sub) == this_draw$n_samples)

# log2(TPM + 1) and batch correction
log2_tpm  <- log2(expr_sub + 1)
corrected <- ComBat(dat = log2_tpm, batch = batch_sub)

# Pearson correlation (genes x genes)
r <- cor(t(corrected))
# genes with no variation in a small draw give NA correlations: treat as no correlation
n_na <- sum(is.na(r))
r[is.na(r)] <- 0
cat("NA correlations replaced:", n_na, "\n")
rm(corrected, log2_tpm, expr_sub); gc()

# mutual rank: geometric mean of the two ranks (lower = stronger)
diag(r) <- NA                                # ignore self-correlation
R <- apply(-r, 2, rank, na.last = "keep", ties.method = "average")
dimnames(R) <- dimnames(r)
mr <- sqrt(R * t(R))
rm(r, R); gc()

# rank-standardise the whole network and divide by the max rank -> 0-1 scale
# (same as rank_std in evaluate_coexpression_networks.R)
rank_std <- function(net) {
  diag(net) <- NA
  rk <- rank(net, na.last = "keep", ties.method = "average")
  net[] <- rk / max(rk, na.rm = TRUE)
  diag(net) <- 0
  net
}
network <- rank_std(-mr)                     # lowest MR -> highest weight
rm(mr); gc()

# EGAD neighbour voting (GO terms with 20-300 genes, 3-fold cross-validation)
go <- read.delim("plasmo_computed_propogated.tsv", check.names = FALSE)
labels <- make_annotations(go[, 1:2], rownames(network), unique(go[[2]]))

set.seed(1)
neigh <- run_GBA(network, labels, min = 20, max = 300, nfold = 3)
nv_auroc          <- neigh[[3]]
degree_null_auroc <- mean(neigh[[1]][, "degree_null_auc"], na.rm = TRUE)

# complex-pair AUPRC (prior = proportion of positive pairs, the random expectation)
pairs <- read.delim("plasmo_complex_pairs.tsv")
pairs <- pairs[pairs$Gene1 %in% rownames(network) & pairs$Gene2 %in% rownames(network), ]
scores <- network[cbind(pairs$Gene1, pairs$Gene2)]
pred   <- prediction(scores, pairs$Label)
auprc  <- performance(pred, "aucpr")@y.values[[1]]
prior  <- mean(pairs$Label)

# one row of results for this draw
out <- data.frame(
  draw_id           = draw_id,
  n_datasets        = this_draw$n_datasets,
  rep               = this_draw$rep,
  n_samples         = this_draw$n_samples,
  nv_auroc          = nv_auroc,
  degree_null_auroc = degree_null_auroc,
  auprc             = auprc,
  prior             = prior,
  n_go_terms        = nrow(neigh[[1]]),
  n_na_corr         = n_na
)
print(out)
dir.create("results_mrrank_full_plas", showWarnings = FALSE)
write.table(out, sprintf("results_mrrank_full_plas/draw_%02d.tsv", draw_id),
            sep = "\t", row.names = FALSE, quote = FALSE)
