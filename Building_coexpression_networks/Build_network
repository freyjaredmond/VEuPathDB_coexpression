# Building co-expression networks
# Builds Pearson, mutual rank (MR) and WGCNA co-expression networks for four organisms
# (P. berghei, C. auris, A. nidulans, T. vaginalis) and extracts the top 50 co-expressed
# genes for each gene of interest.
#
# Input (per organism, set in the `organisms` list):
#   - expr:  combined RNA-seq TPM table (.tsv), one row per gene, first column "Gene ID";
#            non-sample columns are removed with `drop`
#   - genes: list of genes of interest (.txt) (every gene) with a "Gene ID" column
#   - batch: batch (dataset) label for each sample, in column order
#
# Output (per organism):
#   - {org}_expression_matrix.tsv     batch-corrected log2(TPM + 1) matrix
#   - {org}_soft_power.tsv            soft-thresholding power used
#   - {org}_blockwiseTOM-block.1.RData  WGCNA topological overlap matrix
#   - {org}_gene_modules.tsv          WGCNA module assignment for each gene
#   - {org}_pearson_top50.tsv         top 50 genes by Pearson correlation
#   - {org}_mr_top50.tsv              top 50 genes by mutual rank
#   - {org}_wgcna_top50.tsv           top 50 genes by TOM (whole network)
#   - {org}_modules_top50.tsv         top 50 genes by TOM within the gene's module

library(dplyr)
library(sva)
library(WGCNA)
library(limma)

while (sink.number() > 0) sink()

# get all the info needed for each organism
organisms <- list(
  berg  = list(expr = "Berg_combined.tsv", drop = c(2, 3),
               genes = "Berg_genes.txt",
               batch = rep(c("Berg_ExErythro", "Berg_AP2KO", "Berg_LiverStage", "Berg_GamSex",
                             "Berg_DietRestrict", "Berg_AsexGametes", "Berg_AsexSex"),
                           times = c(7, 20, 14, 3, 4, 12, 9))),
  auris = list(expr = "Auris_combined.tsv", drop = c(2, 3, 4, 5),
               genes = "Auris_genes.txt",
               batch = rep(c("Auris_FluconazoleEvol", "Auris_ClinIsol", "Auris_Isolates", "Auris_DonorTimecourse",
                             "Auris_CarbonSource", "Auris_NCPFBiofilm"),
                           times = c(2, 4, 5, 18, 5, 8))),
  nidulan = list(expr = "Nidulan_combined.tsv", drop = c(2, 3, 4, 5),
                 genes = "Nidulan_genes.txt",
                 batch = rep(c("Nid_Development", "Nid_mtfA_veA", "Nid_Nitrogen",
                               "Nid_Conidia", "Nid_AminoAcid", "Nid_Cellulase"),
                             times = c(3, 4, 5, 4, 3, 3))),
  trich = list(expr = "Trich_combined.tsv", drop = c(2, 3, 4, 5),
               genes = "trich_genes.txt",
               batch = rep(c("Trich_MTZ", "Trich_Mycoplasma", "Trich_Apicidin"),
                           times = c(2, 4, 2)))
)

# top 50 by module (TOM within the GOI's module)
get_top50 <- function(genes_of_interest, TOM, data_type, module_colors) {
  results <- data.frame()
  for (gene in genes_of_interest) {
    # genes that were filter out prior to wcgna
    if (!gene %in% names(module_colors)) {
      results <- rbind(results, data.frame(dataset = data_type, gene = gene, correlated_gene = NA,
                                           TOM_similarity = NA,
                                           module = "not in network", stringsAsFactors = FALSE))
      next
    }
    gene_module <- module_colors[gene]
    # genes that were filtered out by wgcna
    if (gene_module == "grey") {
      results <- rbind(results, data.frame(dataset = data_type, gene = gene, correlated_gene = NA,
                                           TOM_similarity = NA,
                                           module = "grey", stringsAsFactors = FALSE))
      next
    }
    # get genes in the same TOM as GOI
    same_module_genes <- names(module_colors[module_colors == gene_module])
    same_module_genes <- same_module_genes[same_module_genes != gene]
    # get the top genes by TOM for GOI
    gene_similarities <- TOM[gene, same_module_genes]
    gene_similarities <- sort(gene_similarities, decreasing = TRUE)
    top50 <- head(gene_similarities, 50)

    results <- rbind(results, data.frame(
      dataset = data_type,
      gene = gene,
      correlated_gene = names(top50),
      TOM_similarity = as.numeric(top50),
      module = gene_module,
      stringsAsFactors = FALSE
    ))
  }
  return(results)
}

# top 50 across the whole network
# decreasing parameter added for mr vs pearson
get_50 <- function(genes_of_interest, TOM, data_type, decreasing = TRUE) {
  results <- data.frame()
  for (gene in genes_of_interest) {
    # genes not in the network
    if (!gene %in% rownames(TOM)) {
      results <- rbind(results, data.frame(dataset = data_type, gene = gene,
                                           correlated_gene = NA, TOM_similarity = NA,
                                           stringsAsFactors = FALSE))
      next
    }
    # rank every other gene to the GOI
    other_genes <- setdiff(rownames(TOM), gene)
    gene_similarities <- sort(TOM[gene, other_genes], decreasing = decreasing)
    top50 <- head(gene_similarities, 50)

    results <- rbind(results, data.frame(
      dataset = data_type,
      gene = gene,
      correlated_gene = names(top50),
      TOM_similarity = as.numeric(top50),
      stringsAsFactors = FALSE
    ))
  }
  rownames(results) <- NULL
  return(results)
}


for (org in names(organisms)) {
  s <- organisms[[org]]

  # Load and prepare data
  data <- readr::read_delim(s$expr, delim = "\t", na = c("", "NA", "N/A"))
  data <- data %>% distinct(`Gene ID`, .keep_all = TRUE)

  data <- data[, -s$drop]
  data <- na.omit(data)
  gene_ids <- data[[1]]
  expr_matrix <- as.matrix(data[, -1])

  #filter out lowly expressed genes
  rownames(expr_matrix) <- gene_ids
  keep <- rowSums(expr_matrix >= 1) >= 0.1 * ncol(expr_matrix)  # >= 1 TPM in at least 10% of samples
  expr_matrix <- expr_matrix[keep, ]
  gene_ids    <- gene_ids[keep]

  # Pearson / MR input: batch correct only
  batch_vector <- s$batch
  log2_tpm <- log2(expr_matrix + 1)
  corrected <- ComBat(dat = log2_tpm, batch = batch_vector)

  rownames(corrected) <- gene_ids
  wgcna_data <- t(corrected)
  write.table(corrected, paste0(org, "_expression_matrix.tsv"), sep = "\t",
              row.names = TRUE, col.names = NA, quote = FALSE)

  # Pearson
  r <- cor(wgcna_data)
  pearsons <- 0.5 * (1 + r)

  # mutual rank: geometric mean of the two ranks (lower = stronger)
  r_mr <- r
  diag(r_mr) <- NA
  R <- apply(-r_mr, 2, rank, na.last = "keep", ties.method = "average")
  dimnames(R) <- dimnames(r_mr)
  mr <- sqrt(R * t(R))
  rm(R, r_mr)

  # WGCNA: quantile normalise then batch correct
  log2_tpm <- normalizeBetweenArrays(log2_tpm, method = "quantile")
  wgcna_data <- t(ComBat(dat = log2_tpm, batch = batch_vector))

  # Soft threshold: WGCNA's estimate (lowest power with scale-free fit R^2 > 0.85)
  powers <- c(1:20)
  sft <- pickSoftThreshold(wgcna_data,
                           powerVector = powers,
                           networkType = "signed",
                           RsquaredCut = 0.85,
                           verbose = 5)
  sp <- sft$powerEstimate
  # no power reached R^2 0.85: WGCNA FAQ default for signed networks, by number of samples
  if (is.na(sp)) {
    n_samples <- nrow(wgcna_data)
    sp <- if (n_samples < 20) 18 else if (n_samples <= 30) 16 else if (n_samples <= 40) 14 else 12
  }
  write.table(data.frame(organism = org, soft_power = sp), paste0(org, "_soft_power.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  # Run WGCNA
  temp_cor <- cor
  cor <- WGCNA::cor
  net <- blockwiseModules(
    wgcna_data,
    power = sp,
    networkType = "signed",
    deepSplit = 2,
    pamRespectsDendro = FALSE,
    minModuleSize = 30,
    maxBlockSize = 15000,
    reassignThreshold = 0,
    mergeCutHeight = 0.25,
    saveTOMs = TRUE,
    saveTOMFileBase = paste0(org, "_blockwiseTOM"),
    numericLabels = FALSE,
    verbose = 3
  )
  cor <- temp_cor
  module_colors <- net$colors

  # Load genes of interest
  data_type <- org
  genes_of_interest <- readr::read_delim(s$genes, delim = "\t") %>%
    pull("Gene ID") %>% unique()

  # Load TOM and assign gene names
  load(net$TOMFiles[1])
  TOM <- as.matrix(TOM)
  # only want gene names that wgcna didn't filter out
  network_genes <- colnames(wgcna_data)[net$goodGenes]
  rownames(TOM) <- network_genes
  colnames(TOM) <- network_genes

  # get module assignments
  mods <- module_colors[net$goodGenes]
  gene_modules <- data.frame(Gene = network_genes, Module = mods)
  write.table(gene_modules, paste0(org, "_gene_modules.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

  # top 50s
  write.table(get_top50(genes_of_interest, TOM, data_type, module_colors),
              paste0(org, "_modules_top50.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(get_50(genes_of_interest, TOM, data_type, decreasing = TRUE),
              paste0(org, "_wgcna_top50.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(get_50(genes_of_interest, pearsons, data_type, decreasing = TRUE),
              paste0(org, "_pearson_top50.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(get_50(genes_of_interest, mr, data_type, decreasing = FALSE),   # lowest MR first
              paste0(org, "_mr_top50.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

  rm(data, expr_matrix, log2_tpm, corrected, wgcna_data, r, pearsons, mr, TOM, net)
  gc()
}
