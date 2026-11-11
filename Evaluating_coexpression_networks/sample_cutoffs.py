# Random dataset draws for subsampling
# Randomly selects combinations of P. falciparum RNA-seq datasets at increasing dataset counts,
# to test whether network performance depends on the number of datasets or samples.
# Each draw is later used to build a network and calculate neighbour-voting AUROC and complex AUPRC.
#
# Input:
#   - dataset names and sample counts below (same order and counts as the ComBat batch vector)
#
# Output:
#   - draws_plas.tsv   draw_id, n_datasets, rep, n_samples, datasets (comma-separated dataset names)

import pandas as pd
import random

random.seed(28)

# number of datasets to draw
cutoffs = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26]

# P. falciparum datasets and number of samples in each
names = ["Plas_1", "Plas_10", "Plas_11", "Plas_12", "Plas_13", "Plas_14", "Plas_15", "Plas_16", "Plas_17",
         "Plas_18", "Plas_19", "Plas_2", "Plas_20", "Plas_21", "Plas_22", "Plas_23", "Plas_24", "Plas_25",
         "Plas_26", "Plas_27", "Plas_3", "Plas_4", "Plas_5", "Plas_6", "Plas_7", "Plas_8", "Plas_9"]
times = [3, 2, 7, 43, 3, 8, 8, 6, 44, 25, 27, 8, 3, 4, 7, 5, 3, 2, 7, 24, 8, 21, 2, 2, 4, 6, 3]
batch_dict = dict(zip(names, times))

# 5 random draws at each cutoff, because sample numbers vary a lot between datasets
rows = []
for k in cutoffs:
    for rep in range(1, 6):
        chosen = random.sample(names, k)
        n_samples = sum(batch_dict[d] for d in chosen)
        rows.append({"n_datasets": k, "rep": rep, "n_samples": n_samples,
                     "datasets": ",".join(chosen)})

# all datasets, once
rows.append({"n_datasets": len(names), "rep": 1, "n_samples": sum(times),
             "datasets": ",".join(names)})

plas_draws = pd.DataFrame(rows)
plas_draws.insert(0, "draw_id", range(1, len(plas_draws) + 1))  # draw ID as the first column
plas_draws.to_csv("draws_plas.tsv", sep="\t", index=False)
print(plas_draws[["draw_id", "n_datasets", "rep", "n_samples"]])
