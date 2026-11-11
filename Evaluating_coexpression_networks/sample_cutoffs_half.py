# Half-sample dataset draws for subsampling
# For each P. falciparum dataset draw in draws_plas.tsv, keeps the same datasets but randomly
# selects half of the samples from each one. Comparing full and half draws separates the effect
# of sample number from the effect of dataset number.
#
# Input:
#   - draws_plas.tsv   full dataset draws (from dataset_draws.py)
#   - dataset names and sample counts below (same order as the batch vector in the R script,
#     so column positions match the combined expression table)
#
# Output:
#   - draws_plas_half.tsv   draw_id, n_datasets, rep, fraction, n_samples, n_samples_full,
#                           datasets, columns (sample columns to keep, numbered from 1 for R)


import pandas as pd
import random
import math

random.seed(28)

# P. falciparum datasets and number of samples in each (same order as the R batch vector)
names = ["Plas_1", "Plas_10", "Plas_11", "Plas_12", "Plas_13", "Plas_14", "Plas_15", "Plas_16", "Plas_17",
         "Plas_18", "Plas_19", "Plas_2", "Plas_20", "Plas_21", "Plas_22", "Plas_23", "Plas_24", "Plas_25",
         "Plas_26", "Plas_27", "Plas_3", "Plas_4", "Plas_5", "Plas_6", "Plas_7", "Plas_8", "Plas_9"]
times = [3, 2, 7, 43, 3, 8, 8, 6, 44, 25, 27, 8, 3, 4, 7, 5, 3, 2, 7, 24, 8, 21, 2, 2, 4, 6, 3]

# proportion of samples to keep from each dataset
fraction = 0.5

# column positions of every sample, grouped by dataset
positions = {}
start = 0
for name, n in zip(names, times):
    positions[name] = list(range(start, start + n))
    start += n

# reuse the full draws, so each half draw has the same datasets as its full draw
full_draws = pd.read_csv("draws_plas.tsv", sep="\t")

rows = []
for idx, draw in full_draws.iterrows():
    chosen = draw["datasets"].split(",")
    cols = []
    for d in chosen:
        # keep half of each dataset (rounded up), but never fewer than 2 (ComBat needs 2 per batch)
        n_keep = max(2, math.ceil(len(positions[d]) * fraction))
        cols += random.sample(positions[d], n_keep)
    cols.sort()
    rows.append({"draw_id": draw["draw_id"], "n_datasets": draw["n_datasets"], "rep": draw["rep"],
                 "fraction": fraction, "n_samples": len(cols), "n_samples_full": draw["n_samples"],
                 "datasets": draw["datasets"],
                 "columns": ",".join(str(c + 1) for c in cols)})   # +1 because R counts columns from 1

plas_half_draws = pd.DataFrame(rows)
plas_half_draws.to_csv("draws_plas_half.tsv", sep="\t", index=False)
print(plas_half_draws[["draw_id", "n_datasets", "rep", "n_samples_full", "n_samples"]])
