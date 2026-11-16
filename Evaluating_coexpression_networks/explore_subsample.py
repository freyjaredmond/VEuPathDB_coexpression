# Define the minimum sample size
# Bins the subsampled rank-standardised MR networks by number of samples and, for each genome,
# counts how many networks in each bin have a neighbour-voting AUROC at or below 0.7.
# A bin passes the cutoff only if every network in it is above 0.7.
#
# Input (one per genome, from subsample_egad_mrrank_full.R):
#   - subsample_results_toxo__mrrank_full.tsv
#   - subsample_results_plas__mrrank_full.tsv
#   - subsample_results_anoph__mrrank_full.tsv
#   - subsample_results_asper__mrrank_full.tsv
#
# Output:
#   - cutoff_bins_all_genomes.tsv   for each genome and sample bin: number of networks, number at or
#                                   below 0.7, mean samples, mean datasets, mean AUPRC, pass/fail
#
# Packages: pandas

import pandas as pd

# results file for each genome
genomes = {"Toxoplasma gondii": "subsample_results_toxo__mrrank_full.tsv",
           "Plasmodium falciparum": "subsample_results_plas__mrrank_full.tsv",
           "Anopheles gambiae": "subsample_results_anoph__mrrank_full.tsv",
           "Aspergillus fumigatus": "subsample_results_asper__mrrank_full.tsv"}

# sample-size bins: (0,10], (10,20], ... (60,70], 70+
bins = [0, 10, 20, 30, 40, 50, 60, 70, float("inf")]
labels = ["0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70+"]

all_results = []
for genome, path in genomes.items():
    df = pd.read_csv(path, sep="\t")
    df["sample_bin"] = pd.cut(df["n_samples"], bins=bins, labels=labels)

    results = []
    for bin_label, data in df.groupby("sample_bin", observed=True):
        # networks in this bin with neighbour-voting AUROC at or below 0.7
        n_below = (data["nv_auroc"] <= 0.7).sum()
        results.append({"Genome": genome,
                        "Bin": bin_label,
                        "n_draws": len(data),
                        "n_below": n_below,
                        # mean samples, datasets and complex AUPRC across the networks in this bin
                        "mean_samples": round(data["n_samples"].mean(), 1),
                        "mean_datasets": round(data["n_datasets"].mean(), 1),
                        "mean_auprc": round(data["auprc"].mean(), 3),
                        # bin passes only if no network in it is at or below 0.7
                        "Cutoff": "Yes" if n_below == 0 else "No"})

    results = pd.DataFrame(results)
    print(f"\n{genome}")
    print(results.drop(columns="Genome").to_string(index=False))
    all_results.append(results)

# all genomes in one table
all_results = pd.concat(all_results, ignore_index=True)
all_results.to_csv("cutoff_bins_all_genomes.tsv", sep="\t", index=False)
