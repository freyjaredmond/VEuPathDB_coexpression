# Plot sample-size analysis
# Plots neighbour-voting AUROC and complex AUPRC against number of samples for the
# P. falciparum MR subsampling draws, comparing full draws with half-sample draws.
# Dot size shows the number of datasets in each draw.
#
# Input:
#   - subsample_results_plas_mr.tsv         
#   - subsample_results_plasmo_half_mr.tsv  
#
# Output:
#   - plas_mr_nv_dot_half.png   neighbour-voting and degree-null AUROC vs number of samples
#   - plas_mr_pr_dot_half.png   complex AUPRC vs number of samples
#
# Packages: pandas, matplotlib

import pandas as pd
import matplotlib.pyplot as plt

# colour-blind safe colours (Okabe-Ito palette)
FULL = "#0072B2"         # blue: all samples
HALF = "#E69F00"         # orange: half samples
NULL = "#CC79A7"         # reddish purple: degree null
THRESHOLD = "#D55E00"    # vermillion: 0.7 threshold line (also dotted, so it doesn't rely on colour)

plas_full = pd.read_csv("subsample_results_plas_mr.tsv", sep="\t")
plas_half = pd.read_csv("subsample_results_plasmo_half_mr.tsv", sep="\t")


### neighbour-voting AUROC
fig, ax = plt.subplots()
sc = ax.scatter(plas_full["n_samples"], plas_full["nv_auroc"], s=plas_full["n_datasets"] * 8,
                color=FULL, alpha=0.7, label="Neighbour voting")
ax.scatter(plas_full["n_samples"], plas_full["degree_null_auroc"], s=plas_full["n_datasets"] * 8,
           color=NULL, alpha=0.7, label="Degree null")
ax.scatter(plas_half["n_samples"], plas_half["nv_auroc"], s=plas_half["n_datasets"] * 8,
           color=HALF, alpha=0.7, label="Neighbour voting (half)")
ax.axhline(0.5, linestyle="--", color="black", label="Random")
ax.axhline(0.7, linestyle=":", linewidth=2, color=THRESHOLD, label="Threshold")
ax.set_xlabel("Number of samples")
ax.set_ylabel("AUROC")
ax.set_xticks(range(0, 301, 50))

# legend 1: colours (method) and lines, placed below the size legend
colour_legend = ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 0.55))
ax.add_artist(colour_legend)          # keep it when the second legend is added

# legend 2: dot sizes (number of datasets); divide by 8 to undo the size scaling above
handles, labels = sc.legend_elements(prop="sizes", num=[5, 10, 15, 20, 25],
                                     color="grey", alpha=0.6, func=lambda s: s / 8)
ax.legend(handles, labels, title="Datasets", frameon=False,
          loc="upper left", bbox_to_anchor=(1.02, 1))

# save before showing, otherwise some backends save a blank figure
fig.savefig("plas_mr_nv_dot_half.png", dpi=300, bbox_inches="tight",
            bbox_extra_artists=[colour_legend])
plt.show()


### complex AUPRC
fig, ax = plt.subplots()
sc = ax.scatter(plas_full["n_samples"], plas_full["auprc"], s=plas_full["n_datasets"] * 8,
                color=FULL, alpha=0.7, label="AUPRC")
ax.scatter(plas_half["n_samples"], plas_half["auprc"], s=plas_half["n_datasets"] * 8,
           color=HALF, alpha=0.7, label="AUPRC (half)")
# random = proportion of positive pairs
ax.axhline(plas_full["prior"].iloc[0], linestyle="--", color="black", label="Random")
ax.set_xlabel("Number of samples")
ax.set_ylabel("AUPRC")
ax.set_xticks(range(0, 301, 50))

# legend 1: colours (method) and lines
colour_legend = ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 0.55))
ax.add_artist(colour_legend)

# legend 2: dot sizes (number of datasets)
handles, labels = sc.legend_elements(prop="sizes", num=[5, 10, 15, 20, 25],
                                     color="grey", alpha=0.6, func=lambda s: s / 8)
ax.legend(handles, labels, title="Datasets", frameon=False,
          loc="upper left", bbox_to_anchor=(1.02, 1))

fig.savefig("plas_mr_pr_dot_half.png", dpi=300, bbox_inches="tight",
            bbox_extra_artists=[colour_legend])
plt.show()
