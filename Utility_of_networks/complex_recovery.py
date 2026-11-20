# Complex partner recovery
# For each subunit of five P. falciparum complexes, calculates the percentage of the other subunits
# found in its top 10, 20 and 50 co-expressed partners (Pearson, MR, WGCNA TOM, WGCNA TOM within module),
# and averages across subunits to give a score per complex. Subunits missing from the network are excluded.
# Random baseline: for each subunit, the same number of randomly selected genes, averaged over 1,000 draws.
#
# Input:
#   - plas_pearson_top50.tsv, plas_mr_top50.tsv, plas_wgcna_top50.tsv, plas_modules_top50.tsv
#                                                top 50 partners per gene (see build_network.R)
#   - plas_complexes dictionary below see complex_pairs.py
#
# Output:
#   - complex_partner_recovery_plas.tsv          mean, SD and median recovery per complex, method and cutoff
#   - complex_partner_recovery_per_gene_plas.tsv recovery for every subunit
#   - complex_partner_recovery.png               mean recovery per complex

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os
import random
import seaborn as sns

# P. falciparum complex subunits
plas_complexes = {
    "Proteasome": ["PF3D7_0108000", "PF3D7_0317000", "PF3D7_0518300", "PF3D7_0608500", "PF3D7_0727400",
                   "PF3D7_0803800", "PF3D7_0807500", "PF3D7_0931800", "PF3D7_1011400", "PF3D7_1328100",
                   "PF3D7_1353800", "PF3D7_1353900", "PF3D7_1470900", "PF3D7_1474800"],
    "Polymerase": ["PF3D7_0110400", "PF3D7_0205500", "PF3D7_0215700", "PF3D7_0303300", "PF3D7_0318200",
                   "PF3D7_0708100", "PF3D7_0923000", "PF3D7_1027400", "PF3D7_1213700", "PF3D7_1304900",
                   "PF3D7_1342700", "PF3D7_1364800"],
    "CCT": ["PF3D7_0214000", "PF3D7_0306800", "PF3D7_0308200", "PF3D7_0320300",
            "PF3D7_0608700", "PF3D7_1132200", "PF3D7_1229500", "PF3D7_1357800"],
    "MCM": ["PF3D7_1417800", "PF3D7_0527000", "PF3D7_1317100",
            "PF3D7_1211700", "PF3D7_1355100", "PF3D7_0705400"],
    "ATP": ["PF3D7_1311900", "PF3D7_0406100", "PF3D7_0106100", "PF3D7_1341900", "PF3D7_0934500",
            "PF3D7_1140100", "PF3D7_1323200", "PF3D7_1306600", "PF3D7_0806800", "PF3D7_0519200",
            "PF3D7_1354400", "PF3D7_1464700", "PF3D7_0721900"],
}

### top 50 files for each network
datasets=["plas_pearson_top50.tsv",
          "plas_modules_top50.tsv",
          "plas_wgcna_top50.tsv",
          "plas_mr_top50.tsv"]

cutoffs=[10,20,50]

network_results = []
gene_rows = []
for ds in datasets:
    ds_name = os.path.basename(ds).replace("plas_", "").replace("_top50.tsv", "")
    df = pd.read_csv(ds, sep="\t")
    # genes with at least one partner in this network
    network_present = set(df.loc[df["correlated_gene"].notna(), "gene"])
    for cut in cutoffs:
        for complex_name, complex_list in plas_complexes.items():
            # drop members missing from the network so they don't score 0 or count in the denominator
            complex_genes = set(complex_list) & network_present
            # need at least 2 members in the network to have a partner to recover
            if len(complex_genes) < 2:
                continue
            pct_match_all = []
            for c in sorted(complex_genes):
                df_c = df[df["gene"] == c].head(cut)
                network_genes = set(df_c["correlated_gene"])
                matches = (complex_genes - {c}) & network_genes
                pct = len(matches) / (len(complex_genes) - 1) * 100
                pct_match_all.append(pct)
                gene_rows.append({"Dataset": ds_name, "Cutoff": cut,
                                  "Complex": complex_name, "Gene": c, "Percent": pct})
            pct_match_all = pd.Series(pct_match_all)
            network_results.append({"Dataset": ds_name,
                                    "Cutoff": cut,
                                    "Complex": complex_name,
                                    "Members_in_network": len(complex_genes),
                                    "Mean Percentage Match": pct_match_all.mean(),
                                    "Percentage Match Standard Deviation": pct_match_all.std(),
                                    "Median": pct_match_all.median()})

### generate the random control
# all networks contain the same genes, so one random control covers every network
random.seed(28)
all_genes = list(df.loc[df["correlated_gene"].notna(), "gene"].unique())
n_random = 1000

for cut in cutoffs:
    for complex_name, complex_list in plas_complexes.items():
        # same filtering as the networks, keeping list order so the seed stays reproducible
        complex_present = [g for g in complex_list if g in network_present]
        complex_genes = set(complex_present)
        if len(complex_genes) < 2:
            continue
        pct_match_all = []
        for c in complex_present:
            other_genes = [g for g in all_genes if g != c]
            # average over many random partner sets, a single draw is almost always 0
            random_pct = []
            for i in range(n_random):
                random_genes = set(random.sample(other_genes, cut))
                matches = (complex_genes - {c}) & random_genes
                random_pct.append(len(matches) / (len(complex_genes) - 1) * 100)
            random_mean = sum(random_pct) / len(random_pct)
            pct_match_all.append(random_mean)
            gene_rows.append({"Dataset": "random", "Cutoff": cut,
                              "Complex": complex_name, "Gene": c, "Percent": random_mean})
        pct_match_all = pd.Series(pct_match_all)
        network_results.append({"Dataset": "random",
                                "Cutoff": cut,
                                "Complex": complex_name,
                                "Members_in_network": len(complex_genes),
                                "Mean Percentage Match": pct_match_all.mean(),
                                "Percentage Match Standard Deviation": pct_match_all.std(),
                                "Median": pct_match_all.median()})

network_results_df = pd.DataFrame(network_results)
print(network_results_df.to_string(index=False))
network_results_df.to_csv("complex_partner_recovery_plas.tsv", sep="\t", index=False)

gene_df = pd.DataFrame(gene_rows)
gene_df.to_csv("complex_partner_recovery_per_gene_plas.tsv", sep="\t", index=False)

### plot
# colour-blind safe (Okabe-Ito) palette, same colours as the GO matching plots
palette = {"wgcna": "#0072B2",     # blue
           "modules": "#E69F00",   # orange
           "pearson": "#CC79A7",   # reddish purple
           "mr": "#009E73",        # bluish green
           "random": "#999999"}    # grey
hue_order = ["wgcna", "modules", "pearson", "mr", "random"]
method_labels = {"wgcna": "WGCNA", "modules": "Modules", "pearson": "Pearson", "mr": "MR", "random": "Random"}

g = sns.catplot(data=gene_df, x="Cutoff", y="Percent", hue="Dataset", hue_order=hue_order,
                col="Complex", kind="point", palette=palette, height=3.2, aspect=0.8,
                estimator="mean", errorbar=None, dodge=0.4, legend=False)
for cx, ax in g.axes_dict.items():
    ax.set_title(cx, fontsize=15)

# x label under every panel and legend on the right, as in the GO match plots
g.set_axis_labels("Top n coexpressed partners", "Complex partners recovered (%)")
g.tight_layout()
# legend added after the layout and attached to the last panel, so there's no gap left for seaborn's own legend
last_ax = g.axes.flat[-1]
handles = [Line2D([0], [0], color=palette[m], marker="o", linewidth=2.5, markersize=7) for m in hue_order]
last_ax.legend(handles, [method_labels[m] for m in hue_order], title="Method", frameon=False,
               loc="upper left", bbox_to_anchor=(1.02, 1))
g.savefig("complex_partner_recovery.png", dpi=300, bbox_inches="tight")
plt.show()
