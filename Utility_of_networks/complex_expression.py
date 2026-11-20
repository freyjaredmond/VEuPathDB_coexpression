# Complex expression in P. falciparum
# Plots mean-centred expression of the five complexes (20S proteasome, RNA polymerase II, CCT, MCM2-7, V-ATPase)
# across all samples used to build the P. falciparum co-expression network, summarises the spread (SD) of each
# complex, and clusters the average expression profiles of the complexes.
#
# Input:
#   - plas_expression_matrix.tsv                       batch-corrected log2(TPM + 1) matrix (genes x samples);
#                                                       sample columns are named "<dataset> | <sample>"
#
# Output:
#   - plas_complex_variance.tsv                    per complex: number of genes, mean expression, mean/min/max SD
#   - plas_dataset_lookup.tsv                      plot number -> dataset file -> number and names of samples
#   - plas_complex_expression.png                  heatmap, datasets in original order
#   - plas_complex_expression_idc.png              heatmap, intraerythrocytic developmental cycle (IDC) datasets first
#   - plas_complex_expression_idc_nosymbols.png    as above, without subunit symbols or title
#   - plas_complex_profiles_clustered.png          mean profile of each complex, complexes clustered
#
# Subunits are labelled with S. cerevisiae nomenclature (best BLASTp hit; see blast.py).


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import groupby

df = pd.read_csv("plas_expression_matrix.tsv", sep="\t", index_col=0)


complexes = {
    "Proteasome": ["PF3D7_0807500", "PF3D7_0608500", "PF3D7_1353800", "PF3D7_1353900", "PF3D7_0727400",
                   "PF3D7_1474800", "PF3D7_0317000", "PF3D7_0931800", "PF3D7_1328100", "PF3D7_0108000",
                   "PF3D7_1470900", "PF3D7_1011400", "PF3D7_0518300", "PF3D7_0803800"],
    "Polymerase": ["PF3D7_0318200", "PF3D7_0215700", "PF3D7_0923000", "PF3D7_0205500",
                   "PF3D7_1364800", "PF3D7_0303300", "PF3D7_1027400", "PF3D7_1213700",
                   "PF3D7_0110400", "PF3D7_0708100", "PF3D7_1304900", "PF3D7_1342700"],
    "CCT": ["PF3D7_1132200", "PF3D7_0306800", "PF3D7_1229500", "PF3D7_1357800",
            "PF3D7_0320300", "PF3D7_0608700", "PF3D7_0308200", "PF3D7_0214000"],
    "MCM": ["PF3D7_1417800", "PF3D7_0527000", "PF3D7_1317100", "PF3D7_1211700", "PF3D7_1355100", "PF3D7_0705400"],
    "ATP": ["PF3D7_1311900", "PF3D7_0406100", "PF3D7_0106100", "PF3D7_1341900", "PF3D7_0934500",
            "PF3D7_1140100", "PF3D7_1323200", "PF3D7_1306600",
            "PF3D7_0806800", "PF3D7_0519200", "PF3D7_1354400", "PF3D7_1464700", "PF3D7_0721900"],
}
symbols_pfal = {
    # Proteasome
    "PF3D7_0807500": "α1", "PF3D7_0608500": "α2", "PF3D7_1353800": "α3", "PF3D7_1353900": "α4",
    "PF3D7_0727400": "α5", "PF3D7_1474800": "α6", "PF3D7_0317000": "α7",
    "PF3D7_0931800": "β1", "PF3D7_1328100": "β2", "PF3D7_0108000": "β3", "PF3D7_1470900": "β4",
    "PF3D7_1011400": "β5", "PF3D7_0518300": "β6", "PF3D7_0803800": "β7",
    # Pol II
    "PF3D7_0318200": "RPB1", "PF3D7_0215700": "RPB2", "PF3D7_0923000": "RPB3", "PF3D7_0205500": "RPB4",
    "PF3D7_1364800": "RPB5", "PF3D7_0303300": "RPB6", "PF3D7_1027400": "RPB7", "PF3D7_1213700": "RPB8",
    "PF3D7_0110400": "RPB9", "PF3D7_0708100": "RPB10", "PF3D7_1304900": "RPB11", "PF3D7_1342700": "RPB12",
    # CCT
    "PF3D7_1132200": "CCT1", "PF3D7_0306800": "CCT2", "PF3D7_1229500": "CCT3", "PF3D7_1357800": "CCT4",
    "PF3D7_0320300": "CCT5", "PF3D7_0608700": "CCT6", "PF3D7_0308200": "CCT7", "PF3D7_0214000": "CCT8",
    # MCM
    "PF3D7_1417800": "MCM2", "PF3D7_0527000": "MCM3", "PF3D7_1317100": "MCM4",
    "PF3D7_1211700": "MCM5", "PF3D7_1355100": "MCM6", "PF3D7_0705400": "MCM7",
    # V-ATPase
    "PF3D7_1311900": "V1-A", "PF3D7_0406100": "V1-B", "PF3D7_0106100": "V1-C", "PF3D7_1341900": "V1-D",
    "PF3D7_0934500": "V1-E", "PF3D7_1140100": "V1-F", "PF3D7_1323200": "V1-G", "PF3D7_1306600": "V1-H",
    "PF3D7_0806800": "V0-a", "PF3D7_0519200": "V0-c′", "PF3D7_1354400": "V0-c″", "PF3D7_1464700": "V0-d",
    "PF3D7_0721900": "V0-e",
}
summary=[]
for name, gene_list in complexes.items():
    found = [g for g in gene_list if g in df.index]
    sd = df.loc[found].std(axis=1)        # each gene's spread across samples

    summary.append({
        "Complex": name,
        "Genes": len(found),
        "Mean expression": df.loc[found].values.mean(),
        "Mean SD": sd.mean(),
        "Min SD": sd.min(),
        "Max SD": sd.max(),
    })

summary = pd.DataFrame(summary).round(2)
summary.to_csv("plas_complex_variance.tsv", sep="\t", index=False)
# get where the complexes start in the list for heatmap splitting
genes = []
sizes = []
starts = []

for name, gene_list in complexes.items():
    found = [g for g in gene_list if g in df.index]
    starts.append(len(genes))
    sizes.append(len(found))
    genes.extend(found)
    print(name, len(found), "of", len(gene_list))

# subtract each gene's mean so the heatmap shows change relative to that gene's average
complex_expr = df.loc[genes]
gene_labels = [symbols_pfal.get(g, g) for g in genes]                 # subunit symbols for the y axis
complex_expr = complex_expr.sub(complex_expr.mean(axis=1), axis=0)

# get dataset assignments: "Plas_12 | hour 2" -> "Plas_12"
def get_dataset(col):
    return col.split(" | ")[0]

col_datasets = [get_dataset(c) for c in complex_expr.columns]

ds_names, ds_starts, ds_sizes = [], [], []
lookup_rows = []
pos = 0
number = 1

for ds, group in groupby(col_datasets):
    n = len(list(group))
    ds_names.append(number)
    ds_starts.append(pos)
    ds_sizes.append(n)

    # the same row in the lookup table: plot number -> file -> its samples
    samples = [c.split(" | ", 1)[1] for c in complex_expr.columns[pos:pos + n]]
    lookup_rows.append({"plot_no": number,
                        "file": ds,
                        "n_samples": n,
                        "samples": "; ".join(samples),
                        "dataset_title": ""})
    pos += n
    number += 1

lookup = pd.DataFrame(lookup_rows)
lookup.to_csv("plas_dataset_lookup.tsv", sep="\t", index=False)
print(lookup[["plot_no", "file", "n_samples"]].to_string(index=False))
print("Repeated files:", lookup["file"][lookup["file"].duplicated()].tolist())   # should be []
fig, ax = plt.subplots(figsize=(12, 8))
sns.heatmap(complex_expr, ax=ax, cmap="RdBu_r", center=0, vmin=-3, vmax=3,
            xticklabels=False, yticklabels=gene_labels,
            cbar_kws={"label": "log2 expression (centred per gene)"})

for name, start, n in zip(complexes, starts, sizes):
    ax.axhline(start, color="black", lw=1.5)
    ax.text(-0.07, start + n / 2, name, ha="right", va="center",
            transform=ax.get_yaxis_transform())                        # complex name left of the symbols

for start in ds_starts[1:]:
    ax.axvline(start, color="black", lw=1.5)
middles = []
for i in range(len(ds_starts)):
    start = ds_starts[i]
    size = ds_sizes[i]
    middle = start + size / 2
    middles.append(middle)

ax.set_xticks(middles)
ax.set_xticklabels(ds_names, rotation=90, fontsize=7)

ax.set_xlabel("Samples (dataset order)")
ax.tick_params(axis="y", labelsize=6, labelrotation=0)
ax.set_ylabel("")
ax.set_title("Plasmodium falciparum complex expression")
plt.tight_layout()
plt.show()
fig.savefig("plas_complex_expression.png", dpi=300)

# heatmap with idc first 
# plot numbers (from plas_dataset_lookup.tsv) of the IDC time-course datasets
idc_ids = [6, 7, 8, 10, 12, 15, 19, 20, 21, 22, 25]
# reorder the lookup rows: IDC first, original order kept within each group
lookup["is_idc"] = lookup["plot_no"].isin(idc_ids)
lookup_idc = lookup.sort_values(["is_idc", "plot_no"], ascending=[False, True])

# columns in the same order as the table rows - get the columns for each dataset and add in order
order = []
for f in lookup_idc["file"]:
    order.extend([c for c in complex_expr.columns if get_dataset(c) == f])
expr_idc = complex_expr[order]

print(lookup_idc)
# labels, sizes and starts come straight from the table
idc_names = lookup_idc["plot_no"].tolist()
idc_sizes = lookup_idc["n_samples"].tolist()
idc_starts = []
pos = 0
for n in idc_sizes:
    idc_starts.append(pos)
    pos += n
fig, ax = plt.subplots(figsize=(12, 8))
sns.heatmap(expr_idc, ax=ax, cmap="RdBu_r", center=0, vmin=-3, vmax=3,
            xticklabels=False, yticklabels=gene_labels,
            cbar_kws={"label": "log2 expression (centred per gene)"})

complex_name_labels = []
for name, start, n in zip(complexes, starts, sizes):
    ax.axhline(start, color="black", lw=1.5)
    complex_name_labels.append(ax.text(-0.07, start + n / 2, name, ha="right", va="center",
                                       transform=ax.get_yaxis_transform()))   # complex name left of the symbols

for start in idc_starts[1:]:
    ax.axvline(start, color="black", lw=1.5)
middles = []
for i in range(len(idc_starts)):
    start = idc_starts[i]
    size = idc_sizes[i]
    middle = start + size / 2
    middles.append(middle)

# IDC boundary: thick black line with a white outline so it stands out on both red and blue (no colour needed)
import matplotlib.patheffects as pe
from matplotlib.patches import Rectangle
from matplotlib.transforms import blended_transform_factory
idc_end = idc_starts[len(idc_ids)]
ax.axvline(idc_end, color="black", lw=3,
           path_effects=[pe.Stroke(linewidth=6, foreground="white"), pe.Normal()])

# labelled bar above the heatmap showing which datasets are IDC and which aren't
trans = blended_transform_factory(ax.transData, ax.transAxes)        # x in samples, y relative to the axes
n_cols = expr_idc.shape[1]
for x0, x1, label, shade, text_colour in [(0, idc_end, "IDC datasets", "#333333", "white"),
                                          (idc_end, n_cols, "Other datasets", "#BBBBBB", "black")]:
    ax.add_patch(Rectangle((x0, 1.01), x1 - x0, 0.045, transform=trans, facecolor=shade,
                           edgecolor="white", linewidth=1.5, clip_on=False))
    ax.text((x0 + x1) / 2, 1.0325, label, transform=trans, ha="center", va="center",
            fontsize=9, fontweight="bold", color=text_colour)
ax.set_xticks(middles)
ax.set_xticklabels(idc_names, rotation=90, fontsize=7)


ax.set_xlabel("Samples (dataset order)")
ax.tick_params(axis="y", labelsize=6, labelrotation=0)
ax.set_ylabel("")
ax.set_title("Plasmodium falciparum complex expression ordered", pad=36)   # room for the IDC bar
plt.tight_layout()
plt.show()
fig.savefig("plas_complex_expression_idc.png", dpi=300)

# same IDC heatmap without the gene symbols and without the title
ax.set_yticks([])                                                    # removes the subunit symbols and their ticks
for label in complex_name_labels:
    label.set_x(-0.01)                                               # complex names move in next to the heatmap
ax.set_title("")
fig.savefig("plas_complex_expression_idc_nosymbols.png", dpi=300, bbox_inches="tight")


# complex clustermap
# get the mean expression across complex members for each sample

cluster = {}
for name, gene_list in complexes.items():
    found = [g for g in gene_list if g in complex_expr.index]
    cluster[name] = complex_expr.loc[found].mean(axis=0)

complex_profiles = pd.DataFrame(cluster).T

# correlation between the complexes' mean profiles
print(complex_profiles.T.corr().round(2))

### cluster the complexes, keep samples in dataset order
g = sns.clustermap(complex_profiles,
                   row_cluster=True,
                   col_cluster=False,
                   metric="correlation",
                   method="average",
                   cmap="RdBu_r", center=0, vmin=-2, vmax=2,
                   xticklabels=False,
                   figsize=(14, 4),
                   dendrogram_ratio=(0.1, 0.05),
                   cbar_pos=(1.01, 0.3, 0.015, 0.4),
                   cbar_kws={"label": "mean centred log2 expression"})

ax = g.ax_heatmap
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)

# dataset boundaries and labels
for start in ds_starts[1:]:
    ax.axvline(start, color="grey", lw=1.5)

middles = []
for i in range(len(ds_starts)):
    middles.append(ds_starts[i] + ds_sizes[i] / 2)
ax.set_xticks(middles)
ax.set_xticklabels(ds_names, rotation=90, fontsize=7)
ax.set_xlabel("Dataset")
ax.set_ylabel("")

g.fig.suptitle("Plasmodium falciparum average complex expression, clustered", y=1.05)
g.savefig("plas_complex_profiles_clustered.png", dpi=300, bbox_inches="tight")
plt.show()
