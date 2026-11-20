# Subunit x subunit MR rank matrices for P. falciparum complexes
# For each of the five complexes (20S proteasome, RNA polymerase II, CCT, MCM2-7, V-ATPase), plots the MR rank of
# every other subunit (columns) within each subunit's co-expression list (rows), binned into 1-10, 11-20, 21-50,
# 51-100, 101-200 and >200. Genes not in the network are shown separately.
#
# Input:
#   - pfal_mr_ranks_complexes.tsv     MR rank of each complex partner in each subunit's co-expression list
#                                     (columns include Gene, Partner, MR_rank)
#
# Output:
#   - complexes_mr_rank_matrix.png    all five complexes in one figure
#   - complexes_mr_rank_matrix.svg    vector version of the same figure
#
# Subunits are labelled with S. cerevisiae nomenclature (best BLASTp hit; see blast.py).


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch

mr=pd.read_csv("pfal_mr_ranks_complexes.tsv", sep="\t")

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

# subunit symbols (S. cerevisiae nomenclature)
symbols = {
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

mr_genes=set(mr["Gene"].to_list())            # genes in the rank file

# rank bins
bins=[0, 10, 20, 50, 100, 200, np.inf]
bin_labels=["1-10", "11-20", "21-50", "51-100", "101-200", ">200"]

colours=["#00441b", "#1b7837", "#5aae61", "#b2dd18", "#db8d17", "#db4517", "salmon"]

# sized to a 16:9 slide (13.33 in wide) so fonts appear at their real size when placed full-width
fig, axes=plt.subplots(1, len(complexes), figsize=(13.33, 3.1), layout="constrained")   # all five in one figure, spacing recalculated so labels don't overlap
fig.get_layout_engine().set(wspace=0.02, h_pad=0.15)                   # h_pad: extra gap between tick labels and the shared x label

for ax, (complex, genes) in zip(axes, complexes.items()):

    sub=mr[mr["Gene"].isin(genes) & mr["Partner"].isin(genes)]
    ranks=sub.pivot(index="Gene", columns="Partner", values="MR_rank").reindex(index=genes, columns=genes)

    mat=ranks.apply(lambda col: pd.cut(col, bins=bins, labels=False))   # rank -> bin number (0-5)

    missing=[g for g in genes if g not in mr_genes]                     # genes not in the mr at all
    mat.loc[missing, :]=len(bin_labels)                                 # salmon
    mat.loc[:, missing]=len(bin_labels)

    mat=mat.mask(np.eye(len(genes), dtype=bool))                        # diagonal shows as grey background
    mat = mat.rename(index=symbols, columns=symbols)
    sns.heatmap(mat, cmap=colours, vmin=-0.5, vmax=len(colours) - 0.5, cbar=False,
                linewidths=0.5, linecolor="lightgrey", square=True, ax=ax,
                xticklabels=True, yticklabels=True)                     # label every gene, not every other one
    ax.set_facecolor("lightgrey")
    ax.set_title(f"{complex}", fontsize=12, fontweight="bold")
    ax.set_xlabel("")                                                   # shared axis labels added once below
    ax.set_ylabel("")
    ax.tick_params(labelsize=8, length=0)
    ax.tick_params(axis="y", labelrotation=0)                           # horizontal gene labels on every panel
    ax.tick_params(axis="x", labelrotation=90)                          # vertical partner labels on every panel

legend=[Patch(facecolor=c, edgecolor="lightgrey", label=l) for c, l in zip(colours, bin_labels + ["Not in network"])]
fig.supxlabel("Complex partner", fontsize=11)
fig.supylabel("Gene", fontsize=11)
fig.legend(handles=legend, title="MR rank", loc="outside right center", ncol=1, frameon=False,
           fontsize=10, title_fontsize=11)

fig.savefig("complexes_mr_rank_matrix.png", dpi=300)
             
plt.show()
