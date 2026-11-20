# Expression variability of complexes in P. falciparum
# Plots the standard deviation of expression across all samples for the members of each complex
# (20S proteasome, RNA polymerase II, CCT, MCM2-7, V-ATPase), with all genes for comparison.
#
# Input:
#   - plas_expression_matrix.tsv     batch-corrected log2(TPM + 1) matrix (genes x samples)
#
# Output:
#   - complex_std_dis.png            boxplot of per-gene expression SD per complex, plus all genes

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# complex members
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

expr = pd.read_csv("plas_expression_matrix.tsv", sep="\t", index_col=0)
sd = expr.std(axis=1)                                          # each gene's SD across all samples

rows = []
for v in sd:
    rows.append({"Complex": "All genes", "SD": v})
for cx, genes in complexes.items():
    rows += [{"Complex": cx, "SD": sd[g]} for g in genes if g in sd.index]
sd_df = pd.DataFrame(rows)
fig, ax = plt.subplots(figsize=(6, 4))
sns.boxplot(data=sd_df, x="Complex", y="SD", color="lightgrey", showfliers=False, ax=ax)
plt.ylabel("Expression SD across samples")
ax.set_xlabel("")
ax.tick_params(axis="x", labelrotation=30)
plt.setp(ax.get_xticklabels(), ha="right", rotation_mode="anchor")                     
plt.tight_layout()
plt.show()
fig.savefig("complex_std_dis.png", dpi=300)
