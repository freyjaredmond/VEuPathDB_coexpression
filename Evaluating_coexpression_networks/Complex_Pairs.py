# Complex pairs
# Builds every pairwise combination of complex subunits for the complex-recovery AUPRC.
# Pairs of subunits from the same complex are positives (Label = 1) and pairs from
# different complexes are negatives (Label = 0).
#
# Input:
#   - complexes dictionary below: P. falciparum subunits of five conserved complexes,
#     identified from annotation and BLASTp to S. cerevisiae (see blast.py)
#
# Output:
#   - plasmo_complex_pairs.tsv   Gene1, Gene2, Complex1, Complex2, Label

import pandas as pd
from itertools import combinations

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

# look up which complex each gene belongs to
gene_to_complex = {}
for complex_name, members in plas_complexes.items():
    for gene in members:
        gene_to_complex[gene] = complex_name

# every unique pair of complex genes, labelled 1 if both are in the same complex, otherwise 0
plas_pairs = []
for gene1, gene2 in combinations(sorted(gene_to_complex), 2):
    complex1 = gene_to_complex[gene1]
    complex2 = gene_to_complex[gene2]
    plas_pairs.append({"Gene1": gene1,
                       "Gene2": gene2,
                       "Complex1": complex1,
                       "Complex2": complex2,
                       "Label": 1 if complex1 == complex2 else 0})

plas_pairs = pd.DataFrame(plas_pairs)
plas_pairs.to_csv("plasmo_complex_pairs.tsv", sep="\t", index=False)
