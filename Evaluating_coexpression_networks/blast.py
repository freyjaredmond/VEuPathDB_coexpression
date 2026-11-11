# BLAST complex subunits against S. cerevisiae
# BLASTp of P. falciparum complex subunits against the S. cerevisiae S288C proteome,
# keeping the top hit (highest bitscore) per gene to confirm subunit identity.
#
# Input:
#   - plasmo_sequences.txt   PlasmoDB protein sequences 
#   - scer_db                BLAST database of the S. cerevisiae S288C proteome built with:
#                            makeblastdb -in FungiDB-71_ScerevisiaeS288C_AnnotatedProteins.fasta -dbtype prot -out scer_db
#
# Output:
#   - plasmo_complex_genes.fasta          protein sequences of the complex subunits
#   - plasmo_vs_scer.tsv                  all BLASTp hits (E-value < 1e-5, up to 5 per gene)
#   - plasmo_complex_top_blast_hits.tsv   top yeast hit for each subunit
#

import pandas as pd
import subprocess

# path to blastp
BLASTP = "blastp"

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

# P. falciparum protein sequences 
seq_df = pd.read_csv("plasmo_sequences.txt", sep="\t")
seq_dict = dict(zip(seq_df["Gene ID"], seq_df["Predicted Protein Sequence"].str.rstrip("*")))

# one row per subunit with its sequence
table = []
for complex_name, genes in plas_complexes.items():
    for g in genes:
        table.append({"Complex": complex_name, "Gene": g, "Sequence": seq_dict.get(g)})
plas_df = pd.DataFrame(table)

# write the sequences to a FASTA file (genes without a sequence are skipped)
with open("plasmo_complex_genes.fasta", "w") as f:
    for idx, row in plas_df.dropna(subset=["Sequence"]).iterrows():
        f.write(f">{row['Gene']}\n{row['Sequence']}\n")

# BLASTp against the yeast proteome
subprocess.run([BLASTP,
                "-query", "plasmo_complex_genes.fasta",
                "-db", "scer_db",
                "-evalue", "1e-5",
                "-max_target_seqs", "5",
                "-outfmt", "6 qseqid sseqid pident length evalue bitscore stitle",
                "-out", "plasmo_vs_scer.tsv"], check=True)

# top hit per gene (highest bitscore)
cols = ["Gene", "Yeast_hit", "Pct_identity", "Align_length", "Evalue", "Bitscore", "Yeast_description"]
hits = pd.read_csv("plasmo_vs_scer.tsv", sep="\t", names=cols)
top = hits.sort_values("Bitscore", ascending=False).drop_duplicates("Gene")

# genes with no hit are kept, with empty BLAST columns
plas_df = plas_df.drop(columns="Sequence").merge(top, on="Gene", how="left")
print(plas_df.to_string(index=False))
plas_df.to_csv("plasmo_complex_top_blast_hits.tsv", sep="\t", index=False)
