# Propagate GO annotations
# Propagates GO terms to all ancestor terms (is_a and part_of relationships) using GOATOOLS,
# removing the three root terms (BP, MF, CC).
# Input:
#   - go-basic.obo:          GO ontology 
#   - Plasmo_computed.txt:   computed GO annotations from PlasmoDB 
#   - Plasmo_curated.txt:    curated GO annotations from PlasmoDB 
# Output:
#   For EGAD neighbour voting 
#   - plasmo_computed_propogated.tsv        
#
#   For GO term matching 
#   - plasmo_combined_propagated.tsv         
#   - plasmo_combined_propagated_filtered.tsv 

import pandas as pd
from goatools.obo_parser import GODag
from goatools.gosubdag.gosubdag import GoSubDag
 
#parse the go term tree using goatools
godag = GODag("go-basic.obo", optional_attrs={"relationship"})
# looking for is_a and part_of relationships
optional_relationships = {"part_of"}
# remove the root terms (BP, MF and CC)
roots = {"GO:0008150", "GO:0003674", "GO:0005575"}
 
# get the go terms for each gene id
df = pd.read_csv("Plasmo_computed.txt", sep="\t")
df = df.drop_duplicates(subset="Gene ID")
go_cols = [col for col in df.columns if "GO" in col]
 
gene_terms = {}
for _, row in df.iterrows():
    terms = set()
    for col in go_cols:
        if isinstance(row[col], str):
            for t in row[col].split(";"):
                t = t.strip()
                # this gets the updated GO in the godag
                if t.startswith("GO:") and t in godag:
                    terms.add(godag[t].id)
    gene_terms[row["Gene ID"]] = terms
# get all possible go terms
all_terms = set()
for terms in gene_terms.values():
    all_terms.update(terms)
all_terms = sorted(all_terms)
 
# use goatools to get the ancestors
subdag = GoSubDag(all_terms, godag, relationships=optional_relationships, prt=None)
go2anc = subdag.rcntobj.go2ancestors
 
# build long table for EGAD
df_long = []
for gene, terms in gene_terms.items():
    full = set(terms)
    for t in terms:
        full |= set(go2anc.get(t, set()))
    for go in full - roots:
        df_long.append({"Gene ID": gene, "GO": go})
 
df_long = pd.DataFrame(df_long)
df_long.to_csv("plasmo_computed_propogated.tsv", sep="\t", index=False)
 
# propagate computed + curated GO for GO term matching
 
computed_df=pd.read_csv("Plasmo_computed.txt", sep="\t")
curated_df=pd.read_csv("Plasmo_curated.txt", sep="\t")
 
combined_df=computed_df.merge(curated_df, on=["Gene ID","source_id","Product Description"], how="outer")
combined_df = combined_df.drop_duplicates(subset="Gene ID")
 
go_cols_comb = [col for col in combined_df.columns if "GO" in col]
 
go_terms_combined = {}
for _, row in combined_df.iterrows():
    terms = set()
    for col in go_cols_comb:
        if isinstance(row[col], str):
            for t in row[col].split(";"):
                t = t.strip()
                # this gets the updated GO in the godag
                if t.startswith("GO:") and t in godag:
                    terms.add(godag[t].id)
    go_terms_combined[row["Gene ID"]] = terms
combined_terms = set()
for terms in go_terms_combined.values():
    combined_terms.update(terms)
combined_terms = sorted(combined_terms)
 
# use goatools to get the ancestors
subdag = GoSubDag(combined_terms, godag, relationships=optional_relationships, prt=None)
go2anc = subdag.rcntobj.go2ancestors
 
# extract terms
def get_ids(cell):
    ids = []
    if isinstance(cell, str):          # skips empty (NaN) cells
        for t in cell.split(";"):
            t = t.strip()
            if t.startswith("GO:") and t in godag:
                ids.append(godag[t].id)
    return ids
 
# get ancestors
def propagate(terms):
    full = set(terms)
    for t in terms:
        full |= set(go2anc.get(t, set()))
    return full - roots
 
# build wide table for GO matching
df_wide=[]
for idx, rows in combined_df.iterrows():
    gene=rows["Gene ID"]
    mf_gos=[]
    mf_gos += get_ids(rows["Computed GO Function IDs"])
    mf_gos += get_ids(rows["Curated GO Function IDs"])
    mf_gos = propagate(mf_gos)
    cc_gos=[]
    cc_gos += get_ids(rows["Computed GO Component IDs"])
    cc_gos += get_ids(rows["Curated GO Component IDs"])
    cc_gos = propagate(cc_gos)
    bp_gos=[]
    bp_gos += get_ids(rows["Computed GO Process IDs"])
    bp_gos += get_ids(rows["Curated GO Process IDs"])
    bp_gos = propagate(bp_gos)
    df_wide.append({"Gene":gene,
    "MF terms": ";".join(sorted(mf_gos)),
    "CC terms": ";".join(sorted(cc_gos)),
    "BP terms": ";".join(sorted(bp_gos))})
 
df_wide = pd.DataFrame(df_wide)
df_wide.to_csv("plasmo_combined_propagated.tsv", sep="\t", index=False)
 
### cut off
df_wide_filtered = df_wide.copy()
for term in ["MF","CC","BP"]:
    # split the joined strings back into one term per row before counting
    vc=df_wide[f"{term} terms"].str.split(";").explode().value_counts()
    vc=vc.reset_index()
    vc.columns=["GO","Count"]
    vc_filtered=vc[vc["Count"]<=300]
    vc_keep=set(vc_filtered["GO"].to_list())
 
    df_wide_filtered[f"{term} terms"] = df_wide[f"{term} terms"].apply(
        lambda s: ";".join(t for t in s.split(";") if t in vc_keep))
df_wide_filtered.to_csv("plasmo_combined_propagated_filtered.tsv", sep="\t", index=False)
