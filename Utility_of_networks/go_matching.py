# GO term matching of top co-expressed partners
# For every annotated P. falciparum gene, calculates the percentage of its top 10, 20 and 50
# co-expressed partners (Pearson, MR, WGCNA TOM, WGCNA TOM within module) that share at least
# one propagated GO term, separately for BP, MF and CC, and compares against a random baseline.
# Partners are only counted if both genes are annotated in that ontology; grey-module genes
# score 0% for the module network.
# Methods are compared with paired t-tests (one-sided vs random, two-sided between methods,
# Benjamini-Hochberg corrected) and plotted as mean bars with interquartile ranges.
#
# Input:
#   - plasmo_combined_propagated_filtered.tsv   propagated computed + curated GO (see propagate.py)
#   - plas_pearson_top50.tsv, plas_mr_top50.tsv, plas_wgcna_top50.tsv, plas_modules_top50.tsv
#                                               top 50 partners per gene (see build_coexpression_networks.R)
#
# Output:
#   - plasmo_methods_go_match_long_propogated.tsv   per-gene percentage match for every method,
#                                                   cutoff and ontology (including random)
#   - methods_ttest_paired_propogated.tsv           paired t-test results with 95% CI of the mean difference
#   - mean_bar_iqr_methods.png                      mean percentage match with interquartile range


import pandas as pd
import os
import math
import random
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

### build the GO dictionaries once, before looping over the networks

## propagated GO (computed + curated), one dictionary per ontology
go_prop=pd.read_csv("plasmo_combined_propagated_filtered.tsv", sep="\t")

mf_dict={}
cc_dict={}
bp_dict={}

##mf
for idx, row in go_prop.iterrows():
    mf_terms=[]
    if pd.notna(row["MF terms"]):
        mf_terms += [mf.strip() for mf in row["MF terms"].split(";")]
    mf_dict[row["Gene"]]=list(set(mf_terms))

##cc
for idx, row in go_prop.iterrows():
    cc_terms=[]
    if pd.notna(row["CC terms"]):
        cc_terms += [cc.strip() for cc in row["CC terms"].split(";")]
    cc_dict[row["Gene"]]=list(set(cc_terms))

##bp
for idx, row in go_prop.iterrows():
    bp_terms=[]
    if pd.notna(row["BP terms"]):
        bp_terms += [bp.strip() for bp in row["BP terms"].split(";")]
    bp_dict[row["Gene"]]=list(set(bp_terms))

go_type="propogated"

networks=["plas_pearson_top50.tsv",
          "plas_wgcna_top50.tsv",
          "plas_mr_top50.tsv",
          "plas_modules_top50.tsv"]

# every results table from every network, stacked at the end for the comparisons and plot
all_results=[]

for ds in networks:
    # network name from the file name, e.g. plas_pearson_top50.tsv -> pearson
    dataset=os.path.basename(ds).split("_")[1]
    tom=pd.read_csv(ds, sep="\t")

    ### extract the GO terms for the network
    ## add the cut offs (top 10, 20 and 50 will be taken but excluding from percentage match if any of them don't have go terms)
    # genes not in the network have no partners (na for coexpressed)
    if "module" in tom.columns:
        grey_genes=tom.loc[tom["module"]=="grey", "gene"].unique()
    else:
        grey_genes=[]
    tom=tom[tom["correlated_gene"].notna()]

    tom_grouped=tom.groupby("gene")

    for cutoff in [10, 20, 50]:
        cc_results=[]
        mf_results=[]
        bp_results=[]

        for gene, group in tom_grouped:
            goi_go_cc=set(cc_dict.get(gene, []))
            goi_go_mf=set(mf_dict.get(gene, []))
            goi_go_bp=set(bp_dict.get(gene, []))

            cc_matches=[]
            cc_no_match=[]
            mf_matches=[]
            mf_no_match=[]
            bp_matches=[]
            bp_no_match=[]

            # partners are sorted strongest to weakest (highest r/TOM first, lowest MR first), so head gives the top n
            for coexp in group["correlated_gene"].head(cutoff):
                coexp_cc=set(cc_dict.get(coexp, []))
                coexp_mf=set(mf_dict.get(coexp, []))
                coexp_bp=set(bp_dict.get(coexp, []))

                # only compare when both the GOI and the partner have terms in that ontology
                if goi_go_cc and coexp_cc:
                    if goi_go_cc & coexp_cc:
                        cc_matches.append(coexp)
                    else:
                        cc_no_match.append(coexp)

                if goi_go_mf and coexp_mf:
                    if goi_go_mf & coexp_mf:
                        mf_matches.append(coexp)
                    else:
                        mf_no_match.append(coexp)

                if goi_go_bp and coexp_bp:
                    if goi_go_bp & coexp_bp:
                        bp_matches.append(coexp)
                    else:
                        bp_no_match.append(coexp)

            cc_total=len(cc_matches)+len(cc_no_match)
            #handle cases where the goi had no corresponding terms - or none of the coexpressed partners do
            if cc_total > 0:
                cc_results.append({"gene": gene, "ontology": "CC", "dataset": dataset, "cutoff": cutoff, "total": cc_total, "n_match": len(cc_matches),
                                   "percent_match": round(len(cc_matches) / cc_total * 100, 1)})

            mf_total=len(mf_matches)+len(mf_no_match)
            if mf_total > 0:
                mf_results.append({"gene": gene, "ontology": "MF", "dataset": dataset, "cutoff": cutoff, "total": mf_total, "n_match": len(mf_matches),
                                   "percent_match": round(len(mf_matches) / mf_total * 100, 1)})

            bp_total=len(bp_matches)+len(bp_no_match)
            if bp_total > 0:
                bp_results.append({"gene": gene, "ontology": "BP", "dataset": dataset, "cutoff": cutoff, "total": bp_total, "n_match": len(bp_matches),
                                   "percent_match": round(len(bp_matches) / bp_total * 100, 1)})

        # penalise grey genes with 0% in each ontology the gene has terms for (total 0 marks them as penalties)
        for gene in grey_genes:
            if cc_dict.get(gene):
                cc_results.append({"gene": gene, "ontology": "CC", "dataset": dataset, "cutoff": cutoff, "total": 0, "n_match": 0, "percent_match": 0.0})
            if mf_dict.get(gene):
                mf_results.append({"gene": gene, "ontology": "MF", "dataset": dataset, "cutoff": cutoff, "total": 0, "n_match": 0, "percent_match": 0.0})
            if bp_dict.get(gene):
                bp_results.append({"gene": gene, "ontology": "BP", "dataset": dataset, "cutoff": cutoff, "total": 0, "n_match": 0, "percent_match": 0.0})

        cc_results=pd.DataFrame(cc_results)
        mf_results=pd.DataFrame(mf_results)
        bp_results=pd.DataFrame(bp_results)
        all_results += [cc_results, mf_results, bp_results]

### generate the random control
# pearson, wgcna and mr contain the same genes with the same number of partners, so one random control serves all three
random.seed(28)
tom=pd.read_csv(networks[0], sep="\t")
tom=tom[tom["correlated_gene"].notna()]
tom_grouped=tom.groupby("gene")

# all genes in the network with partners
tom_genes=list(tom_grouped.groups.keys())
for cutoff in [10, 20, 50]:
    cc_results_random_gene=[]
    mf_results_random_gene=[]
    bp_results_random_gene=[]

    for gene, group in tom_grouped:
        goi_go_cc=set(cc_dict.get(gene, []))
        goi_go_mf=set(mf_dict.get(gene, []))
        goi_go_bp=set(bp_dict.get(gene, []))
        cc_matches_all=[]
        cc_no_match_all=[]
        mf_matches_all=[]
        mf_no_match_all=[]
        bp_matches_all=[]
        bp_no_match_all=[]
        # all other network genes
        all_genes=[g for g in tom_genes if g != gene]
        # sample the same number of partners the gene actually has at this cutoff
        n=len(group["correlated_gene"].head(cutoff))
        random_genes_all=random.sample(all_genes, n)

        for coexp in random_genes_all:
            coexp_cc=set(cc_dict.get(coexp, []))
            coexp_mf=set(mf_dict.get(coexp, []))
            coexp_bp=set(bp_dict.get(coexp, []))

            if goi_go_cc and coexp_cc:
                if goi_go_cc & coexp_cc:
                    cc_matches_all.append(coexp)
                else:
                    cc_no_match_all.append(coexp)

            if goi_go_mf and coexp_mf:
                if goi_go_mf & coexp_mf:
                    mf_matches_all.append(coexp)
                else:
                    mf_no_match_all.append(coexp)

            if goi_go_bp and coexp_bp:
                if goi_go_bp & coexp_bp:
                    bp_matches_all.append(coexp)
                else:
                    bp_no_match_all.append(coexp)

        ## random within network results
        cc_total_all=len(cc_matches_all)+len(cc_no_match_all)
        if cc_total_all > 0:
            cc_results_random_gene.append({"gene": gene, "ontology": "CC", "dataset": "random", "cutoff": cutoff, "total": cc_total_all, "n_match": len(cc_matches_all),
                                           "percent_match": round(len(cc_matches_all) / cc_total_all * 100, 1)})

        mf_total_all=len(mf_matches_all)+len(mf_no_match_all)
        if mf_total_all > 0:
            mf_results_random_gene.append({"gene": gene, "ontology": "MF", "dataset": "random", "cutoff": cutoff, "total": mf_total_all, "n_match": len(mf_matches_all),
                                           "percent_match": round(len(mf_matches_all) / mf_total_all * 100, 1)})

        bp_total_all=len(bp_matches_all)+len(bp_no_match_all)
        if bp_total_all > 0:
            bp_results_random_gene.append({"gene": gene, "ontology": "BP", "dataset": "random", "cutoff": cutoff, "total": bp_total_all, "n_match": len(bp_matches_all),
                                           "percent_match": round(len(bp_matches_all) / bp_total_all * 100, 1)})

    cc_results_random_gene=pd.DataFrame(cc_results_random_gene)
    mf_results_random_gene=pd.DataFrame(mf_results_random_gene)
    bp_results_random_gene=pd.DataFrame(bp_results_random_gene)
    all_results += [cc_results_random_gene, mf_results_random_gene, bp_results_random_gene]

### combine for plotting (all networks together)
go_match_long=pd.concat(all_results, ignore_index=True)
go_match_long.to_csv(f"plasmo_methods_go_match_long_{go_type}.tsv", sep="\t", index=False)

### paired t-test: tests whether the mean per-gene difference is 0
# each method vs random (one-sided: method > random), and each pair of methods (two-sided)
comparisons=[("pearson", "random", "greater"),
             ("wgcna", "random", "greater"),
             ("mr", "random", "greater"),
             ("modules", "random", "greater"),
             ("pearson", "wgcna", "two-sided"),
             ("pearson", "mr", "two-sided"),
             ("wgcna", "mr", "two-sided"),
             ("wgcna", "modules", "two-sided"),
             ("modules", "pearson", "two-sided"),
             ("modules", "mr", "two-sided")]
ttest_rows=[]
for ontology in ["BP", "MF", "CC"]:
    for cutoff in [10, 20, 50]:
        subset=go_match_long[(go_match_long["ontology"]==ontology) & (go_match_long["cutoff"]==cutoff)]
        wide=subset.pivot(index="gene", columns="dataset", values="percent_match")

        for group1, group2, alternative in comparisons:
            pair=wide[[group1, group2]].dropna()          # genes scored in both, so the values are paired
            diff=pair[group1]-pair[group2]
            n=len(pair)
            res=stats.ttest_rel(pair[group1], pair[group2], alternative=alternative)
            # two-sided 95% CI of the mean difference (reported for every comparison so they're comparable)
            se=diff.std(ddof=1)/n**0.5
            t_crit=stats.t.ppf(0.975, n-1)
            ttest_rows.append({"ontology": ontology, "cutoff": cutoff, "group1": group1, "group2": group2,
                               "alternative": alternative, "n_genes": n,
                               "mean_group1": round(pair[group1].mean(), 2),
                               "mean_group2": round(pair[group2].mean(), 2),
                               "mean_diff": round(diff.mean(), 2),
                               "ci95_low": round(diff.mean()-t_crit*se, 2),
                               "ci95_high": round(diff.mean()+t_crit*se, 2),
                               "t": round(res.statistic, 2), "df": n-1, "p": res.pvalue})

ttest_results=pd.DataFrame(ttest_rows)
# correct for all 90 tests (Benjamini-Hochberg)
ttest_results["p_adj"]=stats.false_discovery_control(ttest_results["p"])
ttest_results.to_csv(f"methods_ttest_paired_{go_type}.tsv", sep="\t", index=False)

### plot
# colour-blind safe (Okabe-Ito) palette
palette={"wgcna": "#0072B2",     # blue
         "modules": "#E69F00",   # orange
         "pearson": "#CC79A7",   # reddish purple
         "mr": "#009E73",        # bluish green
         "random": "#999999"}    # grey
hue_order=["wgcna", "modules", "pearson", "mr", "random"]
# display names for the legends (the data keeps the lower-case names)
method_labels={"wgcna": "WGCNA", "modules": "Modules", "pearson": "Pearson", "mr": "MR", "random": "Random"}
def tidy_legend(leg):
    leg.set_title("Method")
    for text in leg.get_texts():
        text.set_text(method_labels.get(text.get_text(), text.get_text()))

# mean as a bar, with the interquartile range (25th-75th percentile) of per-gene values as the error bar
# paired t-test results (methods_ttest_paired_*.tsv) are reported in the caption/supplement
def iqr(x):
    return (x.quantile(0.25), x.quantile(0.75))

# y axis starts at 0 but stops just above the tallest error bar (rounded up to the next 10), not at 100
q75_max = go_match_long.groupby(["ontology", "cutoff", "dataset"])["percent_match"].quantile(0.75).max()
y_top = math.ceil((q75_max + 5) / 10) * 10

fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
for ax, ontology in zip(axes, ["BP", "MF", "CC"]):
    sub = go_match_long[go_match_long["ontology"] == ontology]
    sns.barplot(data=sub, x="cutoff", y="percent_match", hue="dataset",
                hue_order=hue_order, palette=palette, order=[10, 20, 50],
                estimator="mean", errorbar=iqr, capsize=0.15, saturation=1,   # saturation=1: exact palette colours
                err_kws={"color": "black", "linewidth": 1},
                legend=(ontology == "CC"), ax=ax)
    ax.set_title(ontology, fontsize=15)
    ax.set_xlabel("Top n coexpressed partners")
    ax.set_ylim(0, y_top)
axes[0].set_ylabel("Mean partners sharing a GO term (%)")
sns.move_legend(axes[2], "upper left", bbox_to_anchor=(1.02, 1), frameon=False, title="Method")
tidy_legend(axes[2].get_legend())
fig.savefig("mean_bar_iqr_methods.png", dpi=300, bbox_inches="tight")
plt.show()
