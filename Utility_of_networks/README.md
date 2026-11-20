# Exploring the utility of co-expression networks

The top 10, 20 and 50 co-expressed partners of each gene were extracted from four networks: PCC, MR, WGCNA TOM and WGCNA TOM restricted to the gene's module. Networks were built with [Build_network.R](../Evaluating_coexpress_networks/Build_network.R) and are deposited at X.

## GO matching

For every gene with at least one GO annotation in the ontology of interest, the percentage of partners sharing at least one propagated GO term (computed and curated combined; GO terms with ≤300 genes) was calculated for each cutoff and method using [go_matching.py](go_matching.py). For WGCNA restricted to modules, genes in the grey (unassigned) module were scored as a negative match.

At all sample sizes, GO matching remained above random, with a general trend of MR and PCC performing better than WGCNA TOM and WGCNA modules.

## Complex recovery

For each subunit of the five complexes, the percentage of the other subunits in that complex found within its top 10, 20 and 50 co-expressed genes was calculated and averaged across subunits to give a score for each complex using [complex_recovery.py](complex_recovery.py).

PCC and MR generally performed better than WGCNA. Recovery declined as sample size decreased, and RNA polymerase II and the V-ATPase were poorly recovered in all genomes.

![GO matching and complex recovery](fig2.png)

## Why were some complexes poorly recovered?

- **Expression patterns:** mean-centred, batch-corrected log2(TPM + 1) expression was plotted for each complex member using [Complex expression](complex_expression.py). Clear bands of expression were seen for the proteasome, CCT and MCM2–7. RNA polymerase II also showed bands, but these mirrored CCT, while the V-ATPase showed no clear bands.
- **Expression variability:** the standard deviation of expression was plotted for each complex using [standard deviation](std.py). RNA polymerase II and the V-ATPase consistently had low SD, although the *A. fumigatus* proteasome had the lowest SD despite high recovery.
- **Subunit × subunit MR ranks:** a matrix of each subunit's MR rank for every other member of its complex was produced using [MR matrix](mr_rank_matrix.py). In RNA polymerase II, RPB1 and RPB2 consistently ranked each other highly, but no other clear pattern was seen for RNA polymerase II or the V-ATPase.

![Complex expression, variability and MR ranks](fig3.png)
