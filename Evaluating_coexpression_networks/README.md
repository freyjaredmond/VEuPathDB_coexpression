# Evaluating the co-expression networks
Neighbour voting and complex-recovery AUPRC was run on seven organisms of varying RNA-seq coverage using [Evaluating.R](Evaluating.R). PCC, PCC ranled, WGCNA, wgcna
binary module memebership and MR was explored

Organism 

Strain 

Datasets (n) 

Samples (n) 

Aspergillus fumigatus 

Af293 

15 

85 

Aspergillus nidulans 

FGSC A4 

6 

22 

Anopheles gambiae 

PEST 

17 

101 

Candida auris 

B8441 

6 

42 

Plasmodium berghei 

ANKA 

7 

69 

Plasmodium falciparum 

3D7 

27 

285 

Toxoplasma gondii 

ME49 

22 

Neighbour voting requires propagated GO terms that were generated with [Propagate.py](Propogate.py) and complex-recovery requires complex members
which were identifed using cell compartment GO terms to identify candidates and where annotation was missing or unclear, were identified using S. cervisiae BLASR
matching [blast.py](blast.py). Positive and negtaive pairs for AUPRC were generated using [Complex_Pairs.py](Complex_Pairs.py).

The results of Neighbour voting and AUPRC are given in each organisms folder i.e [AF](AF).
This includes:
- Neighbour voting results [egad_results](AF/egad_results_asper.tsv)
- Complex AUPRC results    [complex_results](AF/complex_results_asper.tsv)
- WGCNA soft power         [soft power](AF/asper_soft_power.tsv)
- BLAST hits     [blast hits](AF/afum_complex_top_blast_hits.tsv)
- Complex pairs [complex pais](AF/asper_complex_pairs.tsv)
- Propagated GO terms [Propagated GO terms](AF/asper_computed_propagted)

Overall Neighbour voting was very similar acorss methods, but MR outperformed others in complex recovery- WGCNA was poor here
# Subsampling

Subsampling was run in increasing intervals of 2 datasets for 5 random combinations to determine the optimal sample/dataset size.
This was also run with half the number of samples.
To generate the draws used, the following scripts were used [full](sample_cutoffs.py) and [half](sample_cutoffs_half.py) which were then 
run using [subsample.R](subsample.R). Results were explored with [explore](explore_subsample.py) and plotted with [plot](plot_subsample.py)
[results](Subsampling/cutoff_bins_all_genomes.tsv).

Overall sample size was determined to be the major driver as same samples sizes with more datasets did not generally improve the results.
Asper fum was found to be best performing network and T. gondii the worst.

Using the dataset size corresponding with the sample cut offs we get a lenient cut off [lenient](!Subsampling/lenient_rna_organism.tsv] and a strict [strict](!Subsampling/strict_rna_organism.tsv)


122 

