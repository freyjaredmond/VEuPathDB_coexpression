# Preprocessing comparisons

An initial comparison was run on *A. fumigatus* and *P. falciparum* to select the best-performing preprocessing for Pearson's correlation coefficient (PCC) and WGCNA. Mutual rank is derived directly from PCC, so it was not assessed separately.

Four preprocessing combinations were explored:
- [P1](Preprocessing_comparison/P1): batch correction and quantile normalisation
- [P2](Preprocessing_comparison/P2): batch correction only
- [P3](Preprocessing_comparison/P3): quantile normalisation only
- [P4](Preprocessing_comparison/P4): neither

The metrics assessed are described in [Evaluating networks](Evaluating_networks).

| Run | Organism | Method | Batch correction | Quantile normalisation | GO AUROC | Degree-null AUROC | Complex AUPRC |
|---|---|---|---|---|---|---|---|
| P1 | *A. fumigatus* | Pearson | Y | Y | 0.75 | 0.60 | 0.69 |
| P1 | *A. fumigatus* | WGCNA | Y | Y | 0.73 | 0.60 | 0.38 |
| P1 | *P. falciparum* | Pearson | Y | Y | 0.72 | 0.54 | 0.70 |
| P1 | *P. falciparum* | WGCNA | Y | Y | 0.72 | 0.56 | 0.65 |
| P2 | *A. fumigatus* | Pearson | Y | N | 0.78 | 0.60 | 0.79 |
| P2 | *A. fumigatus* | WGCNA | Y | N | 0.72 | 0.62 | 0.44 |
| P2 | *P. falciparum* | Pearson | Y | N | 0.72 | 0.58 | 0.70 |
| P2 | *P. falciparum* | WGCNA | Y | N | 0.70 | 0.59 | 0.52 |
| P3 | *A. fumigatus* | Pearson | N | Y | 0.73 | 0.56 | 0.61 |
| P3 | *A. fumigatus* | WGCNA | N | Y | 0.70 | 0.57 | 0.40 |
| P3 | *P. falciparum* | Pearson | N | Y | 0.72 | 0.49 | 0.57 |
| P3 | *P. falciparum* | WGCNA | N | Y | 0.71 | 0.57 | 0.40 |
| P4 | *A. fumigatus* | Pearson | N | N | 0.75 | 0.59 | 0.57 |
| P4 | *A. fumigatus* | WGCNA | N | N | 0.72 | 0.61 | 0.38 |
| P4 | *P. falciparum* | Pearson | N | N | 0.71 | 0.57 | 0.48 |
| P4 | *P. falciparum* | WGCNA | N | N | 0.66 | 0.58 | 0.35 |

Overall, batch correction improved complex precision-recall for both methods. Once batch-corrected, quantile normalisation offered no further advantage for Pearson but improved WGCNA.

Overall, batch correction improved complex precision-recall for both methods. Once batch-corrected, quantile normalisation offered no further advantage for Pearson but improved WGCNA.

