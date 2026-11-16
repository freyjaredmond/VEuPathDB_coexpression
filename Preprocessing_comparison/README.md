# Preprocessing comparisons

An initial comparison was run on *A. fumigatus*, *A. gambiae* ,*P. falciparum* and *T.gondii* to select the best-performing preprocessing for Pearson's correlation coefficient (PCC) and WGCNA.

Four preprocessing combinations were explored:
- [P1](P1): batch correction and quantile normalisation
- [P2](P2): batch correction only
- [P3](P3): quantile normalisation only
- [P4](P4): neither

The metrics assessed are described in [Evaluating networks](Evaluating_networks).

| Run | Organism | Method | Batch correction | Quantile normalisation | GO AUROC | Degree-null AUROC | Complex AUPRC |
|---|---|---|---|---|---|---|---|
| P1 | *A. fumigatus* | Pearson | Y | Y | 0.75 | 0.60 | 0.69 |
| P1 | *A. fumigatus* | WGCNA | Y | Y | 0.73 | 0.60 | 0.38 |
| P1 | *A. gambiae* | Pearson | Y | Y | 0.77 | 0.52 | 0.70 |
| P1 | *A. gambiae* | WGCNA | Y | Y | 0.76 | 0.53 | 0.61 |
| P1 | *P. falciparum* | Pearson | Y | Y | 0.72 | 0.54 | 0.71 |
| P1 | *P. falciparum* | WGCNA | Y | Y | 0.72 | 0.56 | 0.65 |
| P1 | *T. gondii* | Pearson | Y | Y | 0.73 | 0.57 | 0.47 |
| P1 | *T. gondii* | WGCNA | Y | Y | 0.72 | 0.57 | 0.37 |
| P2 | *A. fumigatus* | Pearson | Y | N | 0.78 | 0.60 | 0.78 |
| P2 | *A. fumigatus* | WGCNA | Y | N | 0.72 | 0.62 | 0.44 |
| P2 | *A. gambiae* | Pearson | Y | N | 0.77 | 0.57 | 0.66 |
| P2 | *A. gambiae* | WGCNA | Y | N | 0.74 | 0.55 | 0.33 |
| P2 | *P. falciparum* | Pearson | Y | N | 0.72 | 0.58 | 0.70 |
| P2 | *P. falciparum* | WGCNA | Y | N | 0.70 | 0.59 | 0.52 |
| P2 | *T. gondii* | Pearson | Y | N | 0.73 | 0.59 | 0.51 |
| P2 | *T. gondii* | WGCNA | Y | N | 0.69 | 0.56 | 0.33 |
| P3 | *A. fumigatus* | Pearson | N | Y | 0.73 | 0.56 | 0.61 |
| P3 | *A. fumigatus* | WGCNA | N | Y | 0.71 | 0.57 | 0.40 |
| P3 | *A. gambiae* | Pearson | N | Y | 0.76 | 0.52 | 0.69 |
| P3 | *A. gambiae* | WGCNA | N | Y | 0.76 | 0.50 | 0.59 |
| P3 | *P. falciparum* | Pearson | N | Y | 0.72 | 0.49 | 0.57 |
| P3 | *P. falciparum* | WGCNA | N | Y | 0.71 | 0.57 | 0.40 |
| P3 | *T. gondii* | Pearson | N | Y | 0.71 | 0.57 | 0.36 |
| P3 | *T. gondii* | WGCNA | N | Y | 0.70 | 0.48 | 0.28 |
| P4 | *A. fumigatus* | Pearson | N | N | 0.75 | 0.59 | 0.57 |
| P4 | *A. fumigatus* | WGCNA | N | N | 0.72 | 0.61 | 0.38 |
| P4 | *A. gambiae* | Pearson | N | N | 0.76 | 0.58 | 0.65 |
| P4 | *A. gambiae* | WGCNA | N | N | 0.48 | 0.57 | 0.20 |
| P4 | *P. falciparum* | Pearson | N | N | 0.71 | 0.57 | 0.48 |
| P4 | *P. falciparum* | WGCNA | N | N | 0.66 | 0.58 | 0.35 |
| P4 | *T. gondii* | Pearson | N | N | 0.71 | 0.62 | 0.39 |
| P4 | *T. gondii* | WGCNA | N | N | 0.39 | 0.61 | 0.20 |

Overall, batch correction improved complex precision-recall for both methods. Once batch-corrected, quantile normalisation offered no further advantage for Pearson but improved WGCNA.
