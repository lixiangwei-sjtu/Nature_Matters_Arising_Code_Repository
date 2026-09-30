# Adult human amino-acid intake and the Drosophila Lsp2 longevity signal

Code repository accompanying the Nature Matters Arising manuscript:

**Adult human amino-acid intake reveals a translational boundary for the Drosophila Lsp2 longevity signal**

## Scope

This repository contains reproducibility code for the derived amino-acid exposures, the albumin-informed compositional metrics, survey-model specifications, and the main figure. It does **not** redistribute participant-level HRS data or any other restricted-use data.

The human analyses address a translational question distinct from the original Drosophila developmental mechanism. HRS captures habitual diet in later adulthood and uses DunedinPACE; NHANES uses a 24-hour dietary recall and DunedinPoAm. These measures are related but not interchangeable.

## Repository structure

```
code/
  01_build_amino_acid_exposures.py
  02_albumin_profile_similarity.py
  03_analysis_models.R
  04_figure1.py
results/
  locked_main_results.csv
README.md
METHODS_NOTES.md
requirements.txt
LICENSE
CITATION.cff
```

## Data availability

- **HRS**: participant-level data are not redistributed. Authorized users should obtain the Health and Retirement Study / Health Care and Nutrition Study / Venous Blood Study data under HRS data-use conditions.
- **NHANES**: public-use 2001–2002 dietary and DNA-methylation resources are available from the US Centers for Disease Control and Prevention.
- **USDA food composition**: FNDDS 1.0 and USDA Standard Reference 16-1 are public resources used to reconstruct amino-acid intake.

No participant-level data are included in this repository.

## Main derived exposures

### Phe+Tyr
Phenylalanine and tyrosine are summed after participant-level reconstruction and then standardized for regression. Energy adjustment should follow the analytic specification used in the manuscript.

### Albumin-profile similarity
The locked definition is cosine similarity between:

1. each participant's 18-amino-acid dietary **mole-fraction** profile, and
2. the corresponding 18-bin amino-acid residue-frequency profile of mature human serum albumin (585 residues).

Dietary g/day values are converted to molar residue equivalents. Dietary cystine contributes two cysteine-residue equivalents. The resulting similarity is standardized to 1 SD before regression.

### Albumin-enriched score
The exploratory score combines Cys, Lys, Phe, Tyr, Ala, Glu and His, weighted by their relative enrichment in mature human serum albumin versus average human-proteome amino-acid frequencies. Because the manuscript analysis used a locked set of enrichment weights, users should provide those weights in the input/configuration used with `01_build_amino_acid_exposures.py`.

## Statistical models

The repository includes model specifications rather than restricted data. The HRS specification uses survey-weighted linear regression and supports Rubin pooling across multiple imputations. The NHANES specification uses the DNAm subsample weight with strata and primary sampling units.

The primary cross-cohort comparison is deliberately interpretive rather than a strict replication because HRS uses DunedinPACE whereas NHANES uses DunedinPoAm.

## Reproducing Figure 1

`code/04_figure1.py` reproduces the manuscript forest plot from the locked, non-disclosive summary estimates in `results/locked_main_results.csv`.

## Important note on the NHANES albumin-informed confidence intervals

For the current manuscript figure, the two NHANES albumin-informed 95% confidence intervals were reconstructed from the locked beta and P values using the survey inference degrees of freedom (df=15). They should be replaced by the exact model-exported intervals if the final participant-level model output is rerun before publication.

## Environment

Python dependencies are listed in `requirements.txt`. The survey-model script requires R with the `survey` and `mitools` packages.

## Contact

Xiangwei Li, Shanghai Jiao Tong University School of Medicine.
