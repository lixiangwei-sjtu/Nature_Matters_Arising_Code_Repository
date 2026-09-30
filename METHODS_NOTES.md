# Methods notes

## Biological framing

The Drosophila experiments manipulate amino-acid availability during larval development and identify an Lsp2-linked developmental mechanism affecting later-life translation and lifespan. The human analyses test whether a corresponding amino-acid pattern is detectable in adult dietary exposure and DNA-methylation measures of ageing. They do not directly test the fly pathway.

## HRS

- Dietary exposure: HRS 2013 Health Care and Nutrition Study.
- DNA methylation: 2016 Venous Blood Study.
- Primary pace-of-ageing outcome: DunedinPACE.
- Core covariates: age, sex, race/ethnicity, education, smoking, alcohol consumption, physical activity, marital status, socioeconomic position and total energy intake.
- Survey-weighted linear models.
- Component-specific Phe and Tyr estimates were pooled across 50 completed imputations using Rubin's rules.

## NHANES

- Dietary exposure: 2001–2002 individual 24-hour-recall food records.
- Food composition linkage: FNDDS 1.0 and USDA Standard Reference 16-1.
- Main pace-of-ageing bridge: DunedinPoAm.
- Survey design: DNAm subsample weight, strata and primary sampling units.
- The direct reconstructed Phe+Tyr analysis included n=1,085.

## Albumin-profile similarity

18 dietary amino acids:
Trp, Thr, Ile, Leu, Lys, Met, Cys/cystine, Phe, Tyr, Val, Arg, His, Ala, Asp, Glu, Gly, Pro, Ser.

Mature human serum albumin 18-bin residue counts (total=585):

- Trp 1
- Thr 28
- Ile 8
- Leu 61
- Lys 59
- Met 6
- Cys 35
- Phe 31
- Tyr 18
- Val 41
- Arg 24
- His 16
- Ala 62
- Asp+Asn 53
- Glu+Gln 82
- Gly 12
- Pro 24
- Ser 24

Procedure:
1. Convert g/day to mol/day using amino-acid molecular weights.
2. Convert dietary cystine to two cysteine-residue equivalents.
3. Normalize each participant's 18-component molar profile to sum to 1.
4. Normalize albumin residue counts to sum to 1.
5. Calculate cosine similarity.
6. Standardize the similarity score to 1 SD in the analytic sample.

## Interpretation constraints

- Developmental intervention in Drosophila vs habitual adult diet in humans.
- Lifespan/translation pathway in flies vs DNAm biomarkers in humans.
- DunedinPACE and DunedinPoAm are related but distinct biomarkers.
- Phe/Tyr intake is highly correlated with total dietary protein, limiting amino-acid-specific inference in free-living human diets.
