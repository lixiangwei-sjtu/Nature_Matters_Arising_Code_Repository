"""
Reproducibility definition for the albumin-profile similarity used in the
Nature Matters Arising HRS/NHANES analysis.

This file defines the exposure only; participant-level HRS/NHANES data are
not redistributed here.

Locked definition
-----------------
1. Use 18 dietary amino-acid measurements:
   Trp, Thr, Ile, Leu, Lys, Met, Cys/cystine, Phe, Tyr, Val,
   Arg, His, Ala, Asp, Glu, Gly, Pro, Ser.
2. Convert g/day to mol/day.
3. For dietary cystine, multiply molar cystine by 2 to obtain cysteine-residue
   equivalents.
4. Normalize the 18 molar residue-equivalent amounts so they sum to 1.
5. Build the mature human serum albumin 18-bin residue profile:
   W=1, T=28, I=8, L=61, K=59, M=6, C=35, F=31, Y=18,
   V=41, R=24, H=16, A=62, D+N=53, E+Q=82, G=12, P=24, S=24.
   The counts sum to 585 residues; normalize them to sum to 1.
6. Calculate cosine similarity between the participant profile and albumin
   profile.
7. Standardize the resulting similarity score before regression.

Reverse-engineering QC against the locked HRS results:
- correlation with energy-adjusted total protein: r = 0.317879 (~0.32)
- GrimAge acceleration beta: -0.065578 (reported -0.0656)
- GrimAge acceleration beta after total-protein adjustment: -0.054177
  (reported -0.0542)

These fingerprints distinguish this definition from mass-fraction,
Pearson-correlation, Euclidean-distance and other candidate similarity metrics.
"""

import numpy as np
import pandas as pd

AA_ORDER = [
    "TRP", "THR", "ILE", "LEU", "LYS", "MET", "CYS", "PHE", "TYR",
    "VAL", "ARG", "HIS", "ALA", "ASP", "GLU", "GLY", "PRO", "SER"
]

# Free amino-acid / cystine molecular weights used to convert dietary g/day
# to molar amounts. Cystine is subsequently multiplied by two residues.
MOLECULAR_WEIGHT = {
    "TRP": 204.23, "THR": 119.12, "ILE": 131.17, "LEU": 131.17,
    "LYS": 146.19, "MET": 149.21, "CYS": 240.30, "PHE": 165.19,
    "TYR": 181.19, "VAL": 117.15, "ARG": 174.20, "HIS": 155.16,
    "ALA": 89.09, "ASP": 133.10, "GLU": 147.13, "GLY": 75.07,
    "PRO": 115.13, "SER": 105.09
}

# Mature human serum albumin, 585 residues.
# ASP bin = Asp + Asn; GLU bin = Glu + Gln, matching hydrolytic dietary
# amino-acid measurements.
ALBUMIN_COUNTS = {
    "TRP": 1, "THR": 28, "ILE": 8, "LEU": 61, "LYS": 59, "MET": 6,
    "CYS": 35, "PHE": 31, "TYR": 18, "VAL": 41, "ARG": 24, "HIS": 16,
    "ALA": 62, "ASP": 53, "GLU": 82, "GLY": 12, "PRO": 24, "SER": 24
}

assert sum(ALBUMIN_COUNTS.values()) == 585


def albumin_profile_similarity(df: pd.DataFrame, column_map: dict) -> pd.Series:
    """
    Parameters
    ----------
    df : DataFrame
        Participant-level dietary data.
    column_map : dict
        Maps each AA_ORDER key to its g/day column in df.

    Returns
    -------
    Series
        Cosine similarity (unstandardized). Standardize this Series before
        regression in the analytic sample.
    """
    grams = np.column_stack([
        pd.to_numeric(df[column_map[aa]], errors="coerce").to_numpy()
        for aa in AA_ORDER
    ])

    mw = np.array([MOLECULAR_WEIGHT[aa] for aa in AA_ORDER], dtype=float)
    mol = grams / mw

    # Dietary variable is cystine; one cystine molecule yields two cysteine
    # residue equivalents.
    cys_idx = AA_ORDER.index("CYS")
    mol[:, cys_idx] *= 2.0

    row_sum = np.nansum(mol, axis=1)
    profile = mol / row_sum[:, None]

    albumin = np.array([ALBUMIN_COUNTS[aa] for aa in AA_ORDER], dtype=float)
    albumin /= albumin.sum()

    numerator = np.sum(profile * albumin[None, :], axis=1)
    denominator = np.linalg.norm(profile, axis=1) * np.linalg.norm(albumin)
    similarity = numerator / denominator

    similarity[~np.isfinite(similarity)] = np.nan
    return pd.Series(similarity, index=df.index, name="albumin_profile_similarity")


def zscore(x: pd.Series) -> pd.Series:
    """Sample-SD standardization used before regression."""
    return (x - x.mean()) / x.std(ddof=1)
