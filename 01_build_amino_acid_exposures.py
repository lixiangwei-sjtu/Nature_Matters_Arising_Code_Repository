"""Build manuscript amino-acid exposures from participant-level dietary data.

This script intentionally contains no participant-level data and no restricted
HRS variable names. Supply a column map/configuration for your authorized data.
"""

from __future__ import annotations
import numpy as np
import pandas as pd


def zscore(x: pd.Series) -> pd.Series:
    return (x - x.mean()) / x.std(ddof=1)


def combine_phe_tyr(df: pd.DataFrame, phe_col: str, tyr_col: str,
                    out_col: str = "phe_tyr_g_day") -> pd.DataFrame:
    out = df.copy()
    out[out_col] = pd.to_numeric(out[phe_col], errors="coerce") + pd.to_numeric(out[tyr_col], errors="coerce")
    return out


def residual_energy_adjust(exposure: pd.Series, energy: pd.Series) -> pd.Series:
    """Residual-method energy adjustment with the exposure mean added back.

    This generic implementation is provided for reproducibility. Confirm that
    it matches the exact locked analytic specification before generating final
    manuscript estimates.
    """
    x = pd.to_numeric(exposure, errors="coerce")
    e = pd.to_numeric(energy, errors="coerce")
    ok = x.notna() & e.notna()
    result = pd.Series(np.nan, index=x.index, dtype=float)
    if ok.sum() < 3:
        return result
    X = np.column_stack([np.ones(ok.sum()), e.loc[ok].to_numpy()])
    beta = np.linalg.lstsq(X, x.loc[ok].to_numpy(), rcond=None)[0]
    fitted = X @ beta
    result.loc[ok] = x.loc[ok].to_numpy() - fitted + x.loc[ok].mean()
    return result


def albumin_enriched_score(df: pd.DataFrame, column_map: dict[str, str],
                           weights: dict[str, float]) -> pd.Series:
    """Weighted score for CYS, LYS, PHE, TYR, ALA, GLU and HIS.

    `weights` must contain the locked relative-enrichment weights used in the
    final manuscript analysis. The code deliberately does not invent them.
    """
    aas = ["CYS", "LYS", "PHE", "TYR", "ALA", "GLU", "HIS"]
    missing = [aa for aa in aas if aa not in column_map or aa not in weights]
    if missing:
        raise ValueError(f"Missing column mapping or locked weight for: {missing}")
    pieces = []
    for aa in aas:
        x = pd.to_numeric(df[column_map[aa]], errors="coerce")
        pieces.append(x * float(weights[aa]))
    return pd.concat(pieces, axis=1).sum(axis=1, min_count=len(aas)).rename("albumin_enriched_score")


if __name__ == "__main__":
    print("Import this module from an analysis script and supply authorized data + locked mappings.")
