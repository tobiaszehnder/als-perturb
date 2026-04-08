import anndata as ad
import numpy as np
import pandas as pd


def get_sample_size(
    adata: ad.AnnData,
    groups: list[str],
    celltypes: list[str],
    max_n: int = 150,
    celltype_col: str = "Cellstates_LVL3",
    group_col: str = "Group",
) -> int:
    """Return the largest n that fits in every (group, celltype) stratum, capped at max_n."""
    min_cells = min(
        ((adata.obs[group_col] == group) & (adata.obs[celltype_col] == ct)).sum()
        for group in groups
        for ct in celltypes
    )
    return min(min_cells, max_n)


def sample_cells(
    adata: ad.AnnData,
    groups: list[str],
    celltypes: list[str],
    n: int,
    celltype_col: str = "Cellstates_LVL3",
    group_col: str = "Group",
    seed: int = 42,
) -> ad.AnnData:
    """Sample n cells per (group, celltype) combination and return as a single AnnData.

    Prints a summary of cells sampled vs available for each stratum.
    """
    rng = np.random.default_rng(seed)
    idx = []
    for group in groups:
        for ct in celltypes:
            mask = (adata.obs[group_col] == group) & (adata.obs[celltype_col] == ct)
            available = np.where(mask)[0]
            chosen = rng.choice(available, size=min(n, len(available)), replace=False)
            idx.append(chosen)
            print(f"  {group} / {ct}: {len(chosen):,} sampled (of {len(available):,} available)")
    return adata[np.concatenate(idx)].copy()


def get_als_cell_mask(adata: ad.AnnData, group: str, celltypes: list[str]) -> pd.Series:
    """Boolean mask selecting sALS cells from the relevant excitatory subtypes."""
    return (
        (adata.obs["Group"] == group) &
        (adata.obs["Cellstates_LVL3"].isin(celltypes))
    )
