"""In-silico perturbation functions for knock-up / knock-down simulation."""
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def perturb_expression(
    adata: ad.AnnData,
    genes: list[str],
    direction: str,
    cell_mask: pd.Series | None = None,
) -> ad.AnnData:
    """
    Return a copy of adata with expression of `genes` perturbed.

    Parameters
    ----------
    adata      : AnnData object (.X used as expression matrix)
    genes      : gene names to perturb
    direction  : "down" → set to 0; "up" → set to per-cell max + 1
    cell_mask  : boolean Series aligned to adata.obs; restricts perturbation
                 to selected cells (e.g. only sALS cells)

    Notes
    -----
    GeneFormer tokenizes cells by rank-ordering genes internally, so modifying
    .X before embedding is sufficient, the rank shift follows automatically.
    - Knock-down: setting expression to 0 moves the gene to the bottom of the
      rank ordering, effectively silencing it.
    - Knock-up: setting expression to per-cell max + 1 moves the gene to the
      top of the rank ordering, simulating maximal expression.
    """
    assert direction in ("down", "up"), "direction must be 'down' or 'up'"

    missing = [g for g in genes if g not in adata.var_names]
    if missing:
        raise ValueError(f"Genes not found in dataset: {missing}")

    tag = f"[{genes[0]} {direction}]"
    gene_idx = [adata.var_names.get_loc(g) for g in genes]
    cell_idx = np.where(cell_mask.values)[0] if cell_mask is not None \
               else np.arange(adata.n_obs)

    print(f"{tag} tocsc copy...")
    X = adata.X.tocsc().copy() if sp.issparse(adata.X) else adata.X.copy()

    if direction == "up":
        print(f"{tag} computing row max...")
        row_max = X[cell_idx].toarray().max(axis=1)

    print(f"{tag} modifying columns...")
    for gi in gene_idx:
        col = X[:, gi].toarray().ravel()
        if direction == "down":
            col[cell_idx] = 0
        else:
            col[cell_idx] = row_max + 1
        X[:, gi] = col.reshape(-1, 1)

    print(f"{tag} building AnnData...")
    pert = ad.AnnData(
        X=X.tocsr() if sp.issparse(X) else X,
        obs=adata.obs.copy(),
        var=adata.var.copy(),
    )
    pert.uns["perturbation"] = {
        "genes": genes,
        "direction": direction,
        "n_cells_perturbed": len(cell_idx),
    }
    return pert
