"""Constants for the ALS perturbation analysis."""
# celltypes with the largest transcriptomic divergence from PN according to Pineda et al. Fig 4C
RELEVANT_CELLTYPES = [
    "Ex.L5.VAT1L_THSD4",
    "Ex.L5.VAT1L_EYA4",
    "Ex.L3_L5.SCN4B_NEFH",
]

# Gene selection based exclusively on SALS vs PN (MCX) findings in Pineda et al., 2024
# Targeted for Ex.L5.VAT1L and Ex.L3_L5.SCN4B populations
SALS_GENES = [
    # Predicted Master Regulators of L5 VAT1L+ DEGs
    "MYT1L",    # Top predicted master regulator of neuronal identity
    "REST",     # Top predicted master regulator 
    "SREBF2",   # Master regulator of cholesterol biosynthesis

    # Top Transcriptome-wide Transcriptional Divergence (TxD) Drivers
    "LYNX1",    # Top TxD-associated gene; modulates nAChRs to prevent hyperexcitation
    "TBR1",     # Broadly dysregulated TxD gene; normally restricts corticospinal tract origins

    # Axonal/Structural Genes Highly Correlated with TxD
    "KIF5A",    # ALS-linked gene strongly correlated with TxD
    "DCTN1",    # Axonal transport gene strongly correlated with TxD
    "DYNC1H1",  # Axonal transport gene strongly correlated with TxD
    "TUBA4A",   # Explicitly noted as an ALS-linked LMN gene correlated with TxD

    # Highly Specific Differentially Expressed Genes (DEGs) in MCX
    "NEFH",     # Selectively upregulated specifically in MCX L5 VAT1L+ neurons in ALS
    "TTC14",    # Among the most downregulated genes in the MCX
    "HSP90AA1", # Highly upregulated DEG; facilitates toxic aggregate accumulation
]