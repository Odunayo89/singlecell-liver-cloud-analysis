import scanpy as sc
import numpy as np

# Load dataset
adata = sc.read_h5ad("data/feca90bb-00df-4623-8398-1e3e6a90971d.h5ad")

print("Dataset loaded successfully")
print(adata)

# Use existing metadata-aware columns if present
print("Available obs columns:")
print(list(adata.obs.columns))

# Subsample for memory-safe cloud analysis
sc.pp.subsample(adata, n_obs=20000, random_state=42)

print("After subsampling:")
print(adata)

# QC
adata.var["mt"] = adata.var_names.str.startswith("MT-")

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True
)

# Filter cells
adata = adata[
    (adata.obs["n_genes_by_counts"] > 200) &
    (adata.obs["pct_counts_mt"] < 20),
    :
].copy()

print("After QC filtering:")
print(adata)

# Normalize and log transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Highly variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var["highly_variable"]].copy()

# PCA / neighbors / UMAP / clustering
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)

# Save UMAP clustering
sc.pl.umap(
    adata,
    color=["leiden", "cell_type", "disease"],
    save="_liver_cloud_analysis.png",
    show=False
)

# UMAP colored by cell type
sc.pl.umap(
    adata,
    color="cell_type",
    save="_celltypes.png",
    show=False
)

# UMAP colored by disease
sc.pl.umap(
    adata,
    color="disease",
    save="_disease.png",
    show=False
)

# QC violin plots
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    save="_qc_metrics.png",
    show=False
)

# Marker genes by cluster
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")

markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.to_csv("results/cluster_marker_genes.csv", index=False)

# Save processed object
adata.write("results/liver_processed_subset.h5ad")

print("Single-cell liver analysis completed successfully")
