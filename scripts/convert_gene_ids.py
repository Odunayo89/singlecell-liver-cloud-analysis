import scanpy as sc
import pandas as pd

# Load processed object
adata = sc.read_h5ad("results/liver_processed_subset.h5ad")

# Load marker genes
markers = pd.read_csv("results/cluster_marker_genes.csv")

# Create ENSG -> gene symbol mapping
gene_map = dict(zip(adata.var_names, adata.var["feature_name"]))

# Convert IDs to symbols
markers["gene_symbol"] = markers["names"].map(gene_map)

# Save updated file
markers.to_csv("results/cluster_marker_genes_annotated.csv", index=False)

print(markers.head())
print("Annotated marker genes saved successfully")
