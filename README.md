# Cloud-Based Single-Cell RNA-seq Analysis of Human Liver Tissue Using AWS EC2 and Scanpy

## Overview
This project demonstrates a scalable cloud-based single-cell RNA-seq analysis workflow using AWS EC2, AWS S3, and Scanpy.

## Objective
To perform preprocessing, clustering, visualization, and marker gene analysis on human liver single-cell transcriptomic data in a cloud computing environment.

## Dataset
Healthy and IFALD pediatric and adult human liver tissue dataset (81,001 cells)

Source:
CELLxGENE Discover

## Workflow
1. Upload dataset to AWS S3
2. Pull dataset into AWS EC2
3. Create Python virtual environment
4. Install Scanpy ecosystem
5. Perform quality control filtering
6. Normalize and log-transform data
7. Identify highly variable genes
8. Perform PCA and UMAP dimensionality reduction
9. Perform Leiden clustering
10. Generate visualization figures
11. Export marker genes and processed AnnData object
12. Store outputs in AWS S3

## Tools & Technologies
- AWS EC2
- AWS S3
- Linux (bash)
- Python
- Scanpy
- AnnData
- pandas
- matplotlib

## Key Features
- Cloud-native scRNA-seq workflow
- Large dataset handling
- Memory-aware optimization using subsampling
- Reproducible project structure
- Automated visualization generation

## Outputs
- UMAP clustering visualization
- Cell type visualization
- Disease-state visualization
- QC violin plots
- Marker gene tables
- Processed AnnData object

## Skills Demonstrated
- Cloud computing for bioinformatics
- Single-cell RNA-seq analysis
- Linux infrastructure management
- Scanpy workflows
- scRNA-seq preprocessing
- Clustering and visualization
- Memory optimization for large datasets
- Reproducible computational biology workflows

## Author
Adekunle Ajiboye
