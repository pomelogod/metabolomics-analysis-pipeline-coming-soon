# %%
"""
Genes_ranked.py
language: Python3
authors: David Montefusco <David.Montefusco@vcuhealth.org>, L. Xie <xiel4@vcu.edu>
"""
import os
import pandas as pd
import numpy as np
import scipy.stats as stats
import statsmodels.stats.multitest as multitest
import json

# Load config
with open('config.json', 'r') as f:
    config = json.load(f)
merged_proteins_path = config['global']['merged_proteins_path']

# Get input path from user
enrichment_results_path = input("Please enter the path to your enrichment analysis results CSV file: ")
# Strip quotes if present
enrichment_results_path = enrichment_results_path.strip('"\'')

# Normalize path separators for cross-platform compatibility
enrichment_results_path = os.path.normpath(enrichment_results_path)

# Generate output path in the same directory as input
output_dir = os.path.dirname(enrichment_results_path)
output_filename = "genes_ranked.csv"
output_path = os.path.join(output_dir, output_filename)

# Load input files with error handling
try:
    enrichment_results_df = pd.read_csv(enrichment_results_path)
    merged_proteins_df = pd.read_csv(merged_proteins_path, low_memory=False, dtype={"Gene Name": str})
    print("Successfully loaded all input files.")
except FileNotFoundError as e:
    print(f"ERROR: Missing file - {e}")
    exit()

# Standardize column names
merged_proteins_df.columns = merged_proteins_df.columns.str.strip()
enrichment_results_df.columns = enrichment_results_df.columns.str.strip()

#Rename columns for consistency
merged_proteins_df.rename(columns={"Pathway Name": "Pathway", "Gene Name": "Gene"}, inplace=True)

# Standardize pathway names for comparison
merged_proteins_df["Pathway"] = merged_proteins_df["Pathway"].str.strip().str.lower()
enrichment_results_df["Pathway"] = enrichment_results_df["Pathway"].str.strip().str.lower()

#Step 1: Filter Out 'Drug Action' Pathways
merged_proteins_df = merged_proteins_df[merged_proteins_df["Pathway Subject"] != "Drug Action"]

# Step 2: Create Pathway-to-Gene Mapping
gene_to_pathways = {}
for _, row in merged_proteins_df.iterrows():
    gene = row["Gene"]
    pathway = row["Pathway"]

    if gene not in gene_to_pathways:
        gene_to_pathways[gene] = set()
    gene_to_pathways[gene].add(pathway)

#Step 3: Create Pathway-to-Compound and Compound Name Mapping
pathway_to_compounds = {}
pathway_to_names = {}

for _, row in enrichment_results_df.iterrows():
    pathway = row["Pathway"]
    compounds = str(row["Hits Detail (HMDB)"]).split(", ")  # Split compounds if multiple exist
    compound_names = str(row["Hits Detail (Name)"]).split(", ")  # Split names if multiple exist

    if pathway not in pathway_to_compounds:
        pathway_to_compounds[pathway] = set()
        pathway_to_names[pathway] = set()

    pathway_to_compounds[pathway].update(compounds)
    pathway_to_names[pathway].update(compound_names)

#Convert input pathways to set for matching
input_pathways = set(enrichment_results_df["Pathway"])  # Unique input pathways
N = len(set(merged_proteins_df["Pathway"]))  # Total unique pathways
n = len(input_pathways)  # Total input pathways being tested

#Step 4: Compute Enrichment per Gene (Exclude zero-overlap genes)
gene_enrichment_results = []

for gene, pathways in gene_to_pathways.items():
    K = len(pathways)  #Total pathways for this gene
    overlapping_pathways = pathways & input_pathways  #Overlapping pathways
    k = len(overlapping_pathways)  #Overlapping pathway count

    if k == 0:
        continue  #Skip genes with zero overlap

    #Collect compounds from the matched pathways
    related_compounds = set()
    related_compound_names = set()

    for pathway in overlapping_pathways:
        related_compounds.update(pathway_to_compounds.get(pathway, []))
        related_compound_names.update(pathway_to_names.get(pathway, []))

    compound_list = ", ".join(related_compounds) if related_compounds else "None"
    compound_name_list = ", ".join(related_compound_names) if related_compound_names else "None"

    #Compute hypergeometric p-value
    hypergeo_p_value = 1 - stats.hypergeom.cdf(k - 1, N, K, n) if k > 0 else 1.0

    #Compute Fisher's Exact Test for weighting
    background_gene_count = len(gene_to_pathways)  # Total genes
    genes_with_overlap = sum(1 for g, p in gene_to_pathways.items() if len(p & input_pathways) > 0)  # Genes with overlap

    contingency_table = [
        [k, K - k],  # Genes overlapping vs. not overlapping
        [genes_with_overlap - k, background_gene_count - genes_with_overlap]  # Non-overlapping genes
    ]
    
    fisher_p_value = stats.fisher_exact(contingency_table, alternative="greater")[1]

    #Combine p-values: Weighted Fisher's + Hypergeometric (simple average)
    combined_p_value = (hypergeo_p_value + fisher_p_value) / 2

    gene_enrichment_results.append([
        gene, K, k, hypergeo_p_value, fisher_p_value, combined_p_value, compound_list, compound_name_list
    ])

#Convert to DataFrame
enrichment_df = pd.DataFrame(gene_enrichment_results, columns=[
    "Gene", "Total Pathways", "Overlap", "Hypergeometric P-Value", "Fisher P-Value",
    "Combined P-Value", "Related Compounds", "Compound Names"
])

#Step 5: Apply FDR Correction
if not enrichment_df.empty:
    _, fdr_p_values, _, _ = multitest.multipletests(enrichment_df["Combined P-Value"], method="fdr_bh")
    enrichment_df["FDR-Corrected P-Value"] = fdr_p_values

#Save Output
enrichment_df.to_csv(output_path, index=False)
print(f"\nGene ranking saved to: {output_path}")
