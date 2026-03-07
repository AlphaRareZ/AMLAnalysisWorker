import pandas as pd
import numpy as np
import logging
from .utils import ensure_dir_for_file
from mygene import MyGeneInfo

logger = logging.getLogger(__name__)

def map_exons_to_genes(expression_file, mapping_file, output_file):
    logger.info("--------------------Mapping Exons To Genes--------------------")
    
    # Load expression data (using float32 to cut RAM usage in half)
    expr = pd.read_csv(expression_file, sep="\t", index_col=0)

    # Load and format the mapping file
    mapping = (
        pd.read_csv(mapping_file, sep="\t")[["id", "gene"]]
        .dropna()
        .drop_duplicates()
        .set_index("id") # Keep "id" as the index to match expr
    )

    # The Fix: Use an inner join instead of adding a column.
    # This automatically filters to the intersection AND adds the 'gene' column 
    # cleanly without fragmenting the memory.
    expr = expr.join(mapping["gene"], how="inner")

    # Group by the gene and calculate the mean
    gene_expr = expr.groupby("gene").mean()

    ensure_dir_for_file(output_file) 
    gene_expr.to_csv(output_file)

    logger.info("Gene-level expression saved to:", output_file)
    logger.info("Final shape:", gene_expr.shape)

def filter_protein_coding(input_file, output_file, batch_size=1000):
    logger.info("--------------------Filtering Protein Coding--------------------")
    mg = MyGeneInfo()
    gene_expr = pd.read_csv(input_file, index_col=0)

    all_genes = gene_expr.index.tolist()
    query_results = []
    for i in range(0, len(all_genes), batch_size):
        batch = all_genes[i : i + batch_size]
        try:
            query_results.extend(
                mg.querymany(
                    batch, scopes="symbol", fields="type_of_gene", species="human"
                )
            )
        except Exception as e:
            logger.info(f"MyGene query batch failed: {e}")

    gene_types = {
        q["query"]: q.get("type_of_gene", None) for q in query_results if "query" in q
    }
    keep_genes = [g for g, t in gene_types.items() if t == "protein-coding"]

    filtered = gene_expr.loc[gene_expr.index.intersection(keep_genes)]

    ensure_dir_for_file(output_file)  # Create dir if needed
    filtered.to_csv(output_file)

    logger.info(f"Filtered {gene_expr.shape[0]} -> {filtered.shape[0]} (protein-coding only)")
    logger.info(f"Saved to: {output_file}")


def normalize_log(expr):
    return np.log2(expr + 1.0)


def select_hvgs(expr_file, out_hvgs, out_stats, n_hvgs=2000, min_mean=0.5):
    logger.info("--------------------Selecting HVGS--------------------")
    expr = pd.read_csv(expr_file, index_col=0)
    expr_norm = normalize_log(expr)

    mean_expr = expr_norm.mean(axis=1)
    expr_f = expr_norm.loc[mean_expr >= min_mean]
    logger.info("After mean filter:", expr_f.shape)

    var = expr_f.var(axis=1)
    hvgs = var.sort_values(ascending=False).head(n_hvgs)
    hvgs_df = expr_f.loc[hvgs.index]

    # Save HVGs
    ensure_dir_for_file(out_hvgs)  # Create dir if needed
    hvgs_df.to_csv(out_hvgs)

    # Save variance table
    stats_df = pd.DataFrame(
        {"mean": mean_expr.loc[hvgs.index], "variance": var.loc[hvgs.index]}
    )
    ensure_dir_for_file(out_stats)  # Create dir if needed
    stats_df.to_csv(out_stats)

    logger.info(f"Saved top {n_hvgs} HVGs to {out_hvgs}")
    logger.info(f"Saved HVG stats to {out_stats}")
