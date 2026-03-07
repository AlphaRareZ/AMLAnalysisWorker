import pandas as pd
import numpy as np
import logging
import gc
from .utils import ensure_dir_for_file
from mygene import MyGeneInfo

logger = logging.getLogger(__name__)

def map_exons_to_genes(expression_file, mapping_file, output_file):
    logger.info("--------------------Mapping Exons To Genes--------------------")
    
    expr = pd.read_csv(expression_file, sep="\t", index_col=0)

    mapping = (
        pd.read_csv(mapping_file, sep="\t")[["id", "gene"]]
        .dropna()
        .drop_duplicates()
        .set_index("id") 
    )

    expr = expr.join(mapping["gene"], how="inner")
    
    # Clean up mapping before grouping
    del mapping
    gc.collect()

    gene_expr = expr.groupby("gene").mean()

    ensure_dir_for_file(output_file) 
    gene_expr.to_csv(output_file)

    logger.info(f"Gene-level expression saved to: {output_file}")
    logger.info(f"Final shape: {gene_expr.shape}")
    
    # Final cleanup
    del expr, gene_expr
    gc.collect()

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

    ensure_dir_for_file(output_file) 
    filtered.to_csv(output_file)

    logger.info(f"Filtered {gene_expr.shape[0]} -> {filtered.shape[0]} (protein-coding only)")
    logger.info(f"Saved to: {output_file}")
    
    # Final cleanup
    del gene_expr, filtered, all_genes, query_results, gene_types, keep_genes
    gc.collect()


def normalize_log(expr):
    return np.log2(expr + 1.0)


def select_hvgs(expr_file, out_hvgs, out_stats, n_hvgs=2000, min_mean=0.5):
    logger.info("--------------------Selecting HVGS--------------------")
    expr = pd.read_csv(expr_file, index_col=0)
    expr_norm = normalize_log(expr)

    mean_expr = expr_norm.mean(axis=1)
    expr_f = expr_norm.loc[mean_expr >= min_mean]
    logger.info(f"After mean filter: {expr_f.shape}")

    var = expr_f.var(axis=1)
    hvgs = var.sort_values(ascending=False).head(n_hvgs)
    hvgs_df = expr_f.loc[hvgs.index]

    ensure_dir_for_file(out_hvgs) 
    hvgs_df.to_csv(out_hvgs)

    stats_df = pd.DataFrame(
        {"mean": mean_expr.loc[hvgs.index], "variance": var.loc[hvgs.index]}
    )
    ensure_dir_for_file(out_stats) 
    stats_df.to_csv(out_stats)

    logger.info(f"Saved top {n_hvgs} HVGs to {out_hvgs}")
    logger.info(f"Saved HVG stats to {out_stats}")
    
    # Final cleanup
    del expr, expr_norm, mean_expr, expr_f, var, hvgs, hvgs_df, stats_df
    gc.collect()