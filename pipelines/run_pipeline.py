import json
import gc
import pandas as pd
from logging import getLogger

from .utils import ensure_dir_for_file
from .data_preprocessing import map_exons_to_genes, filter_protein_coding, select_hvgs
from .network_analysis import build_adjacency, module_detection, module_eigengenes, intramodular_connectivity, rank_and_annotate
from  .building_structures import fetch_alphafold, render_proteins, combine_images
from .eda_analysis import run_simple_eda, run_advanced_eda
logger = getLogger(__name__)

def main(expression_file, mapping_file, config_path="config.json"):

    # --- 1. Load Configuration ---
    global LOG_FILE_PATH
    try:
        with open(config_path, "r") as f:
            config = json.load(f)
    except FileNotFoundError:
        logger.info(f"FATAL: Configuration file not found at {config_path}")
        return
    except json.JSONDecodeError:
        logger.info(f"FATAL: Configuration file {config_path} is not valid JSON.")
        return

    LOG_FILE_PATH = config["files"]["log_file"]
    ensure_dir_for_file(LOG_FILE_PATH)

    logger.info("Starting a New Process")
    logger.info(f"Loaded configuration from {config_path}")

    files_cfg = config["files"]
    dirs_cfg = config["directories"]

    gene_expr_raw = files_cfg["gene_expression_raw"]
    gene_expr_coding = files_cfg["gene_expression_coding"]
    hvgs_file = files_cfg["hvgs"]
    hvgs_stats_file = files_cfg["hvgs_stats"]

    # --- 2. Pre-processing ---
    map_exons_to_genes(expression_file, mapping_file, gene_expr_raw)
    gc.collect() # Safety net
    
    filter_protein_coding(gene_expr_raw, gene_expr_coding)
    gc.collect()
    
    select_hvgs(gene_expr_coding, hvgs_file, hvgs_stats_file, n_hvgs=2000, min_mean=0.5)
    gc.collect()

    # --- 3. Coexpression and Modules ---
    logger.info("--------------------Coexpression and Modules--------------------")
    try:
        expr = pd.read_csv(hvgs_file, index_col=0)
        logger.info(f"HVG matrix shape: {expr.shape}")

        adj = build_adjacency(expr, power=6)
        modules, G = module_detection(adj)
        
        # Fixing the f-string issue you mentioned previously
        modules_count = {k: len(v) for k, v in modules.items()}
        logger.info(f"Detected modules: {modules_count}")

        eig_df = module_eigengenes(expr, modules)
        me_file = files_cfg["module_eigengenes"]
        ensure_dir_for_file(me_file)
        eig_df.to_csv(me_file)

        k_within = intramodular_connectivity(adj, modules)
        k_file = files_cfg["intramodular_connectivity"]
        ensure_dir_for_file(k_file)
        k_within.to_csv(k_file, header=["kWithin"])

        mod_json_file = files_cfg["modules_json"]
        ensure_dir_for_file(mod_json_file)
        with open(mod_json_file, "w") as f:
            json.dump(modules, f, indent=2)

        logger.info("Saved modules, eigengenes, and intramodular connectivity.")
        
        # Cleanup
        del expr, adj, modules, G, eig_df, k_within
        gc.collect()
        
    except FileNotFoundError:
        logger.info(f"Failed to read HVG file {hvgs_file}. Skipping module detection.")
        return  
    except Exception as e:
        logger.info(f"Error during module detection: {e}")
        return
        
    # --- 4. EDA For the filtered Genes ---
    run_simple_eda(gene_expr_coding, files_cfg)
    gc.collect()
    
    # --- 5. Rank Hubs and Annotate (WGCNA-like) ---
    ranked_cfg = files_cfg["ranked_biomarkers"]
    final = rank_and_annotate(
        hvgs_file,
        files_cfg["intramodular_connectivity"],
        files_cfg["modules_json"],
        gene_expr_coding,
        top_n=20,
        out_file=ranked_cfg["top20_annotated"],
        all_ranked_file=ranked_cfg["all_ranked"],
    )
    if final is not None and not final.empty:
        logger.info("Top ranked biomarkers (WGCNA-like score):")
        # Fixing the f-string parsing here as well
        cols_to_log = ["score", "kWithin", "variance", "mean", "type_of_gene", "Entry", "targetability"]
        logger.info(f"\n{final[cols_to_log]}")
    else:
        logger.info("Ranking and annotation step failed or returned empty.")

    del final
    gc.collect()

    # --- 6. Advanced Biomarkers EDA ---
    run_advanced_eda(
        gene_expr_raw, 
        ranked_cfg, 
        files_cfg
    )
    gc.collect()
    
    # --- 7. Fetching & Visualizing 3D Structures ---
    protein_cfg = files_cfg["protein"]
    csv_input = ranked_cfg["top20_annotated"]

    report_df = fetch_alphafold(
        csv_input, dirs_cfg["protein_structures"], protein_cfg["fetch_report"]
    )
    
    del report_df
    gc.collect()

    render_proteins(protein_cfg["fetch_report"])

    combine_images(
        dirs_cfg["protein_images"],  
        protein_cfg["combined_image"]
    )
    
    logger.info("Pipeline complete.")
    
    # Final worker wipe
    gc.collect()