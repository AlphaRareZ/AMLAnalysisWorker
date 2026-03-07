import json
import pandas as pd
from logging import getLogger

# بنستورد الـ Functions من الفايلات بتاعتنا
from utils import ensure_dir_for_file
from data_preprocessing import map_exons_to_genes, filter_protein_coding, select_hvgs
from network_analysis import build_adjacency, module_detection, module_eigengenes, intramodular_connectivity, rank_and_annotate
from  building_structures import fetch_alphafold, render_proteins, combine_images
from eda_analysis import run_simple_eda, run_advanced_eda
logger = getLogger(__name__)

# Main Pipeline
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

    # Set log file path from config and ensure its directory exists
    LOG_FILE_PATH = config["files"]["log_file"]
    ensure_dir_for_file(LOG_FILE_PATH)

    logger.info("Starting a New Process")
    logger.info(f"Loaded configuration from {config_path}")

    # Define file path variables from config
    files_cfg = config["files"]
    dirs_cfg = config["directories"]

    gene_expr_raw = files_cfg["gene_expression_raw"]
    gene_expr_coding = files_cfg["gene_expression_coding"]
    hvgs_file = files_cfg["hvgs"]
    hvgs_stats_file = files_cfg["hvgs_stats"]

    # --- 2. Pre-processing ---
    map_exons_to_genes(expression_file, mapping_file, gene_expr_raw)
    filter_protein_coding(gene_expr_raw, gene_expr_coding)
    select_hvgs(gene_expr_coding, hvgs_file, hvgs_stats_file, n_hvgs=2000, min_mean=0.5)

    # --- 3. Coexpression and Modules ---
    logger.info("--------------------Coexpression and Modules--------------------")
    try:
        expr = pd.read_csv(hvgs_file, index_col=0)
        logger.info("HVG matrix shape:", expr.shape)

        adj = build_adjacency(expr, power=6)
        modules, G = module_detection(adj)
        logger.info("Detected modules:", {k: len(v) for k, v in modules.items()})

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
    except FileNotFoundError:
        logger.info(f"Failed to read HVG file {hvgs_file}. Skipping module detection.")
        return  # Stop pipeline if this critical step fails
    except Exception as e:
        logger.info(f"Error during module detection: {e}")
        return
    # --- 4. EDA For the filtered Genes ---
    run_simple_eda(gene_expr_coding, files_cfg)
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
        logger.info(
            final[
                [
                    "score",
                    "kWithin",
                    "variance",
                    "mean",
                    "type_of_gene",
                    "Entry",
                    "targetability",
                ]
            ]
        )
    else:
        logger.info("Ranking and annotation step failed or returned empty.")
        # We might want to stop here if downstream steps depend on this
        # return
    # --- 6. Advanced Biomarkers EDA ---
    run_advanced_eda(
        gene_expr_raw, 
        ranked_cfg, 
        files_cfg
    )
    # --- 7. Fetching & Visualizing 3D Structures ---
    protein_cfg = files_cfg["protein"]
    csv_input = ranked_cfg["top20_annotated"]

    # Step 1: Fetch
    report_df = fetch_alphafold(
        csv_input, dirs_cfg["protein_structures"], protein_cfg["fetch_report"]
    )

    # Step 2: Render (Skipped)
    render_proteins(protein_cfg["fetch_report"])

    # Step 3: Combine (This function is for a different set of images)
    # The script seems to have *two* image combination/collage steps.
    # This one combines images from a different dir name.
    # We'll call it as intended.
    combine_images(
        dirs_cfg["protein_images"],  # "results/proteins_3d/images_Point_Cloud"
        protein_cfg[
            "combined_image"
        ],  # "results/proteins_3d/all_proteins_combined.png"
    )
    logger.info("Pipeline complete.")