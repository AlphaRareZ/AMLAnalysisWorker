import os
import math
import json
import requests
import pandas as pd
import numpy as np

# ---------------------------------------------------------
# HEADLESS CONFIGURATION: Must happen before pyplot/seaborn
# ---------------------------------------------------------
import matplotlib

matplotlib.use("Agg")

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, minmax_scale
import umap.umap_ as umap
import networkx as nx
from networkx.algorithms import community
from scipy.cluster.hierarchy import linkage, dendrogram
import mygene
from biotite.structure.io.pdb import PDBFile
from biotite.structure import filter_amino_acids
from PIL import Image, ImageDraw, ImageFont
from datetime import datetime

# import pyvista as pv
# # ---------------------------------------------------------
# # PYVISTA HEADLESS CONFIGURATION
# # ---------------------------------------------------------
# pv.OFF_SCREEN = True
# try:
#     pv.start_xvfb()  # Starts a virtual framebuffer for 3D rendering
# except Exception as e:
#     print(
#         f"Note: Could not start xvfb programmatically. Ensure xvfb-run is used or ignore if rendering works. Error: {e}"
#     )


def ensure_dir_for_file(file_path):
    """Ensures the directory for a given file path exists."""
    directory = os.path.dirname(file_path)
    if directory and not os.path.exists(directory):
        os.makedirs(directory, exist_ok=True)


def log(*args, sep=" ", end="\n", flush=False):
    """Writes a timestamped message to the global log file."""
    message = sep.join(str(arg) for arg in args)
    timestamp = datetime.now().strftime("[%Y-%m-%d %H:%M:%S] ")

    # Note: Directory for log file is created in main()
    try:
        with open(LOG_FILE_PATH, "a", encoding="utf-8") as log_file:
            log_file.write(timestamp + message + end)
            if flush:
                log_file.flush()
    except Exception as e:
        print(f"Failed to write to log file {LOG_FILE_PATH}: {e}")


# def map_exons_to_genes(expression_file, mapping_file, output_file):
#     log("--------------------Mapping Exons To Genes--------------------")
#     expr = pd.read_csv(expression_file, sep="\t", index_col=0)

#     mapping = (
#         pd.read_csv(mapping_file, sep="\t")[["id", "gene"]].dropna().drop_duplicates()
#     )

#     mapping = mapping.rename(columns={"id": "exon"}).set_index("exon")
#     expr = expr.loc[expr.index.intersection(mapping.index)]
#     expr["gene"] = mapping.loc[expr.index, "gene"]

#     gene_expr = expr.groupby("gene").mean()

#     ensure_dir_for_file(output_file)  # Create dir if needed
#     gene_expr.to_csv(output_file)

#     log("Gene-level expression saved to:", output_file)
#     log("Final shape:", gene_expr.shape)

def map_exons_to_genes(expression_file, mapping_file, output_file):
    log("--------------------Mapping Exons To Genes--------------------")
    
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

    log("Gene-level expression saved to:", output_file)
    log("Final shape:", gene_expr.shape)

def filter_protein_coding(input_file, output_file, batch_size=1000):
    log("--------------------Filtering Protein Coding--------------------")
    mg = mygene.MyGeneInfo()
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
            log(f"MyGene query batch failed: {e}")

    gene_types = {
        q["query"]: q.get("type_of_gene", None) for q in query_results if "query" in q
    }
    keep_genes = [g for g, t in gene_types.items() if t == "protein-coding"]

    filtered = gene_expr.loc[gene_expr.index.intersection(keep_genes)]

    ensure_dir_for_file(output_file)  # Create dir if needed
    filtered.to_csv(output_file)

    log(f"Filtered {gene_expr.shape[0]} -> {filtered.shape[0]} (protein-coding only)")
    log(f"Saved to: {output_file}")


def normalize_log(expr):
    return np.log2(expr + 1.0)


def select_hvgs(expr_file, out_hvgs, out_stats, n_hvgs=2000, min_mean=0.5):
    log("--------------------Selecting HVGS--------------------")
    expr = pd.read_csv(expr_file, index_col=0)
    expr_norm = normalize_log(expr)

    mean_expr = expr_norm.mean(axis=1)
    expr_f = expr_norm.loc[mean_expr >= min_mean]
    log("After mean filter:", expr_f.shape)

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

    log(f"Saved top {n_hvgs} HVGs to {out_hvgs}")
    log(f"Saved HVG stats to {out_stats}")


# def build_adjacency(expr_df, power=6):
#     corr = expr_df.T.corr()
#     adj = corr.abs() ** power
#     np.fill_diagonal(adj.values, 0.0)
#     return adj


def build_adjacency(expr_df, power=6):
    corr = expr_df.T.corr()
    adj = corr.abs() ** power
    np.fill_diagonal(adj.to_numpy().copy(), 0.0)
    return adj


def module_detection(adj_matrix):
    G = nx.from_pandas_adjacency(adj_matrix)
    communities_generator = community.greedy_modularity_communities(G, weight="weight")
    communities = list(communities_generator)
    modules = {f"module_{i+1}": list(c) for i, c in enumerate(communities)}
    return modules, G


def module_eigengenes(expr_df, modules):
    eig = {}
    for mname, genes in modules.items():
        if len(genes) < 2:
            eig[mname] = pd.Series([0] * expr_df.shape[1], index=expr_df.columns)
            continue
        pca = PCA(n_components=1)
        comp = pca.fit_transform(expr_df.loc[genes].T).flatten()
        eig[mname] = pd.Series(comp, index=expr_df.columns)
    eig_df = pd.DataFrame(eig)
    return eig_df


def intramodular_connectivity(adj, modules):
    k_within = {}
    for mname, genes in modules.items():
        if not genes:
            continue
        sub = adj.loc[genes, genes]
        s = sub.sum(axis=1)
        for g in genes:
            k_within[g] = s.loc[g]
    return pd.Series(k_within)


def rank_and_annotate(
    hvgs_file,
    intramod_file,
    modules_json,
    gene_expr_file,
    out_file,
    all_ranked_file,
    top_n=20,
):
    log("--------------------Rank Hubs and Annotate--------------------")
    # --- Load files ---
    try:
        hvgs = pd.read_csv(hvgs_file, index_col=0)
        k_within = pd.read_csv(intramod_file, index_col=0)
        expr = pd.read_csv(gene_expr_file, index_col=0)
        with open(modules_json, "r") as f:
            modules_data = json.load(f)
        modules = pd.DataFrame(
            dict([(k, pd.Series(v)) for k, v in modules_data.items()])
        )
    except Exception as e:
        log(f"️Error loading input files: {e}")
        return

    if k_within.shape[1] == 1:
        k_within = k_within.iloc[:, 0]

    # --- Expression stats ---
    mean_expr = expr.mean(axis=1)
    var_expr = expr.var(axis=1)

    # --- Combine metrics ---
    genes = hvgs.index.intersection(expr.index)
    df = pd.DataFrame(index=genes)
    df["mean"] = mean_expr.reindex(df.index)
    df["variance"] = var_expr.reindex(df.index)
    df["kWithin"] = k_within.reindex(df.index).fillna(0.0)

    # --- Normalize metrics ---
    df["mean_s"] = minmax_scale(df["mean"].fillna(0))
    df["var_s"] = minmax_scale(df["variance"].fillna(0))
    df["k_s"] = minmax_scale(df["kWithin"].fillna(0))

    # --- Composite score ---
    df["score"] = 0.5 * df["k_s"] + 0.3 * df["var_s"] + 0.2 * df["mean_s"]
    df = df.sort_values("score", ascending=False)

    # --- Annotate top genes ---
    mg = mygene.MyGeneInfo()
    top = df.head(max(top_n, 200)).index.tolist()
    try:
        query = mg.querymany(
            top,
            scopes="symbol",
            fields="type_of_gene,uniprot,entrezgene,go,alias",
            species="human",
            verbose=False,
        )
    except Exception as e:
        log(f"MyGene query failed: {e}")
        query = []

    ann = {}
    for q in query:
        if "notfound" in q and q["notfound"]:
            continue
        name = q.get("query")
        if not name:
            continue

        uni_raw = q.get("uniprot")
        accession = None
        if isinstance(uni_raw, dict):
            accession = uni_raw.get("Swiss-Prot") or uni_raw.get("TrEMBL")
            if isinstance(accession, list):
                accession = accession[0] if accession else None
        elif isinstance(uni_raw, list):
            accession = uni_raw[0] if uni_raw else None
        elif isinstance(uni_raw, str):
            accession = uni_raw

        ann[name] = {
            "type_of_gene": q.get("type_of_gene"),
            "uniprot": str(uni_raw),  # Ensure serializable
            "Entry": accession,
            "entrez": q.get("entrezgene"),
            "go": q.get("go"),
            "alias": q.get("alias"),
        }

    ann_df = pd.DataFrame.from_dict(ann, orient="index")
    merged = df.merge(ann_df, left_index=True, right_index=True, how="left")

    # --- Targetability flag ---
    def flag_targetable(go):
        if not isinstance(go, dict):
            return ""
        cc = go.get("CC")
        terms = ""
        if isinstance(cc, list):
            terms = " ".join(
                [t.get("term", "").lower() for t in cc if isinstance(t, dict)]
            )
        elif isinstance(cc, dict):
            terms = cc.get("term", "").lower()
        if any(
            x in terms
            for x in ["membrane", "plasma membrane", "extracellular", "secreted"]
        ):
            return "surface/secreted"
        return ""

    merged["targetability"] = merged["go"].apply(flag_targetable)

    # --- Save outputs ---
    final = merged.head(top_n)

    ensure_dir_for_file(all_ranked_file)  # Create dir
    merged.to_csv(all_ranked_file)

    ensure_dir_for_file(out_file)  # Create dir
    final.to_csv(out_file)

    log(f" Saved top {top_n} biomarkers -> {out_file}")
    log(f" Full ranking saved -> {all_ranked_file}")
    return final


def fetch_alphafold(csv_file, out_dir, out_csv_report):
    log("--------------------Fetching Structures--------------------")
    # Ensure output directories exist
    os.makedirs(out_dir, exist_ok=True)
    ensure_dir_for_file(out_csv_report)

    try:
        df = pd.read_csv(csv_file)
    except FileNotFoundError:
        log(f"[ERROR] Could not find input file: {csv_file}")
        return pd.DataFrame()  # Return empty df

    if "Entry" not in df.columns:
        log("[ERROR] CSV must contain 'Entry' column with UniProt accessions")
        return pd.DataFrame()

    results = []
    for _, row in df.iterrows():
        accession = str(row["Entry"]).strip()
        gene_col = (
            "Unnamed: 0"
            if "Unnamed: 0" in row
            else (df.columns[0] if len(df.columns) > 0 else "gene")
        )
        gene = str(row.get(gene_col, accession)).strip()

        if not accession or accession == "nan":
            log(f"[SKIP] No UniProt Entry for gene {gene}")
            results.append(
                {
                    "gene": gene,
                    "accession": "",
                    "pdb_file": "",
                    "status": "missing_accession",
                }
            )
            continue

        name = gene if gene != accession else accession
        url = f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v6.pdb"
        save_path = os.path.join(out_dir, f"{name}.pdb")

        try:
            r = requests.get(url)
            if r.status_code == 200:
                with open(save_path, "wb") as f:
                    f.write(r.content)
                log(f"[OK] Downloaded {gene} ({accession}) -> {save_path}")
                results.append(
                    {
                        "gene": gene,
                        "accession": accession,
                        "pdb_file": save_path,
                        "status": "downloaded",
                    }
                )
            else:
                log(
                    f"[MISS] No structure for {gene} ({accession}) (Code: {r.status_code})"
                )
                results.append(
                    {
                        "gene": gene,
                        "accession": accession,
                        "pdb_file": "",
                        "status": "missing_on_server",
                    }
                )
        except requests.exceptions.RequestException as e:
            log(f"[ERROR] Request failed for {gene} ({accession}): {e}")
            results.append(
                {
                    "gene": gene,
                    "accession": accession,
                    "pdb_file": "",
                    "status": "request_failed",
                }
            )

    results_df = pd.DataFrame(results)
    results_df.to_csv(out_csv_report, index=False)
    log(f" Structure fetch complete -> {out_csv_report}")
    return results_df


def render_proteins(*args, **kwargs):
    log(" Skipping rendering — py3Dmol not available in script environment.")
    return


def combine_images(image_dir, out_file):
    log("--------------------Combining Images--------------------")
    if not os.path.exists(image_dir):
        log(f" Image directory not found, skipping: {image_dir}")
        return

    images = [
        os.path.join(image_dir, f) for f in os.listdir(image_dir) if f.endswith(".png")
    ]
    images.sort()
    n = len(images)

    if n == 0:
        log(" No .png images found to combine.")
        return

    cols = int(math.ceil(math.sqrt(n)))
    rows = math.ceil(n / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 4, rows * 4), squeeze=False)
    axes = axes.flatten()  # Flatten to 1D array for easy iteration

    for i, img_path in enumerate(images):
        ax = axes[i]
        try:
            img = mpimg.imread(img_path)
            ax.imshow(img)
            ax.axis("off")
            name = os.path.splitext(os.path.basename(img_path))[0]
            ax.set_title(name, fontsize=9)
        except Exception as e:
            log(f"Failed to read image {img_path}: {e}")
            ax.axis("off")

    # Hide any unused subplots
    for j in range(i + 1, len(axes)):
        axes[j].axis("off")

    plt.tight_layout()

    ensure_dir_for_file(out_file)  # Create dir
    plt.savefig(out_file, dpi=200)
    plt.close()
    log(f" Combined image saved -> {out_file}")


# Main Pipeline
def main(expression_file, mapping_file, config_path="config.json"):

    # --- 1. Load Configuration ---
    global LOG_FILE_PATH
    try:
        with open(config_path, "r") as f:
            config = json.load(f)
    except FileNotFoundError:
        log(f"FATAL: Configuration file not found at {config_path}")
        return
    except json.JSONDecodeError:
        log(f"FATAL: Configuration file {config_path} is not valid JSON.")
        return

    # Set log file path from config and ensure its directory exists
    LOG_FILE_PATH = config["files"]["log_file"]
    ensure_dir_for_file(LOG_FILE_PATH)

    log("Starting a New Process")
    log(f"Loaded configuration from {config_path}")

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
    log("--------------------Coexpression and Modules--------------------")
    try:
        expr = pd.read_csv(hvgs_file, index_col=0)
        log("HVG matrix shape:", expr.shape)

        adj = build_adjacency(expr, power=6)
        modules, G = module_detection(adj)
        log("Detected modules:", {k: len(v) for k, v in modules.items()})

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

        log("Saved modules, eigengenes, and intramodular connectivity.")
    except FileNotFoundError:
        log(f"Failed to read HVG file {hvgs_file}. Skipping module detection.")
        return  # Stop pipeline if this critical step fails
    except Exception as e:
        log(f"Error during module detection: {e}")
        return

    # --- 4. EDA For the filtered Genes (Simple Variance-based) ---
    log("--------------------EDA For the filtered Genes--------------------")
    try:
        expr = pd.read_csv(gene_expr_coding, index_col=0)
        gene_var = expr.var(axis=1).sort_values(ascending=False)
        filtered_genes = gene_var[
            ~gene_var.index.str.contains("HB|globin", case=False, regex=True)
        ]
        top20 = filtered_genes.head(20)
        top20_df = expr.loc[top20.index]

        simple_cfg = files_cfg["simple_biomarkers"]

        ensure_dir_for_file(simple_cfg["top20_csv"])
        top20_df.to_csv(simple_cfg["top20_csv"])
        log("\nTop 20 biomarker genes (variance-based, hemoglobins excluded):")
        log(top20)
        log(f"\nSaved to: {simple_cfg['top20_csv']}")

        summary_stats = top20_df.describe()
        ensure_dir_for_file(simple_cfg["summary"])
        summary_stats.to_csv(simple_cfg["summary"])

        plt.figure(figsize=(8, 5))
        sns.histplot(top20_df.to_numpy().copy().flatten(), bins=50, kde=True)
        plt.title("Expression Distribution (Top 20 Biomarkers)")
        ensure_dir_for_file(simple_cfg["distribution_png"])
        plt.savefig(simple_cfg["distribution_png"])
        plt.close()

        clust = sns.clustermap(
            top20_df,
            cmap="plasma",
            metric="correlation",
            method="average",
            standard_scale=1,
            figsize=(12, 8),
        )
        clust.fig.suptitle("Clustered Heatmap (Top 20 Biomarkers)", y=1.02)
        ensure_dir_for_file(simple_cfg["heatmap_png"])
        clust.savefig(simple_cfg["heatmap_png"])
        plt.close()

        pca = PCA(n_components=2)
        pca_res = pca.fit_transform(top20_df.T)
        pca_df = pd.DataFrame(pca_res, columns=["PC1", "PC2"])
        ensure_dir_for_file(simple_cfg["pca_csv"])
        pca_df.to_csv(simple_cfg["pca_csv"])

        plt.figure(figsize=(7, 6))
        sns.scatterplot(data=pca_df, x="PC1", y="PC2")
        plt.title("PCA of Samples (Top 20 Biomarkers)")
        ensure_dir_for_file(simple_cfg["pca_png"])
        plt.savefig(simple_cfg["pca_png"])
        plt.close()

        reducer = umap.UMAP(random_state=42)
        umap_res = reducer.fit_transform(top20_df.T)
        umap_df = pd.DataFrame(umap_res, columns=["UMAP1", "UMAP2"])
        ensure_dir_for_file(simple_cfg["umap_csv"])
        umap_df.to_csv(simple_cfg["umap_csv"])

        plt.figure(figsize=(7, 6))
        sns.scatterplot(data=umap_df, x="UMAP1", y="UMAP2")
        plt.title("UMAP of Samples (Top 20 Biomarkers)")
        ensure_dir_for_file(simple_cfg["umap_png"])
        plt.savefig(simple_cfg["umap_png"])
        plt.close()

        corr = top20_df.T.corr()
        G = nx.Graph()
        for g1 in corr.columns:
            for g2 in corr.columns:
                if g1 < g2 and corr.loc[g1, g2] > 0.7:
                    G.add_edge(g1, g2, weight=corr.loc[g1, g2])
        plt.figure(figsize=(10, 8))
        nx.draw_networkx(G, node_size=300, font_size=8)
        plt.title("Gene Co-expression Network (Top 20 Biomarkers)")
        ensure_dir_for_file(simple_cfg["network_png"])
        plt.savefig(simple_cfg["network_png"])
        plt.close()

    except Exception as e:
        log(f"Error during simple EDA: {e}")

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
        log("Top ranked biomarkers (WGCNA-like score):")
        log(
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
        log("Ranking and annotation step failed or returned empty.")
        # We might want to stop here if downstream steps depend on this
        # return

    # --- 6. Advanced Biomarkers EDA (based on WGCNA-ranked genes) ---
    log("--------------------Advanced Biomarkers EDA--------------------")
    try:
        expr = pd.read_csv(gene_expr_raw, index_col=0)  # full matrix
        biomarkers = pd.read_csv(ranked_cfg["top20_annotated"], index_col=0)
        genes = biomarkers.index.tolist()

        if not genes:
            raise ValueError("No biomarker genes found in annotated file.")

        adv_cfg = files_cfg["advanced_biomarkers"]

        # 1. Heatmap
        sns.clustermap(
            expr.loc[genes],
            cmap="RdBu_r",
            z_score=0,
            metric="correlation",
            method="average",
            figsize=(12, 8),
            xticklabels=False,
        )
        plt.title("Top 20 Ranked Biomarkers Heatmap", y=1.02)
        ensure_dir_for_file(adv_cfg["heatmap"])
        plt.savefig(adv_cfg["heatmap"])
        plt.close()

        # 2. PCA
        pca = PCA(n_components=2)
        scaled = StandardScaler().fit_transform(expr.loc[genes].T)
        pca_res = pca.fit_transform(scaled)
        pca_df = pd.DataFrame(pca_res, columns=["PC1", "PC2"], index=expr.columns)
        plt.figure(figsize=(7, 6))
        sns.scatterplot(x="PC1", y="PC2", data=pca_df, s=80)
        plt.title("PCA using Top 20 Ranked Biomarkers")
        ensure_dir_for_file(adv_cfg["pca"])
        plt.savefig(adv_cfg["pca"])
        plt.close()

        # 3. Correlation Heatmap
        corr = expr.loc[genes].T.corr()
        plt.figure(figsize=(10, 8))
        sns.heatmap(
            corr,
            cmap="coolwarm",
            center=0,
            annot=True,
            fmt=".2f",
            annot_kws={"size": 8},
        )
        plt.title("Correlation Among Top 20 Ranked Biomarkers")
        plt.tight_layout()
        ensure_dir_for_file(adv_cfg["correlation_heatmap"])
        plt.savefig(adv_cfg["correlation_heatmap"])
        plt.close()

        # 4. Biomarker Network
        G = nx.Graph()
        for g1 in corr.columns:
            for g2 in corr.columns:
                if g1 < g2 and abs(corr.loc[g1, g2]) > 0.7:
                    G.add_edge(g1, g2, weight=corr.loc[g1, g2])
        plt.figure(figsize=(8, 6))
        nx.draw_networkx(
            G, with_labels=True, node_size=800, node_color="skyblue", font_size=8
        )
        plt.title("Top 20 Biomarker Network (Corr > 0.7)")
        ensure_dir_for_file(adv_cfg["network"])
        plt.savefig(adv_cfg["network"])
        plt.close()

        # ... (Other plots) ...

        # 6. Boxplots
        plt.figure(figsize=(14, 6))
        expr.loc[genes].T.boxplot(rot=90)
        plt.ylabel("Expression")
        plt.title("Boxplots of Top 20 Biomarkers Across Samples")
        plt.tight_layout()
        ensure_dir_for_file(adv_cfg["boxplots"])
        plt.savefig(adv_cfg["boxplots"])
        plt.close()

        # 7. Violin plots
        plt.figure(figsize=(14, 6))
        sns.violinplot(data=expr.loc[genes].T)
        plt.xticks(rotation=90)
        plt.ylabel("Expression")
        plt.title("Violin Plots of Top 20 Biomarkers")
        plt.tight_layout()
        ensure_dir_for_file(adv_cfg["violins"])
        plt.savefig(adv_cfg["violins"])
        plt.close()

        # 9. Hierarchical Dendrogram (Samples)
        Z = linkage(expr.loc[genes].T, method="average", metric="euclidean")
        plt.figure(figsize=(10, 6))
        dendrogram(Z, labels=expr.columns, leaf_rotation=90, leaf_font_size=8)
        plt.title("Hierarchical Clustering of Samples (Top 20 Biomarkers)")
        plt.tight_layout()
        ensure_dir_for_file(adv_cfg["dendrogram"])
        plt.savefig(adv_cfg["dendrogram"])
        plt.close()

        # 10. UMAP Projection
        reducer = umap.UMAP(random_state=42)
        umap_res = reducer.fit_transform(expr.loc[genes].T)
        umap_df = pd.DataFrame(umap_res, columns=["UMAP1", "UMAP2"], index=expr.columns)
        plt.figure(figsize=(7, 6))
        sns.scatterplot(x="UMAP1", y="UMAP2", data=umap_df, s=80)
        plt.title("UMAP of Samples (Top 20 Biomarkers)")
        ensure_dir_for_file(adv_cfg["umap"])
        plt.savefig(adv_cfg["umap"])
        plt.close()

        log("All advanced biomarker plots generated in Biomarkers_results_graphs/")

    except FileNotFoundError:
        log(
            f"Skipping advanced EDA: Could not find {ranked_cfg['top20_annotated']} or {gene_expr_raw}"
        )
    except Exception as e:
        log(f"Error during advanced biomarker EDA: {e}")

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

    # log(
    #     "--------------------Viewing & Exporting Proteins (PyVista)--------------------"
    # )

    # csv_file = protein_cfg["fetch_report"]
    # pdb_dir = dirs_cfg["protein_structures"]
    # img_dir = dirs_cfg["protein_images"]  # This is the *output* dir for this step
    # collage_file = protein_cfg["collage"]

    # os.makedirs(img_dir, exist_ok=True)  # Ensure output dir exists

    # try:
    #     df = pd.read_csv(csv_file)
    # except FileNotFoundError:
    #     log(f"Cannot generate protein images, report file not found: {csv_file}")
    #     df = pd.DataFrame()  # empty

    # saved_images = []

    # for _, row in df.iterrows():
    #     if row.get("status") != "downloaded":
    #         log(f" Skipping {row.get('gene')} (status: {row.get('status')})")
    #         continue

    #     protein_id = str(row.get("gene", row.get("accession", "unknown")))
    #     # Use the PDB file path from the report
    #     pdb_file = row.get("pdb_file")

    #     if not pdb_file or not os.path.exists(pdb_file):
    #         log(f" Skipping {protein_id} (PDB file not found at {pdb_file})")
    #         continue

    #     try:
    #         pdb = PDBFile.read(pdb_file)
    #         atoms = pdb.get_structure()[0]
    #         atoms = atoms[filter_amino_acids(atoms)]
    #         coords = atoms.coord

    #         if "bfactor" in atoms.get_annotation_categories():
    #             scalars = atoms.bfactor
    #             cmap = "plasma"
    #             clim = [0, 100]
    #         else:
    #             scalars = np.arange(len(coords))
    #             cmap = "rainbow"
    #             clim = None

    #         plotter = pv.Plotter(off_screen=True)
    #         plotter.add_points(
    #             coords,
    #             scalars=scalars,
    #             render_points_as_spheres=True,
    #             point_size=12,
    #             cmap=cmap,
    #             clim=clim,
    #         )

    #         # This is the output path for the individual image
    #         out_path = os.path.join(img_dir, f"{protein_id}.png")
    #         # We don't need ensure_dir_for_file; os.makedirs(img_dir) was called

    #         plotter.screenshot(out_path)
    #         plotter.close()

    #         img = Image.open(out_path)
    #         draw = ImageDraw.Draw(img)
    #         try:
    #             font = ImageFont.truetype("arial.ttf", 20)
    #         except IOError:
    #             font = ImageFont.load_default()

    #         text = protein_id
    #         try:
    #             bbox = draw.textbbox((0, 0), text, font=font)
    #             text_w, text_h = bbox[2] - bbox[0], bbox[3] - bbox[1]
    #         except AttributeError:
    #             # Fallback for older PIL versions
    #             text_w, text_h = draw.textsize(text, font=font)

    #         x = (img.width - text_w) // 2
    #         y = img.height - text_h - 5
    #         draw.text((x, y), text, fill=(0, 0, 0), font=font)
    #         img.save(out_path)

    #         saved_images.append(out_path)
    #         log(f" ✅ Saved image for {protein_id}")

    #     except Exception as e:
    #         log(f" ❌ Failed to render {protein_id}: {e}")

    # # --- Create Collage ---
    # if saved_images:
    #     imgs = [Image.open(p) for p in saved_images]
    #     n_cols = 5
    #     n_rows = (len(imgs) + n_cols - 1) // n_cols

    #     cell_w, cell_h = 300, 380
    #     collage = Image.new("RGB", (n_cols * cell_w, n_rows * cell_h), (255, 255, 255))
    #     draw = ImageDraw.Draw(collage)

    #     try:
    #         font = ImageFont.truetype("arial.ttf", 24)
    #     except IOError:
    #         font = ImageFont.load_default()

    #     for idx, img in enumerate(imgs):
    #         img = img.resize((300, 300))
    #         x = (idx % n_cols) * cell_w
    #         y = (idx // n_cols) * cell_h
    #         collage.paste(img, (x, y))

    #         protein_name = os.path.splitext(os.path.basename(saved_images[idx]))[0]
    #         try:
    #             bbox = draw.textbbox((0, 0), protein_name, font=font)
    #             text_w, text_h = bbox[2] - bbox[0], bbox[3] - bbox[1]
    #         except AttributeError:
    #             text_w, text_h = draw.textsize(protein_name, font=font)

    #         text_x = x + (cell_w - text_w) // 2
    #         text_y = y + 310
    #         draw.text((text_x, text_y), protein_name, fill=(0, 0, 0), font=font)

    #     ensure_dir_for_file(collage_file)  # Create dir
    #     collage.save(collage_file)
    #     log(f" 🖼️ Collage saved -> {collage_file}")
    # else:
    #     log(" No protein images generated, skipping collage.")

    # log("Process Ended ✅")
