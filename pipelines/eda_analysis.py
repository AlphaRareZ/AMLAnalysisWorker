# eda_analysis.py

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import umap
import networkx as nx
from scipy.cluster.hierarchy import linkage, dendrogram
from .utils import ensure_dir_for_file
from logging import getLogger

logger = getLogger(__name__)

def run_simple_eda(gene_expr_coding, files_cfg):
    # --- 4. EDA For the filtered Genes (Simple Variance-based) ---
    logger.info("--------------------EDA For the filtered Genes--------------------")
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
        logger.info("\nTop 20 biomarker genes (variance-based, hemoglobins excluded):")
        logger.info(top20)
        logger.info(f"\nSaved to: {simple_cfg['top20_csv']}")

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
        logger.info(f"Error during simple EDA: {e}")



def run_advanced_eda(gene_expr_raw, ranked_cfg, files_cfg):
    logger.info("--------------------Advanced Biomarkers EDA--------------------")
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

        logger.info("All advanced biomarker plots generated in Biomarkers_results_graphs/")

    except FileNotFoundError:
        logger.info(
            f"Skipping advanced EDA: Could not find {ranked_cfg['top20_annotated']} or {gene_expr_raw}"
        )
    except Exception as e:
        logger.info(f"Error during advanced biomarker EDA: {e}")
