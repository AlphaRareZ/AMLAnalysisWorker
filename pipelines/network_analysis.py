import numpy as np
import pandas as pd
import networkx as nx
from networkx.algorithms import community
from sklearn.decomposition import PCA
from sklearn.preprocessing import minmax_scale
from logging import getLogger
import json
from .utils import ensure_dir_for_file
from mygene import MyGeneInfo

logger = getLogger(__name__)
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
    logger.info("--------------------Rank Hubs and Annotate--------------------")
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
        logger.info(f"️Error loading input files: {e}")
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
    mg = MyGeneInfo()
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
        logger.info(f"MyGene query failed: {e}")
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

    logger.info(f" Saved top {top_n} biomarkers -> {out_file}")
    logger.info(f" Full ranking saved -> {all_ranked_file}")
    return final

