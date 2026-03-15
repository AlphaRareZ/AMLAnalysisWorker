import os
import math
import requests
import pandas as pd
import logging
import gc
import matplotlib
from .utils import ensure_dir_for_file
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.image as mpimg

logger = logging.getLogger(__name__)

def fetch_alphafold(csv_file, out_dir, out_csv_report):
    logger.info("--------------------Fetching Structures--------------------")
    os.makedirs(out_dir, exist_ok=True)
    ensure_dir_for_file(out_csv_report)

    try:
        df = pd.read_csv(csv_file)
    except FileNotFoundError:
        logger.info(f"[ERROR] Could not find input file: {csv_file}")
        return pd.DataFrame() 

    if "Entry" not in df.columns:
        logger.info("[ERROR] CSV must contain 'Entry' column with UniProt accessions")
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
            logger.info(f"[SKIP] No UniProt Entry for gene {gene}")
            results.append(
                {
                    "gene": gene,
                    "accession": "",
                    "pdb_file": "",
                    "status": "missing_accession",
                }
            )
            continue

        name = gene # if gene != accession else accession
        url = f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v6.pdb"
        save_path = os.path.join(out_dir, f"{name}.pdb")

        try:
            r = requests.get(url)
            if r.status_code == 200:
                with open(save_path, "wb") as f:
                    f.write(r.content)
                logger.info(f"[OK] Downloaded {gene} ({accession}) -> {save_path}")
                results.append(
                    {
                        "gene": gene,
                        "accession": accession,
                        "pdb_file": save_path,
                        "status": "downloaded",
                    }
                )
            else:
                logger.info(
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
            logger.info(f"[ERROR] Request failed for {gene} ({accession}): {e}")
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
    logger.info(f" Structure fetch complete -> {out_csv_report}")
    
    # Cleanup input df
    del df
    gc.collect()
    
    return results_df


def render_proteins(*args, **kwargs):
    logger.info(" Skipping rendering — py3Dmol not available in script environment.")
    return


def combine_images(image_dir, out_file):
    logger.info("--------------------Combining Images--------------------")
    if not os.path.exists(image_dir):
        logger.info(f" Image directory not found, skipping: {image_dir}")
        return

    images = [
        os.path.join(image_dir, f) for f in os.listdir(image_dir) if f.endswith(".png")
    ]
    images.sort()
    n = len(images)

    if n == 0:
        logger.info(" No .png images found to combine.")
        return

    cols = int(math.ceil(math.sqrt(n)))
    rows = math.ceil(n / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 4, rows * 4), squeeze=False)
    axes = axes.flatten() 

    for i, img_path in enumerate(images):
        ax = axes[i]
        try:
            img = mpimg.imread(img_path)
            ax.imshow(img)
            ax.axis("off")
            name = os.path.splitext(os.path.basename(img_path))[0]
            ax.set_title(name, fontsize=9)
        except Exception as e:
            logger.info(f"Failed to read image {img_path}: {e}")
            ax.axis("off")

    for j in range(i + 1, len(axes)):
        axes[j].axis("off")

    plt.tight_layout()

    ensure_dir_for_file(out_file) 
    plt.savefig(out_file, dpi=200)
    
    # Extensive cleanup to prevent matplotlib memory leaks
    plt.close('all')
    del fig, axes, images
    gc.collect()
    
    logger.info(f" Combined image saved -> {out_file}")