import os
import csv
from typing import Dict, List, Union, Optional


def get_top_10_rows_from_output(
    output_dir: str = "Output",
    dest: Optional[Dict[str, Union[List[dict], dict]]] = None,
) -> Dict[str, Union[List[dict], dict]]:
    """Traverse `output_dir`, find CSV files, and add their top 10 rows to
    the provided `dest` dictionary (or a new dict if none provided).

    Each file produces a key named "<file_name>_top_10_rows" whose value is a
    list of up to 10 row dictionaries. On read error the value will be
    `{"error": "..."}`.

    Args:
            output_dir: Path to folder containing CSV files.
            dest: Optional dict to update; function returns this dict.

    Returns:
            The `dest` dictionary containing added entries.
    """
    if dest is None:
        dest = {}

    if not os.path.isdir(output_dir):
        return dest

    for root, _dirs, files in os.walk(output_dir):
        for fname in files:
            if not fname.lower().endswith(".csv"):
                continue

            path = os.path.join(root, fname)
            key = f"{os.path.splitext(fname)[0]}_top_10_rows"
            rows: List[dict] = []

            try:
                with open(path, newline="", encoding="utf-8") as f:
                    reader = csv.DictReader(f)
                    for i, row in enumerate(reader):
                        rows.append(row)
                        if i >= 9:
                            break
            except Exception as e:
                dest[key] = {"error": str(e)}
            else:
                dest[key] = rows

    return dest
