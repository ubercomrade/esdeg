"""Expression annotation for AUC-based enrichment results."""

from __future__ import annotations

import logging
import re

import numpy as np
import pandas as pd

from esdeg.parsers import read_table

logger = logging.getLogger(__name__)
_ATTRIBUTE = re.compile(r"(?:^|;)\s*([^\s;]+)\s+(?:\"([^\"]*)\"|([^;\s]+))")


def _read_numeric_table(path, required, label) -> pd.DataFrame:
    table = read_table(path)
    missing = sorted(set(required) - set(table.columns))
    if missing:
        raise ValueError(f"{label} table is missing required columns: {', '.join(missing)}.")
    return table


def create_gname_converter(gtf_path):
    """Read gene_id and gene_name regardless of GTF attribute order."""
    gtf = pd.read_csv(gtf_path, sep="\t", comment="#", header=None, dtype=str)
    if gtf.shape[1] < 9:
        raise ValueError("GTF must contain at least nine tab-separated columns.")
    converter = {}
    for attributes in gtf.loc[gtf[2].eq("gene"), 8].dropna():
        parsed = {}
        for match in _ATTRIBUTE.finditer(str(attributes)):
            parsed[match.group(1)] = match.group(2) or match.group(3) or ""
        if parsed.get("gene_id") and parsed.get("gene_name"):
            converter[parsed["gene_id"]] = parsed["gene_name"].upper()
    return converter


def _split_metadata(value) -> list[str]:
    return [item.strip() for item in str(value).split("::") if item.strip()]


def _expand_dimers(table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in table.iterrows():
        values = {
            column: _split_metadata(row[column]) for column in ("tf_name", "tf_class", "tf_family")
        }
        target = max(map(len, values.values()))
        if target == 0:
            continue
        for column, parts in values.items():
            if len(parts) == 1 and target > 1:
                values[column] = parts * target
            elif len(parts) not in {1, target}:
                raise ValueError(
                    f"Incompatible dimer metadata lengths in {column}: {len(parts)} vs {target}."
                )
        for index in range(target):
            expanded = row.copy()
            for column in values:
                expanded[column] = (
                    values[column][index].upper() if column == "tf_name" else values[column][index]
                )
            rows.append(expanded)
    return pd.DataFrame(rows, columns=table.columns).reset_index(drop=True)


def read_esdeg_table(esdeg_path):
    table = read_table(esdeg_path)
    required = {"motif_id", "tf_name", "tf_class", "tf_family", "auc_roc", "p_value_roc_adj"}
    missing = sorted(required - set(table.columns))
    if missing:
        raise ValueError(f"Enrichment table is missing required columns: {', '.join(missing)}.")
    for column in ("auc_roc", "p_value_roc_adj"):
        numeric = pd.to_numeric(table[column], errors="coerce")
        if numeric.isna().any():
            raise ValueError(f"Enrichment column {column!r} must be numeric.")
        table[column] = numeric
    table["me_padj"] = table["p_value_roc_adj"]
    table = table[table["tf_name"].astype(str).str.upper().ne("NA")].copy()
    table = table[table["tf_family"].astype(str).str.upper().ne("NA")]
    return _expand_dimers(table)


def _expression_maps(deg_path, counts_path, id_to_name):
    deg = _read_numeric_table(deg_path, ("id", "log2FoldChange", "padj"), "DEG")
    counts = _read_numeric_table(counts_path, ("id", "counts"), "Counts")
    for table, columns, label in (
        (deg, ("log2FoldChange", "padj"), "DEG"),
        (counts, ("counts",), "Counts"),
    ):
        for column in columns:
            numeric = pd.to_numeric(table[column], errors="coerce")
            invalid = table[column].notna() & numeric.isna()
            if invalid.any():
                raise ValueError(f"{label} column {column!r} contains non-numeric values.")
            table[column] = numeric

    counts = counts[counts["id"].isin(id_to_name)].copy()
    counts["tf_name"] = counts["id"].map(id_to_name)
    counts_map = counts.groupby("tf_name")["counts"].max().to_dict()

    deg = deg[deg["id"].isin(id_to_name)].copy()
    deg["tf_name"] = deg["id"].map(id_to_name)
    deg["abs_lfc"] = deg["log2FoldChange"].abs()
    deg = deg.sort_values(
        ["tf_name", "padj", "abs_lfc"], ascending=[True, True, False], kind="mergesort"
    )
    deg_map = (
        deg.drop_duplicates("tf_name")
        .set_index("tf_name")[["log2FoldChange", "padj"]]
        .to_dict("index")
    )
    return counts_map, deg_map


def annotation(
    deg_path,
    counts_path,
    esdeg_path,
    gtf_path,
    filter_flag=False,
    best_flag=False,
    me_padj_thr=0.05,
    de_padj_thr=0.05,
    lfc_thr=1.0,
    counts_filter=5.0,
    min_auc=0.5,
):
    """Add expression columns and optionally filter AUC enrichment results."""
    id_to_name = create_gname_converter(gtf_path)
    table = read_esdeg_table(esdeg_path)
    counts_map, deg_map = _expression_maps(deg_path, counts_path, id_to_name)

    table["counts"] = table["tf_name"].map(counts_map).astype(float)
    table["lfc"] = table["tf_name"].map(
        lambda name: deg_map.get(name, {}).get("log2FoldChange", np.nan)
    )
    table["de_padj"] = table["tf_name"].map(lambda name: deg_map.get(name, {}).get("padj", np.nan))

    if filter_flag:
        motif_mask = (table["p_value_roc_adj"] < me_padj_thr) & (table["auc_roc"] >= min_auc)
        expression_mask = ((table["de_padj"] < de_padj_thr) & (table["lfc"].abs() >= lfc_thr)) | (
            table["counts"] >= counts_filter
        )
        table = table[motif_mask & expression_mask].copy()

    if best_flag:
        table = table.sort_values(
            ["tf_name", "p_value_roc_adj", "auc_roc", "motif_id"],
            ascending=[True, True, False, True],
            kind="mergesort",
        ).drop_duplicates("tf_name", keep="first")

    return table.reset_index(drop=True)
