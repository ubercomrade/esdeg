"""Motif enrichment, composition-matched background, and Mann-Whitney statistics."""

from __future__ import annotations

import logging
from importlib import resources

import mimosa
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from esdeg.parsers import (
    promoters_parser,
    read_gene_set,
    read_model_records,
    read_motifs_from_db,
    read_table,
)

logger = logging.getLogger(__name__)
DEFAULT_BACKGROUND_LFC = np.log2(5 / 4)
MAX_SMD = 0.1


def calculate_gc(sequences) -> np.ndarray:
    """Return per-sequence GC fractions from Mimosa codes or DNA strings."""
    if hasattr(sequences, "data") and hasattr(sequences, "offsets"):
        data = np.asarray(sequences.data)
        offsets = np.asarray(sequences.offsets)
        lengths = np.diff(offsets)
        if np.any(lengths == 0):
            raise ValueError("GC content cannot be calculated for an empty sequence.")
        gc = np.add.reduceat((data == 1) | (data == 2), offsets[:-1])
        return gc.astype(np.float64) / lengths

    strings = list(sequences)
    if not strings or any(not sequence for sequence in strings):
        raise ValueError("GC content requires at least one non-empty sequence.")
    return np.asarray(
        [sum(base in "CG" for base in sequence.upper()) / len(sequence) for sequence in strings],
        dtype=np.float64,
    )


def calculate_cpg_log_oe(sequences) -> np.ndarray:
    """Return per-sequence log2 CpG observed/expected ratios with a pseudocount of one."""
    data = np.asarray(sequences.data)
    offsets = np.asarray(sequences.offsets)
    lengths = np.diff(offsets)
    if np.any(lengths == 0):
        raise ValueError("CpG content cannot be calculated for an empty sequence.")
    is_c = data == 1
    is_g = data == 2
    n_c = np.add.reduceat(is_c, offsets[:-1]).astype(np.float64)
    n_g = np.add.reduceat(is_g, offsets[:-1]).astype(np.float64)
    cg = np.append(is_c[:-1] & is_g[1:], False)
    cg[offsets[1:-1] - 1] = False  # a dinucleotide must not span two sequences
    n_cg = np.add.reduceat(cg, offsets[:-1]).astype(np.float64)
    return np.log2((n_cg + 1.0) / (n_c * n_g / lengths + 1.0))


def _validate_features(values, name: str) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    if values.ndim != 2 or values.shape[0] == 0 or values.shape[1] == 0:
        raise ValueError(f"{name} must be a non-empty two-dimensional array.")
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{name} must contain finite values.")
    return values


def _standardized_difference(foreground: np.ndarray, selected: np.ndarray) -> np.ndarray:
    spread = np.sqrt((foreground.var(axis=0) + selected.var(axis=0)) / 2)
    difference = foreground.mean(axis=0) - selected.mean(axis=0)
    return np.divide(difference, spread, out=np.zeros_like(difference), where=spread > 0)


def select_matched_background(
    foreground_features: np.ndarray,
    background_features: np.ndarray,
    match_ratio: int = 3,
    random_state: int | None = None,
    feature_names: tuple[str, ...] = (),
) -> np.ndarray:
    """Select a unique background pool by nearest-neighbour matching on promoter features.

    Distances are Mahalanobis distances under the covariance of the foreground and the
    candidate pool. Matching runs in rounds: every foreground promoter takes its nearest
    unused candidate before any promoter takes a second one. If the standardized mean
    difference of any feature exceeds MAX_SMD, the ratio is lowered by one and matching repeats.
    """
    foreground = _validate_features(foreground_features, "foreground_features")
    pool = _validate_features(background_features, "background_features")
    if foreground.shape[1] != pool.shape[1]:
        raise ValueError("foreground and background features must have the same columns.")
    if isinstance(match_ratio, bool) or not isinstance(match_ratio, (int, np.integer)):
        raise ValueError("match_ratio must be a positive integer.")
    if match_ratio <= 0:
        raise ValueError("match_ratio must be a positive integer.")

    target = min(pool.shape[0], int(match_ratio) * foreground.shape[0])
    combined = np.vstack((foreground, pool))
    covariance = np.atleast_2d(np.cov(combined, rowvar=False))
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    keep = eigenvalues > 1e-12 * max(float(eigenvalues.max()), 1e-300)
    whitening = eigenvectors[:, keep] / np.sqrt(eigenvalues[keep])
    center = combined.mean(axis=0)
    foreground_white = (foreground - center) @ whitening
    pool_white = (pool - center) @ whitening

    rng = np.random.default_rng(random_state)
    available = np.ones(pool.shape[0], dtype=bool)
    selected: list[int] = []
    while len(selected) < target:
        for row in rng.permutation(foreground.shape[0]):
            if len(selected) == target:
                break
            distance = np.sum((pool_white - foreground_white[row]) ** 2, axis=1)
            distance[~available] = np.inf
            best = int(np.argmin(distance))
            available[best] = False
            selected.append(best)

    result = np.asarray(selected, dtype=np.int32)
    if result.size != target or np.unique(result).size != result.size:
        raise ValueError("Could not select a unique matched background pool.")
    names = feature_names or tuple(f"feature{index}" for index in range(foreground.shape[1]))
    smd = _standardized_difference(foreground, pool[result])
    if int(match_ratio) > 1 and np.any(np.abs(smd) > MAX_SMD):
        logger.info(
            "Background with match_ratio=%d is unbalanced (max |SMD|=%.3f); retrying with %d.",
            match_ratio,
            np.abs(smd).max(),
            match_ratio - 1,
        )
        return select_matched_background(
            foreground, pool, int(match_ratio) - 1, random_state, feature_names
        )
    logger.info(
        "Matched background: n=%d, %s",
        result.size,
        ", ".join(
            f"{name}: foreground={foreground[:, index].mean():.4f} "
            f"pool={pool[:, index].mean():.4f} background={pool[result, index].mean():.4f} "
            f"smd={smd[index]:+.3f}"
            for index, name in enumerate(names)
        ),
    )
    for name, value in zip(names, smd, strict=True):
        if abs(value) > MAX_SMD:
            logger.warning(
                "Background is poorly matched on %s (|SMD|=%.3f > %.1f); "
                "p-values may be affected by promoter composition.",
                name,
                abs(value),
                MAX_SMD,
            )
    return result


def mann_whitney_enrichment(
    foreground_scores: np.ndarray, background_scores: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Return ROC AUC and one-sided p-values for each score column.

    The AUC is the Mann-Whitney statistic U / (n1 * n2). The p-value is the upper tail of its
    normal approximation with tie and continuity corrections, used for every group size.
    """
    foreground_scores = np.asarray(foreground_scores, dtype=np.float64)
    background_scores = np.asarray(background_scores, dtype=np.float64)
    if foreground_scores.ndim != 2 or background_scores.ndim != 2:
        raise ValueError("score matrices must be two-dimensional.")
    if foreground_scores.shape[1] != background_scores.shape[1]:
        raise ValueError("score matrices must have the same number of columns.")
    if not foreground_scores.shape[0] or not background_scores.shape[0]:
        raise ValueError("foreground and background groups must both be non-empty.")
    if not np.all(np.isfinite(foreground_scores)) or not np.all(np.isfinite(background_scores)):
        raise ValueError("scores must be finite.")
    result = mannwhitneyu(
        foreground_scores,
        background_scores,
        alternative="greater",
        use_continuity=True,
        method="asymptotic",
        axis=0,
    )
    auc = result.statistic / (foreground_scores.shape[0] * background_scores.shape[0])
    # Columns with all scores equal have zero variance and no evidence of enrichment.
    pvalues = np.where(np.isnan(result.pvalue), 1.0, result.pvalue)
    return auc, pvalues


def _validate_deg_table(df: pd.DataFrame) -> pd.DataFrame:
    required = {"id", "log2FoldChange", "padj"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"DEG table is missing required columns: {', '.join(missing)}.")
    result = df.copy()
    for column in ("log2FoldChange", "padj"):
        numeric = pd.to_numeric(result[column], errors="coerce")
        invalid = result[column].notna() & numeric.isna()
        if invalid.any():
            raise ValueError(f"DEG column {column!r} contains non-numeric values.")
        result[column] = numeric
    dropped = int(result["padj"].isna().sum())
    if dropped:
        logger.warning("Ignoring %d DEG rows with NaN padj.", dropped)
    return result[result["padj"].notna()].copy()


def get_deg_gene_ids(df, condition="all", padj_threshold=0.05, log2fc_threshold=1.0):
    if condition not in {"all", "up", "down"}:
        raise ValueError("condition must be one of: all, up, down.")
    log2fc_threshold = abs(float(log2fc_threshold))
    significant = df[df["padj"] <= padj_threshold]
    if condition == "up":
        significant = significant[significant["log2FoldChange"] >= log2fc_threshold]
    elif condition == "down":
        significant = significant[significant["log2FoldChange"] <= -log2fc_threshold]
    else:
        significant = significant[significant["log2FoldChange"].abs() >= log2fc_threshold]
    return pd.unique(significant.loc[significant["id"].notna(), "id"].astype(str)).astype(object)


def get_background_gene_ids(df, padj_threshold=0.05, log2fc_threshold=DEFAULT_BACKGROUND_LFC):
    threshold = abs(float(log2fc_threshold))
    background = df[
        (df["log2FoldChange"] >= -threshold)
        & (df["log2FoldChange"] <= threshold)
        & (df["padj"] > padj_threshold)
    ]
    return pd.unique(background.loc[background["id"].notna(), "id"].astype(str)).astype(object)


def get_indexes(all_ids: np.ndarray, sub_ids: np.ndarray) -> np.ndarray:
    wanted = set(map(str, sub_ids))
    return np.asarray(
        [index for index, value in enumerate(all_ids) if str(value) in wanted], dtype=np.int32
    )


def get_motif_to_cluster(cluster_path: str) -> dict[str, str]:
    clusters = pd.read_csv(cluster_path, sep="\t")
    if "motif_ids" not in clusters.columns:
        raise ValueError("Cluster table must contain a motif_ids column.")
    mapping = {}
    for _, row in clusters.iterrows():
        for motif_id in str(row["motif_ids"]).split(","):
            mapping[motif_id.strip()] = str(row.iloc[0])
    return mapping


def _adjust_pvalues(pvalues) -> np.ndarray:
    pvalues = np.asarray(pvalues, dtype=np.float64)
    if not pvalues.size:
        return pvalues
    order = np.argsort(pvalues, kind="stable")
    adjusted = np.minimum.accumulate(
        (pvalues[order] * pvalues.size / np.arange(1, pvalues.size + 1))[::-1]
    )[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def process_all_models(sequences, motifs) -> np.ndarray:
    """Scan all motif records and return one maximum score per sequence."""
    scores = np.empty((len(sequences), len(motifs)), dtype=np.float32)
    for motif_index, record in enumerate(motifs):
        track = mimosa.scan(record["model"], sequences, strands="best")
        row_lengths = np.diff(track.offsets)
        if np.any(row_lengths == 0):
            raise ValueError(
                f"Mimosa returned an empty score row for motif_id {record['motif_id']}."
            )
        scores[:, motif_index] = np.maximum.reduceat(track.data, track.offsets[:-1])
    return scores


def esdeg(
    motif_db: str = "hocomoco",
    taxon: str = "human",
    path_to_promoters: str | None = None,
    path_to_data: str | None = None,
    type_of_data: str = "deg",
    condition: str = "all",
    log2fc_thr_deg: float = 1.0,
    log2fc_thr_background: float = np.log2(5 / 4),
    padj_thr: float = 0.05,
    match_ratio: int = 3,
    seed: int = 0,
    model_paths=None,
    background_path: str | None = None,
) -> pd.DataFrame:
    """Run motif enrichment on a DEG table or a gene set."""
    if path_to_promoters is None or path_to_data is None:
        raise ValueError("path_to_promoters and path_to_data are required.")
    if not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("seed must be a non-negative integer.")
    if type_of_data not in {"deg", "set"}:
        raise ValueError("type_of_data must be 'deg' or 'set'.")
    if condition not in {"all", "up", "down"}:
        raise ValueError("condition must be one of: all, up, down.")

    sequences, ids = promoters_parser(path_to_promoters)
    motifs = (
        read_model_records(model_paths) if model_paths else read_motifs_from_db(motif_db, taxon)
    )
    if not motifs:
        raise ValueError("No motif models were loaded.")
    max_window = max(mimosa.window_size(record["model"]) for record in motifs)
    lengths = np.diff(sequences.offsets)
    if np.any(lengths < max_window):
        shortest = int(lengths.min())
        raise ValueError(
            f"FASTA contains a sequence shorter than the largest motif window "
            f"({shortest} < {max_window})."
        )

    if type_of_data == "set":
        foreground_ids = read_gene_set(path_to_data)
        foreground_set = set(foreground_ids)
        background_ids = np.asarray(
            [item for item in ids if item not in foreground_set], dtype=object
        )
    else:
        deg = _validate_deg_table(read_table(path_to_data))
        foreground_ids = get_deg_gene_ids(deg, condition, padj_thr, log2fc_thr_deg)
        background_ids = get_background_gene_ids(deg, padj_thr, log2fc_thr_background)

    foreground_idx = get_indexes(ids, foreground_ids)
    background_pool_idx = get_indexes(ids, background_ids)
    if foreground_idx.size == 0:
        raise ValueError("Foreground IDs do not intersect the FASTA promoter IDs.")
    if background_pool_idx.size == 0:
        raise ValueError("Background pool does not intersect the FASTA promoter IDs.")
    logger.info(
        "Input: seed=%d, match_ratio=%d, models=%d, "
        "foreground=%d, background_pool=%d",
        seed,
        match_ratio,
        len(motifs),
        foreground_idx.size,
        background_pool_idx.size,
    )

    features = np.column_stack((calculate_gc(sequences), calculate_cpg_log_oe(sequences)))
    selected_bg_idx = select_matched_background(
        features[foreground_idx],
        features[background_pool_idx],
        match_ratio=match_ratio,
        random_state=int(seed),
        feature_names=("GC", "log2 CpG o/e"),
    )
    background_idx = background_pool_idx[selected_bg_idx]
    if background_path is not None:
        pd.DataFrame({"id": np.asarray(ids)[background_idx]}).to_csv(
            background_path, sep="\t", index=False
        )
    selected_idx = np.concatenate((foreground_idx, background_idx))
    selected_sequences = mimosa.EncodedSequences.from_rows(
        sequences[index] for index in selected_idx
    )
    scores = process_all_models(selected_sequences, motifs)
    auc_roc, p_value_roc = mann_whitney_enrichment(
        scores[: foreground_idx.size], scores[foreground_idx.size :]
    )

    result = pd.DataFrame(
        {
            "motif_id": [motif["motif_id"] for motif in motifs],
            "tf_name": [motif["tf_name"] for motif in motifs],
            "tf_class": [motif["tf_class"] for motif in motifs],
            "tf_family": [motif["tf_family"] for motif in motifs],
            "auc_roc": auc_roc,
            "p_value_roc": p_value_roc,
        }
    )
    result["p_value_roc_adj"] = _adjust_pvalues(result["p_value_roc"])
    if motif_db == "jaspar" and not model_paths:
        cluster_path = resources.files("esdeg").joinpath(f"clusters/{taxon}.tsv")
        result["jaspar_cluster"] = result["motif_id"].map(get_motif_to_cluster(str(cluster_path)))
    return result.sort_values(
        ["auc_roc", "motif_id"], ascending=[False, True], kind="mergesort"
    ).reset_index(drop=True)
