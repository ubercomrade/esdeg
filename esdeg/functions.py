"""Motif enrichment, GC matching, and permutation statistics."""

from __future__ import annotations

import logging
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from importlib import resources

import mimosa
import numpy as np
import pandas as pd
from numba import njit

from esdeg.parsers import (
    promoters_parser,
    read_gene_set,
    read_model_records,
    read_motifs_from_db,
    read_table,
)

logger = logging.getLogger(__name__)
DEFAULT_BACKGROUND_LFC = np.log2(5 / 4)


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


def _validate_gc(values, name: str) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    if values.ndim != 1 or values.size == 0:
        raise ValueError(f"{name} must be a non-empty one-dimensional array.")
    if not np.all(np.isfinite(values)) or np.any((values < 0) | (values > 1)):
        raise ValueError(f"{name} must contain finite values in [0, 1].")
    return values


def _largest_remainder(counts: np.ndarray, total: int) -> np.ndarray:
    ideal = total * counts / counts.sum()
    quota = np.floor(ideal).astype(np.int64)
    remainder = total - int(quota.sum())
    if remainder:
        order = np.argsort(-(ideal - quota), kind="stable")
        quota[order[:remainder]] += 1
    return quota


def select_gc_matched_background(
    foreground_gc: np.ndarray,
    background_gc: np.ndarray,
    match_ratio: int = 5,
    random_state: int | None = None,
) -> np.ndarray:
    """Select a unique background pool with a robust foreground GC match."""
    foreground_gc = _validate_gc(foreground_gc, "foreground_gc")
    background_gc = _validate_gc(background_gc, "background_gc")
    if isinstance(match_ratio, bool) or not isinstance(match_ratio, (int, np.integer)):
        raise ValueError("match_ratio must be a positive integer.")
    if match_ratio <= 0:
        raise ValueError("match_ratio must be a positive integer.")

    target = min(background_gc.size, int(match_ratio) * foreground_gc.size)
    n_bins = max(1, min(10, foreground_gc.size // 50))
    quantiles = np.quantile(foreground_gc, np.linspace(0, 1, n_bins + 1))
    inner_edges = np.unique(quantiles[1:-1])
    foreground_bins = np.searchsorted(inner_edges, foreground_gc, side="right")
    background_bins = np.searchsorted(inner_edges, background_gc, side="right")
    n_bins = len(inner_edges) + 1

    counts = np.bincount(foreground_bins, minlength=n_bins)
    quota = _largest_remainder(counts, target)
    rng = np.random.default_rng(random_state)
    selected: list[int] = []
    available = np.ones(background_gc.size, dtype=bool)

    for bin_index, needed in enumerate(quota):
        if needed == 0:
            continue
        local = np.flatnonzero(available & (background_bins == bin_index))
        take = min(int(needed), local.size)
        if take:
            chosen = rng.choice(local, size=take, replace=False)
            selected.extend(int(index) for index in chosen)
            available[chosen] = False

        deficit = int(needed) - take
        if deficit:
            foreground_in_bin = foreground_gc[foreground_bins == bin_index]
            center = (
                float(foreground_in_bin.mean())
                if foreground_in_bin.size
                else float(background_gc.mean())
            )
            candidates = np.flatnonzero(available)
            order = np.argsort(np.abs(background_gc[candidates] - center), kind="stable")
            chosen = candidates[order[:deficit]]
            selected.extend(int(index) for index in chosen)
            available[chosen] = False

    result = np.asarray(selected, dtype=np.int32)
    if result.size != target or np.unique(result).size != result.size:
        raise ValueError("Could not select a unique GC-matched background pool.")
    logger.info(
        "GC-matched background: n=%d, foreground_gc=%.4f, background_gc=%.4f, bins=%d",
        result.size,
        foreground_gc.mean(),
        background_gc[result].mean(),
        n_bins,
    )
    return result


@dataclass(frozen=True)
class PermutationResult:
    auc_roc: float
    auc_prc: float
    p_value_roc: float
    p_value_prc: float


@njit(cache=True)
def precision_recall_curve(classification, scores):
    """Prevalence-adjusted precision-recall curve used by ESDEG."""
    if len(scores) == 0:
        return np.array([1.0]), np.array([0.0]), np.array([np.inf])
    indexes = np.argsort(scores)[::-1]
    sorted_scores = scores[indexes]
    sorted_classification = classification[indexes]
    max_size = np.unique(scores).shape[0] + 1
    precision = np.zeros(max_size)
    recall = np.zeros(max_size)
    thresholds = np.zeros(max_size)
    precision[0] = 1.0
    thresholds[0] = np.inf

    true_count = np.sum(classification == 1)
    false_count = np.sum(classification == 0)
    true_false_ratio = 1.0 if false_count == 0 else true_count / false_count
    true_positive = 0
    false_positive = 0
    position = 1
    score = sorted_scores[0]
    for index in range(len(scores)):
        if sorted_classification[index] == 1:
            true_positive += 1
        else:
            false_positive += 1
        if index == len(scores) - 1 or score != sorted_scores[index + 1]:
            thresholds[position] = sorted_scores[index]
            precision[position] = true_positive / (
                true_positive + true_false_ratio * false_positive
            )
            recall[position] = 0.0 if true_count == 0 else true_positive / true_count
            position += 1
            if index < len(scores) - 1:
                score = sorted_scores[index + 1]
    return precision[:position], recall[:position], thresholds[:position]


@njit(cache=True)
def roc_curve(classification, scores):
    """ROC curve with score ties aggregated before integration."""
    if len(scores) == 0:
        return np.array([0.0]), np.array([0.0]), np.array([np.inf])
    indexes = np.argsort(scores)[::-1]
    sorted_scores = scores[indexes]
    sorted_classification = classification[indexes]
    max_size = np.unique(scores).shape[0] + 1
    tpr = np.zeros(max_size)
    fpr = np.zeros(max_size)
    thresholds = np.zeros(max_size)
    thresholds[0] = np.inf
    true_count = np.sum(classification == 1)
    false_count = np.sum(classification == 0)
    true_positive = 0
    false_positive = 0
    position = 1
    score = sorted_scores[0]
    for index in range(len(scores)):
        if sorted_classification[index] == 1:
            true_positive += 1
        else:
            false_positive += 1
        if index == len(scores) - 1 or score != sorted_scores[index + 1]:
            thresholds[position] = sorted_scores[index]
            tpr[position] = 0.0 if true_count == 0 else true_positive / true_count
            fpr[position] = 0.0 if false_count == 0 else false_positive / false_count
            position += 1
            if index < len(scores) - 1:
                score = sorted_scores[index + 1]
    return tpr[:position], fpr[:position], thresholds[:position]


@njit(cache=True)
def compute_aucs(classification, scores):
    """Return ROC AUC and ESDEG's prevalence-adjusted PR AUC."""
    tpr, fpr, _ = roc_curve(classification, scores)
    auc_roc = np.trapezoid(tpr, x=fpr)
    precision, recall, _ = precision_recall_curve(classification, scores)
    auc_prc = np.trapezoid(precision, x=recall)
    return auc_roc, auc_prc


@njit(cache=True)
def permutation_loop(all_labels, all_scores, n_permutations, seed):
    """Run an upper-tail permutation loop with a private Numba RNG seed."""
    np.random.seed(seed)
    n_samples = len(all_labels)
    permuted_auc_roc = np.empty(n_permutations, dtype=np.float32)
    permuted_auc_prc = np.empty(n_permutations, dtype=np.float32)
    shuffled_labels = all_labels.copy()
    for permutation in range(n_permutations):
        for index in range(n_samples - 1, 0, -1):
            swap = np.random.randint(0, index + 1)
            shuffled_labels[index], shuffled_labels[swap] = (
                shuffled_labels[swap],
                shuffled_labels[index],
            )
        auc_roc, auc_prc = compute_aucs(shuffled_labels, all_scores)
        permuted_auc_roc[permutation] = auc_roc
        permuted_auc_prc[permutation] = auc_prc
    return permuted_auc_roc, permuted_auc_prc


def permutation_test(
    foreground_scores: np.ndarray,
    background_scores: np.ndarray,
    n_permutations: int = 10000,
    seed: int = 0,
    random_state: int | None = None,
) -> PermutationResult:
    """Calculate AUCs and one-sided upper-tail permutation p-values."""
    if random_state is not None:
        seed = random_state
    if isinstance(n_permutations, bool) or n_permutations < 1:
        raise ValueError("n_permutations must be at least 1.")
    foreground_scores = np.asarray(foreground_scores, dtype=np.float64)
    background_scores = np.asarray(background_scores, dtype=np.float64)
    if foreground_scores.ndim != 1 or background_scores.ndim != 1:
        raise ValueError("foreground_scores and background_scores must be one-dimensional.")
    if not foreground_scores.size or not background_scores.size:
        raise ValueError("foreground and background groups must both be non-empty.")
    if not np.all(np.isfinite(foreground_scores)) or not np.all(np.isfinite(background_scores)):
        raise ValueError("scores must be finite.")
    if not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("seed must be a non-negative integer.")

    all_scores = np.concatenate((foreground_scores, background_scores))
    all_labels = np.concatenate(
        (
            np.ones(foreground_scores.size, dtype=np.float64),
            np.zeros(background_scores.size, dtype=np.float64),
        )
    )
    observed_auc_roc, observed_auc_prc = compute_aucs(all_labels, all_scores)
    permuted_auc_roc, permuted_auc_prc = permutation_loop(
        all_labels, all_scores, int(n_permutations), int(seed) % (2**32)
    )
    return PermutationResult(
        auc_roc=float(observed_auc_roc),
        auc_prc=float(observed_auc_prc),
        p_value_roc=float(
            (np.count_nonzero(permuted_auc_roc >= observed_auc_roc) + 1) / (n_permutations + 1)
        ),
        p_value_prc=float(
            (np.count_nonzero(permuted_auc_prc >= observed_auc_prc) + 1) / (n_permutations + 1)
        ),
    )


def process_single_motif(args) -> dict:
    motif_data, foreground, background, n_permutations, seed = args
    motif_id = motif_data["motif_id"]
    try:
        result = permutation_test(foreground, background, n_permutations, seed=seed)
    except Exception as exc:
        raise RuntimeError(f"Permutation failed for motif_id {motif_id}: {exc}") from exc
    return {
        "motif_id": motif_id,
        "tf_name": motif_data["tf_name"],
        "tf_class": motif_data["tf_class"],
        "tf_family": motif_data["tf_family"],
        "auc_roc": result.auc_roc,
        "auc_prc": result.auc_prc,
        "p_value_roc": result.p_value_roc,
        "p_value_prc": result.p_value_prc,
    }


def parallel_permutation_runner(
    motifs: list[dict],
    foreground_scores: np.ndarray,
    background_scores: np.ndarray,
    n_permutations: int = 10000,
    n_workers: int = 1,
    seeds: list[int] | None = None,
) -> list[dict]:
    """Run one independently seeded permutation stream per motif."""
    if isinstance(n_workers, bool) or n_workers < 1:
        raise ValueError("nproc must be at least 1.")
    if n_permutations < 1:
        raise ValueError("n_permutations must be at least 1.")
    if foreground_scores.ndim != 2 or background_scores.ndim != 2:
        raise ValueError("score matrices must be two-dimensional.")
    if foreground_scores.shape[1] != len(motifs) or background_scores.shape[1] != len(motifs):
        raise ValueError("score matrix columns must match motif count.")
    if seeds is None:
        seeds = list(range(len(motifs)))
    if len(seeds) != len(motifs):
        raise ValueError("one RNG seed is required per motif.")
    if not motifs:
        return []

    args = [
        (
            motif,
            foreground_scores[:, index],
            background_scores[:, index],
            n_permutations,
            seeds[index],
        )
        for index, motif in enumerate(motifs)
    ]
    workers = min(int(n_workers), os.cpu_count() or 1, len(args))
    if workers == 1:
        return [process_single_motif(item) for item in args]

    results = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(process_single_motif, item): item[0]["motif_id"] for item in args
        }
        try:
            for future in as_completed(futures):
                results.append(future.result())
        except Exception as exc:
            for future in futures:
                future.cancel()
            motif_id = futures.get(future, "unknown")
            if "motif_id" in str(exc):
                raise
            raise RuntimeError(f"Permutation failed for motif_id {motif_id}: {exc}") from exc
    return results


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
    nproc: int = 4,
    type_of_data: str = "deg",
    condition: str = "all",
    log2fc_thr_deg: float = 1.0,
    log2fc_thr_background: float = np.log2(5 / 4),
    padj_thr: float = 0.05,
    n_permutations: int = 10000,
    match_ratio: int = 5,
    seed: int = 0,
    model_paths=None,
) -> pd.DataFrame:
    """Run motif enrichment on a DEG table or a gene set."""
    if path_to_promoters is None or path_to_data is None:
        raise ValueError("path_to_promoters and path_to_data are required.")
    if nproc < 1:
        raise ValueError("nproc must be at least 1.")
    if n_permutations < 1:
        raise ValueError("n_permutations must be at least 1.")
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
        "Input: seed=%d, n_permutations=%d, match_ratio=%d, models=%d, "
        "foreground=%d, background_pool=%d",
        seed,
        n_permutations,
        match_ratio,
        len(motifs),
        foreground_idx.size,
        background_pool_idx.size,
    )

    scores = process_all_models(sequences, motifs)
    gc_content = calculate_gc(sequences)
    foreground_gc = gc_content[foreground_idx]
    background_pool_gc = gc_content[background_pool_idx]
    streams = np.random.SeedSequence(int(seed)).spawn(len(motifs) + 1)
    gc_seed = int(streams[0].generate_state(1, dtype=np.uint32)[0])
    selected_bg_idx = select_gc_matched_background(
        foreground_gc, background_pool_gc, match_ratio=match_ratio, random_state=gc_seed
    )
    foreground_scores = scores[foreground_idx]
    background_scores = scores[background_pool_idx[selected_bg_idx]]
    motif_seeds = [int(stream.generate_state(1, dtype=np.uint32)[0]) for stream in streams[1:]]
    results = parallel_permutation_runner(
        motifs,
        foreground_scores,
        background_scores,
        n_permutations=n_permutations,
        n_workers=nproc,
        seeds=motif_seeds,
    )

    result = pd.DataFrame(results)
    result["p_value_roc_adj"] = _adjust_pvalues(result["p_value_roc"])
    result["p_value_prc_adj"] = _adjust_pvalues(result["p_value_prc"])
    if motif_db == "jaspar" and not model_paths:
        cluster_path = resources.files("esdeg").joinpath(f"clusters/{taxon}.tsv")
        result["jaspar_cluster"] = result["motif_id"].map(get_motif_to_cluster(str(cluster_path)))
    return result.sort_values(
        ["auc_roc", "motif_id"], ascending=[False, True], kind="mergesort"
    ).reset_index(drop=True)
