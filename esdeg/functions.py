"""
Motif enrichment analysis with GC-stratified permutation testing.

This module provides functions for:
1. GC-matched background selection from promoter sequences
2. Stratified permutation testing for AUC significance
3. Parallel processing of multiple motifs
"""

import os
import sys
import collections
from typing import List, Dict, Tuple, Optional, Any
from dataclasses import dataclass

import numpy as np
import pandas as pd
import scipy.stats as st
import torch
import torch.nn.functional as F
from collections import defaultdict
from sklearn.metrics import auc, precision_recall_curve, roc_auc_score
from concurrent.futures import ProcessPoolExecutor, as_completed
from numba import njit
from tqdm import tqdm
from importlib import resources
from esdeg.parsers import promoters_parser, read_motifs_from_db


# ============================================================================
# Core utility functions
# ============================================================================

def calculate_gc(sequences: List[str]) -> np.ndarray:
    """Calculate GC content for each sequence.

    Args:
        sequences: List of DNA sequences (strings).

    Returns:
        Array of GC content values (0.0 to 1.0).
    """
    gc = []
    length = len(sequences[0])
    for seq in sequences:
        counter = collections.Counter(seq)
        gc.append((counter['C'] + counter['G']) / length)
    return np.array(gc)


def seq_to_int(sequences: List[str], device: str = 'cuda') -> torch.Tensor:
    """Convert DNA sequences to integer tensor.

    Args:
        sequences: List of DNA sequences.
        device: Target device ('cuda' or 'cpu').

    Returns:
        Integer tensor of shape [num_sequences, sequence_length].
    """
    converter = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    return torch.stack([
        torch.tensor([converter[c] for c in seq], dtype=torch.int64, device=device)
        for seq in sequences
    ])


# ============================================================================
# PWM scoring functions
# ============================================================================

def calculate_scores(
    numeric_sequences: torch.Tensor,
    models: torch.Tensor,
    chunk_size: int = 64
) -> torch.Tensor:
    """Calculate PWM scores for sequences using convolution.

    Args:
        numeric_sequences: Integer tensor [batch, length].
        models: PWM models [num_models, 5, motif_len].
        chunk_size: Batch size for processing.

    Returns:
        Score tensor [batch, num_models].
    """
    device = numeric_sequences.device
    B, L = numeric_sequences.shape
    M, C, motif_len = models.shape

    # Precompute reverse complement models
    rc_models = models[:, [3, 2, 1, 0, 4]].flip(-1)
    all_models = torch.cat([models, rc_models], dim=0)

    result = torch.zeros((B, M), device=device)

    # Process in chunks to manage memory
    for i in range(0, B, chunk_size):
        chunk_end = min(i + chunk_size, B)
        chunk = numeric_sequences[i:chunk_end]

        # One-hot encode
        one_hot = F.one_hot(chunk, 5).float().permute(0, 2, 1)

        # Convolve with all models (forward + reverse)
        all_scores = F.conv1d(one_hot, all_models)

        # Split and take best orientation
        forward_scores, rc_scores = all_scores.chunk(2, dim=1)
        best_scores = torch.maximum(forward_scores, rc_scores)
        chunk_result = torch.max(best_scores, dim=-1)[0]

        result[i:chunk_end] = chunk_result

    return result


def process_all_models(
    numeric_sequences: torch.Tensor,
    sorted_models: List[Dict[str, Any]],
    chunk_size: int = 64
) -> torch.Tensor:
    """Process all PWM models grouped by length.

    Args:
        numeric_sequences: Sequences as integer tensor [batch, length].
        sorted_models: List of model dictionaries with 'length' and 'pwm'.
        chunk_size: Batch size for scoring.

    Returns:
        Score tensor [batch, num_models].
    """
    device = numeric_sequences.device
    B = numeric_sequences.shape[0]

    # Group models by motif length for efficient processing
    length_groups = defaultdict(list)
    for orig_idx, model in enumerate(sorted_models):
        length = model['length']
        length_groups[length].append((orig_idx, model['pwm']))

    total_models = len(sorted_models)
    all_scores = torch.zeros((B, total_models), device=device)

    # Process each length group
    for length, models_group in length_groups.items():
        orig_indices = [idx for idx, _ in models_group]
        pwms = [pwm for _, pwm in models_group]

        models_tensor = torch.stack([
            torch.from_numpy(pwm).to(dtype=torch.float32)
            for pwm in pwms
        ]).to(device)

        group_scores = calculate_scores(numeric_sequences, models_tensor, chunk_size)
        all_scores[:, orig_indices] = group_scores

    return all_scores


# ============================================================================
# GC-matched background selection
# ============================================================================

def select_gc_matched_background(
    foreground_gc: np.ndarray,
    background_gc: np.ndarray,
    match_ratio: int = 10,
    n_quantiles: int = 10,
    random_state: Optional[int] = None
) -> np.ndarray:
    """Select background samples matched to foreground GC distribution.

    Args:
        foreground_gc: GC content array for foreground [n_foreground].
        background_gc: GC content array for background pool [n_background].
        match_ratio: Ratio of background to foreground samples.
        n_quantiles: Number of GC quantiles for stratification.
        random_state: Random seed for reproducibility.

    Returns:
        Indices of selected background samples.
    """
    if random_state is not None:
        np.random.seed(random_state)

    n_foreground = len(foreground_gc)
    n_quantiles = min(n_quantiles, n_foreground // 50)  # At least 50 per quantile

    # Define GC quantiles from foreground
    fg_quantiles = np.percentile(
        foreground_gc,
        np.linspace(0, 100, n_quantiles + 1)
    )

    selected_indices = []

    for i in range(n_quantiles):
        gc_low = fg_quantiles[i]
        gc_high = fg_quantiles[i + 1]

        # Count foreground in this quantile
        fg_mask = (foreground_gc >= gc_low) & (foreground_gc < gc_high)
        fg_count = np.sum(fg_mask)

        if fg_count == 0:
            continue

        # Find background in same GC range
        bg_mask = (background_gc >= gc_low) & (background_gc < gc_high)
        bg_indices_in_range = np.where(bg_mask)[0]

        bg_count_needed = fg_count * match_ratio

        # Handle edge cases
        if len(bg_indices_in_range) == 0:
            # No background in range, use nearest
            distances = np.abs(background_gc - np.mean([gc_low, gc_high]))
            bg_indices_in_range = np.argsort(distances)[:bg_count_needed]

        # Sample with or without replacement
        if len(bg_indices_in_range) >= bg_count_needed:
            selected = np.random.choice(
                bg_indices_in_range,
                size=bg_count_needed,
                replace=False
            )
        else:
            selected = np.random.choice(
                bg_indices_in_range,
                size=bg_count_needed,
                replace=True
            )

        selected_indices.extend(selected)

    return np.array(selected_indices, dtype=np.int32)


# ============================================================================
# Stratified permutation test
# ============================================================================


@dataclass
class PermutationResult:
    """Results from permutation test."""
    auc_roc: float
    auc_prc: float
    p_value_roc: float
    p_value_prc: float


@njit
def precision_recall_curve(classification, scores):
    """Compute precision-recall curve (JIT-compiled)."""
    if len(scores) == 0:
        return np.array([1.0]), np.array([0.0]), np.array([np.inf])

    # Получаем индексы сортировки оценок по убыванию
    indexes = np.argsort(scores)[::-1]
    sorted_scores = scores[indexes]
    sorted_classification = classification[indexes]

    # Инициализируем массивы (с запасом +1 для начальной точки)
    number_of_uniq_scores = np.unique(scores).shape[0]
    max_size = number_of_uniq_scores + 1

    precision = np.zeros(max_size)
    recall = np.zeros(max_size)
    uniq_scores = np.zeros(max_size)

    # Начальная точка: (recall=0, precision=1, threshold=inf)
    precision[0] = 1.0
    recall[0] = 0.0
    uniq_scores[0] = np.inf  # ← ИСПРАВЛЕНИЕ: было sorted_scores[0]

    TP, FP = 0, 0
    number_of_true = np.sum(classification == 1)
    number_of_false = np.sum(classification == 0)

    if number_of_false == 0:
        true_false_ratio = 1.0
    else:
        true_false_ratio = number_of_true / number_of_false

    position = 1
    score = sorted_scores[0]

    for i in range(len(scores)):
        _score = sorted_scores[i]
        _flag = sorted_classification[i]

        # Обновляем TP и FP
        if _flag == 1:
            TP += 1
        else:
            FP += 1

        # Проверяем, изменилась ли оценка
        if i == len(scores) - 1 or score != sorted_scores[i + 1]:
            uniq_scores[position] = _score

            if TP + FP > 0:
                precision[position] = TP / (TP + true_false_ratio * FP)
            else:
                precision[position] = 1.0

            if number_of_true > 0:
                recall[position] = TP / number_of_true
            else:
                recall[position] = 0.0

            position += 1
            if i < len(scores) - 1:
                score = sorted_scores[i + 1]

    return precision[:position], recall[:position], uniq_scores[:position]


@njit
def roc_curve(classification, scores):
    """Compute ROC curve (JIT-compiled)."""
    if len(scores) == 0:
        return np.array([0.0]), np.array([0.0]), np.array([np.inf])

    # Получаем индексы сортировки оценок по убыванию
    indexes = np.argsort(scores)[::-1]
    sorted_scores = scores[indexes]
    sorted_classification = classification[indexes]

    # Инициализируем массивы
    number_of_uniq_scores = np.unique(scores).shape[0]
    max_size = number_of_uniq_scores + 1

    tpr = np.zeros(max_size)
    fpr = np.zeros(max_size)
    uniq_scores = np.zeros(max_size)

    # Начальная точка: (fpr=0, tpr=0, threshold=inf)
    tpr[0] = 0.0
    fpr[0] = 0.0
    uniq_scores[0] = np.inf

    TP, FP = 0, 0
    number_of_true = np.sum(classification == 1)
    number_of_false = np.sum(classification == 0)
    position = 1
    score = sorted_scores[0]

    for i in range(len(scores)):
        _score = sorted_scores[i]
        _flag = sorted_classification[i]

        # Обновляем TP и FP
        if _flag == 1:
            TP += 1
        else:
            FP += 1

        # Проверяем, изменилась ли оценка
        if i == len(scores) - 1 or score != sorted_scores[i + 1]:
            uniq_scores[position] = _score

            if number_of_true > 0:
                tpr[position] = TP / number_of_true
            else:
                tpr[position] = 0.0

            if number_of_false > 0:
                fpr[position] = FP / number_of_false
            else:
                fpr[position] = 0.0

            position += 1
            if i < len(scores) - 1:
                score = sorted_scores[i + 1]

    return tpr[:position], fpr[:position], uniq_scores[:position]


@njit
def compute_aucs(classification, scores):
    """Compute both ROC AUC and PRC AUC for given labels and scores.

    Args:
        classification: Binary labels (0 or 1).
        scores: Predicted scores.

    Returns:
        Tuple of (auc_roc, auc_prc).
    """
    # Compute ROC curve and AUC
    tpr, fpr, _ = roc_curve(classification, scores)
    # Правильный порядок - np.trapz(y, x)
    # Для ROC: y=TPR (ось Y), x=FPR (ось X)
    auc_roc = np.trapz(tpr, x=fpr)

    # Compute PR curve and AUC
    precision, recall, _ = precision_recall_curve(classification, scores)
    # Правильный порядок - np.trapz(y, x)
    # Для PRC: y=Precision (ось Y), x=Recall (ось X)
    auc_prc = np.trapz(precision, x=recall)

    return auc_roc, auc_prc


@njit
def permutation_loop(
    all_labels: np.ndarray,
    all_scores: np.ndarray,
    n_permutations: int
) -> tuple:
    """JIT-compiled permutation loop for maximum speed.

    Args:
        all_labels: Binary labels for all samples.
        all_scores: Scores for all samples.
        n_permutations: Number of permutations to perform.

    Returns:
        Tuple of (permuted_auc_roc, permuted_auc_prc) arrays.
    """
    n_samples = len(all_labels)
    permuted_auc_roc = np.empty(n_permutations, dtype=np.float32)
    permuted_auc_prc = np.empty(n_permutations, dtype=np.float32)

    # Copy labels for shuffling
    shuffled_labels = all_labels.copy()

    for i in range(n_permutations):
        # In-place shuffle using Fisher-Yates algorithm
        for j in range(n_samples - 1, 0, -1):
            k = np.random.randint(0, j + 1)
            shuffled_labels[j], shuffled_labels[k] = shuffled_labels[k], shuffled_labels[j]

        # Compute both AUCs
        auc_roc, auc_prc = compute_aucs(shuffled_labels, all_scores)
        permuted_auc_roc[i] = auc_roc
        permuted_auc_prc[i] = auc_prc

    return permuted_auc_roc, permuted_auc_prc


def permutation_test(
    foreground_scores: np.ndarray,
    background_scores: np.ndarray,
    n_permutations: int = 5000,
    random_state: Optional[int] = None
) -> PermutationResult:
    """Perform permutation test computing both ROC AUC and PRC AUC.

    Optimized for speed using Numba JIT compilation. Use this function
    when background is already GC-matched to foreground (no stratification needed).

    Args:
        foreground_scores: Motif scores for foreground samples (positives).
        background_scores: Motif scores for background samples (negatives).
        n_permutations: Number of permutations for testing (default 5000).
        random_state: Random seed for reproducibility.

    Returns:
        PermutationResult with observed AUCs, p-values, and permutation means.

    Example:
        >>> fg = np.array([0.8, 0.9, 0.7, 0.85])
        >>> bg = np.array([0.3, 0.4, 0.5, 0.35, 0.45])
        >>> result = permutation_test(fg, bg, n_permutations=1000)
        >>> print(f"AUC ROC: {result.auc_roc:.3f}, p={result.p_value_roc:.4f}")
    """
    if random_state is not None:
        np.random.seed(random_state)

    # Combine and prepare data
    all_scores = np.concatenate([foreground_scores, background_scores])
    all_labels = np.concatenate([
        np.ones(len(foreground_scores), dtype=np.float64),
        np.zeros(len(background_scores), dtype=np.float64)
    ])

    # Compute observed AUCs
    observed_auc_roc, observed_auc_prc = compute_aucs(all_labels, all_scores)

    # Run permutation test (JIT-compiled loop)
    permuted_auc_roc, permuted_auc_prc = permutation_loop(
        all_labels, all_scores, n_permutations
    )

    # Calculate p-values
    count_roc = np.sum(permuted_auc_roc >= observed_auc_roc)
    p_value_roc = (count_roc + 1) / (n_permutations + 1)

    count_prc = np.sum(permuted_auc_prc >= observed_auc_prc)
    p_value_prc = (count_prc + 1) / (n_permutations + 1)

    return PermutationResult(
        auc_roc=float(observed_auc_roc),
        auc_prc=float(observed_auc_prc),
        p_value_roc=float(p_value_roc),
        p_value_prc=float(p_value_prc),
    )



# ============================================================================
# Parallel processing
# ============================================================================

def process_single_motif(args: Tuple) -> Dict[str, Any]:
    """Process a single motif with permutation test.

    Args:
        args: Tuple of (motif_data, foreground, background, fg_gc, bg_gc, n_perm).

    Returns:
        Dictionary with motif info and test results.
    """
    motif_data, foreground, background, fg_gc, bg_gc, n_perm = args


    result = permutation_test(
        foreground_scores=foreground,
        background_scores=background,
        n_permutations=n_perm
    )


    return {
        'motif_id': motif_data['motif_id'],
        'tf_name': motif_data['tf_name'],
        'tf_class': motif_data['tf_class'],
        'tf_family': motif_data['tf_family'],
        'auc_roc': result.auc_roc,
        'auc_prc': result.auc_prc,
        'p_value_roc': result.p_value_roc,
        'p_value_prc': result.p_value_prc,
    }


def parallel_permutation_runner(
    motifs: List[Dict],
    foreground_scores: np.ndarray,
    background_scores: np.ndarray,
    foreground_gc: np.ndarray,
    background_gc: np.ndarray,
    n_permutations: int = 10000,
    n_workers: int = 4
) -> List[Dict]:
    """Run permutation tests in parallel for multiple motifs.

    Args:
        motifs: List of motif dictionaries.
        foreground_scores: Score matrix [n_foreground, n_motifs].
        background_scores: Score matrix [n_background, n_motifs].
        foreground_gc: GC content for foreground.
        background_gc: GC content for background.
        n_permutations: Number of permutations per motif.
        n_workers: Number of parallel workers.

    Returns:
        List of result dictionaries.
    """
    if not motifs:
        return []

    n_workers = min(n_workers, os.cpu_count() or 4)

    # Prepare arguments for each motif
    args_list = [
        (
            motif,
            foreground_scores[:, idx],
            background_scores[:, idx],
            foreground_gc,
            background_gc,
            n_permutations
        )
        for idx, motif in enumerate(motifs)
    ]

    results = []
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = [executor.submit(process_single_motif, arg) for arg in args_list]

        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing motifs"):
            result = future.result()
            if result is not None:
                results.append(result)

    return results


# ============================================================================
# Gene ID utilities
# ============================================================================

def get_deg_gene_ids(
    df: pd.DataFrame,
    condition: str = 'all',
    padj_threshold: float = 0.05,
    log2fc_threshold: float = 1.0
) -> np.ndarray:
    """Extract differentially expressed gene IDs.

    Args:
        df: DESeq2 results DataFrame.
        condition: 'all', 'up', or 'down'.
        padj_threshold: Adjusted p-value cutoff.
        log2fc_threshold: Log2 fold-change cutoff (absolute value).

    Returns:
        Array of unique gene IDs.
    """
    log2fc_threshold = abs(log2fc_threshold)

    if condition == 'all':
        df = df[
            (df['log2FoldChange'] <= -log2fc_threshold) |
            (df['log2FoldChange'] >= log2fc_threshold)
        ]
    elif condition == 'down':
        df = df[df['log2FoldChange'] <= -log2fc_threshold]
    elif condition == 'up':
        df = df[df['log2FoldChange'] >= log2fc_threshold]

    df = df[df['padj'] <= padj_threshold]
    gene_ids = np.array([i for i in df['id'] if isinstance(i, str)])
    return np.unique(gene_ids)


def get_background_gene_ids(
    df: pd.DataFrame,
    padj_threshold: float = 0.05,
    log2fc_threshold: float = np.log2(5/4)
) -> np.ndarray:
    """Extract non-DEG gene IDs for background.

    Args:
        df: DESeq2 results DataFrame.
        padj_threshold: Adjusted p-value cutoff.
        log2fc_threshold: Log2 fold-change cutoff.

    Returns:
        Array of background gene IDs.
    """
    log2fc_threshold = abs(log2fc_threshold)
    df = df[
        (df['log2FoldChange'] >= -log2fc_threshold) &
        (df['log2FoldChange'] <= log2fc_threshold) &
        (df['padj'] > padj_threshold)
    ]
    gene_ids = np.array([i for i in df['id'] if isinstance(i, str)])
    return np.unique(gene_ids)


def read_gene_set(path: str) -> np.ndarray:
    """Read gene IDs from file.

    Args:
        path: Path to file with gene IDs (one per line).

    Returns:
        Array of unique gene IDs.
    """
    with open(path) as file:
        gene_ids = [line.strip() for line in file]
    return np.array(list(set(gene_ids)))


def get_indexes(all_ids: np.ndarray, sub_ids: np.ndarray) -> np.ndarray:
    """Get indices of sub_ids in all_ids.

    Args:
        all_ids: Array of all IDs.
        sub_ids: Array of subset IDs.

    Returns:
        Indices where sub_ids appear in all_ids.
    """
    _, indices, _ = np.intersect1d(all_ids, sub_ids, assume_unique=False, return_indices=True)
    return indices


def get_motif_to_cluster(cluster_path: str) -> Dict[str, str]:
    """Load motif-to-cluster mapping.

    Args:
        cluster_path: Path to cluster TSV file.

    Returns:
        Dictionary mapping motif_id to cluster_id.
    """
    clusters = pd.read_csv(cluster_path, sep='\t')
    motif_to_cluster = {}
    for _, row in clusters.iterrows():
        cluster_id = row.iloc[0]
        for motif_id in row['motif_ids'].split(','):
            motif_to_cluster[motif_id] = cluster_id
    return motif_to_cluster


# ============================================================================
# Main pipeline
# ============================================================================

def esdeg(
    motif_db: str,
    taxon: str,
    gc_threshold: float,
    path_to_promoters: str,
    path_to_data: str,
    nproc: int,
    type_of_data: str = 'deg',
    condition: str = 'all',
    log2fc_thr_deg: float = 1.0,
    log2fc_thr_background: float = np.log2(5/4),
    padj_thr: float = 0.05,
    n_permutations: int = 1000,
    match_ratio: int = 5
) -> pd.DataFrame:
    """Main enrichment analysis pipeline with stratified permutation testing.

    Args:
        motif_db: Database name ('jaspar' or 'hocomoco').
        taxon: Organism taxon.
        gc_threshold: GC matching threshold (deprecated, kept for compatibility).
        path_to_promoters: Path to promoter FASTA file.
        path_to_data: Path to DEG table or gene set file.
        nproc: Number of parallel processes.
        type_of_data: 'deg' or 'set'.
        condition: DEG condition ('all', 'up', 'down').
        log2fc_thr_deg: Log2FC threshold for DEGs.
        log2fc_thr_background: Log2FC threshold for background.
        padj_thr: Adjusted p-value threshold.
        n_permutations: Number of permutations for testing.
        match_ratio: Ratio of background to foreground samples.

    Returns:
        DataFrame with enrichment results.
    """
    torch.set_num_threads(2)
    #device = 'cuda' if torch.cuda.is_available() else 'cpu'
    device = 'cpu'

    # Load and filter promoters
    print('-' * 30)
    print('Reading promoters...')
    promoters, ids = promoters_parser(path_to_promoters)

    # Filter by length
    length_counts = collections.Counter(len(seq) for seq in promoters)
    if len(length_counts) > 1:
        print('#' * 10 + ' WARNING! ' + '#' * 10)
        print('FASTA file contains sequences with different lengths')
        for length, count in length_counts.items():
            print(f'  {count} sequences with length {length}')

        target_length, _ = length_counts.most_common(1)[0]
        print(f'Using only sequences with most common length: {target_length}')

        valid_mask = [len(seq) == target_length for seq in promoters]
        promoters = [seq for seq, valid in zip(promoters, valid_mask) if valid]
        ids = np.array([id_ for id_, valid in zip(ids, valid_mask) if valid])
        print('#' * 10 + ' DATA FILTERED! ' + '#' * 10)

    gc_content = calculate_gc(promoters)
    promoters = seq_to_int(promoters, device=device)
    print('-' * 30)

    # Load gene IDs
    if type_of_data == 'set':
        print('Reading gene set...')
        foreground_ids = read_gene_set(path_to_data)
        check = np.sum(np.isin(foreground_ids, ids))
        if check == 0:
            print('ERROR: No common IDs found. Check ID format.')
            sys.exit(1)
        background_ids = np.setdiff1d(ids, foreground_ids)

    elif type_of_data == 'deg':
        print('Reading DEG table...')
        deg_table = pd.read_csv(path_to_data, sep='\t', comment='#')
        deg_table = deg_table[deg_table['padj'] <= 1]

        foreground_ids = get_deg_gene_ids(
            deg_table, condition, padj_thr, log2fc_thr_deg
        )

        if len(foreground_ids) <= 5:
            print(f'WARNING: Number of DEGs ({len(foreground_ids)}) is very low!')

        check = np.sum(np.isin(foreground_ids, ids))
        if check == 0:
            print('ERROR: No common IDs found. Check ID format.')
            sys.exit(1)

        background_ids = get_background_gene_ids(
            deg_table, padj_thr, log2fc_thr_background
        )

    print('-' * 30)

    # Load motifs and scan promoters
    print('Loading motifs...')
    motifs = read_motifs_from_db(motif_db, taxon)

    print('Scanning promoters...')
    scores = process_all_models(promoters, motifs, chunk_size=64)
    print('-' * 30)

    # Prepare data
    foreground_idx = get_indexes(ids, foreground_ids)
    background_pool_idx = get_indexes(ids, background_ids)

    foreground_scores = scores[foreground_idx, :].cpu().numpy()
    background_pool_scores = scores[background_pool_idx, :].cpu().numpy()
    foreground_gc = gc_content[foreground_idx]
    background_pool_gc = gc_content[background_pool_idx]

    # Select GC-matched background
    print('Selecting GC-matched background...')
    selected_bg_idx = select_gc_matched_background(
        foreground_gc=foreground_gc,
        background_gc=background_pool_gc,
        match_ratio=match_ratio,
        n_quantiles=10
    )

    background_scores = background_pool_scores[selected_bg_idx, :]
    background_gc = background_pool_gc[selected_bg_idx]

    print(f'Foreground: {len(foreground_scores)} samples')
    print(f'Background: {len(background_scores)} samples (matched from {len(background_pool_scores)})')
    print(f'GC - Foreground: {foreground_gc.mean():.3f}, Background: {background_gc.mean():.3f}')
    print('-' * 30)

    # Run permutation tests
    print('Running stratified permutation tests...')
    results = parallel_permutation_runner(
        motifs=motifs,
        foreground_scores=foreground_scores,
        background_scores=background_scores,
        foreground_gc=foreground_gc,
        background_gc=background_gc,
        n_permutations=n_permutations,
        n_workers=nproc
    )
    print('-' * 30)

    # Format results
    df = pd.DataFrame(results)
    df['p_value_roc_adj'] = st.false_discovery_control(df['p_value_roc'])
    df['p_value_prc_adj'] = st.false_discovery_control(df['p_value_prc'])

    # Add cluster information for JASPAR
    if motif_db == 'jaspar':
        cluster_path = resources.files('esdeg').joinpath(f'clusters/{taxon}.tsv')
        motif_to_cluster = get_motif_to_cluster(cluster_path)
        df['jaspar_cluster'] = df['motif_id'].map(motif_to_cluster)

        # Handle dimers
        df['tf_name'] = df['tf_name'].str.split('::')
        df['tf_class'] = df['tf_class'].str.split('::')
        df['tf_family'] = df['tf_family'].str.split('::')

        # Fix length mismatches
        for col in ['tf_name', 'tf_class', 'tf_family']:
            max_len = df[col].apply(len).max()
            df[col] = df[col].apply(lambda x: x * max_len if len(x) < max_len else x)

        df = df.explode(column=['tf_name', 'tf_class', 'tf_family'], ignore_index=True)

    df = df.sort_values(by='auc_roc', ascending=False)

    return df
