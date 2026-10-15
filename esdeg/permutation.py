"""auPRC and auROC enrichment with permutation p-values calibrated to the matched design."""

from __future__ import annotations

import numba as nb
import numpy as np
from scipy.special import ndtr
from scipy.stats import rankdata

PAIRED_PERMUTATIONS = 2000
CF_LIMIT = 12.0


@nb.njit(cache=True)
def _statistics(labels_in_rank_order, inverse, out):
    """Write U (foreground-above-background pairs) and the AP sum (= N_foreground * AP)."""
    hits = 0
    pairs = 0
    total = 0.0
    for i in range(labels_in_rank_order.shape[0]):
        value = np.int64(labels_in_rank_order[i])
        hits += value
        pairs += (1 - value) * hits
        total += value * hits * inverse[i]
    out[0] = pairs
    out[1] = total


@nb.njit(parallel=True, cache=True)
def _observed(order, labels, inverse, out):
    n_models, n = order.shape
    for j in nb.prange(n_models):
        ranked = np.empty(n, dtype=np.uint8)
        for i in range(n):
            ranked[i] = labels[order[j, i]]
        _statistics(ranked, inverse, out[j])


@nb.njit(cache=True)
def _random_statistics(n, n_foreground, seed, inverse, out):
    """U and AP sum of one uniformly random arrangement, by selection sampling in one pass.

    Position i is foreground with probability (foreground still needed) / (positions left), which
    yields exactly n_foreground foreground positions. Randomness is a splitmix64 stream seeded by
    the permutation index, so results do not depend on the number of threads.
    """
    state = np.uint64(seed)
    need = n_foreground
    hits = 0
    pairs = 0
    total = 0.0
    for i in range(n):
        state += np.uint64(0x9E3779B97F4A7C15)
        z = (state ^ (state >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
        z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
        z = z ^ (z >> np.uint64(31))
        if ((z >> np.uint64(32)) * np.uint64(n - i)) >> np.uint64(32) < np.uint64(need):
            need -= 1
            hits += 1
            total += hits * inverse[i]
        else:
            pairs += hits
    out[0] = pairs
    out[1] = total


@nb.njit(parallel=True, cache=True)
def _shared_null(n, n_foreground, n_permutations, seed, inverse, out):
    for b in nb.prange(n_permutations):
        _random_statistics(n, n_foreground, seed + b, inverse, out[b])


@nb.njit(parallel=True, cache=True)
def _paired_labels(start, items, seed, labels):
    """Each matched group (a foreground promoter and its backgrounds) gets one random foreground."""
    for b in nb.prange(labels.shape[0]):
        np.random.seed(seed + b)
        labels[b, :] = 0
        for g in range(start.shape[0] - 1):
            labels[b, items[start[g] + np.random.randint(0, start[g + 1] - start[g])]] = 1


@nb.njit(parallel=True, cache=True)
def _paired_moments(order, labels, inverse, mean, std):
    n_models, n = order.shape
    n_permutations = labels.shape[0]
    for j in nb.prange(n_models):
        total = np.zeros(2)
        total_sq = np.zeros(2)
        ranked = np.empty(n, dtype=np.uint8)
        value = np.empty(2)
        for b in range(n_permutations):
            for i in range(n):
                ranked[i] = labels[b, order[j, i]]
            _statistics(ranked, inverse, value)
            total += value
            total_sq += value * value
        for k in range(2):
            mean[j, k] = total[k] / n_permutations
            variance = (total_sq[k] - total[k] * total[k] / n_permutations) / (n_permutations - 1)
            std[j, k] = np.sqrt(max(variance, 0.0))


def _upper_tail(sample, values):
    """P(X >= value) from a Cornish-Fisher fit of the sample's mean, SD, skewness and kurtosis.

    The expansion w(z) = z + g1 (z^2 - 1) / 6 + g2 (z^3 - 3 z) / 24 - g1^2 (2 z^3 - 5 z) / 36 maps a
    standard normal z to the standardized statistic; it is inverted on a grid over the monotone run
    that contains z = 0. Beyond that run z continues with unit slope in w.
    """
    mean, sd = sample.mean(), sample.std(ddof=1)
    standardized = (sample - mean) / sd
    g1 = np.mean(standardized**3)
    g2 = np.mean(standardized**4) - 3.0
    z = np.linspace(-CF_LIMIT, CF_LIMIT, 4801)
    w = z + g1 * (z**2 - 1) / 6 + g2 * (z**3 - 3 * z) / 24 - g1**2 * (2 * z**3 - 5 * z) / 36
    rising = np.diff(w) > 0
    centre = z.size // 2
    up, down = rising[centre:], rising[:centre][::-1]
    top = centre + (up.size if up.all() else int(np.argmin(up)))
    bottom = centre - (down.size if down.all() else int(np.argmin(down)))
    z, w = z[bottom : top + 1], w[bottom : top + 1]
    scaled = (np.asarray(values, dtype=np.float64) - mean) / sd
    z_values = np.interp(scaled, w, z)
    z_values = np.where(scaled > w[-1], z[-1] + (scaled - w[-1]), z_values)
    z_values = np.where(scaled < w[0], z[0] + (scaled - w[0]), z_values)
    return ndtr(-z_values)


def _validate(foreground_scores, background_scores, matched_to):
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
    matched_to = np.asarray(matched_to, dtype=np.int64)
    if matched_to.shape != (background_scores.shape[0],):
        raise ValueError("matched_to must hold one foreground index per background item.")
    if matched_to.min() < 0 or matched_to.max() >= foreground_scores.shape[0]:
        raise ValueError("matched_to must index foreground rows.")
    return foreground_scores, background_scores, matched_to


def permutation_enrichment(
    foreground_scores: np.ndarray,
    background_scores: np.ndarray,
    matched_to: np.ndarray,
    n_permutations: int = 100_000,
    seed: int = 0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return auPRC, auROC and their one-sided p-values for each score column.

    auPRC is the average precision of the pooled foreground and background ranked by score; its
    chance level is N_foreground / N_total. auROC is U / (N_foreground * N_background) with
    midranks for ties.

    Without tied scores the labels of a column, read in descending score order, are a uniformly
    random arrangement under H0. That null depends only on the group sizes, so one sample of
    n_permutations arrangements serves every column; its upper tail is taken from a Cornish-Fisher
    fit to the sample moments, which reaches far below 1 / n_permutations. The background was chosen
    by matching, so labels are exchangeable only inside a matched group (a foreground promoter and
    the backgrounds assigned to it). The true null of a column has the mean and standard deviation
    of PAIRED_PERMUTATIONS permutations that choose one foreground per group. The observed
    statistic is mapped from that moment scale to the scale of the shared sample before the tail is
    read.

    Ties in a column are broken by one random row order, so the result is reproducible by seed. A
    column with all scores equal has no information and gets p = 1.
    """
    foreground_scores, background_scores, matched_to = _validate(
        foreground_scores, background_scores, matched_to
    )
    if n_permutations < 1:
        raise ValueError("n_permutations must be positive.")
    n_foreground, n_background = foreground_scores.shape[0], background_scores.shape[0]
    n = n_foreground + n_background
    scores = np.vstack((foreground_scores, background_scores))
    n_models = scores.shape[1]
    rng = np.random.default_rng(seed)
    row_order = rng.permutation(n)
    shared_seed, paired_seed = (int(v) for v in rng.integers(0, 2**31, size=2))

    order = np.ascontiguousarray(
        np.argsort(-scores[row_order], axis=0, kind="stable").T, dtype=np.int32
    )
    observed_labels = (row_order < n_foreground).astype(np.uint8)
    inverse = 1.0 / np.arange(1, n + 1)
    observed = np.empty((n_models, 2))
    _observed(order, observed_labels, inverse, observed)

    group_of_item = np.concatenate((np.arange(n_foreground), matched_to))[row_order]
    items = np.argsort(group_of_item, kind="stable").astype(np.int64)
    start = np.concatenate(([0], np.cumsum(np.bincount(group_of_item, minlength=n_foreground))))
    paired = np.empty((PAIRED_PERMUTATIONS, n), dtype=np.uint8)
    _paired_labels(start, items, paired_seed, paired)
    paired_mean = np.empty((n_models, 2))
    paired_std = np.empty((n_models, 2))
    _paired_moments(order, paired, inverse, paired_mean, paired_std)

    null = np.empty((n_permutations, 2))
    _shared_null(n, n_foreground, n_permutations, shared_seed, inverse, null)

    pvalues = np.empty((n_models, 2))
    for k in range(2):
        sample = null[:, k]
        usable = paired_std[:, k] > 0
        scale = np.where(usable, sample.std(ddof=1) / np.where(usable, paired_std[:, k], 1.0), 1.0)
        mapped = np.where(
            usable, sample.mean() + (observed[:, k] - paired_mean[:, k]) * scale, observed[:, k]
        )
        pvalues[:, k] = _upper_tail(sample, mapped)

    auprc = observed[:, 1] / n_foreground
    ranks = rankdata(scores, axis=0)
    auc = (ranks[:n_foreground].sum(axis=0) - n_foreground * (n_foreground + 1) / 2) / (
        n_foreground * n_background
    )
    constant = np.ptp(scores, axis=0) == 0
    auprc = np.where(constant, n_foreground / n, auprc)
    p_pr = np.where(constant, 1.0, pvalues[:, 1])
    p_roc = np.where(constant, 1.0, pvalues[:, 0])
    return auprc, auc, p_pr, p_roc
