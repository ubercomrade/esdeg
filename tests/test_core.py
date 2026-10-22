import mimosa
import numpy as np
import pandas as pd
import pytest

from esdeg.annotation import _expand_dimers
from esdeg.functions import (
    calculate_cpg_log_oe,
    calculate_gc,
    esdeg,
    process_all_models,
    select_matched_background,
)
from esdeg.parsers import read_table
from esdeg.permutation import _shared_null, _upper_tail, permutation_enrichment


def test_matching_is_unique_reproducible_and_close():
    rng = np.random.default_rng(0)
    pool = rng.normal(size=(500, 2))
    foreground = rng.normal(loc=1.0, scale=0.3, size=(20, 2))
    first, first_owner = select_matched_background(foreground, pool, 5, random_state=7)
    second, second_owner = select_matched_background(foreground, pool, 5, random_state=7)
    assert np.array_equal(first, second) and np.array_equal(first_owner, second_owner)
    assert len(first) == len(np.unique(first))
    assert np.array_equal(np.bincount(first_owner, minlength=20), np.full(20, len(first) // 20))
    assert len(first) % 20 == 0 and 20 <= len(first) <= 100
    assert np.all(pool[first].mean(axis=0) > pool.mean(axis=0))
    partial, partial_owner = select_matched_background(foreground, pool[:10], 1, random_state=7)
    assert len(partial) == len(partial_owner) == 10


def test_matching_lowers_ratio_when_background_is_unbalanced():
    rng = np.random.default_rng(0)
    pool = np.vstack((rng.normal(size=(400, 2)), rng.normal(loc=3.0, size=(40, 2))))
    foreground = rng.normal(loc=3.0, size=(20, 2))
    result, _ = select_matched_background(foreground, pool, 5, random_state=1)
    assert len(result) % 20 == 0 and len(result) < 100
    assert np.abs(pool[result].mean(axis=0) - foreground.mean(axis=0)).max() < 0.5


def test_cpg_log_oe_does_not_cross_sequence_boundaries():
    rows = ["ACGT", "CCGG", "CATG", "GCAC"]
    sequences = mimosa.EncodedSequences.from_rows(mimosa.encode_sequence(row) for row in rows)
    expected = [
        np.log2((n_cg + 1) / (row.count("C") * row.count("G") / len(row) + 1))
        for row, n_cg in zip(rows, [1, 1, 0, 0], strict=True)
    ]
    np.testing.assert_allclose(calculate_cpg_log_oe(sequences), expected)


def _reference_ap(foreground, background):
    """Average precision of one column by the textbook definitions (no ties)."""
    scores = np.concatenate((foreground, background))
    labels = np.concatenate((np.ones(len(foreground)), np.zeros(len(background))))
    ranked = labels[np.argsort(-scores)]
    hits = np.cumsum(ranked)
    ap = np.sum(ranked * hits / np.arange(1, len(ranked) + 1)) / len(foreground)
    return ap


def test_ap_enrichment_matches_definitions_and_is_reproducible():
    rng = np.random.default_rng(3)
    foreground = rng.normal(size=(40, 3)) + np.array([0.0, 0.8, -0.8])
    background = rng.normal(size=(80, 3))
    matched_to = np.tile(np.arange(40), 2)
    args = (foreground, background, matched_to)
    first = permutation_enrichment(*args, n_permutations=20_000, seed=5)
    second = permutation_enrichment(*args, n_permutations=20_000, seed=5)
    for a, b in zip(first, second, strict=True):
        np.testing.assert_array_equal(a, b)
    auprc, pvalue = first
    for column in range(3):
        ap = _reference_ap(foreground[:, column], background[:, column])
        assert auprc[column] == pytest.approx(ap)
    assert pvalue[1] < 0.01
    assert pvalue[2] > 0.5
    assert pvalue.min() >= 1 / 20_001


def test_ap_enrichment_handles_perfect_constant_and_tied_columns():
    foreground = np.column_stack((np.arange(30, 40.0), np.ones(10), np.repeat([1.0, 2.0], 5)))
    background = np.column_stack((np.arange(0, 20.0), np.ones(20), np.tile([1.0, 2.0], 10)))
    matched_to = np.tile(np.arange(10), 2)
    auprc, pvalue = permutation_enrichment(foreground, background, matched_to, 5_000, seed=1)
    assert auprc[0] == 1.0
    assert pvalue[0] < 1e-3
    assert auprc[1] == pytest.approx(10 / 30)
    assert pvalue[1] == 1.0
    assert 0.05 < pvalue[2] <= 1.0
    with pytest.raises(ValueError, match="matched_to"):
        permutation_enrichment(foreground, background, np.zeros(19, dtype=int), 100)


def test_paired_calibration_keeps_the_false_positive_rate_under_matched_null():
    rng = np.random.default_rng(11)
    n_groups, n_columns = 40, 300
    group_effect = rng.normal(scale=3.0, size=(n_groups, n_columns))
    first = group_effect + rng.normal(size=group_effect.shape)
    second = group_effect + rng.normal(size=group_effect.shape)
    swap = rng.random((n_groups, n_columns)) < 0.5  # which member of a pair is the foreground
    foreground = np.where(swap, first, second)
    background = np.where(swap, second, first)
    _, pvalue = permutation_enrichment(foreground, background, np.arange(n_groups), 20_000, seed=2)
    assert 0.02 < np.mean(pvalue < 0.05) < 0.09
    assert 0.4 < pvalue.mean() < 0.6


def test_enrichment_scans_only_selected_promoters_and_checks_all_lengths(monkeypatch):
    rows = ["AAAAAAAC", "ACACACAC", "CCCCCCCC", "GGGGGGGG", "TTTTTTTT", "AACCGGTT", "TGTGTGTG"]
    sequences = mimosa.EncodedSequences.from_rows(mimosa.encode_sequence(row) for row in rows)
    ids = np.asarray(list("abcdefg"), dtype=object)
    pfm = np.array([[0.7] * 6, [0.1] * 6, [0.1] * 6, [0.1] * 6], dtype=np.float32)
    motif = {"motif_id": "M1", "tf_name": "M1", "tf_class": "NA", "tf_family": "NA"}
    motif["model"] = mimosa.pwm_from_pfm(pfm, background=0.25, name="M1")
    monkeypatch.setattr("esdeg.functions.promoters_parser", lambda path: (sequences, ids))
    monkeypatch.setattr("esdeg.functions.read_motifs_from_db", lambda db, taxon: [motif])
    monkeypatch.setattr("esdeg.functions.read_gene_set", lambda path: np.asarray(["a", "b"]))
    deg = pd.DataFrame(
        {
            "id": list("abcde"),
            "log2FoldChange": [2, 2, 0, 0, 0],
            "padj": [0.01, 0.01, 0.5, 0.5, 0.5],
        }
    )
    monkeypatch.setattr("esdeg.functions.read_table", lambda path: deg)
    scanned = []

    def capture(batch, motifs):
        scanned.append(batch)
        return process_all_models(batch, motifs)

    monkeypatch.setattr("esdeg.functions.process_all_models", capture)
    full_scores = process_all_models(sequences, [motif])
    for mode, ratio, pool in (("deg", 5, np.arange(2, 5)), ("set", 1, np.arange(2, 7))):
        result = esdeg(
            path_to_promoters="unused",
            path_to_data="unused",
            type_of_data=mode,
            match_ratio=ratio,
            seed=9,
        )
        features = np.column_stack((calculate_gc(sequences), calculate_cpg_log_oe(sequences)))
        selected, _ = select_matched_background(features[:2], features[pool], ratio, 9)
        expected_idx = np.concatenate((np.arange(2), pool[selected]))
        np.testing.assert_array_equal(
            process_all_models(scanned[-1], [motif]), full_scores[expected_idx]
        )
        assert len(scanned[-1]) == len(expected_idx)
        assert result["motif_id"].tolist() == ["M1"]

    short = mimosa.EncodedSequences.from_rows(
        mimosa.encode_sequence(row) for row in rows[:-1] + ["AAA"]
    )
    monkeypatch.setattr("esdeg.functions.promoters_parser", lambda path: (short, ids))
    with pytest.raises(ValueError, match="shorter than the largest motif window"):
        esdeg(path_to_promoters="unused", path_to_data="unused")


def test_dimer_metadata_broadcasts_singletons():
    table = pd.DataFrame(
        [
            {
                "motif_id": "M1",
                "tf_name": "TF1::TF2",
                "tf_class": "class",
                "tf_family": "family1::family2",
            }
        ]
    )
    result = _expand_dimers(table)
    assert result["tf_name"].tolist() == ["TF1", "TF2"]
    assert result["tf_class"].tolist() == ["class", "class"]


def test_common_reader_accepts_csv_and_tsv(tmp_path):
    frame = pd.DataFrame({"id": ["a"], "value": [1]})
    csv_path = tmp_path / "data.csv"
    tsv_path = tmp_path / "data.tsv"
    frame.to_csv(csv_path, index=False)
    frame.to_csv(tsv_path, sep="\t", index=False)
    assert read_table(csv_path).equals(frame)
    assert read_table(tsv_path).equals(frame)


def test_cornish_fisher_tail_matches_a_large_permutation_sample():
    n_foreground, n_background = 60, 180
    n = n_foreground + n_background
    inverse = 1.0 / np.arange(1, n + 1)
    large = np.empty((2_000_000, 2))
    _shared_null(n, n_foreground, len(large), 1, inverse, large)
    for column in range(2):
        truth = np.sort(large[:, column])
        small = large[:100_000, column]
        for probability in (1e-1, 1e-2, 1e-3, 1e-4):
            value = truth[len(truth) - round(probability * len(truth))]
            exact = (len(truth) - np.searchsorted(truth, value)) / len(truth)
            assert _upper_tail(small, [value])[0] == pytest.approx(exact, rel=0.5)
