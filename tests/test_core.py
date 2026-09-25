import mimosa
import numpy as np
import pandas as pd
import pytest
from numba import njit

from esdeg.annotation import _expand_dimers
from esdeg.functions import (
    compute_aucs,
    esdeg,
    permutation_loop,
    permutation_test,
    process_all_models,
    select_gc_matched_background,
)
from esdeg.parsers import read_table


def test_gc_matching_is_unique_and_reproducible():
    foreground = np.array([0.1, 0.1, 0.9])
    background = np.array([0.1, 0.2, 0.8, 0.9, 0.5])
    first = select_gc_matched_background(foreground, background, 1, random_state=7)
    second = select_gc_matched_background(foreground, background, 1, random_state=7)
    assert np.array_equal(first, second)
    assert len(first) == len(np.unique(first)) == 3


def test_auc_and_permutation_are_one_sided_and_reproducible():
    labels = np.array([1.0, 1.0, 0.0, 0.0])
    assert compute_aucs(labels, np.array([4.0, 3.0, 2.0, 1.0]))[0] == 1.0
    first = permutation_test(np.array([4.0, 3.0]), np.array([2.0, 1.0]), 10, seed=3)
    second = permutation_test(np.array([4.0, 3.0]), np.array([2.0, 1.0]), 10, seed=3)
    assert first == second
    assert first.p_value_roc >= 1 / 11


@njit
def _reference_permutations(labels, scores, count, seed):
    np.random.seed(seed)
    shuffled = labels.copy()
    result = np.empty((count, 2), dtype=np.float32)
    for permutation in range(count):
        for index in range(len(labels) - 1, 0, -1):
            swap = np.random.randint(0, index + 1)
            shuffled[index], shuffled[swap] = shuffled[swap], shuffled[index]
        result[permutation] = compute_aucs(shuffled, scores)
    return result


def test_permutation_loop_matches_resorting_reference():
    cases = [
        ([4, 3], [2, 1]),
        ([3, 3, 1], [3, 2, 1, 1]),
        ([1, 2, 3], [4, 5]),
        ([1, 1], [1, 1]),
    ]
    for foreground, background in cases:
        scores = np.asarray(foreground + background, dtype=np.float64)
        labels = np.asarray([1] * len(foreground) + [0] * len(background), dtype=np.float64)
        observed = compute_aucs(labels, scores)
        for seed, count in ((0, 1), (3, 101), (19, 101)):
            expected = _reference_permutations(labels, scores, count, seed)
            actual = np.column_stack(permutation_loop(labels, scores, count, seed))
            np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-7)
            np.testing.assert_array_equal(
                np.count_nonzero(actual >= observed, axis=0),
                np.count_nonzero(expected >= observed, axis=0),
            )


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
            n_permutations=5,
            nproc=1,
            seed=9,
        )
        gc = np.array([sum(base in "CG" for base in row) / len(row) for row in rows])
        gc_seed = int(np.random.SeedSequence(9).spawn(2)[0].generate_state(1)[0])
        selected = select_gc_matched_background(gc[:2], gc[pool], ratio, gc_seed)
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
        esdeg(path_to_promoters="unused", path_to_data="unused", n_permutations=1)


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
