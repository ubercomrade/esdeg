import mimosa
import numpy as np
import pandas as pd
import pytest
from scipy.stats import mannwhitneyu

from esdeg.annotation import _expand_dimers
from esdeg.functions import (
    calculate_cpg_log_oe,
    calculate_gc,
    esdeg,
    mann_whitney_enrichment,
    process_all_models,
    select_matched_background,
)
from esdeg.parsers import read_table


def test_matching_is_unique_reproducible_and_close():
    rng = np.random.default_rng(0)
    pool = rng.normal(size=(500, 2))
    foreground = rng.normal(loc=1.0, scale=0.3, size=(20, 2))
    first = select_matched_background(foreground, pool, 5, random_state=7)
    second = select_matched_background(foreground, pool, 5, random_state=7)
    assert np.array_equal(first, second)
    assert len(first) == len(np.unique(first))
    assert len(first) % 20 == 0 and 20 <= len(first) <= 100
    assert np.all(pool[first].mean(axis=0) > pool.mean(axis=0))
    assert len(select_matched_background(foreground, pool[:10], 1, random_state=7)) == 10


def test_matching_lowers_ratio_when_background_is_unbalanced():
    rng = np.random.default_rng(0)
    pool = np.vstack((rng.normal(size=(400, 2)), rng.normal(loc=3.0, size=(40, 2))))
    foreground = rng.normal(loc=3.0, size=(20, 2))
    result = select_matched_background(foreground, pool, 5, random_state=1)
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


def test_mann_whitney_matches_scipy_with_ties_and_is_one_sided():
    rng = np.random.default_rng(1)
    foreground = np.round(rng.normal(0.4, 1.0, size=(40, 3)), 1)
    background = np.round(rng.normal(0.0, 1.0, size=(200, 3)), 1)
    auc, pvalue = mann_whitney_enrichment(foreground, background)
    for column in range(3):
        expected = mannwhitneyu(
            foreground[:, column], background[:, column], alternative="greater"
        )
        assert auc[column] == pytest.approx(expected.statistic / (40 * 200))
        assert pvalue[column] == pytest.approx(expected.pvalue)
    _, depleted_p = mann_whitney_enrichment(background[:40, :1], foreground[:, :1])
    assert depleted_p[0] > 0.5
    perfect_auc, perfect_p = mann_whitney_enrichment(
        np.array([[4.0], [3.0]]), np.array([[2.0], [1.0]])
    )
    assert perfect_auc[0] == 1.0 and perfect_p[0] < 0.5
    constant_auc, constant_p = mann_whitney_enrichment(np.ones((3, 1)), np.ones((6, 1)))
    assert constant_auc[0] == 0.5 and constant_p[0] == 1.0


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
        selected = select_matched_background(features[:2], features[pool], ratio, 9)
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
