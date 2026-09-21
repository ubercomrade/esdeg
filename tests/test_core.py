import numpy as np
import pandas as pd

from esdeg.annotation import _expand_dimers
from esdeg.functions import compute_aucs, permutation_test, select_gc_matched_background
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
