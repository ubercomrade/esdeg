"""Input and motif database readers."""

from __future__ import annotations

import json
import logging
from importlib import resources
from pathlib import Path

import mimosa
import numpy as np
import pandas as pd
from pyjaspar import jaspardb

logger = logging.getLogger(__name__)
JASPAR_RELEASE = "JASPAR2024"
DNA_ALPHABET = "ACGT"


def promoters_parser(path: str | Path):
    """Read FASTA through Mimosa and normalize IDs before the first colon."""
    sequences, names = mimosa.read_fasta(path)
    if len(sequences) == 0:
        raise ValueError(f"FASTA file {path!s} is empty.")

    ids = [str(name).split(":", 1)[0].strip() for name in names]
    if any(not promoter_id for promoter_id in ids):
        raise ValueError("FASTA contains an empty promoter ID.")
    if len(ids) != len(set(ids)):
        raise ValueError("FASTA contains duplicate promoter IDs.")
    lengths = np.diff(sequences.offsets)
    if np.any(lengths == 0):
        raise ValueError("FASTA contains an empty sequence.")
    return sequences, np.asarray(ids, dtype=object)


def read_table(path: str | Path) -> pd.DataFrame:
    """Read either CSV or TSV input using the same parser everywhere."""
    return pd.read_csv(path, sep=None, engine="python", comment="#")


def read_gene_set(path: str | Path) -> np.ndarray:
    """Read a deterministic, non-empty set of IDs."""
    with open(path, encoding="utf-8") as handle:
        ids = list(dict.fromkeys(line.strip() for line in handle if line.strip()))
    if not ids:
        raise ValueError(f"Gene set {path!s} is empty.")
    return np.asarray(ids, dtype=object)


def _annotation(value) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float) and np.isnan(value):
        return "NA"
    if isinstance(value, (list, tuple, np.ndarray)):
        values = [str(item).strip() for item in value if str(item).strip()]
        return "NA" if not values else "::".join(values)
    text = str(value).strip()
    return text or "NA"


def _record(model, motif_id: str, tf_name, tf_class, tf_family) -> dict:
    return {
        "model": model,
        "motif_id": motif_id,
        "tf_name": _annotation(tf_name).upper(),
        "tf_class": _annotation(tf_class),
        "tf_family": _annotation(tf_family),
    }


def _hocomoco_records() -> list[dict]:
    path = resources.files("esdeg").joinpath("hocomoco/H12CORE_annotation.jsonl")
    records = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            try:
                item = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f"Invalid HOCOMOCO JSONL at line {line_number}.") from exc

            required = {"name", "pcm", "masterlist_info"}
            missing = sorted(required - item.keys())
            if missing:
                raise ValueError(
                    f"HOCOMOCO motif at line {line_number} misses: {', '.join(missing)}."
                )
            raw_pcm = np.asarray(item["pcm"], dtype=np.float32)
            if raw_pcm.ndim != 2:
                raise ValueError(f"HOCOMOCO motif {item['name']} PCM must be 2-D.")
            if raw_pcm.shape[1] == 4:
                pcm = raw_pcm.T
            elif raw_pcm.shape[0] == 4:
                pcm = raw_pcm
            else:
                raise ValueError(
                    f"HOCOMOCO motif {item['name']} PCM must have shape (length, 4) or (4, length)."
                )
            if pcm.shape[0] != 4 or pcm.shape[1] < 1:
                raise ValueError(f"HOCOMOCO motif {item['name']} PCM has invalid shape.")
            if int(item.get("length", pcm.shape[1])) != pcm.shape[1]:
                raise ValueError(f"HOCOMOCO motif {item['name']} length disagrees with PCM.")
            if pcm.shape[1] < 6:
                continue

            info = item["masterlist_info"]
            species = info.get("species", {})
            species_info = species.get("HUMAN") or species.get("MOUSE") or {}
            pfm = mimosa.pcm_to_pfm(pcm, pseudocount=0.25)
            model = mimosa.pwm_from_pfm(pfm, background=0.25, name=item["name"])
            records.append(
                _record(
                    model,
                    item["name"],
                    species_info.get("gene_symbol", "NA"),
                    info.get("tfclass_class", "NA"),
                    info.get("tfclass_family", "NA"),
                )
            )
    return records


def _jaspar_pcm(counts, motif_id: str) -> np.ndarray:
    try:
        rows = [counts[base] if base in counts else counts[base.lower()] for base in DNA_ALPHABET]
    except (KeyError, TypeError) as exc:
        raise ValueError(f"JASPAR motif {motif_id} has no complete A/C/G/T mapping.") from exc
    pcm = np.asarray(rows, dtype=np.float32)
    if pcm.ndim != 2 or pcm.shape[0] != 4 or pcm.shape[1] < 1:
        raise ValueError(f"JASPAR motif {motif_id} PCM must have shape (4, length).")
    return pcm


def _jaspar_records(taxon: str) -> list[dict]:
    try:
        database = jaspardb(release=JASPAR_RELEASE)
        motifs = database.fetch_motifs(collection="CORE", tax_group=[taxon], min_length=6)
    except Exception as exc:
        raise RuntimeError(f"Failed to load {JASPAR_RELEASE} for taxon {taxon}: {exc}") from exc

    records = []
    for motif in motifs:
        motif_id = str(motif.matrix_id)
        pcm = _jaspar_pcm(motif.counts, motif_id)
        pfm = mimosa.pcm_to_pfm(pcm, pseudocount=0.25)
        model = mimosa.pwm_from_pfm(pfm, background=0.25, name=motif_id)
        records.append(_record(model, motif_id, motif.name, motif.tf_class, motif.tf_family))
    return records


def read_motifs_from_db(motif_db: str, taxon: str) -> list[dict]:
    """Load database motifs using the single record contract."""
    if motif_db not in {"jaspar", "hocomoco"}:
        raise ValueError(f"Unknown motif database: {motif_db!r}.")
    records = _hocomoco_records() if motif_db == "hocomoco" else _jaspar_records(taxon)
    logger.info("Loaded %d motifs from %s.", len(records), motif_db)
    return records


def read_model_records(paths) -> list[dict]:
    """Load one Mimosa model per path for the optional file-model interface."""
    records = []
    for path in paths:
        model = mimosa.read_model(path, format="auto")
        records.append(_record(model, model.name, "NA", "NA", "NA"))
    return records
