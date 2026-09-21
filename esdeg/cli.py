"""Command-line interface for enrichment and annotation."""

from __future__ import annotations

import argparse
import logging
import sys

from esdeg.annotation import annotation
from esdeg.functions import esdeg as run
from esdeg.writers import (
    create_picture,
    write_report,
    write_table,
    write_table_ann,
    write_xlsx,
    write_xlsx_ann,
)


def _common_outputs(parser):
    parser.add_argument("-v", "--visualization", default=None, help="Path to an HTML plot.")
    parser.add_argument("-x", "--xlsx", default=None, help="Path to an XLSX table.")


def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="subparser_name", required=True)

    enrichment = subparsers.add_parser("enrichment", help="Run motif enrichment.")
    enrichment.add_argument("promoters", help="Promoter FASTA path.")
    enrichment.add_argument("input", help="DEG CSV/TSV or gene-set path.")
    enrichment.add_argument("output", help="Output enrichment TSV path.")
    enrichment.add_argument("-f", "--format", choices=["deg", "set"], default="deg")
    motifs = enrichment.add_mutually_exclusive_group()
    motifs.add_argument("-m", "--motifs", choices=["jaspar", "hocomoco"], default=None)
    motifs.add_argument("--model", action="append", dest="model_paths", metavar="PATH")
    enrichment.add_argument(
        "-t",
        "--taxon",
        choices=["plants", "vertebrates", "insects", "urochordates", "nematodes", "fungi"],
        default="vertebrates",
    )
    enrichment.add_argument("-r", "--regulated", choices=["all", "up", "down"], default="all")
    enrichment.add_argument("-P", "--pvalue", type=float, default=0.05)
    enrichment.add_argument("-l", "--log2fc-deg", "--log2fc_deg", type=float, default=1.0)
    enrichment.add_argument(
        "-L", "--log2fc-back", "--log2fc_back", type=float, default=0.32192809488736235
    )
    enrichment.add_argument("--match-ratio", type=int, default=5)
    enrichment.add_argument("--n-permutations", type=int, default=10000)
    enrichment.add_argument("--seed", type=int, default=0)
    enrichment.add_argument("-p", "--nproc", type=int, default=4)
    _common_outputs(enrichment)
    enrichment.add_argument("-R", "--report", default=None, help="Path to an HTML table report.")

    annotated = subparsers.add_parser("annotation", help="Add expression annotation.")
    annotated.add_argument("esdeg", help="Enrichment TSV path.")
    annotated.add_argument("deg", help="DEG CSV/TSV path.")
    annotated.add_argument("counts", help="Counts CSV/TSV path.")
    annotated.add_argument("gtf", help="GTF path.")
    annotated.add_argument("output", help="Annotated TSV path.")
    annotated.add_argument("-x", "--xlsx", default=None)
    annotated.add_argument("-b", "--best", action="store_true")
    annotated.add_argument("-f", "--filter", action="store_true", dest="filter_flag")
    annotated.add_argument("-m", "--me-padj", "--me_padj", type=float, default=0.05)
    annotated.add_argument("--min-auc", type=float, default=0.5)
    annotated.add_argument("-d", "--de-padj", "--de_padj", type=float, default=0.05)
    annotated.add_argument("-L", "--log-fc", "--log_fc", type=float, default=1.0)
    annotated.add_argument("-n", "--ncounts", type=float, default=5.0)
    return parser.parse_args(argv)


def main_cli(argv=None):
    logging.basicConfig(stream=sys.stderr, level=logging.INFO, format="%(message)s")
    args = parse_args(argv)
    if args.subparser_name == "enrichment":
        motif_db = args.motifs or "hocomoco"
        taxon = "human" if motif_db == "hocomoco" else args.taxon
        result = run(
            motif_db,
            taxon,
            args.promoters,
            args.input,
            args.nproc,
            type_of_data=args.format,
            condition=args.regulated,
            log2fc_thr_deg=args.log2fc_deg,
            log2fc_thr_background=args.log2fc_back,
            padj_thr=args.pvalue,
            n_permutations=args.n_permutations,
            match_ratio=args.match_ratio,
            seed=args.seed,
            model_paths=args.model_paths,
        )
        if args.visualization is not None:
            create_picture(result, args.visualization)
        if args.report is not None:
            write_report(result, taxon, args.report)
        if args.xlsx is not None:
            write_xlsx(result, taxon, args.xlsx)
        write_table(result, args.output)
    else:
        result = annotation(
            args.deg,
            args.counts,
            args.esdeg,
            args.gtf,
            filter_flag=args.filter_flag,
            best_flag=args.best,
            me_padj_thr=args.me_padj,
            de_padj_thr=args.de_padj,
            lfc_thr=args.log_fc,
            counts_filter=args.ncounts,
            min_auc=args.min_auc,
        )
        if args.xlsx is not None:
            write_xlsx_ann(result, "vertebrates", args.xlsx)
        write_table_ann(result, args.output)


if __name__ == "__main__":
    main_cli()
