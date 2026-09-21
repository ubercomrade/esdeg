"""Output writers and the optional visualization."""

from __future__ import annotations

from contextlib import ExitStack
from importlib import resources

import numpy as np
import pandas as pd
import panel as pn
import plotly.express as px


def write_table(df: pd.DataFrame, path_to_output) -> None:
    """Write all columns in the supplied frame to TSV."""
    df.to_csv(path_to_output, sep="\t", index=False)


write_table_ann = write_table


def _find_logo(motif_id: str, taxon: str):
    root = resources.files("esdeg").joinpath("logos")
    preferred = root.joinpath(str(taxon), f"{motif_id}.png")
    if preferred.is_file():
        return preferred
    for directory in root.iterdir():
        candidate = directory.joinpath(f"{motif_id}.png")
        if candidate.is_file():
            return candidate
    return None


def write_xlsx(df: pd.DataFrame, taxon: str, path_to_output) -> None:
    """Write a complete table and available motif logos to XLSX."""
    with pd.ExcelWriter(path_to_output, engine="xlsxwriter") as writer:
        df.to_excel(writer, sheet_name="ESDEG", index=False)
        worksheet = writer.sheets["ESDEG"]
        workbook = writer.book
        header = workbook.add_format(
            {"bold": True, "border": 1, "align": "center", "valign": "top"}
        )
        centered = workbook.add_format({"align": "center", "valign": "vcenter"})
        worksheet.set_row(0, 18, header)
        if len(df.columns):
            worksheet.set_column(0, len(df.columns) - 1, 14, centered)

        image_col = len(df.columns)
        worksheet.write(0, image_col, "logo", header)
        worksheet.set_column(image_col, image_col, 48)
        with ExitStack() as stack:
            for row, motif_id in enumerate(df["motif_id"].astype(str), 1):
                logo = _find_logo(motif_id, taxon)
                if logo is None:
                    continue
                local_logo = stack.enter_context(resources.as_file(logo))
                worksheet.insert_image(
                    row,
                    image_col,
                    str(local_logo),
                    {"x_scale": 1.2, "y_scale": 1.2, "x_offset": 5, "y_offset": 5},
                )
                worksheet.set_row(row, 27)


write_xlsx_ann = write_xlsx


def write_report(df: pd.DataFrame, taxon: str, path_to_output) -> None:
    columns = [
        column
        for column in df.columns
        if column
        in {
            "motif_id",
            "tf_name",
            "tf_class",
            "tf_family",
            "auc_roc",
            "p_value_roc_adj",
            "p_value_prc_adj",
        }
    ]
    report = df.loc[:, columns].copy()
    report["logo"] = report["motif_id"].map(
        lambda motif_id: (
            f"https://raw.githubusercontent.com/ubercomrade/esdeg/main/esdeg/logos/{taxon}/{motif_id}.png"
        )
    )
    table = pn.widgets.Tabulator(
        report, formatters={"logo": {"type": "image"}}, pagination=None, text_align="center"
    )
    table.save(path_to_output)


def create_picture(df: pd.DataFrame, path_to_output, threshold: float = 0.05):
    """Write an AUC versus adjusted-p-value report, including the empty case."""
    if threshold <= 0:
        raise ValueError("threshold must be positive.")
    plot = df.copy()
    pvalues = np.clip(plot["p_value_roc_adj"].astype(float), np.finfo(float).tiny, 1.0)
    plot["-log10(p_value_roc_adj)"] = -np.log10(pvalues)
    significant = plot[plot["p_value_roc_adj"] < threshold]
    if significant.empty:
        significant = pd.DataFrame(
            columns=["-log10(p_value_roc_adj)", "auc_roc", "tf_class", "tf_name"]
        )
    figure = px.scatter(
        significant,
        y="auc_roc",
        x="-log10(p_value_roc_adj)",
        color="tf_class",
        symbol="tf_name",
        range_y=[0, 1],
    )
    figure.add_vline(x=-np.log10(threshold), line_dash="dash", line_color="green")
    if plot[plot["p_value_roc_adj"] < threshold].empty:
        figure.add_annotation(
            text="No significant motifs", xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False
        )
    figure.update_layout(
        font_family="Courier New",
        font_size=16,
        legend_title="TF class, TF name",
        yaxis_title="ROC AUC",
        xaxis_title="-log10(adjusted p-value)",
    )
    figure.write_html(path_to_output)
    return figure
