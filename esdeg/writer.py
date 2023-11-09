import plotly.express as px
import numpy as np
import pandas as pd


def write_table(df, path_to_output):
    #df = df[['motif_id', 'tf_name', 'tf_class', 'log(or)', 'pval', 'adj.pval', 'genes']]
    df = df[['motif_id', 'tf_name', 'tf_class', 'jaspar_cluster', 'log2(or)', 'log10(pval)', 'log10(adj.pval)', 'adj.pval', 'genes_low_thr', 'genes_high_thr']]
    df = df.sort_values(by='adj.pval')
    df.to_csv(path_to_output, sep='\t', index=False)
    print('All done. Exit')
    pass


def create_picture(df, path_to_output):
    df['-log10(adj.pval)'] = -df['log10(adj.pval)']
    df = df[df['adj.pval'] < 0.05]
    fig = px.scatter(df, y="log(or)", x="-log10(adj.pval)", color="tf_class", symbol="tf_name",
                    range_x=[0, np.max(df["-log10(adj.pval)"]) + 1],
                    range_y=[0, np.max(df["log2(or)"]) + .02])
    fig.add_vline(x=-np.log10(0.05), line_width=2.5, line_dash="dash", line_color="green",
                 annotation_text="<b>-log10(0.05)</b>", annotation_position="top left",
                 annotation=dict(font_color="green"))
    fig.update_layout(
        font_family="Courier New",
        font_size=16,
        legend_title="TF class, TF name",
        yaxis_title="Motif relative abundance, DEG vs non-DEG - log2(OR)",
        xaxis_title="-log10(adj.p-value)",
    )
    fig.update_layout(
            annotations=[
            dict(
                text="<b>-log10(0.05)</b>",
                textangle=270,
                font=dict(
                    color="green",
                    size=14
                ))]
    )
    fig.write_html(path_to_output)
    pass
