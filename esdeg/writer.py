import os
import plotly.express as px
import xlsxwriter
import numpy as np
import pandas as pd
import panel as pn
from importlib import resources

def write_table(df, path_to_output):
    if 'jaspar_cluster' in df.columns:
        df = df[['motif_id', 'tf_name', 'tf_class', 'tf_family', 'jaspar_cluster', 'log2(or)', 'log10(pval)', 'log10(adj.pval)', 'adj.pval', 'genes']]
    else:
        df = df[['motif_id', 'tf_name', 'tf_class', 'tf_family', 'log2(or)', 'log10(pval)', 'log10(adj.pval)', 'adj.pval', 'genes']]
    df = df.sort_values(by='log10(adj.pval)')
    df.to_csv(path_to_output, sep='\t', index=False)
    print('All done. Exit')
    pass


def write_table_ann(df, path_to_output):
    if 'jaspar_cluster' in df.columns:
        df = df[['motif_id', 'tf_name', 'tf_class', 'tf_family', 'jaspar_cluster', 'log2(or)', 'me_padj', 'counts', 'lfc', 'de_padj', 'genes']]
    else:
        df = df[['motif_id', 'tf_name', 'tf_class', 'tf_family', 'log2(or)', 'me_padj', 'counts', 'lfc', 'de_padj', 'genes']]
    df = df.sort_values(by='me_padj')
    df.to_csv(path_to_output, sep='\t', index=False)
    print('All done. Exit')
    pass



def write_xlsx(df, taxon, path_to_output):

    if 'jaspar_cluster' in df.columns:
        df = df[['motif_id', 'tf_name', 'tf_class', 'tf_family', 'jaspar_cluster', 'log2(or)', 'log10(pval)', 'log10(adj.pval)', 'adj.pval']]
    else:
        df = df[['motif_id', 'tf_name', 'tf_class', 'tf_family', 'log2(or)', 'log10(pval)', 'log10(adj.pval)', 'adj.pval']]
    df = df.sort_values(by='log10(adj.pval)')

    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writer = pd.ExcelWriter(path_to_output, engine='xlsxwriter')

    # Convert the dataframe to an XlsxWriter Excel object.
    df.to_excel(writer, sheet_name='ESDEG', index=False)

    # Get the xlsxwriter objects from the dataframe writer object.
    workbook  = writer.book
    worksheet = writer.sheets['ESDEG']
    worksheet.set_default_row(27)

    image_row = 1
    image_col = len(df.columns)
    images = list(df.motif_id.map(lambda id: resources.files('esdeg').joinpath(f'logos/{taxon}/{id}.png')))
    for image in images:
        worksheet.insert_image(image_row,
                               image_col,
                               image,
                               {'x_scale': 1.2, 'y_scale': 1.2,
                                'x_offset': 5, 'y_offset': 5,
                                'positioning': 1})
        # positioning = 1 allows move and size with cells (may not always perform as expected)
        image_row += 1
    cell_format = workbook.add_format()
    cell_format.set_bold(True)
    cell_format.set_border(True)
    cell_format.set_align('center')
    cell_format.set_align('top')

    cell_format_2 = workbook.add_format()
    cell_format_2.set_align('center')
    cell_format_2.set_align('vcenter')


    worksheet.set_column(0, 0, 10, cell_format_2)
    worksheet.set_column(1, 1, 14, cell_format_2)
    worksheet.set_column(2, 2, 30, cell_format_2)
    worksheet.set_column(3, 7, 14, cell_format_2)
    # worksheet.set_column(4, 4, 13, cell_format_2)
    # worksheet.set_column(5, 5, 11, cell_format_2)
    # worksheet.set_column(6, 6, 14, cell_format_2)
    worksheet.set_column(image_col, image_col, 48)
    worksheet.set_row_pixels(0, 18, cell_format)
    #worksheet.autofit()
    worksheet.write(0, image_col, 'logo', cell_format)
    writer.close()
    return 0


def write_xlsx_ann(df, path_to_output):

    if 'jaspar_cluster' in df.columns:
        df = df[['motif_id', 'tf_name', 'tf_class', 'tf_family', 'jaspar_cluster', 'log2(or)', 'me_padj', 'counts', 'lfc', 'de_padj']]
    else:
        df = df[['motif_id', 'tf_name', 'tf_class', 'tf_family', 'log2(or)', 'me_padj', 'counts', 'lfc', 'de_padj']]


    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writer = pd.ExcelWriter(path_to_output, engine='xlsxwriter')

    # Convert the dataframe to an XlsxWriter Excel object.
    df.to_excel(writer, sheet_name='ESDEG', index=False)

    # Get the xlsxwriter objects from the dataframe writer object.
    workbook  = writer.book
    worksheet = writer.sheets['ESDEG']
    worksheet.set_default_row(27)

    image_row = 1
    image_col = len(df.columns)

    taxons = os.listdir(resources.files('esdeg').joinpath('logos'))
    number_of_taxons = len(taxons)
    images = []
    for i in df['motif_id']:
        img_flag = True
        for index, taxon in enumerate(taxons):
            img_path = resources.files('esdeg').joinpath(f'logos/{taxon}/{i}.png')
            if os.path.exists(img_path):
                images.append(img_path)
                img_flag = False
                break
        if img_flag:
            print(f'Probably logo of motif {i} doesnt exist')
            images.append('NA')

    for image in images:
        worksheet.insert_image(image_row,
                               image_col,
                               image,
                               {'x_scale': 1.2, 'y_scale': 1.2,
                                'x_offset': 5, 'y_offset': 5,
                                'positioning': 1})
        # positioning = 1 allows move and size with cells (may not always perform as expected)
        image_row += 1
    cell_format = workbook.add_format()
    cell_format.set_bold(True)
    cell_format.set_border(True)
    cell_format.set_align('center')
    cell_format.set_align('top')

    cell_format_2 = workbook.add_format()
    cell_format_2.set_align('center')
    cell_format_2.set_align('vcenter')


    worksheet.set_column(0, 0, 10, cell_format_2)
    worksheet.set_column(1, 1, 14, cell_format_2)
    worksheet.set_column(2, 2, 30, cell_format_2)
    worksheet.set_column(3, 2, 30, cell_format_2)
    worksheet.set_column(4, 7, 14, cell_format_2)
    # worksheet.set_column(4, 4, 13, cell_format_2)
    # worksheet.set_column(5, 5, 11, cell_format_2)
    # worksheet.set_column(6, 6, 14, cell_format_2)
    worksheet.set_column(image_col, image_col, 48)
    worksheet.set_row_pixels(0, 18, cell_format)
    #worksheet.autofit()
    worksheet.write(0, image_col, 'logo', cell_format)
    writer.close()
    return 0


def write_report(df, taxon, path_to_output):
    if 'jaspar_cluster' in df.columns:
        df = df[['motif_id', 'tf_name', 'tf_class', 'tf_family', 'jaspar_cluster', 'log2(or)', 'log10(pval)', 'log10(adj.pval)', 'adj.pval', 'genes']]
    else:
        df = df[['motif_id', 'tf_name', 'tf_class', 'tf_family', 'log2(or)', 'log10(pval)', 'log10(adj.pval)', 'adj.pval', 'genes']]
    df = df.sort_values(by='log10(adj.pval)')
    df = df[df.columns[:-1]]
    df['logo'] = df.motif_id.map(lambda id: f'https://raw.githubusercontent.com/ubercomrade/esdeg/main/esdeg/logos/{taxon}/{id}.png')
    table = pn.widgets.Tabulator(df, formatters={'logo': {'type': 'image'}}, pagination=None, text_align='center')
    table.save(path_to_output)
    pass


def create_picture(df, path_to_output):
    df['-log10(adj.pval)'] = -df['log10(adj.pval)']
    df = df[df['adj.pval'] < 0.05]
    fig = px.scatter(df, y="log2(or)", x="-log10(adj.pval)", color="tf_class", symbol="tf_name",
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
