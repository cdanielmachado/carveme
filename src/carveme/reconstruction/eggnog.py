import pandas as pd


def split_and_expand(df, col, sep):
    split_col = df[col].str.split(sep).apply(pd.Series, 1).stack()
    split_col.index = split_col.index.droplevel(-1)
    split_col.name = col
    df = df.drop(col, axis=1).join(split_col)
    df.reset_index(drop=True, inplace=True)
    return df


def load_eggnog_data(filename, drop_unannotated=True, drop_unused_cols=True):
    """ Load and parse an eggnog results for new eggnog-mapper version.

    Args:
        filename (str): input file
        drop_unannotated (bool): remove entries without BiGG annotation (default: True)
        drop_unused_cols (bool): remove columns not used for model carving (default: True)

    Returns:
        pandas.DataFrame: eggnog data

    """
    columns = ['query_gene', 'seed_eggNOG_ortholog', 'evalue', 'score',
               'predicted_gene_name', 'GO_terms', 'KEGG_pathways', 'Annotation_tax_scope', 'OGs',
               'bestOG_evalue_score', 'COG_cat', 'eggNOG_annotation', 'BiGG_gene']

    data = pd.read_csv(filename, comment='#', sep='\t', names=columns)

    if drop_unannotated:
        data.dropna(subset=['BiGG_gene'], inplace=True)

    if drop_unused_cols:
        data = data[['query_gene', 'BiGG_gene', 'score']]

    data = split_and_expand(data, 'BiGG_gene', ',')

    return data



