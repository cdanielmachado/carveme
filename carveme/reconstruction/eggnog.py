import pandas as pd


def split_and_expand(df, col, sep):
    split_col = df[col].str.split(sep).apply(pd.Series, 1).stack()
    split_col.index = split_col.index.droplevel(-1)
    split_col.name = col
    df = df.drop(col, axis=1).join(split_col)
    df.reset_index(drop=True, inplace=True)
    return df


def load_eggnog_data(filename):
    """ Load and parse an eggnog results from eggnog-mapper annotations.

    Args:
        filename (str): input file

    Returns:
        pandas.DataFrame: eggnog data

    """

    columns = ['query_gene', 'score', 'BiGG_gene']
    data = pd.read_csv(filename, comment='#', sep='\t', usecols=[0, 3, 19], names=columns)
    data = data[['query_gene', 'BiGG_gene', 'score']].query("BiGG_gene != '-'")
    data['BiGG_gene'] = data['BiGG_gene'].apply(lambda x: x.split(','))
    data = data.explode('BiGG_gene', ignore_index=True)

    return data



