import pandas as pd

from framed import Environment


def load_media_db(filename, sep='\t', medium_col='medium', compound_col='compound'):

    data = pd.read_csv(filename, sep=sep)
    media_db = data[[medium_col, compound_col]].groupby(medium_col).agg(lambda x: list(x))

    return media_db[compound_col].to_dict()


def medium_to_constraints(model, compounds, max_uptake=10, inplace=False, verbose=False):

    env = Environment.from_compounds(compounds, max_uptake=max_uptake)
    return env.apply(model, inplace=inplace, warning=verbose)