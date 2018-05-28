from builtins import str
from builtins import zip
from collections import OrderedDict
from warnings import warn

import pandas as pd
from framed import Environment
from framed import CBReaction
from framed.experimental.elements import molecular_weight
import numpy as np


def create_exchange_reactions(model, default_lb=0, default_ub=None):
    ex_rxns = []

    for m_id, met in model.metabolites.items():
        if m_id.endswith('_e'):
            lb = str(default_lb) if default_lb is not None else ''
            ub = str(default_ub) if default_ub is not None else ''
            rxn_str = 'R_EX_{}: {} <-> [{},{}]'.format(m_id[2:], m_id, lb, ub)
            r_id = model.add_reaction_from_str(rxn_str, clear_tmp=False)
            model.reactions[r_id].is_exchange = True
            ex_rxns.append(r_id)

    model._clear_temp()

    return ex_rxns


def set_exchange_bounds(model, lb, ub):
    for r_id in model.reactions:
        if r_id.startswith('R_EX_'):
            model.set_flux_bounds(r_id, lb, ub)


def create_sink_reactions(model, metabolites):
    for m_id in metabolites:
        if m_id in model.metabolites:
            rxn_str = 'R_sink_{}: {} --> '.format(m_id[2:], m_id)
            r_id = model.add_reaction_from_str(rxn_str, clear_tmp=False)
            model.reactions[r_id].is_sink = True
    model._clear_temp()


def add_maintenance_atp(model, lb=0, ub=1000):
    rxn_str = 'R_ATPM: M_atp_c + M_h2o_c --> M_adp_c + M_h_c + M_pi_c [{}, {}]'.format(lb, ub)
    model.add_reaction_from_str(rxn_str)


def tab2fasta(inputfile, outputfile, filter_by_model=None):

    data = pd.read_csv(inputfile, sep='\t', header=0)
    data.dropna(subset=['SEQUENCE'], inplace=True)

    if filter_by_model:
        data = data.query('MODEL == "{}"'.format(filter_by_model))

    with open(outputfile, 'w') as f:
        for _, row in data.iterrows():
            f.write('>{}.{}\n{}\n'.format(row['MODEL'], row['ID'], row['SEQUENCE']))


def load_media_db(filename, sep='\t', medium_col='medium', compound_col='compound'):

    data = pd.read_csv(filename, sep=sep)
    media_db = data[[medium_col, compound_col]].groupby(medium_col).agg(lambda x: list(x))

    return media_db[compound_col].to_dict()


def medium_to_constraints(model, compounds, max_uptake=10, inplace=False, verbose=False, exchange_format=None):

    if not exchange_format:
        exchange_format = "'R_EX_{}_e'"

    env = Environment.from_compounds(compounds, max_uptake=max_uptake, exchange_format=exchange_format)
    return env.apply(model, inplace=inplace, warning=verbose)


def load_biomass_db(filename, sep='\t', normalize_weight=False, model=None):
    data = pd.read_csv(filename, sep=sep)
    data.dropna(subset=['bigg_id', 'comp'], inplace=True)
    data.sort_values(by=['bigg_id', 'comp'], inplace=True)

    biomass_db = {}
    for column in data.columns:
        if column.startswith('@'):
            col_slice = data[['bigg_id', 'comp', column]].dropna().values
            biomass_db[column[1:]] = OrderedDict(('M_{}_{}'.format(x, y), z) for x, y, z in col_slice)

    if normalize_weight:
        if model is None:
            raise RuntimeError('To normalize the biomass weight please provide a model with metabolite formulas.')

        for biomass_id, coeffs in biomass_db.items():
            normalize_coeffs(biomass_id, coeffs, model)

    return biomass_db


def biomass_weight(biomass_id, coeffs, model):
    bio_weight = 0
    for m_id, coeff in coeffs.items():
        metabolite = model.metabolites[m_id]
        if 'FORMULA' in metabolite.metadata:
            formulae = metabolite.metadata['FORMULA'].split(';')
            met_weight = np.mean([molecular_weight(formula) for formula in formulae])
            contribution = -coeff * met_weight
            bio_weight += contribution
#            print '\t'.join([biomass_id, m_id, str(met_weight), str(coeff), str(contribution)])
        else:
            warn('Unable to normalize {} due to missing formula for {}:'.format(biomass_id, m_id))
            break

    return bio_weight


def normalize_coeffs(biomass_id, coeffs, model):
    bio_weight = biomass_weight(biomass_id, coeffs, model)

    if bio_weight > 0:
        for x, val in coeffs.items():
            coeffs[x] = val * 1000.0 / bio_weight


def add_biomass_equation(model, stoichiometry, label=None):
    r_id = 'Growth_' + label if label else 'Growth'
    name = 'Biomass reaction'
    reaction = CBReaction(r_id, name=name, reversible=False, stoichiometry=stoichiometry, objective=1.0)
    model.add_reaction(reaction)


def load_soft_constraints(filename):
    df = pd.read_csv(filename, sep='\t', header=None)
    return dict(zip(df[0], df[1]))


def load_hard_constraints(filename):
    df = pd.read_csv(filename, sep='\t', header=None)
    return dict(zip(df[0], zip(df[1], df[2])))

