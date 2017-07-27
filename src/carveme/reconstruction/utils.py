import pandas as pd

from carveme.reconstruction.media import medium_to_constraints
from framed import FBA


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


def add_biomass_equation(model, equation, label=None, objective=1.0):
    r_id = 'Growth_' + label if label else 'Growth'
    biomass_rxn = '{}: {} [0, 1000] @{}'.format(r_id, equation, objective)
    model.add_reaction_from_str(biomass_rxn)

    return r_id


def add_maintenance_atp(model, lb=0, ub=1000):
    rxn_str = 'R_ATPM: M_atp_c + M_h2o_c --> M_adp_c + M_h_c + M_pi_c [{}, {}]'.format(lb, ub)
    model.add_reaction_from_str(rxn_str)


def simulate_biomass_vs_medium(model, biomass_dict, media_db, constraints=None):

    result = {}

    for organism, equation in biomass_dict.items():
        result[organism] = {}
        r_id = add_biomass_equation(model, equation, organism)

        for medium_id, compounds in media_db.items():
            constr = medium_to_constraints(model, compounds)
            if constraints:
                constr.update(constraints)
            sol = FBA(model, constraints=constr, objective={r_id: 1})
            result[organism][medium_id] = sol.fobj

        model.remove_reaction(r_id)

    return pd.DataFrame(result).T


def tab2fasta(inputfile, outputfile, filter_by_model=None):

    data = pd.read_csv(inputfile, sep='\t', header=0)
    data.dropna(subset=['SEQUENCE'], inplace=True)

    if filter_by_model:
        data = data.query('MODEL == "{}"'.format(filter_by_model))

    with open(outputfile, 'w') as f:
        for _, row in data.iterrows():
            f.write('>{}.{}\n{}\n'.format(row['MODEL'], row['ID'], row['SEQUENCE']))
