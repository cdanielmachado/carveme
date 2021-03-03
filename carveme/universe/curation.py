from carveme.reconstruction.utils import create_exchange_reactions, add_biomass_equation, add_maintenance_atp
from reframed import save_cbmodel, simplify
from reframed.core.transformation import disconnected_metabolites
from reframed.core.elements import parse_formula
from .annotate import compute_missing_formulae, annotate_charges

import pandas as pd


def remove_compartments(model, taxa):

    print(f'Removing compartments that do not belong to {taxa}..')

    if taxa == 'cyanobacteria':
        compartments = {'C_c', 'C_p', 'C_u', 'C_e'}
    else:
        compartments = {'C_c', 'C_p', 'C_e'}

    to_remove = set(model.compartments.keys()) - compartments
    model.remove_compartments(to_remove)

    print(f'Current model size: {len(model.metabolites)} x {len(model.reactions)}')


def filter_reactions_by_kingdom(model, taxa, model_specific_data, bigg_models):

    print(f'Removing reactions that do not belong to {taxa}..')

    df = pd.merge(model_specific_data, bigg_models)
    df = df.groupby('kingdom').agg({'reaction': set})

    if taxa == 'archaea':
        valid_rxns = df.loc['Archaea', 'reaction'] | df.loc['Bacteria', 'reaction']
    else:
        valid_rxns = df.loc['Bacteria', 'reaction']

    rxns_to_remove = set(model.reactions) - valid_rxns
    model.remove_reactions(rxns_to_remove)

    mets_to_remove = disconnected_metabolites(model)
    model.remove_metabolites(mets_to_remove)

    print(f'Current model size: {len(model.metabolites)} x {len(model.reactions)}')


def elemental_balance(formulas, tol=1e-6):
    """ Calculate elemental balance for a given reaction.

    Args:
        formulas (list): formulas of compounds and respective stoichiometric coefficients (as tuples)
        tol (float): absolute tolerance to consider an element to be balanced (default: 1e-6)

    Returns:
        dict: elements and respective balance (only if element is unbalanced)
    """

    balance = {}
    for formula, coeff in formulas:
        elements = parse_formula(formula)
        for element, number in elements.items():
            if element in balance:
                balance[element] += coeff * number
            else:
                balance[element] = coeff * number

    balance_clean = {}
    for element, coeff in balance.items():
        if abs(coeff) > tol:
            balance_clean[element] = coeff

    return balance_clean


def check_elemental_balance(model, r_id, allow_groups=False, proton_relax=True):
    balanced = False
    formulas = []

    for m_id, coeff in model.reactions[r_id].stoichiometry.items():
        metabolite = model.metabolites[m_id]

        if 'FORMULA' not in metabolite.metadata:
            break

        formula = metabolite.metadata['FORMULA']

        if '*' in formula:
            if allow_groups:
                formula = formula.replace('*', '')
            else:
                break

        formulas.append((formula, coeff))

    else:
        balance = elemental_balance(formulas)

        if balance == {}:
            balanced = True
        elif proton_relax:
            balanced = set(balance.keys()) == {'H'}

    return balanced


def remove_unbalanced_reactions(model):

    print('Removing unbalanced reactions..')

    unbalanced = [r_id for r_id in model.reactions if not check_elemental_balance(model, r_id)]
    model.remove_reactions(unbalanced)

    mets_to_remove = disconnected_metabolites(model)
    model.remove_metabolites(mets_to_remove)

    print(f'Current model size: {len(model.metabolites)} x {len(model.reactions)}')


def curate_universe(model, outputfile, model_specific_data, bigg_models, taxa, biomass_eq):

    """ Curate universal reaction database from initial database dump.

    Args:
        model (CBModel): universal model
        outputfile (str): output SBML file (optional)
        model_specific_data (pandas.DataFrame): model specific data downloaded from BiGG
        bigg_models (pandas.DataFrame): Additional information on BiGG models
        taxa (str): filter by taxa (optional)
        biomass_eq (str): default biomass equation

    """

    print(f'Curating {taxa} universe...')
    print(f'Initial model size: {len(model.metabolites)} x {len(model.reactions)}')

    remove_compartments(model, taxa)

    model_specific_data['reaction'] = model_specific_data['reaction'].apply(lambda x: 'R_' + x)
    filter_reactions_by_kingdom(model, taxa, model_specific_data, bigg_models)

    compute_missing_formulae(model)
    annotate_charges(model)

    remove_unbalanced_reactions(model)

    create_exchange_reactions(model, default_lb=-1000, default_ub=1000)
    add_biomass_equation(model, biomass_eq)
    add_maintenance_atp(model)

    print('Removing blocked reactions and dead-end metabolites...')
    simplify(model)
    print(f'Final model size: {len(model.metabolites)} x {len(model.reactions)}')

    save_cbmodel(model, outputfile)